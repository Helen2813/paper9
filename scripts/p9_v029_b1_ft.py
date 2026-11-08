#!/usr/bin/env python3
# Paper 9 | v0.0.29 | B1 (exploratory, rules CHANGELOG v0.2.19): are endometriosis risk associations preferentially
#   enriched in fallopian-tube secretory / ciliated programs relative to matched endometrial epithelial programs?
#   Reference: CELLxGENE 0436a180-cb44-47ba-8ffa-807b7a468469 (Weigert et al. 2025), HECA control donors.
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). ~1.6 GB download, ~30-50 min.

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
try:
    os.nice(10)
except Exception:
    pass

import datetime as dt
import json
import math
import re
import subprocess
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------
# Hardcoded locations
# ----------------------------------------------------------------------------
HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RAW = REPO / "data" / "raw"
INTERIM = REPO / "data" / "interim"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
HECA_GS = PROC / "genesets" / "heca_v014"
CXG_API = "https://api.cellxgene.cziscience.com/curation/v1/collections/380ade76-e561-49a8-afb2-0f10b39c2c72"
DATASET_ID = "0436a180-cb44-47ba-8ffa-807b7a468469"
FT_FILE = RAW / "cellxgene" / ("ft_weigert2025_" + DATASET_ID + ".h5ad")
PB_DIR = PROC / "ft_pseudobulk_v029"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "b1_ft" / ("p9_v029_" + TS)
WORK = INTERIM / "magma_v029"
LOG = REPO / "logs" / ("p9_v029_b1_ft_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.19)
MIN_CELLS_DONOR_TYPE = 10
MIN_DONORS_TYPE = 3
MIN_CELLS_TYPE = 50
FT_SEC_RE = re.compile(r"secret|SE[0-9]|PAX8", re.I)
FT_CIL_RE = re.compile(r"cilia|\bCE\b", re.I)
HECA_EPI_RE = re.compile(r"SOX9|Lumin|Gland|MUC5B|KRT5")
HECA_CIL_RE = re.compile(r"iliated")
TRAITS = ["KOLLER_EUR", "GCST90205183", "ENDO_OVARY", "ENDO_PERITONEUM", "HYPERTENSION", "T2D", "LEIOMYOMA"]
JK_TRAITS = ["KOLLER_EUR", "GCST90205183", "ENDO_OVARY"]
N_BLOCKS = 50

for d in (OUT, WORK, LOG.parent, FT_FILE.parent, PB_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V029_B1", "timestamp_utc": TS, "results": {}}


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def sh(cmd, timeout=None):
    try:
        p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout, env=ENV1)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"


def entrez_unique():
    sym2, dup = {}, set()
    for line in open(GENE_LOC_B37):
        x = line.split()
        if len(x) >= 6 and x[1] in [str(c) for c in range(1, 23)]:
            if x[5] in sym2:
                dup.add(x[5])
            sym2[x[5]] = x[0]
    for d in dup:
        sym2.pop(d, None)
    return sym2


def gene_positions():
    pos = {}
    for line in open(GENE_LOC_B37):
        x = line.split()
        if len(x) >= 6 and x[1] in [str(c) for c in range(1, 23)]:
            pos[x[0]] = (int(x[1]), int(x[2]))
    return pos


def magma_gp(trait, covf, cond, tag):
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    if not raw.exists():
        return None
    outp = WORK / ("%s__%s" % (trait, tag))
    gsa = Path(str(outp) + ".gsa.out")
    o = e = ""
    for key in ("condition-hide=", "condition="):
        if gsa.exists():
            gsa.unlink()
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--gene-covar", covf, "--model", key + cond,
                              "direction=greater", "--out", outp], timeout=3600)
        if rc == 0 and gsa.exists():
            rows = [l for l in open(gsa) if not l.startswith("#") and l.strip()]
            hdr = rows[0].split()
            res = {}
            for l in rows[1:]:
                d = dict(zip(hdr, l.split()))
                name = d.get("FULL_NAME", d.get("VARIABLE"))
                try:
                    res[name] = {"beta": float(d["BETA"]), "beta_std": float(d.get("BETA_STD", "nan")),
                                 "se": float(d["SE"]), "p": float(d["P"])}
                except (KeyError, ValueError):
                    continue
            return res
    say("    [MAGMA_FAIL] %s %s %s" % (trait, tag, (o + e).strip()[-200:]))
    return None


# ----------------------------------------------------------------------------
# 1. Download and read the atlas (h5py, streaming)
# ----------------------------------------------------------------------------
def download():
    if FT_FILE.exists() and FT_FILE.stat().st_size > 100_000_000:
        say("[DL] present %.0f MB" % (FT_FILE.stat().st_size / 1e6))
        return True
    p = subprocess.run(["curl", "-sSL", "--fail", "--max-time", "120", "-A", "paper9", CXG_API], capture_output=True)
    url = None
    if p.returncode == 0:
        rec = json.loads(p.stdout.decode("utf-8", "replace"))
        for ds in rec.get("datasets", []):
            if ds.get("dataset_id") == DATASET_ID:
                for a in ds.get("assets", []):
                    if a.get("filetype", "").upper() == "H5AD":
                        url = a.get("url")
    if not url:
        say("[DL] asset URL not found via API")
        return False
    part = FT_FILE.with_name(FT_FILE.name + ".part")
    t0 = time.time()
    rc, o, e = sh(["curl", "-L", "-sS", "--fail", "--retry", "3", "-C", "-", "-A", "paper9", "-o", part, url], timeout=4 * 3600)
    if rc == 0 and part.exists():
        part.rename(FT_FILE)
        say("[DL] OK %.0f MB in %.0f min" % (FT_FILE.stat().st_size / 1e6, (time.time() - t0) / 60))
        return True
    say("[DL] FAILED rc=%d %s" % (rc, e.strip()[:200]))
    return False


def dec(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def obs_col(f, name):
    g = f["obs"][name]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = np.array([dec(c) for c in g["categories"][:]], dtype=object)
        codes = g["codes"][:]
        out = np.empty(len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        return out
    if "__categories" in f["obs"] and name in f["obs"]["__categories"]:
        cats = np.array([dec(c) for c in f["obs"]["__categories"][name][:]], dtype=object)
        codes = g[:]
        out = np.empty(len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        return out
    v = g[:]
    return np.array([dec(x) for x in v], dtype=object) if v.dtype.kind in "SO" else v


def var_symbols(f, matrix_group):
    # gene symbols for the columns of the chosen matrix
    base = "raw/var" if matrix_group.startswith("raw") and "raw/var" in f else "var"
    v = f[base]
    for col in ("feature_name", "gene_symbols", "gene_name"):
        if col in v:
            g = v[col]
            if isinstance(g, h5py.Group) and "categories" in g:
                cats = np.array([dec(c) for c in g["categories"][:]], dtype=object)
                return cats[g["codes"][:]]
            return np.array([dec(x) for x in g[:]], dtype=object)
    if base == "raw/var" and "feature_name" in f["var"]:
        return var_symbols(f, "X")
    idx = v.attrs.get("_index", "_index")
    return np.array([dec(x) for x in v[idx][:]], dtype=object)


def build_pseudobulk():
    done = PB_DIR / "pseudobulk_counts.npy"
    if done.exists():
        gt = pd.read_csv(PB_DIR / "groups.csv")
        genes = list(pd.read_csv(PB_DIR / "genes.csv")["symbol"].astype(str))
        say("[PB] cached groups=%d genes=%d" % (len(gt), len(genes)))
        return np.load(done), gt, genes
    f = h5py.File(FT_FILE, "r")
    cols = [k for k in f["obs"].keys() if k not in ("_index", "__categories")]
    say("[OBS] columns: " + ", ".join(cols))
    info = {}
    for c in cols:
        try:
            v = obs_col(f, c)
            if v.dtype == object:
                vc = pd.Series(v).value_counts()
                info[c] = vc
                if len(vc) <= 40:
                    say("    [OBS] %-28s n=%d %s" % (c, len(vc), " | ".join("%s=%d" % (k, n) for k, n in vc.head(25).items())[:600]))
        except Exception:
            continue
    # cell-type column: author annotation with most categories, else cell_type
    cand = [c for c in info if re.search(r"author|annotation|celltype|cell_type_|cluster", c, re.I) and c != "cell_type_ontology_term_id"]
    cand = [c for c in cand if 3 <= len(info[c]) <= 80]
    ctype_col = max(cand, key=lambda c: len(info[c])) if cand else "cell_type"
    say("[CELLTYPE_COLUMN] %s (%d labels)" % (ctype_col, len(info.get(ctype_col, []))))
    ctype = obs_col(f, ctype_col)
    donor = obs_col(f, "donor_id") if "donor_id" in cols else obs_col(f, next(c for c in cols if "donor" in c.lower()))
    keep = np.ones(len(ctype), dtype=bool)
    meno_col = None
    for c, vc in info.items():
        if any(re.search(r"pre.?menopaus", str(k), re.I) for k in vc.index):
            meno_col = c
            break
    if meno_col:
        v = obs_col(f, meno_col)
        keep = np.array([bool(re.search(r"pre.?menopaus", str(x), re.I)) for x in v])
        say("[DONORS] pre-menopausal by column %s: cells %d / %d" % (meno_col, int(keep.sum()), len(keep)))
    elif "development_stage" in cols:
        v = obs_col(f, "development_stage")
        age = np.array([float(m.group(1)) if (m := re.search(r"(\d+)-year", str(x))) else np.nan for x in v])
        if np.isfinite(age).any():
            keep = age < 50
            say("[DONORS] pre-menopausal by age < 50: cells %d / %d" % (int(keep.sum()), len(keep)))
        else:
            say("[DONORS] no menopause / age information: all donors kept")
    mg = "raw/X" if "raw" in f and "X" in f["raw"] else "X"
    X = f[mg]
    genes = var_symbols(f, mg)
    say("[MATRIX] %s; genes=%d; encoding=%s" % (mg, len(genes), dict(X.attrs).get("encoding-type") if isinstance(X, h5py.Group) else "dense"))
    if not isinstance(X, h5py.Group):
        raise RuntimeError("dense matrix not supported")
    data, indices, indptr = X["data"], X["indices"], X["indptr"][:]
    sample = np.asarray(data[:2_000_000])
    say("[MATRIX] integer-valued sample: %s (max %.0f)" % (bool(np.all(np.mod(sample, 1) == 0)), sample.max()))
    keys = pd.Series([("%s||%s" % (d, t)) if k else None for d, t, k in zip(donor, ctype, keep)])
    codes, uniq = pd.factorize(keys)
    ng = len(genes)
    acc = np.zeros(len(uniq) * ng, dtype=np.float64)
    n = len(indptr) - 1
    step = 20000
    t0 = time.time()
    for a in range(0, n, step):
        b = min(n, a + step)
        s, e = int(indptr[a]), int(indptr[b])
        ii = np.asarray(indices[s:e], dtype=np.int64)
        xx = np.asarray(data[s:e], dtype=np.float64)
        g = np.repeat(codes[a:b], np.diff(indptr[a:b + 1]))
        m = g >= 0
        acc += np.bincount(g[m] * ng + ii[m], weights=xx[m], minlength=len(uniq) * ng)
    pb = acc.reshape(len(uniq), ng)
    gt = pd.DataFrame({"group": list(uniq)})
    gt["donor"] = gt["group"].str.split("||", regex=False).str[0]
    gt["type"] = gt["group"].str.split("||", regex=False).str[1]
    gt["n_cells"] = pd.Series(codes[codes >= 0]).value_counts().sort_index().values
    np.save(done, pb)
    gt.to_csv(PB_DIR / "groups.csv", index=False)
    pd.Series(genes).to_csv(PB_DIR / "genes.csv", index=False, header=["symbol"])
    say("[PB] built groups=%d genes=%d (%.0f min)" % (len(gt), ng, (time.time() - t0) / 60))
    return pb, gt, list(genes)


# ----------------------------------------------------------------------------
# 2. Specificity and covariates
# ----------------------------------------------------------------------------
def ft_specificity(pb, gt, genes):
    sub = gt[gt["n_cells"] >= MIN_CELLS_DONOR_TYPE]
    agg = sub.groupby("type").agg(n_donors=("donor", "nunique"), n_cells=("n_cells", "sum"))
    ok = agg[(agg["n_donors"] >= MIN_DONORS_TYPE) & (agg["n_cells"] >= MIN_CELLS_TYPE)].index.tolist()
    say("[FT_TYPES] eligible %d / %d" % (len(ok), len(agg)))
    for t, r in agg.sort_values("n_cells", ascending=False).iterrows():
        comp = "FT_SEC" if FT_SEC_RE.search(t) else ("FT_CIL" if FT_CIL_RE.search(t) else "-")
        say("    [FT_TYPE] %-40s donors=%-3d cells=%-7d %-4s %s" % (t[:40], r["n_donors"], r["n_cells"], "OK" if t in ok else "DROP", comp))
    means = {}
    for t in ok:
        rows = sub.index[sub["type"] == t].to_numpy()
        mats = pb[rows]
        lib = mats.sum(axis=1, keepdims=True)
        means[t] = (mats / np.where(lib > 0, lib, 1) * 1e6).mean(axis=0)
    M = pd.DataFrame(means, index=genes)
    M = M.groupby(level=0).sum()
    M = M[(M >= 1.0).any(axis=1)]
    S = M.div(M.sum(axis=1), axis=0)
    return S, M


def joint_covariates(S_ft, M_ft, sym2):
    S_h = pd.read_csv(HECA_GS / "specificity_by_type.csv.gz", index_col=0)
    M_h = pd.read_csv(HECA_GS / "mean_cpm_by_type.csv.gz", index_col=0).loc[S_h.index, S_h.columns]
    epi = [t for t in S_h.columns if HECA_EPI_RE.search(t) and not HECA_CIL_RE.search(t)]
    cil = [t for t in S_h.columns if HECA_CIL_RE.search(t)]
    sec = [t for t in S_ft.columns if FT_SEC_RE.search(t)]
    ftc = [t for t in S_ft.columns if FT_CIL_RE.search(t) and not FT_SEC_RE.search(t)]
    say("[COMPARTMENTS] FT_SEC=%s | FT_CIL=%s" % (", ".join(sec) or "-", ", ".join(ftc) or "-"))
    say("[COMPARTMENTS] HECA_EPI=%s | HECA_CIL=%s" % (", ".join(epi) or "-", ", ".join(cil) or "-"))
    genes = [g for g in S_ft.index.intersection(S_h.index) if g in sym2]
    J = pd.DataFrame(index=genes)
    J["FT_SEC"] = S_ft.loc[genes, sec].sum(axis=1) if sec else np.nan
    J["FT_CIL"] = S_ft.loc[genes, ftc].sum(axis=1) if ftc else np.nan
    J["HECA_EPI"] = S_h.loc[genes, epi].sum(axis=1) if epi else np.nan
    J["HECA_CIL"] = S_h.loc[genes, cil].sum(axis=1) if cil else np.nan
    J["AVG_FT"] = np.log2(M_ft.loc[genes].mean(axis=1) + 1.0)
    J["AVG_H"] = np.log2(M_h.loc[genes].mean(axis=1) + 1.0)
    J = J.dropna(axis=1, how="all").dropna()
    J.insert(0, "GENE", [sym2[g] for g in J.index])
    fj = WORK / "B1_joint.covar.txt"
    J.to_csv(fj, sep="\t", index=False, float_format="%.6g")
    keys = ["K%02d" % (i + 1) for i in range(S_ft.shape[1])]
    lab = dict(zip(keys, S_ft.columns))
    F = pd.DataFrame(S_ft.loc[genes].values, columns=keys, index=genes)
    F["Average"] = J["AVG_FT"]
    F = F.dropna()
    F.insert(0, "GENE", [sym2[g] for g in F.index])
    ff = WORK / "B1_ft_types.covar.txt"
    F.to_csv(ff, sep="\t", index=False, float_format="%.6g")
    say("[COVARIATES] shared genes=%d" % len(J))
    return fj, ff, lab, list(J.columns[1:])


# ----------------------------------------------------------------------------
# 3. Tests
# ----------------------------------------------------------------------------
def jk(full, parts):
    parts = np.asarray([x for x in parts if np.isfinite(x)])
    B = len(parts)
    se = math.sqrt((B - 1) / B * np.sum((parts - parts.mean()) ** 2)) if B > 1 else float("nan")
    return se


def tests(fj, ff, lab, jcols):
    say("")
    say("-" * 72)
    say("[SECTION] B1 tests (MAGMA gene-property, direction=greater; primary KOLLER_EUR)")
    rows = []
    for trait in TRAITS:
        r = magma_gp(trait, ff, "Average", "B1_types")
        if r:
            df = pd.DataFrame([{"type": lab[k], **r[k]} for k in lab if k in r]).sort_values("p")
            say("    [FT_TOP] %-15s %s" % (trait[:15], "; ".join("%s p=%.1e" % (t[:22], p) for t, p in zip(df["type"].head(5), df["p"].head(5)))))
            df.assign(trait=trait).to_csv(OUT / ("B1_types_%s.csv" % trait), index=False)
    comps = [c for c in ("FT_SEC", "FT_CIL", "HECA_EPI", "HECA_CIL") if c in jcols]
    pairs = [("FT_SEC", "HECA_EPI"), ("FT_CIL", "HECA_CIL")]
    for trait in TRAITS:
        base = magma_gp(trait, fj, "AVG_FT,AVG_H", "B1_marg")
        if not base:
            continue
        line = []
        for a, b in pairs:
            if a in comps and b in comps:
                ca = magma_gp(trait, fj, "AVG_FT,AVG_H," + b, "B1_%s_given_%s" % (a, b))
                cb = magma_gp(trait, fj, "AVG_FT,AVG_H," + a, "B1_%s_given_%s" % (b, a))
                pa, pb_ = (ca or {}).get(a, {}).get("p", np.nan), (cb or {}).get(b, {}).get("p", np.nan)
                line.append("%s p=%.3g (|%s %.3g) vs %s p=%.3g (|%s %.3g)" % (a, base[a]["p"], b, pa, b, base[b]["p"], a, pb_))
                rows.append({"trait": trait, "pair": a + "_vs_" + b, "p_a": base[a]["p"], "p_a_cond": pa, "p_b": base[b]["p"], "p_b_cond": pb_,
                             "bstd_a": base[a]["beta_std"], "bstd_b": base[b]["beta_std"]})
        say("    [B1] %-15s %s" % (trait[:15], " || ".join(line)))
    R = pd.DataFrame(rows)
    # jackknife of standardized differences (marginal model)
    pos = gene_positions()
    J = pd.read_csv(fj, sep="\t")
    J["GENE"] = J["GENE"].astype(str)
    g = sorted([x for x in J["GENE"] if x in pos], key=lambda x: pos[x])
    blocks = [set(np.array(g)[i]) for i in np.array_split(np.arange(len(g)), N_BLOCKS)]
    for trait in JK_TRAITS:
        parts = {p: [] for p in pairs}
        for bl in blocks:
            f = WORK / "B1_blk.covar.txt"
            J[~J["GENE"].isin(bl)].to_csv(f, sep="\t", index=False, float_format="%.6g")
            r = magma_gp(trait, f, "AVG_FT,AVG_H", "B1_blk")
            for a, b in pairs:
                if r and a in r and b in r:
                    parts[(a, b)].append(r[a]["beta_std"] - r[b]["beta_std"])
        for a, b in pairs:
            s = R[(R["trait"] == trait) & (R["pair"] == a + "_vs_" + b)]
            if not len(s) or not parts[(a, b)]:
                continue
            full = float(s["bstd_a"].iloc[0] - s["bstd_b"].iloc[0])
            se = jk(full, parts[(a, b)])
            z = full / se if se and se > 0 else np.nan
            p1 = 0.5 * math.erfc(z / math.sqrt(2)) if np.isfinite(z) else np.nan
            say("    [B1_JK] %-15s BETA_STD(%s) - BETA_STD(%s) = %+.4f (jk se %.4f, z %+.2f, one-sided p %.3f)" % (trait[:15], a, b, full, se, z, p1))
            R.loc[s.index, "jk_diff"] = full
            R.loc[s.index, "jk_p"] = p1
    R.to_csv(OUT / "B1_compartments.csv", index=False)
    # Fixed reading (v0.2.19)
    for a, b in pairs:
        k = R[(R["trait"] == "KOLLER_EUR") & (R["pair"] == a + "_vs_" + b)]
        h = R[(R["trait"] == "HYPERTENSION") & (R["pair"] == a + "_vs_" + b)]
        t2 = R[(R["trait"] == "T2D") & (R["pair"] == a + "_vs_" + b)]
        if not len(k):
            continue
        ok = (k["p_a_cond"].iloc[0] < 0.05) and (k.get("jk_p", pd.Series([np.nan])).iloc[0] < 0.05)
        ctrl = (len(h) and h["p_a_cond"].iloc[0] < 0.05) or (len(t2) and t2["p_a_cond"].iloc[0] < 0.05)
        v = "YES" if ok and not ctrl else ("NO (shared with a control trait)" if ok else "NO")
        say("[B1_VERDICT] %s preferential over %s: %s" % (a, b, v))
        SUM["results"]["B1_" + a] = v


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.29 | B1 fallopian-tube programs vs endometrial epithelium (exploratory) | " + TS)
    say("=" * 72)
    t0 = time.time()
    try:
        if download():
            pb, gt, genes = build_pseudobulk()
            S_ft, M_ft = ft_specificity(pb, gt, genes)
            fj, ff, lab, jcols = joint_covariates(S_ft, M_ft, entrez_unique())
            tests(fj, ff, lab, jcols)
    except Exception as e:
        import traceback
        say("[STEP_ERROR] %s: %s" % (type(e).__name__, str(e)[:300]))
        say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v029_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.19 before any fallopian-tube data were read")
    say("[DONE] PAPER9_V029_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
