#!/usr/bin/env python3
# Paper 9 | v0.0.44 | (rules CHANGELOG v0.2.34)
#   A  HECA follow-up of L2c: is the |log2FC| venous top hit carried by the stromal loci? (200 random locus subsets)
#   B  External atlas outside the reproductive tract (GTEx v9 snRNA-seq, Eraslan et al. 2022, CELLxGENE): Koller-style
#      marker families (|log2FC| / up / down, size-matched, continuous specificity) per cell type, applied to six GWAS,
#      with positive controls fixed in advance (leiomyoma -> smooth muscle; hypertension -> vascular cells).
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). Download ~2-4 GB once; ~60-90 min.

import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
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

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RAW = REPO / "data" / "raw"
INTERIM = REPO / "data" / "interim"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
GS_P = PROC / "genesets" / "heca_koller_v031" / "P"
LOCUS_DIR = REPO / "results" / "locus_atlas"
CXG_API = "https://api.cellxgene.cziscience.com/curation/v1/collections/a3ffde6c-7ad2-498a-903c-d58e732f7470"
GTEX_DIR = RAW / "cellxgene" / "gtex_eraslan2022"
GTEX_DE = PROC / "genesets" / "gtex_v044_de"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "external_atlas" / ("p9_v044_" + TS)
WORK = INTERIM / "magma_v044"
LOG = REPO / "logs" / ("p9_v044_external_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.34)
SEED = 20261108
N_PERM = 200
DROP_WIN_BP = 500_000
CAP = 400
MIN_CELLS = 50
MIN_DONORS = 3
LFC, PADJ, MIN_SET, N_MATCH = 1.2, 0.01, 20, 400
FAMILIES = ["ABS", "UP", "DOWN", "ABSN", "UPN", "DOWNN", "SPEC"]
TRAITS = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL", "LEIOMYOMA", "HYPERTENSION", "T2D"]
EXPECT = {"LEIOMYOMA": (r"smooth muscle", "smooth muscle cell"),
          "HYPERTENSION": (r"smooth muscle|endothel|pericyte|vascular", "vascular smooth muscle / endothelial / pericyte")}
EXPLORE = {"T2D": (r"adipocyte", "adipocyte (exploratory)")}

for d in (OUT, WORK, LOG.parent, GTEX_DIR, GTEX_DE):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V044_EXTERNAL", "timestamp_utc": TS, "results": {}}


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


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def gene_loc():
    sym2ent, pos = {}, {}
    dup = set()
    for line in open(GENE_LOC_B37):
        t = line.split()
        if len(t) >= 6 and t[1] in [str(c) for c in range(1, 23)]:
            pos[t[0]] = (int(t[1]), int(t[2]), int(t[3]))
            if t[5] in sym2ent:
                dup.add(t[5])
            sym2ent[t[5]] = t[0]
    for d in dup:
        sym2ent.pop(d, None)
    return sym2ent, pos


def parse_gsa(gsa):
    rows = [l for l in open(gsa) if not l.startswith("#") and l.strip()]
    hdr = rows[0].split()
    res = {}
    for l in rows[1:]:
        d = dict(zip(hdr, l.split()))
        try:
            res[d.get("FULL_NAME", d.get("VARIABLE"))] = float(d["P"])
        except (KeyError, ValueError):
            continue
    return res


def magma_sets(raw, sf, tag):
    outp = WORK / tag
    gsa = Path(str(outp) + ".gsa.out")
    for model in (["--model", "condition-hide=ALLGENES"], []):
        if gsa.exists():
            gsa.unlink()
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--set-annot", sf] + model + ["--out", outp], timeout=3600)
        if gsa.exists():
            return parse_gsa(gsa)
    say("    [MAGMA_FAIL] sets %s" % (o + e).strip()[-200:])
    return None


def magma_covar(raw, cf, tag):
    outp = WORK / tag
    gsa = Path(str(outp) + ".gsa.out")
    for key in ("condition-hide=", "condition="):
        if gsa.exists():
            gsa.unlink()
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--gene-covar", cf, "--model", key + "Average", "direction=greater", "--out", outp],
                      timeout=3600)
        if gsa.exists():
            return parse_gsa(gsa)
    say("    [MAGMA_FAIL] covar %s" % (o + e).strip()[-200:])
    return None


def short(label, n=26):
    return str(label)[:n]


# ----------------------------------------------------------------------------
# A. Stromal loci and the |log2FC| venous hit (HECA, v031 P sets)
# ----------------------------------------------------------------------------
def part_a(sym2ent, pos):
    say("")
    say("-" * 72)
    say("[SECTION] A  stromal loci vs the |log2FC| venous top hit (KOLLER_EUR; v042 loci; %d random subsets)" % N_PERM)
    lf = sorted(LOCUS_DIR.glob("p9_v042_*/l1_loci_compartments.csv"))
    mf = GS_P / "manifest.csv"
    raw = MAGMA_OUT / "KOLLER_EUR_35_10.genes.raw"
    if not lf or not mf.exists() or not raw.exists():
        say("[A] inputs missing (v042 loci / v031 sets); skipped")
        return
    L = pd.read_csv(lf[-1])
    K = L[L["trait"] == "KOLLER_EUR"].reset_index(drop=True)
    man = pd.read_csv(mf)
    lab = dict(zip(man["set_key"], man["label"] + "/" + man["status"]))
    sets = {k: [sym2ent[s] for s in open(GS_P / (k + ".txt")).read().split() if s in sym2ent] for k in man["set_key"]}
    track = [k for k in man["set_key"] if lab[k] in ("Venous/CTL", "Arterial/CTL", "Venous/ENDO") and k.startswith("ABS_")]
    track += [k for k in man["set_key"] if lab[k] in ("Ciliated/ENDO", "Luminal/ENDO") and k.startswith("DOWN_")]
    track += [k for k in man["set_key"] if lab[k] == "dStromal_mid/ENDO" and k.startswith("UP_")]
    sf = WORK / "a_sets.txt"

    def genes_of(ids):
        d = set()
        for _, r in K[K["lead_snp"].isin(ids)].iterrows():
            d |= {e for e, (c, s, en) in pos.items() if c == r["chr"] and en >= r["bp"] - DROP_WIN_BP and s <= r["bp"] + DROP_WIN_BP}
        return d

    def run(drop, tag):
        with open(sf, "w") as o:
            for k, g in sets.items():
                o.write(k + " " + " ".join(g if k == "ALLGENES" else [x for x in g if x not in drop]) + "\n")
        return magma_sets(raw, sf, tag)

    base = run(set(), "a_base")
    if not base:
        return
    stro = list(K.loc[K["compartment"] == "stroma", "lead_snp"])
    endo = list(K.loc[K["compartment"] == "endothelium", "lead_snp"])
    rs, re_ = run(genes_of(stro), "a_stroma"), run(genes_of(endo), "a_endo")
    rng = np.random.default_rng(SEED)
    others = [x for x in K["lead_snp"] if x not in stro]
    null = {k: [] for k in track}
    t0 = time.time()
    for i in range(N_PERM):
        pick = list(rng.choice(others, size=len(stro), replace=False))
        r2 = run(genes_of(pick), "a_rand")
        if r2:
            for k in track:
                if np.isfinite(r2.get(k, np.nan)):
                    null[k].append(math.log10(r2[k]) - math.log10(base[k]))
    say("[A] %d stromal loci (%s); %d endothelial loci; %d random %d-locus subsets (%.0f min)" % (
        len(stro), ",".join(K.loc[K["compartment"] == "stroma", "prioritized_gene"].astype(str)), len(endo), N_PERM, len(stro), (time.time() - t0) / 60))
    rows = []
    for k in track:
        ds = math.log10(rs[k]) - math.log10(base[k])
        de = math.log10(re_[k]) - math.log10(base[k])
        nl = np.array(null[k])
        pp = (1 + int((nl >= ds).sum())) / (1 + len(nl))
        rows.append({"set": lab[k], "family": k.split("_")[0], "p_base": base[k], "p_no_stroma": rs[k], "p_no_endothelium": re_[k],
                     "delta_stroma": ds, "delta_endothelium": de, "null_median": float(np.median(nl)), "null_p95": float(np.percentile(nl, 95)), "p_perm": pp})
        say("    %-4s %-18s base %.1e | without stromal loci %.1e (+%.2f) | without endothelial loci %.1e (+%.2f) | random median +%.2f, 95th +%.2f | perm p %.3f" % (
            k.split("_")[0], lab[k], base[k], rs[k], ds, re_[k], de, np.median(nl), np.percentile(nl, 95), pp))
    pd.DataFrame(rows).to_csv(OUT / "a_stromal_loci_vs_sets.csv", index=False)
    v = next((r for r in rows if r["set"] == "Venous/CTL" and r["family"] == "ABS"), None)
    if v:
        ok = v["p_perm"] < 0.05 and v["delta_stroma"] > v["delta_endothelium"]
        say("[H_STRVEN_VERDICT] |log2FC| Venous/CTL: stromal loci raise p by %.2f (perm p %.3f) vs endothelial loci %.2f -> %s" % (
            v["delta_stroma"], v["p_perm"], v["delta_endothelium"], "SUPPORTED" if ok else "NOT SUPPORTED"))
        SUM["results"]["h_strven"] = {"supported": ok, **v}


# ----------------------------------------------------------------------------
# B. External atlas
# ----------------------------------------------------------------------------
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
    v = g[:]
    return np.array([dec(x) for x in v], dtype=object) if v.dtype.kind in "SO" else v


def var_symbols(f, base):
    v = f[base]
    for col in ("feature_name", "gene_symbols", "gene_name"):
        if col in v:
            g = v[col]
            if isinstance(g, h5py.Group) and "categories" in g:
                cats = np.array([dec(c) for c in g["categories"][:]], dtype=object)
                return cats[g["codes"][:]]
            return np.array([dec(x) for x in g[:]], dtype=object)
    if base == "raw/var":
        return var_symbols(f, "var")
    idx = v.attrs.get("_index", "_index")
    return np.array([dec(x) for x in v[idx][:]], dtype=object)


def download():
    have = sorted(GTEX_DIR.glob("*.h5ad"))
    if have and have[0].stat().st_size > 200_000_000:
        say("[DL] present %s %.0f MB" % (have[0].name, have[0].stat().st_size / 1e6))
        return have[0]
    rc, o, e = sh(["curl", "-sSL", "--fail", "--max-time", "120", "-A", "paper9", CXG_API], timeout=200)
    if rc != 0:
        say("[DL] collection API failed: %s" % (e or o)[:200])
        return None
    rec = json.loads(o)
    ds = rec.get("datasets", [])
    for d in ds:
        say("    [DL] dataset %s cells=%s %s" % (d.get("dataset_id"), d.get("cell_count"), str(d.get("title"))[:80]))
    if not ds:
        return None
    best = max(ds, key=lambda d: d.get("cell_count") or 0)
    url = next((a.get("url") for a in best.get("assets", []) if a.get("filetype", "").upper() == "H5AD"), None)
    if not url:
        say("[DL] no H5AD asset")
        return None
    fn = GTEX_DIR / ("gtex_" + best["dataset_id"] + ".h5ad")
    part = fn.with_name(fn.name + ".part")
    t0 = time.time()
    rc, o, e = sh(["curl", "-L", "-sS", "--fail", "--retry", "3", "-C", "-", "-A", "paper9", "-o", part, url], timeout=6 * 3600)
    if rc == 0 and part.exists():
        part.rename(fn)
        say("[DL] OK %s %.0f MB in %.0f min" % (fn.name, fn.stat().st_size / 1e6, (time.time() - t0) / 60))
        return fn
    say("[DL] FAILED rc=%d %s" % (rc, e.strip()[:200]))
    return None


def gtex_de(fn, sym2ent):
    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp
    done = GTEX_DE / "de_all.csv.gz"
    if done.exists():
        DE = pd.read_csv(done)
        CPM = pd.read_csv(GTEX_DE / "mean_cpm_by_type.csv.gz", index_col=0)
        say("[GTEX] cached DE: %d cell types" % DE["group"].nunique())
        return DE, CPM
    t0 = time.time()
    f = h5py.File(fn, "r")
    ctype = obs_col(f, "cell_type")
    donor = obs_col(f, "donor_id") if "donor_id" in f["obs"] else np.full(len(ctype), "d", dtype=object)
    tissue = obs_col(f, "tissue") if "tissue" in f["obs"] else np.full(len(ctype), "-", dtype=object)
    say("[GTEX] cells %d; tissues: %s" % (len(ctype), ", ".join("%s %d" % (k, v) for k, v in pd.Series(tissue).value_counts().items())))
    base = "raw/X" if "raw" in f and "X" in f["raw"] else "X"
    vbase = "raw/var" if base == "raw/X" and "raw/var" in f else "var"
    X = f[base]
    syms = var_symbols(f, vbase)
    gidx = np.array([i for i, s in enumerate(syms) if s in sym2ent])
    sample_vals = np.asarray(X["data"][:2000])
    say("[GTEX] matrix %s; genes %d (mapped %d); integer counts: %s" % (base, len(syms), len(gidx), bool(np.allclose(sample_vals, np.round(sample_vals)))))
    tab = pd.DataFrame({"ct": ctype, "donor": donor})
    stats = tab.groupby("ct").agg(n=("donor", "size"), donors=("donor", "nunique"))
    keep_types = sorted(stats[(stats["n"] >= MIN_CELLS) & (stats["donors"] >= MIN_DONORS)].index)
    say("[GTEX] cell types kept %d of %d (>= %d cells, >= %d donors): %s" % (len(keep_types), len(stats), MIN_CELLS, MIN_DONORS, ", ".join(keep_types)[:900]))
    rng = np.random.default_rng(SEED)
    picks = []
    for t in keep_types:
        idx = np.where(ctype == t)[0]
        if len(idx) > CAP:
            idx = np.sort(rng.choice(idx, CAP, replace=False))
        picks.append(idx)
    want = np.unique(np.concatenate(picks))
    indptr = X["indptr"][:]
    rows = {}
    n = len(indptr) - 1
    for a in range(0, n, 20000):
        b = min(n, a + 20000)
        lo, hi = np.searchsorted(want, a), np.searchsorted(want, b)
        if hi <= lo:
            continue
        s0, e0 = int(indptr[a]), int(indptr[b])
        ind = np.asarray(X["indices"][s0:e0])
        val = np.asarray(X["data"][s0:e0], dtype=np.float32)
        loc = indptr[a:b + 1] - s0
        for i in want[lo:hi]:
            r = i - a
            rows[int(i)] = (ind[loc[r]:loc[r + 1]], val[loc[r]:loc[r + 1]])
    f.close()
    allidx = want
    data = [rows[int(i)][1] for i in allidx]
    ind = [rows[int(i)][0] for i in allidx]
    ptr = np.concatenate([[0], np.cumsum([len(x) for x in ind])])
    M = sp.csr_matrix((np.concatenate(data), np.concatenate(ind), ptr), shape=(len(allidx), len(syms)))[:, gidx]
    del rows
    lib = np.asarray(M.sum(axis=1)).ravel()
    lib[lib == 0] = 1
    Mc = sp.diags(1e6 / lib) @ M
    types = ctype[allidx]
    CPM = pd.DataFrame({t: np.asarray(Mc[np.asarray(types == t)].mean(axis=0)).ravel() for t in keep_types}, index=list(syms[gidx]))
    CPM = CPM.groupby(level=0).mean()
    A = ad.AnnData(X=M)
    A.var_names = list(syms[gidx])
    A.var_names_make_unique()
    A.obs["type"] = pd.Categorical(types)
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    sc.tl.rank_genes_groups(A, "type", method="wilcoxon", reference="rest", use_raw=False, pts=False)
    parts = []
    for t in keep_types:
        df = sc.get.rank_genes_groups_df(A, group=t)
        parts.append(pd.DataFrame({"group": t, "gene": df["names"].astype(str).values, "score": df["scores"].values,
                                   "lfc": df["logfoldchanges"].values, "padj": df["pvals_adj"].values}))
    DE = pd.concat(parts, ignore_index=True)
    DE.to_csv(done, index=False)
    CPM.to_csv(GTEX_DE / "mean_cpm_by_type.csv.gz")
    say("[GTEX] DE for %d types on %d sampled cells (%.0f min)" % (len(keep_types), len(allidx), (time.time() - t0) / 60))
    return DE, CPM


def families(DE, CPM, sym2ent):
    fams, meta = {}, {}
    for g, d in DE.groupby("group", sort=False):
        sig = d["padj"] < PADJ
        sel = {"ABS": d.loc[(d["lfc"].abs() > LFC) & sig, "gene"], "UP": d.loc[(d["lfc"] > LFC) & sig, "gene"],
               "DOWN": d.loc[(d["lfc"] < -LFC) & sig, "gene"],
               "ABSN": d.reindex(d["score"].abs().sort_values(ascending=False).index)["gene"].head(N_MATCH),
               "UPN": d.sort_values("score", ascending=False)["gene"].head(N_MATCH),
               "DOWNN": d.sort_values("score", ascending=True)["gene"].head(N_MATCH)}
        for fam, genes in sel.items():
            gl = sorted({sym2ent[x] for x in genes.astype(str) if x in sym2ent})
            if len(gl) >= MIN_SET:
                key = ("%s_%s" % (fam, re.sub(r"[^A-Za-z0-9]+", "_", g)))[:60]
                fams[key], meta[key] = gl, (fam, g)
    fams["ALLGENES"] = sorted({sym2ent[x] for x in DE["gene"].astype(str) if x in sym2ent})
    C = CPM[(CPM >= 1.0).any(axis=1)]
    C = C.loc[[g for g in C.index if g in sym2ent]]
    S = C.div(C.sum(axis=1), axis=0)
    names = {"C%03d" % (i + 1): c for i, c in enumerate(S.columns)}
    S.columns = list(names)
    S["Average"] = np.log2(C.mean(axis=1) + 1.0)
    S.index = [sym2ent[g] for g in S.index]
    S = S[~S.index.duplicated()]
    return fams, meta, S, names


def part_b(sym2ent):
    say("")
    say("-" * 72)
    say("[SECTION] B  external atlas: GTEx v9 snRNA-seq (Eraslan et al. 2022; 8 non-reproductive tissues)")
    fn = download()
    if fn is None:
        return
    DE, CPM = gtex_de(fn, sym2ent)
    fams, meta, S, names = families(DE, CPM, sym2ent)
    say("[GTEX] sets: %s; continuous covariates %d" % (", ".join("%s %d" % (f, sum(1 for k in meta if meta[k][0] == f)) for f in FAMILIES if f != "SPEC"), len(names)))
    sf, cf = WORK / "gtex_sets.txt", WORK / "gtex_covar.txt"
    with open(sf, "w") as o:
        for k, g in fams.items():
            o.write(k + " " + " ".join(g) + "\n")
    S.to_csv(cf, sep="\t", index_label="GENE")
    obs, rows = {}, []
    for t in TRAITS:
        raw = MAGMA_OUT / (t + "_35_10.genes.raw")
        if not raw.exists():
            say("[B] %s genes.raw missing" % t)
            continue
        rs, rc = magma_sets(raw, sf, "gtex_sets_" + t), magma_covar(raw, cf, "gtex_cov_" + t)
        obs[t] = {}
        say("  %s" % t)
        for fam in FAMILIES:
            if fam == "SPEC":
                items = [(names[k], p) for k, p in (rc or {}).items() if k in names]
            else:
                items = [(meta[k][1], p) for k, p in (rs or {}).items() if k in meta and meta[k][0] == fam]
            if not items:
                continue
            items.sort(key=lambda x: x[1])
            q = bh(np.array([p for _, p in items]))
            r = [(l, p, qq) for (l, p), qq in zip(items, q)]
            obs[t][fam] = r
            rows += [{"trait": t, "family": fam, "rank": i + 1, "label": l, "p": p, "q": qq} for i, (l, p, qq) in enumerate(r)]
            say("    %-6s %s" % (fam, "; ".join("%s %.1e/%.2f" % (short(l), p, qq) for l, p, qq in r[:3])))
    pd.DataFrame(rows).to_csv(OUT / "b_gtex_family_rankings.csv", index=False)
    say("")
    say("[SECTION] B positive controls (expected group ranked first)")
    for t, (rx, desc) in list(EXPECT.items()) + list(EXPLORE.items()):
        if t not in obs:
            continue
        parts, okf = [], []
        for fam in FAMILIES:
            r = obs[t].get(fam, [])
            if not r:
                continue
            ranks = [i + 1 for i, (l, p, q) in enumerate(r) if re.search(rx, str(l), re.I)]
            best = ranks[0] if ranks else None
            if best == 1:
                okf.append(fam)
            parts.append("%s %s (rank %s, q %.2f)" % (fam, "YES" if best == 1 else "no", best or "-", r[best - 1][2] if best else float("nan")))
        tag = "EXPLORATORY" if t in EXPLORE else "B_VERDICT"
        say("[B] %-12s expected %s: %s" % (t, desc, "; ".join(parts)))
        say("[%s] %s: families ranking the expected group first: %s" % (tag, t, ", ".join(okf) or "none"))
        SUM["results"]["b_" + t] = okf


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.44 | stromal loci vs the venous hit; external non-reproductive atlas | " + TS)
    say("=" * 72)
    t0 = time.time()
    sym2ent, pos = gene_loc()
    for name, fn in (("A", lambda: part_a(sym2ent, pos)), ("B", lambda: part_b(sym2ent))):
        try:
            fn()
        except Exception as e:
            import traceback
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v044_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.34 before this run (A formulated after v042)")
    say("[DONE] PAPER9_V044_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
