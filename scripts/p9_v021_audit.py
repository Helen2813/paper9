#!/usr/bin/env python3
# Paper 9 | v0.0.21 | AUDIT (rules CHANGELOG v0.2.12) + runner for v019 and v020 if they have not finished yet.
#   A0 Koller et al. 2026 Zenodo record: file list only (no download)
#   A1 sumstats rows / non-missing     A2 alleles and SNP mapping     A3 LDSC h2 diagnostics
#   A4 GSE266265 count provenance      A5 GSE266265 leave-one-donor-out (MAGMA, dS / stromal types)
#   A6 G1 sensitivity: no filter, >=50 cells, leave-one-patient-out, patient-label permutation
#   A7 cross-trait difference of G1 betas with rho from the LDSC cross-trait intercept
# At the end: a digest of the key lines of the v019 / v020 logs, so only this log needs to be sent.
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA / LDSC single thread, nice 10).

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
try:
    os.nice(10)
except Exception:
    pass

import datetime as dt
import gzip
import json
import math
import re
import subprocess
import sys
import time
from pathlib import Path

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
EXT = REPO / "data" / "external"
TOOLS = REPO / "tools"
VENV_PY = REPO / ".venv" / "bin" / "python"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
LDSC_PY = TOOLS / "ldsc_venv" / "bin" / "python"
LDSC_DIR = TOOLS / "ldsc_py3"
REF_UNPACKED = RAW / "ldsc_ref" / "unpacked"
MUNGED = PROC / "munged"
HM3_DIR = INTERIM / "sumstats_hm3"
MERGE_ALLELES = INTERIM / "ldsc" / "hm3_merge_alleles.tsv"
LDSC_RESULTS = REPO / "results" / "ldsc"
EXPORT = (EXT / "p8_GSE266265_v232c_export").resolve() / "rna_counts_and_metadata"
PB266 = PROC / "gse266265_pseudobulk"
PB213 = PROC / "gse213216_pseudobulk_clusters"
LOGS = REPO / "logs"
SCRIPTS = REPO / "scripts"
ZENODO_API = "https://zenodo.org/api/records/18983492"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "audit" / ("p9_v021_audit_" + TS)
WORK = INTERIM / "magma_v021"
LOG = LOGS / ("p9_v021_audit_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed settings (CHANGELOG v0.2.12 and earlier rules)
CTL_LABEL = "CTL"
PHASE3 = {"proliferative": "proliferative", "periovulatory": "periovulatory",
          "secretory_early": "secretory", "secretory_mid": "secretory", "secretory_late": "secretory"}
G4_A = ["dS1-myofibroblast", "dS2", "dS3", "dS4"]
MIN_CELLS_PER_DONOR_GROUP = 10
MIN_DONORS_PER_TYPE = 3
MIN_CELLS_PER_TYPE = 50
MIN_CPM_ANY_TYPE = 1.0
LINEAGES = ["stromal", "epithelial", "endothelial", "myeloid", "lymphocyte"]
MIN_PATIENTS = 3
N_PERM = 100
PERM_SEED = 20261106
PERM_TIME_LIMIT_MIN = 45

for d in (OUT, WORK, LOGS):
    d.mkdir(parents=True, exist_ok=True)
_fh = None
SUM = {"task": "PAPER9_V021_AUDIT", "timestamp_utc": TS, "results": {}}


def say(msg=""):
    print(msg, flush=True)
    if _fh:
        _fh.write(msg + "\n")
        _fh.flush()


def sh(cmd, timeout=None, cwd=None):
    try:
        p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout, env=ENV1, cwd=cwd)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"


def holm(p):
    p = np.asarray(p, dtype=float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m)
    run = 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (m - r) * p[i]))
        adj[i] = run
    return adj


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


def magma_gp(trait, covf, cond, tag):
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    if not raw.exists():
        return None
    outp = WORK / ("%s__%s" % (trait, tag))
    gsa = Path(str(outp) + ".gsa.out")
    keys = ["condition-hide=", "condition="] if cond else [None]
    for key in keys:
        if gsa.exists():
            gsa.unlink()
        model = ["direction=greater"] if key is None else [key + cond, "direction=greater"]
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--gene-covar", covf, "--model"] + model + ["--out", outp],
                      timeout=3 * 3600)
        if rc == 0 and gsa.exists():
            rows = [l for l in open(gsa) if not l.startswith("#") and l.strip()]
            hdr = rows[0].split()
            res = {}
            for l in rows[1:]:
                d = dict(zip(hdr, l.split()))
                name = d.get("FULL_NAME", d.get("VARIABLE"))
                try:
                    res[name] = {"beta": float(d["BETA"]), "se": float(d["SE"]), "p": float(d["P"])}
                except (KeyError, ValueError):
                    continue
            return res
    say("    [MAGMA_FAIL] %s %s %s" % (trait, tag, (o + e).strip()[-200:]))
    return None


# ----------------------------------------------------------------------------
# 0. Runner: v019 and v020 first (if not finished)
# ----------------------------------------------------------------------------
def finished(prefix, marker):
    logs = sorted(LOGS.glob(prefix + "*.log"), key=lambda p: p.stat().st_mtime)
    for p in logs:
        if marker in p.read_text(errors="ignore"):
            return "done", p
    if logs and time.time() - logs[-1].stat().st_mtime < 3 * 3600:
        return "running?", logs[-1]
    return "missing", None


def run_previous():
    for prefix, marker, script in [("p9_v019_tissue_g1s_", "PAPER9_V019_COMPLETE", "p9_v019_tissue_g1s.py"),
                                   ("p9_v020_g4c5_limits_", "PAPER9_V020_COMPLETE", "p9_v020_g4c5_limits.py")]:
        st, p = finished(prefix, marker)
        if st == "done":
            print("[RUNNER] %s already complete (%s)" % (script, p.name), flush=True)
            continue
        if st == "running?":
            print("[RUNNER] %s has a recent unfinished log (%s): it may still be running." % (script, p.name), flush=True)
            print("[RUNNER] wait until it finishes (or stop it), then start this script again.", flush=True)
            sys.exit(3)
        print("[RUNNER] running %s first ..." % script, flush=True)
        rc = subprocess.run([str(VENV_PY), str(SCRIPTS / script)], cwd=str(REPO)).returncode
        print("[RUNNER] %s exit %d" % (script, rc), flush=True)


# ----------------------------------------------------------------------------
# A0 Koller 2026 summary statistics record
# ----------------------------------------------------------------------------
def audit_zenodo():
    say("")
    say("-" * 72)
    say("[SECTION] A0 Koller et al. 2026 summary statistics (Zenodo 10.5281/zenodo.18983492), listing only")
    p = subprocess.run(["curl", "-sSL", "--fail", "--max-time", "60", "-A", "paper9", ZENODO_API], capture_output=True)
    if p.returncode != 0:
        say("[A0] Zenodo API not reachable (rc=%d)" % p.returncode)
        return
    rec = json.loads(p.stdout.decode("utf-8", "replace"))
    md = rec.get("metadata", {})
    desc = re.sub(r"<[^>]+>", " ", md.get("description", ""))
    say("[A0] title: %s" % md.get("title", "")[:200])
    say("[A0] access=%s license=%s version=%s" % (md.get("access_right"), (md.get("license") or {}).get("id"), md.get("version")))
    for i in range(0, min(len(desc), 1500), 150):
        say("    [A0_DESC] " + " ".join(desc[i:i + 150].split()))
    for f in rec.get("files", [])[:60]:
        say("    [A0_FILE] %-70s %8.1f MB" % (str(f.get("key"))[:70], f.get("size", 0) / 1e6))


# ----------------------------------------------------------------------------
# A1-A3 genetics audit
# ----------------------------------------------------------------------------
def read_munged(label):
    f = MUNGED / (label + ".sumstats.gz")
    if not f.exists():
        return None
    return pd.read_csv(f, sep="\t", compression="gzip")


def audit_sumstats():
    say("")
    say("-" * 72)
    say("[SECTION] A1 sumstats: rows vs non-missing (munge keeps every merge-allele SNP)")
    nmerge = sum(1 for _ in open(MERGE_ALLELES)) - 1 if MERGE_ALLELES.exists() else -1
    say("[A1] merge-alleles list (HapMap3 x 1000G .bim) = %d SNPs" % nmerge)
    data = {}
    rows = []
    for f in sorted(MUNGED.glob("*.sumstats.gz")):
        label = f.name.replace(".sumstats.gz", "")
        d = pd.read_csv(f, sep="\t", compression="gzip")
        data[label] = d
        nn = d["Z"].notna()
        pre = HM3_DIR / (label + ".tsv")
        npre = sum(1 for _ in open(pre)) - 1 if pre.exists() else -1
        chi = float((d.loc[nn, "Z"] ** 2).mean())
        rows.append({"trait": label, "rows": len(d), "nonmissing": int(nn.sum()), "prefilter_rows": npre,
                     "median_N": float(d.loc[nn, "N"].median()), "mean_chi2": chi})
        say("    [A1] %-34s rows=%d nonmissing=%d prefilter=%d medianN=%.0f mean_chi2=%.3f" % (
            label[:34], len(d), int(nn.sum()), npre, d.loc[nn, "N"].median(), chi))
    pd.DataFrame(rows).to_csv(OUT / "A1_sumstats.csv", index=False)
    return data


def audit_alleles(data):
    say("")
    say("-" * 72)
    say("[SECTION] A2 alleles / mapping")
    g = data.get("GCST90205183")
    if g is None:
        say("[A2] GCST90205183 munged file missing")
        return
    g = g.dropna(subset=["Z"])
    for label, d in data.items():
        if label == "GCST90205183":
            continue
        m = g.merge(d.dropna(subset=["Z"]), on="SNP", suffixes=("_g", "_o"))
        same = (m["A1_g"] == m["A1_o"]).mean()
        r = np.corrcoef(m["Z_g"], m["Z_o"])[0, 1]
        sig = m[(m["Z_g"].abs() > 5.45) & (m["Z_o"].abs() > 2)]
        conc = float((np.sign(sig["Z_g"]) == np.sign(sig["Z_o"])).mean()) if len(sig) else float("nan")
        say("    [A2] GCST vs %-28s shared=%d same A1=%.4f r(Z)=%.3f | GWS(GCST)&|Z|>2: n=%d sign concordance=%.3f" % (
            label[:28], len(m), same, r, len(sig), conc))
    top = g.reindex(g["Z"].abs().sort_values(ascending=False).index).head(8)
    pfile = MAGMA_OUT / "GCST90205183.pval.txt"
    pv = {}
    if pfile.exists():
        want = set(top["SNP"])
        with open(pfile) as f:
            next(f)
            for line in f:
                t = line.split()
                if t[0] in want:
                    pv[t[0]] = t[1]
    say("[A2_TOP] GCST90205183 strongest HapMap3 SNPs: munged Z -> p vs MAGMA input p")
    for _, r in top.iterrows():
        p = math.erfc(abs(r["Z"]) / math.sqrt(2))
        say("    %-12s Z=%7.2f p(Z)=%.2e p(MAGMA file)=%s" % (r["SNP"], r["Z"], p, pv.get(r["SNP"], "absent")))
    say("[A2_MAGMA] SNP mapping lines from gene-analysis logs")
    for lf in sorted(MAGMA_OUT.glob("*_35_10.log")):
        txt = lf.read_text(errors="ignore").splitlines()
        keep = [l.strip() for l in txt if re.search(r"\d", l) and re.search(r"SNP|genes", l, re.I)
                and re.search(r"read|found|mapped|containing|valid|skipp|analys", l, re.I)]
        say("    [%s] %s" % (lf.name.replace("_35_10.log", "")[:16], " | ".join(k[:110] for k in keep[:4])))


def audit_h2():
    say("")
    say("-" * 72)
    say("[SECTION] A3 LDSC h2 diagnostics (from v004 logs)")
    logs = sorted(LDSC_RESULTS.glob("p9_v004*/h2_*.log"))
    if not logs:
        say("[A3] no h2 logs found under " + str(LDSC_RESULTS))
        return
    latest = {}
    for lf in logs:
        latest[lf.name] = lf
    rows = []
    for name, lf in sorted(latest.items()):
        t = lf.read_text(errors="ignore")
        def grab(pat):
            m = re.search(pat, t)
            return float(m.group(1)) if m else float("nan")
        r = {"trait": name[3:-4], "h2_obs": grab(r"Total Observed scale h2:\s*([-\d.e]+)"),
             "intercept": grab(r"Intercept:\s*([-\d.e]+)"), "ratio": grab(r"Ratio:\s*([-\d.e]+)"),
             "mean_chi2": grab(r"Mean Chi\^2:\s*([-\d.e]+)"), "lambda_gc": grab(r"Lambda GC:\s*([-\d.e]+)")}
        if "Ratio < 0" in t:
            r["ratio"] = -1.0
        rows.append(r)
        say("    [A3] %-22s h2_obs=%.4f intercept=%.3f ratio=%.3f mean_chi2=%.3f lambdaGC=%.3f %s" % (
            r["trait"][:22], r["h2_obs"], r["intercept"], r["ratio"], r["mean_chi2"], r["lambda_gc"],
            "FLAG ratio>0.3" if r["ratio"] > 0.3 else ""))
    pd.DataFrame(rows).to_csv(OUT / "A3_h2_diagnostics.csv", index=False)


# ----------------------------------------------------------------------------
# A4-A5 GSE266265
# ----------------------------------------------------------------------------
def file_by_prefix(tag):
    for p in EXPORT.iterdir():
        if tag in p.name:
            return p
    raise FileNotFoundError(tag)


def audit_counts():
    say("")
    say("-" * 72)
    say("[SECTION] A4 GSE266265 count provenance (Paper 8 export, read-only)")
    if not EXPORT.exists():
        say("[A4] export not found: " + str(EXPORT))
        return
    files = sorted(EXPORT.parent.rglob("*"))
    small = [f for f in files if f.is_file() and f.suffix.lower() in (".txt", ".md", ".json", ".yaml", ".csv")
             and f.stat().st_size < 6000 and "meta.data_" not in f.name and "Dimnames" not in f.name]
    say("[A4] export files=%d; small text files: %s" % (len(files), ", ".join(f.name for f in small[:10]) or "none"))
    for f in small[:4]:
        lines = f.read_text(errors="ignore").splitlines()[:12]
        for l in lines:
            say("    [A4_TXT %s] %s" % (f.name[:24], l[:150]))
    x = np.load(file_by_prefix("counts_x.npy"), mmap_mode="r")
    p = np.load(file_by_prefix("counts_p.npy"), mmap_mode="r")
    s = np.asarray(x[:5_000_000], dtype=np.float64)
    frac = s - np.floor(s)
    isint = np.isclose(frac, 0) | np.isclose(frac, 1)
    say("[A4] sample of %d nonzero entries: integer=%.4f min=%.4g median=%.4g max=%.4g" % (
        len(s), isint.mean(), s.min(), np.median(s), s.max()))
    h, _ = np.histogram(frac[~isint], bins=10, range=(0, 1))
    say("[A4] fractional parts of non-integer entries (10 bins 0..1): " + " ".join(str(int(v)) for v in h))
    ncell = min(20000, len(p) - 1)
    sums = np.add.reduceat(np.asarray(x[:int(p[ncell])], dtype=np.float64), np.asarray(p[:ncell], dtype=np.int64))
    say("[A4] per-cell sums (first %d cells): median=%.1f, integer-valued sums=%.3f" % (
        ncell, np.median(sums), np.mean(np.isclose(sums, np.round(sums)))))
    names = json.loads(file_by_prefix("meta.data_names.json").read_text())
    for col in ("nCount_RNA", "nCount_originalexp", "nCount_SCT"):
        if col in names:
            idx = names.index(col)
            f = next((q for q in EXPORT.iterdir() if re.match(r"^\d+_meta\.data_%d__%s\.(npy|json)$" % (idx, re.escape(col)), q.name)), None)
            if f is not None:
                v = np.load(f) if f.suffix == ".npy" else np.array(json.loads(f.read_text()), dtype=float)
                ratio = sums / np.asarray(v[:ncell], dtype=float)
                say("[A4] sums / %s: median=%.4f p5=%.4f p95=%.4f (1.000 = counts slot equals the stored totals)" % (
                    col, np.median(ratio), np.percentile(ratio, 5), np.percentile(ratio, 95)))
    say("[A4] metadata columns: " + ", ".join(names[:60]))


def family_tables(pb, gt, genes, sym2, type_cols, drop_donor=None):
    sub = gt[(gt["status"] == CTL_LABEL) & (gt["n_cells"] >= MIN_CELLS_PER_DONOR_GROUP)].copy()
    if drop_donor is not None:
        sub = sub[sub["sample"] != drop_donor]
    sub["type"] = sub[type_cols].astype(str).apply(lambda r: "@".join(r.values), axis=1)
    agg = sub.groupby("type").agg(n_donors=("sample", "nunique"), n_cells=("n_cells", "sum"))
    ok = agg[(agg["n_donors"] >= MIN_DONORS_PER_TYPE) & (agg["n_cells"] >= MIN_CELLS_PER_TYPE)].index.tolist()
    mapped = np.array([g in sym2 for g in genes])
    means = {}
    for t in ok:
        rows = sub.index[sub["type"] == t].to_numpy()
        mats = np.asarray(pb[rows, :], dtype=np.float64)
        lib = mats.sum(axis=1, keepdims=True)
        means[t] = (mats / np.where(lib > 0, lib, 1) * 1e6).mean(axis=0)
    M = pd.DataFrame(means, index=genes)[mapped]
    M = M[(M >= MIN_CPM_ANY_TYPE).any(axis=1)]
    S = M.div(M.sum(axis=1), axis=0)
    return M, S, sub


def audit_loo(sym2):
    say("")
    say("-" * 72)
    say("[SECTION] A5 GSE266265 leave-one-control-donor-out (MAGMA gene-property, GCST90205183)")
    pb = np.load(PB266 / "pseudobulk_counts.npy", mmap_mode="r")
    gt = pd.read_csv(PB266 / "pseudobulk_groups.csv")
    genes = list(pd.read_csv(PB266 / "genes.csv")["symbol"].astype(str))
    gt["phase3"] = gt["phase5"].map(PHASE3)
    out = []
    for fam, cols, pick in [("A_refined_CTL", ["refined"], lambda t: t in G4_A),
                            ("B_main_phase3_CTL", ["main", "phase3"], lambda t: t.lower().startswith("stromal@"))]:
        M, S, sub = family_tables(pb, gt, genes, sym2, cols)
        test = [t for t in S.columns if pick(t)]
        donors = sorted(sub[sub[cols[0]].astype(str).isin([t.split("@")[0] for t in test])]["sample"].unique())
        say("[A5] %s tested types: %s | control donors: %d" % (fam, ", ".join(test), len(donors)))
        for drop in [None] + donors:
            M, S, _ = family_tables(pb, gt, genes, sym2, cols, drop_donor=drop)
            tt = [t for t in test if t in S.columns]
            keys = ["T%d" % i for i in range(len(tt))]
            cov = pd.DataFrame(S[tt].values, columns=keys)
            cov["Average"] = np.log2(M.mean(axis=1).values + 1.0)
            cov.insert(0, "GENE", [sym2[g] for g in S.index])
            covf = WORK / ("A5_%s.covar.txt" % fam)
            cov.to_csv(covf, sep="\t", index=False, float_format="%.6g")
            res = magma_gp("GCST90205183", covf, "Average", "A5_%s" % fam)
            if not res:
                continue
            ps = [res[k]["p"] for k in keys if k in res]
            hs = holm(ps) if ps else []
            for k, t, h in zip(keys, tt, hs):
                if k in res:
                    out.append({"family": fam, "dropped": drop or "none", "type": t, "beta": res[k]["beta"],
                                "p": res[k]["p"], "holm": h})
    df = pd.DataFrame(out)
    df.to_csv(OUT / "A5_leave_one_donor_out.csv", index=False)
    for (fam, t), g in df.groupby(["family", "type"]):
        full = g[g["dropped"] == "none"]
        loo = g[g["dropped"] != "none"]
        if not len(full) or not len(loo):
            continue
        say("    [A5] %-18s %-20s full p=%.3f (Holm %.3f) | LOO p min=%.3f median=%.3f max=%.3f | Holm<0.05 in %d/%d | beta<0 in %d/%d" % (
            fam[:18], t[:20], full["p"].iloc[0], full["holm"].iloc[0], loo["p"].min(), loo["p"].median(), loo["p"].max(),
            int((loo["holm"] < 0.05).sum()), len(loo), int((loo["beta"] < 0).sum()), len(loo)))


# ----------------------------------------------------------------------------
# A6-A7 G1 sensitivity
# ----------------------------------------------------------------------------
def unit_profiles(M, G, lineage, context, agree_only=True, min_cells=10):
    sub = G[(G["lineage_ref"] == lineage) & (G["context"] == context)]
    if agree_only:
        sub = sub[sub["agree"]]
    prof = {}
    for pid, g in sub.groupby("patient"):
        if int(g["n_cells"].sum()) < min_cells:
            continue
        v = M[g.index.to_numpy()].sum(axis=0)
        prof[(pid, context)] = np.log2(v / max(v.sum(), 1) * 1e6 + 1)
    return prof


def write_contrast(A_list, B_list, genes, mapped_sym, sym2, tag):
    A = np.vstack(A_list).mean(axis=0)
    B = np.vstack(B_list).mean(axis=0)
    expressed = (np.maximum(2 ** A - 1, 2 ** B - 1) >= 1.0) & mapped_sym
    D = (A - B)[expressed]
    AV = ((A + B) / 2)[expressed]
    g = np.asarray(genes)[expressed]
    covf = WORK / ("A6_%s.covar.txt" % tag)
    pd.DataFrame({"GENE": [sym2[x] for x in g], "D": D, "DN": -D, "AVG": AV}).to_csv(covf, sep="\t", index=False, float_format="%.6g")
    return covf, int(expressed.sum())


def g1_fit(trait, covf, tag):
    r = magma_gp(trait, covf, "AVG", tag)
    if not r or "D" not in r or "DN" not in r:
        return None
    return {"beta": r["D"]["beta"], "se": r["D"]["se"], "p_greater": r["D"]["p"], "p_less": r["DN"]["p"]}


def audit_g1(sym2, traits):
    say("")
    say("-" * 72)
    say("[SECTION] A6 G1 sensitivity (peritoneal lesion vs eutopic)")
    M = np.load(PB213 / "pseudobulk_counts.npy", mmap_mode="r")
    G = pd.read_csv(PB213 / "pseudobulk_groups.csv").reset_index(drop=True)
    genes = list(pd.read_csv(PB213 / "genes.csv")["symbol"].astype(str))
    mapped = np.array([g in sym2 for g in genes])
    rows = []
    base = {}
    for setting, agree_only, min_cells in [("primary", True, 10), ("no_filter", False, 10), ("min50cells", True, 50)]:
        for lin in LINEAGES:
            pa = unit_profiles(M, G, lin, "peritoneal_lesion", agree_only, min_cells)
            pb = unit_profiles(M, G, lin, "eutopic", agree_only, min_cells)
            if len(pa) < MIN_PATIENTS or len(pb) < MIN_PATIENTS:
                say("    [A6] %-10s %-11s skipped (patients lesion=%d eutopic=%d)" % (setting, lin, len(pa), len(pb)))
                continue
            covf, ng = write_contrast(list(pa.values()), list(pb.values()), genes, mapped, sym2, "%s_%s" % (setting, lin))
            for trait in traits:
                r = g1_fit(trait, covf, "A6_%s_%s" % (setting, lin))
                if r:
                    rows.append({"setting": setting, "lineage": lin, "trait": trait, "n_les": len(pa), "n_eu": len(pb), "genes": ng, **r})
                    if setting == "primary":
                        base[(lin, trait)] = r
    df = pd.DataFrame(rows)
    if len(df):
        for (setting, trait), g in df.groupby(["setting", "trait"]):
            hl = holm(g["p_less"].values)
            say("    [A6] %-10s %-13s " % (setting, trait[:13]) + " ".join(
                "%s:%+.3f(pL=%.3f,H=%.2f,n=%d/%d)" % (l[:5], b, p, h, a, e)
                for l, b, p, h, a, e in zip(g["lineage"], g["beta"], g["p_less"], hl, g["n_les"], g["n_eu"])))
    df.to_csv(OUT / "A6_g1_settings.csv", index=False)

    # (iii) leave-one-unit-out and (iv) permutation, lymphocyte
    lin = "lymphocyte"
    pa = unit_profiles(M, G, lin, "peritoneal_lesion")
    pb = unit_profiles(M, G, lin, "eutopic")
    if len(pa) >= MIN_PATIENTS and len(pb) >= MIN_PATIENTS:
        units = list(pa.items()) + list(pb.items())
        lab = np.array([1] * len(pa) + [0] * len(pb))
        loo = []
        for i in range(len(units)):
            keep = [j for j in range(len(units)) if j != i]
            A = [units[j][1] for j in keep if lab[j] == 1]
            B = [units[j][1] for j in keep if lab[j] == 0]
            if len(A) < MIN_PATIENTS or len(B) < MIN_PATIENTS:
                continue
            covf, _ = write_contrast(A, B, genes, mapped, sym2, "loo")
            for trait in ["GCST90205183", "HYPERTENSION"]:
                r = g1_fit(trait, covf, "A6_loo")
                if r:
                    loo.append({"dropped": "%s|%s" % units[i][0], "trait": trait, **r})
        L = pd.DataFrame(loo)
        L.to_csv(OUT / "A6_lymphocyte_leave_one_out.csv", index=False)
        for trait, g in L.groupby("trait"):
            say("    [A6_LOO] lymphocyte %-13s beta range %.4f..%.4f | p_less range %.2e..%.2e | p_less<0.05 in %d/%d" % (
                trait[:13], g["beta"].min(), g["beta"].max(), g["p_less"].min(), g["p_less"].max(),
                int((g["p_less"] < 0.05).sum()), len(g)))
        rng = np.random.default_rng(PERM_SEED)
        t0 = time.time()
        perm = {"GCST90205183": [], "HYPERTENSION": []}
        for b in range(N_PERM):
            if (time.time() - t0) / 60 > PERM_TIME_LIMIT_MIN:
                break
            pl = rng.permutation(lab)
            A = [units[j][1] for j in range(len(units)) if pl[j] == 1]
            B = [units[j][1] for j in range(len(units)) if pl[j] == 0]
            covf, _ = write_contrast(A, B, genes, mapped, sym2, "perm")
            for trait in perm:
                r = g1_fit(trait, covf, "A6_perm")
                if r:
                    perm[trait].append(r["beta"] / r["se"])
        for trait, ts in perm.items():
            if (lin, trait) in base and ts:
                tobs = base[(lin, trait)]["beta"] / base[(lin, trait)]["se"]
                ts = np.array(ts)
                pe = (1 + np.sum(ts <= tobs)) / (1 + len(ts))
                say("    [A6_PERM] lymphocyte %-13s t_obs=%.2f | permutation p_less=%.3f (B=%d) vs parametric %.4f | perm t sd=%.2f" % (
                    trait[:13], tobs, pe, len(ts), base[(lin, trait)]["p_less"], ts.std()))
                SUM["results"]["A6_perm_" + trait] = {"t_obs": tobs, "p_perm": pe, "B": len(ts)}
        say("    [A6_PERM] time %.0f min" % ((time.time() - t0) / 60))
    return base


def ld_prefixes():
    def find(dirname):
        d = REF_UNPACKED / dirname
        for p in sorted(d.rglob("*.l2.ldscore.gz")) if d.exists() else []:
            m = re.match(r"^(.*\.)1\.l2\.ldscore\.gz$", p.name)
            if m:
                return str(p.parent / m.group(1))
        return None
    return find("LDscore"), find("1000G_Phase3_weights_hm3_no_MHC")


def audit_crosstrait(base, traits):
    say("")
    say("-" * 72)
    say("[SECTION] A7 cross-trait difference of G1 betas (rho = LDSC gcov_int)")
    ldsc_py = LDSC_DIR / "ldsc.py"
    if not ldsc_py.exists():
        c = list((TOOLS / "ldsc_venv").rglob("ldsc.py"))
        ldsc_py = c[0] if c else None
    ref, w = ld_prefixes()
    others = [t for t in traits if t != "GCST90205183" and (MUNGED / (t + ".sumstats.gz")).exists()]
    rho = {}
    if ldsc_py and ref and w and others:
        outp = WORK / "A7_rg"
        rc, o, e = sh([LDSC_PY, ldsc_py, "--rg", ",".join(str(MUNGED / (t + ".sumstats.gz")) for t in ["GCST90205183"] + others),
                       "--ref-ld-chr", ref, "--w-ld-chr", w, "--out", outp], timeout=7200)
        txt = Path(str(outp) + ".log").read_text(errors="ignore") if Path(str(outp) + ".log").exists() else o + e
        if "Summary of Genetic Correlation Results" in txt:
            block = txt.split("Summary of Genetic Correlation Results")[1].strip().splitlines()
            hdr = block[0].split()
            for line in block[1:]:
                t = line.split()
                if len(t) != len(hdr):
                    break
                d = dict(zip(hdr, t))
                p2 = Path(d["p2"]).name.replace(".sumstats.gz", "")
                rho[p2] = float(d.get("gcov_int", "nan"))
                say("    [A7_RG] GCST90205183 ~ %-13s rg=%s se=%s gcov_int=%s" % (p2, d.get("rg"), d.get("se"), d.get("gcov_int")))
        else:
            say("    [A7_RG] LDSC failed rc=%d %s" % (rc, (o + e).strip()[-200:]))
    for t2 in others:
        r = rho.get(t2, 0.0)
        r = 0.0 if not np.isfinite(r) else r
        for lin in LINEAGES:
            a, b = base.get((lin, "GCST90205183")), base.get((lin, t2))
            if not a or not b:
                continue
            se = math.sqrt(max(a["se"] ** 2 + b["se"] ** 2 - 2 * r * a["se"] * b["se"], 1e-12))
            z = (a["beta"] - b["beta"]) / se
            p = math.erfc(abs(z) / math.sqrt(2))
            say("    [A7] %-11s GCST %+.4f vs %-12s %+.4f | diff %+.4f se %.4f z %+.2f p %.3f (rho %.3f)" % (
                lin, a["beta"], t2[:12], b["beta"], a["beta"] - b["beta"], se, z, p, r))


# ----------------------------------------------------------------------------
# Digest of v019 / v020
# ----------------------------------------------------------------------------
DIGEST_TAGS = ("[G6_VERDICT]", "[G6_POSCTRL]", "[G6_TOP]", "[REPRO]", "[S1]", "[S1_REFIT]", "[S2]", "[S2_GENES]",
               "[S2_SHRINK]", "[S2_VERDICT]", "[S3]", "[S4]", "[C5]", "[C5_VERDICT]", "[MDE_", "[FONSECA", "[STEP_ERROR]",
               "[MAGMA_FAIL]", "[DL]", "[GTEX", "[GNOMAD]", "[MAGMA_GENE]", "[GATE]", "[TIME]", "[dS", "[stromal]",
               "[epithelial]", "[endothelial]", "[myeloid]", "[lymphocyte]", "[FILE]")


def digest():
    say("")
    say("-" * 72)
    say("[SECTION] DIGEST of v019 / v020 logs")
    for prefix in ("p9_v019_tissue_g1s_", "p9_v020_g4c5_limits_"):
        logs = sorted(LOGS.glob(prefix + "*.log"), key=lambda p: p.stat().st_mtime)
        if not logs:
            say("[DIGEST] %s: no log" % prefix)
            continue
        lf = logs[-1]
        say("[DIGEST] " + lf.name)
        n = 0
        for line in lf.read_text(errors="ignore").splitlines():
            s = line.strip()
            if s.startswith(DIGEST_TAGS) or (s[:3].strip().isdigit() and "beta=" in s):
                say("  " + line.rstrip()[:200])
                n += 1
                if n > 260:
                    say("  [DIGEST] truncated")
                    break


# ----------------------------------------------------------------------------
def main():
    global _fh
    run_previous()
    _fh = open(LOG, "w", encoding="utf-8")
    say("=" * 72)
    say("PAPER 9 | v0.0.21 | AUDIT (CHANGELOG v0.2.12) | " + TS)
    say("=" * 72)
    t0 = time.time()
    sym2 = entrez_unique()
    traits = [t for t in ["GCST90205183", "HYPERTENSION", "T2D"] if (MAGMA_OUT / (t + "_35_10.genes.raw")).exists()]
    data = {}
    steps = [("A0", audit_zenodo), ("A1", lambda: data.update(audit_sumstats() or {})), ("A2", lambda: audit_alleles(data)),
             ("A3", audit_h2), ("A4", audit_counts), ("A5", lambda: audit_loo(sym2))]
    for name, fn in steps:
        try:
            fn()
        except Exception as e:
            import traceback
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-800:])
    data.clear()
    base = {}
    try:
        base = audit_g1(sym2, traits)
    except Exception as e:
        import traceback
        say("[STEP_ERROR] A6 %s: %s" % (type(e).__name__, str(e)[:300]))
        say(traceback.format_exc()[-800:])
    try:
        audit_crosstrait(base, traits)
    except Exception as e:
        say("[STEP_ERROR] A7 %s: %s" % (type(e).__name__, str(e)[:300]))
    digest()
    with open(OUT / "p9_v021_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] audit items fixed in CHANGELOG v0.2.12 before this run")
    say("[DONE] PAPER9_V021_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
