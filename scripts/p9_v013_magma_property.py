#!/usr/bin/env python3
# Paper 9 | v0.0.13 | G4 check 2: MAGMA gene-property analysis (continuous cell-type specificity,
#                    conditioned on average expression, one-sided) for GCST90205183, then MAGMA gene
#                    analysis + gene-property for FinnGen traits (low-heat, one at a time).
#                    Also: HECA h5ad structure inspection (no matrix loaded).
# Runs in ~/Desktop/paper9/.venv. Rules: CHANGELOG v0.2.5 (fixed before this run).

import datetime as dt
import gzip
import json
import os
import re
import subprocess
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
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
G1000 = TOOLS / "magma" / "ref_eur" / "g1000_eur"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
ANNOT = MAGMA_OUT / "annot_b37_35_10.genes.annot"
FINNGEN_DIR = RAW / "finngen"
FINNGEN_MANIFEST = FINNGEN_DIR / "finngen_R13_manifest.tsv"
GS_DIR = PROC / "genesets" / "gse266265_v007"
COVAR_DIR = PROC / "magma_covar"
HECA_FILE = RAW / "heca" / "endometriumAtlasV2_cells_with_counts.h5ad"
VENV_PY = REPO / ".venv" / "bin" / "python"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "magma_property" / ("p9_v013_magma_property_" + TS)
LOG = REPO / "logs" / ("p9_v013_magma_property_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

G4_PRIMARY = ["dS1-myofibroblast", "dS2", "dS3", "dS4"]
FAMILIES = ["A_refined_CTL", "B_main_phase3_CTL"]
# FinnGen traits for MAGMA gene analysis, in run order
FG_TRAITS = [
    ("ENDO_OVARY", "N14_ENDOMETRIOSIS_OVARY"),
    ("ENDO_PERITONEUM", "N14_ENDOMETRIOSIS_PELVICPERITONEUM"),
    ("ENDO_DEEP", "N14_ENDOMETRIOSIS_DEEP"),
    ("HYPERTENSION", "I9_HYPTENS"),
    ("ENDO_ALL", "N14_ENDOMETRIOSIS"),
    ("LEIOMYOMA", "CD2_BENIGN_LEIOMYOMA_UTERI"),
]
VASCULAR_EXPECTED = ["VSMC", "Prv-STEAP4", "Prv-MYH11", "Prv_VSMC secretory", "Prv_VSMC secretory_late", "EC-artery"]

for d in (OUT, LOG.parent, COVAR_DIR, MAGMA_OUT):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V013_MAGMA_PROPERTY", "timestamp_utc": TS, "gates": {}, "results": {}}


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def sh(cmd, timeout=None):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=ENV1)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"


def safe(name):
    s = name.replace("$\\Phi$", "Phi").replace("$", "").replace("{", "").replace("}", "").replace("\\", "")
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")


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


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


# ----------------------------------------------------------------------------
# Covariate files: continuous specificity + Average (log2 mean CPM + 1)
# ----------------------------------------------------------------------------
def entrez_map():
    sym2, dup = {}, set()
    for line in open(GENE_LOC_B37):
        t = line.split()
        if len(t) >= 6 and t[1] in [str(c) for c in range(1, 23)]:
            if t[5] in sym2:
                dup.add(t[5])
            sym2[t[5]] = t[0]
    for d in dup:
        sym2.pop(d, None)
    return sym2


def build_covar(family, sym2):
    out = COVAR_DIR / (family + ".covar.txt")
    names = COVAR_DIR / (family + ".names.csv")
    S = pd.read_csv(GS_DIR / family / "specificity_by_type.csv.gz", index_col=0)
    M = pd.read_csv(GS_DIR / family / "mean_cpm_by_type.csv.gz", index_col=0)
    cols = list(S.columns)
    keys = []
    for c in cols:
        k = "C_" + safe(c)
        while k in keys:
            k += "_x"
        keys.append(k)
    df = S.copy()
    df.columns = keys
    df["Average"] = np.log2(M.loc[S.index].mean(axis=1) + 1.0)
    df.insert(0, "GENE", [sym2.get(g) for g in df.index])
    df = df.dropna(subset=["GENE"])
    df.to_csv(out, sep="\t", index=False, float_format="%.6g")
    pd.DataFrame({"key": keys, "label": cols}).to_csv(names, index=False)
    say("[COVAR] %-18s genes=%d cell types=%d" % (family, len(df), len(keys)))
    return out, dict(zip(keys, cols))


def gene_property(trait, family, covar, keymap):
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    outp = OUT / ("%s__%s" % (trait, family))
    rc, o, e = sh(NICE + [str(MAGMA), "--gene-results", str(raw), "--gene-covar", str(covar),
                          "--model", "condition-hide=Average", "direction=greater", "--out", str(outp)], timeout=3600)
    gsa = Path(str(outp) + ".gsa.out")
    if rc != 0 or not gsa.exists():
        say("[PROPERTY] %s %s FAILED rc=%d %s" % (trait, family, rc, (o + e).strip()[-300:]))
        return None
    rows = [l for l in open(gsa) if not l.startswith("#")]
    hdr = rows[0].split()
    data = [dict(zip(hdr, l.split())) for l in rows[1:] if l.strip()]
    df = pd.DataFrame(data)
    for c in ["NGENES", "BETA", "BETA_STD", "SE", "P"]:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["label"] = df["VARIABLE"].map(keymap)
    df = df.dropna(subset=["label"]).copy()
    df["q_bh"] = bh(df["P"].values)
    df["trait"] = trait
    df["family"] = family
    df = df.sort_values("P")
    say("[PROPERTY_TOP] %s | %s" % (trait, family))
    for _, r in df.head(8).iterrows():
        say("    %-30s beta=%8.4f se=%7.4f p=%.2e q=%.3f" % (str(r["label"])[:30], r["BETA"], r["SE"], r["P"], r["q_bh"]))
    return df


def run_gene_analysis(trait, code, counts):
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    if raw.exists():
        return True
    src = next((p for p in FINNGEN_DIR.glob("*" + code + ".gz") if p.name.endswith("_" + code + ".gz")), None)
    if src is None:
        say("[MAGMA_GENE] %s source missing" % trait)
        return False
    pfile = MAGMA_OUT / (trait + ".pval.txt")
    if not pfile.exists():
        n = 0
        with gzip.open(src, "rt") as f, open(pfile, "w") as o:
            hdr = f.readline().rstrip("\n").split("\t")
            ix = {h: i for i, h in enumerate(hdr)}
            o.write("SNP\tP\n")
            for line in f:
                t = line.rstrip("\n").split("\t")
                rs = t[ix["rsids"]].split(",")[0]
                if rs.startswith("rs"):
                    o.write(rs + "\t" + t[ix["pval"]] + "\n")
                    n += 1
        say("[MAGMA_PFILE] %s %d SNPs" % (trait, n))
    ncas, ncon = counts[code]
    t0 = time.time()
    say("[MAGMA_GENE] %s running (N=%d; 20-40 min, quiet)" % (trait, ncas + ncon))
    rc, o, e = sh(NICE + [str(MAGMA), "--bfile", str(G1000), "--pval", str(pfile), "N=%d" % (ncas + ncon),
                          "--gene-annot", str(ANNOT), "--out", str(MAGMA_OUT / (trait + "_35_10"))], timeout=6 * 3600)
    say("[MAGMA_GENE] %s rc=%d %.0f min" % (trait, rc, (time.time() - t0) / 60))
    return raw.exists()


def manifest_counts():
    out = {}
    lines = FINNGEN_MANIFEST.read_text().splitlines()
    ix = {h: i for i, h in enumerate(lines[0].split("\t"))}
    for line in lines[1:]:
        t = line.split("\t")
        out[t[0]] = (int(t[ix["num_cases"]]), int(t[ix["num_controls"]]))
    return out


# ----------------------------------------------------------------------------
# HECA structure (obs only)
# ----------------------------------------------------------------------------
H5AD_INSPECT = r"""
import sys, h5py, numpy as np
f = h5py.File(sys.argv[1], "r")
print("[H5AD] keys:", list(f.keys()))
X = f["X"]
if isinstance(X, h5py.Group):
    print("[H5AD] X:", dict(X.attrs).get("encoding-type"), "shape", X.attrs.get("shape"), "nnz", X["data"].shape[0], "dtype", X["data"].dtype)
    d = X["data"][:2000000]
    print("[H5AD] X integer-valued:", bool(np.all(np.mod(d, 1) == 0)), "max", float(d.max()))
if "layers" in f:
    for k in f["layers"].keys():
        g = f["layers"][k]
        if isinstance(g, h5py.Group) and "data" in g:
            d = g["data"][:2000000]
            print("[H5AD] layer", k, "integer-valued:", bool(np.all(np.mod(d, 1) == 0)), "max", float(d.max()))
if "raw" in f:
    print("[H5AD] raw keys:", list(f["raw"].keys()))
obs = f["obs"]
oldcats = obs["__categories"] if "__categories" in obs else None
cols = [k for k in obs.keys() if k not in ("_index", "__categories")]
print("[H5AD] obs columns:", len(cols))
for k in cols:
    g = obs[k]
    cats = None
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in g["categories"][:]]
        codes = g["codes"][:]
    elif oldcats is not None and k in oldcats:
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in oldcats[k][:]]
        codes = g[:]
    if cats is not None:
        cnt = np.bincount(codes[codes >= 0], minlength=len(cats))
        pairs = sorted(zip(cats, cnt), key=lambda x: -x[1])
        print("[OBS] %-26s n=%-4d %s" % (k, len(cats), " | ".join("%s=%d" % (c, n) for c, n in pairs[:60])[:1500]))
    else:
        ds = g if not isinstance(g, h5py.Group) else None
        print("[OBS] %-26s %s" % (k, ds.dtype if ds is not None else "group"))
var = f["var"]
print("[H5AD] var columns:", list(var.keys())[:20])
idx = var[var.attrs.get("_index", "_index")][:5]
print("[H5AD] var index sample:", [x.decode() if isinstance(x, bytes) else x for x in idx])
"""


def heca_inspect():
    say("")
    say("-" * 72)
    say("[SECTION] HECA STRUCTURE (obs only)")
    if not HECA_FILE.exists():
        say("[HECA] file missing")
        return
    p = subprocess.run([str(VENV_PY), "-c", H5AD_INSPECT, str(HECA_FILE)], capture_output=True, text=True, timeout=3600)
    for line in (p.stdout + p.stderr).strip().splitlines():
        say(line[:1600])


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.13 | MAGMA GENE-PROPERTY (G4 check 2) + FinnGen MAGMA | " + TS)
    say("=" * 72)
    heca_inspect()
    sym2 = entrez_map()
    covars = {f: build_covar(f, sym2) for f in FAMILIES}
    allres = []

    say("")
    say("-" * 72)
    say("[SECTION] G4 CHECK 2: GCST90205183 MAGMA gene-property (family A, Holm over dS types)")
    for fam in FAMILIES:
        df = gene_property("GCST90205183", fam, *covars[fam])
        if df is not None:
            allres.append(df)
    a = [d for d in allres if d["family"].iloc[0] == "A_refined_CTL"]
    if a:
        sub = a[0][a[0]["label"].isin(G4_PRIMARY)].copy()
        sub["holm"] = holm(sub["P"].values)
        for _, r in sub.sort_values("P").iterrows():
            say("[G4_CHECK2] %-20s beta=%8.4f se=%7.4f p=%.2e Holm=%.3f" % (r["label"], r["BETA"], r["SE"], r["P"], r["holm"]))
        rep = bool((sub["holm"] < 0.05).any())
        verdict = ("DECIDUAL ENRICHMENT REPLICATED by MAGMA (check 2 passed)" if rep
                   else "DECIDUAL ENRICHMENT NOT REPLICATED by MAGMA (check 2 failed)")
        say("[G4_CHECK2_VERDICT] " + verdict)
        ranks = a[0].reset_index(drop=True)
        for t in G4_PRIMARY:
            if (ranks["label"] == t).any():
                say("[G4_CHECK2_RANK] %-20s rank %d of %d" % (t, int(ranks.index[ranks["label"] == t][0]) + 1, len(ranks)))
        SUM["results"]["G4_check2"] = {"verdict": verdict, "rows": sub[["label", "BETA", "SE", "P", "holm"]].to_dict("records")}
    pd.concat(allres).to_csv(OUT / "magma_property_results.csv", index=False)

    say("")
    say("-" * 72)
    say("[SECTION] FinnGen MAGMA gene analysis + gene-property (exploratory; positive control = HYPERTENSION)")
    counts = manifest_counts()
    for trait, code in FG_TRAITS:
        if not run_gene_analysis(trait, code, counts):
            continue
        for fam in FAMILIES:
            df = gene_property(trait, fam, *covars[fam])
            if df is not None:
                allres.append(df)
                if trait == "HYPERTENSION" and fam == "A_refined_CTL":
                    top5 = df["label"].head(5).tolist()
                    hit = [x for x in top5 if x in VASCULAR_EXPECTED]
                    say("[POSITIVE_CONTROL] hypertension top-5 vascular types: %d (%s)" % (len(hit), ", ".join(hit)))
                    SUM["gates"]["POSITIVE_CONTROL_VASCULAR_TOP5"] = len(hit)
        pd.concat(allres).to_csv(OUT / "magma_property_results.csv", index=False)
    with open(OUT / "p9_v013_magma_property_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    for k, v in SUM["gates"].items():
        say("[GATE] %-32s %s" % (k, v))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] G4 check 2 rule fixed in CHANGELOG v0.2.5 before this run; FinnGen results exploratory")
    say("[DONE] PAPER9_V013_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
