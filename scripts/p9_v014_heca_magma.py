#!/usr/bin/env python3
# Paper 9 | v0.0.14 | G4 check 3: HECA's own cell types (dStromal early/mid/late) with MAGMA gene-property
#                    (method with a passing positive control). HECA donor-level pseudobulk from the raw
#                    'counts' layer, control donors (pathology C) without hormonal treatment.
#                    Positive controls: hypertension -> vascular/perivascular; leiomyoma -> uSMCs.
# Runs in ~/Desktop/paper9/.venv; reads the h5ad in row chunks (no full matrix in memory).
# Rules: CHANGELOG v0.2.6 (fixed before this run).

import datetime as dt
import json
import os
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
HECA_FILE = RAW / "heca" / "endometriumAtlasV2_cells_with_counts.h5ad"
PB_DIR = PROC / "heca_pseudobulk"
GS_DIR = PROC / "genesets" / "heca_v014"
COVAR_DIR = PROC / "magma_covar"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "magma_property" / ("p9_v014_heca_magma_" + TS)
LOG = REPO / "logs" / ("p9_v014_heca_magma_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Pre-registered rules (CHANGELOG v0.2.6)
DONOR_COL = "genotype"
TYPE_COL = "celltype"
PATHOLOGY_COL = "Endometrial_pathology"
CONTROL_PATHOLOGY = "C"
STAGE_COL = "Binary Stage"
EXCLUDE_STAGE = "Hormones"
MIN_CELLS_PER_DONOR_GROUP = 10
MIN_DONORS_PER_TYPE = 3
MIN_CELLS_PER_TYPE = 50
MIN_CPM_ANY_TYPE = 1.0
G4_CHECK3 = ["dStromal_early", "dStromal_mid", "dStromal_late"]
TRAITS = ["GCST90205183", "ENDO_ALL", "ENDO_OVARY", "ENDO_PERITONEUM", "ENDO_DEEP", "LEIOMYOMA", "HYPERTENSION"]
VASCULAR = ["Arterial", "Venous", "ePV_1a", "ePV_1b", "ePV_2", "mPV", "uSMCs", "Lymphatic"]

for d in (OUT, LOG.parent, PB_DIR, GS_DIR, COVAR_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V014_HECA_MAGMA", "timestamp_utc": TS, "gates": {}, "results": {}}


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


def safe(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


# ----------------------------------------------------------------------------
# HECA readers
# ----------------------------------------------------------------------------
def dec(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def obs_cat(f, name):
    obs = f["obs"]
    g = obs[name]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = np.array([dec(c) for c in g["categories"][:]], dtype=object)
        codes = g["codes"][:]
    else:
        cats = np.array([dec(c) for c in obs["__categories"][name][:]], dtype=object)
        codes = g[:]
    out = np.empty(len(codes), dtype=object)
    ok = codes >= 0
    out[ok] = cats[codes[ok]]
    out[~ok] = None
    return out


def build_pseudobulk():
    say("")
    say("-" * 72)
    say("[SECTION] HECA PSEUDOBULK (donor x celltype, raw counts layer)")
    done = PB_DIR / "pseudobulk_counts.npy"
    with h5py.File(HECA_FILE, "r") as f:
        layer = f["layers"]["counts"]
        enc = dec(layer.attrs.get("encoding-type", "csr_matrix"))
        shape = tuple(int(x) for x in layer.attrs.get("shape"))
        say("[COUNTS] encoding=%s shape=%s nnz=%d" % (enc, shape, layer["data"].shape[0]))
        if "csr" not in enc:
            raise RuntimeError("counts layer is not CSR")
        var = f["var"]
        genes = [dec(x) for x in var[dec(var.attrs.get("_index", "_index"))][:]]
        donor = obs_cat(f, DONOR_COL)
        ctype = obs_cat(f, TYPE_COL)
        path = obs_cat(f, PATHOLOGY_COL)
        stage = obs_cat(f, STAGE_COL)
        group_col = obs_cat(f, "Group")
        keep = (path == CONTROL_PATHOLOGY) & (stage != EXCLUDE_STAGE) & pd.notna(donor) & pd.notna(ctype)
        say("[CELLS] total=%d control(C)=%d control non-hormone=%d" % (len(path), int((path == CONTROL_PATHOLOGY).sum()), int(keep.sum())))
        vc = pd.Series(group_col[keep]).value_counts()
        say("[CELLS_BY_GROUP] " + " | ".join("%s=%d" % (k, v) for k, v in vc.items()))
        vc = pd.Series(stage[keep]).value_counts()
        say("[CELLS_BY_STAGE] " + " | ".join("%s=%d" % (k, v) for k, v in vc.items()))
        say("[DONORS_KEPT] %d" % pd.Series(donor[keep]).nunique())
        keys = pd.Series(donor[keep]).astype(str) + "||" + pd.Series(ctype[keep]).astype(str)
        codes, uniq = pd.factorize(keys)
        gid = np.full(len(path), -1, dtype=np.int64)
        gid[np.where(keep)[0]] = codes
        G, ng = len(uniq), len(genes)
        say("[GROUPS] donor x celltype = %d; genes = %d" % (G, ng))
        if not done.exists():
            indptr = layer["indptr"][:]
            acc = np.zeros(G * ng, dtype=np.float64)
            n = len(path)
            step = 20000
            t0 = time.time()
            for a in range(0, n, step):
                b = min(n, a + step)
                g = gid[a:b]
                if (g < 0).all():
                    continue
                s, e = int(indptr[a]), int(indptr[b])
                idx = layer["indices"][s:e].astype(np.int64)
                val = layer["data"][s:e].astype(np.float64)
                rg = np.repeat(g, np.diff(indptr[a:b + 1]))
                m = rg >= 0
                acc += np.bincount(rg[m] * ng + idx[m], weights=val[m], minlength=G * ng)
                if (a // step) % 4 == 0:
                    say("    [PB] cells %d/%d %.0fs" % (b, n, time.time() - t0))
            np.save(done, acc.reshape(G, ng))
            say("[PSEUDOBULK] built %.0f s" % (time.time() - t0))
    sizes = pd.Series(codes).value_counts().sort_index().values
    gt = pd.DataFrame({"group": list(uniq), "n_cells": sizes})
    gt["donor"] = gt["group"].str.split("||", regex=False).str[0]
    gt["celltype"] = gt["group"].str.split("||", regex=False).str[1]
    gt.to_csv(PB_DIR / "pseudobulk_groups.csv", index=False)
    pd.Series(genes).to_csv(PB_DIR / "genes.csv", index=False, header=["symbol"])
    return np.load(done, mmap_mode="r"), gt, genes


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


def specificity(pb, gt, genes, sym2):
    sub = gt[gt["n_cells"] >= MIN_CELLS_PER_DONOR_GROUP]
    agg = sub.groupby("celltype").agg(n_donors=("donor", "nunique"), n_cells=("n_cells", "sum"))
    ok = agg[(agg["n_donors"] >= MIN_DONORS_PER_TYPE) & (agg["n_cells"] >= MIN_CELLS_PER_TYPE)].index.tolist()
    say("")
    say("[H_celltype_CTL] eligible types %d / %d" % (len(ok), len(agg)))
    for t, r in agg.sort_values("n_cells", ascending=False).iterrows():
        say("    [TYPE] %-28s donors=%-3d cells=%-7d %s" % (t, r["n_donors"], r["n_cells"], "OK" if t in ok else "DROP"))
    mapped = np.array([g in sym2 for g in genes])
    means = {}
    for t in ok:
        rows = sub.index[sub["celltype"] == t].to_numpy()
        mats = np.asarray(pb[rows, :], dtype=np.float64)
        lib = mats.sum(axis=1, keepdims=True)
        means[t] = (mats / np.where(lib > 0, lib, 1) * 1e6).mean(axis=0)
    M = pd.DataFrame(means, index=genes)[mapped]
    M = M[(M >= MIN_CPM_ANY_TYPE).any(axis=1)]
    S = M.div(M.sum(axis=1), axis=0)
    M.to_csv(GS_DIR / "mean_cpm_by_type.csv.gz")
    S.to_csv(GS_DIR / "specificity_by_type.csv.gz")
    say("[H_celltype_CTL] universe=%d genes" % len(S))
    keys = ["C_" + safe(c) for c in S.columns]
    cov = S.copy()
    cov.columns = keys
    cov["Average"] = np.log2(M.mean(axis=1) + 1.0)
    cov.insert(0, "GENE", [sym2[g] for g in cov.index])
    covf = COVAR_DIR / "H_celltype_CTL.covar.txt"
    cov.to_csv(covf, sep="\t", index=False, float_format="%.6g")
    return covf, dict(zip(keys, S.columns)), agg


def gene_property(trait, covar, keymap):
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    if not raw.exists():
        say("[PROPERTY] %s genes.raw missing" % trait)
        return None
    outp = OUT / ("%s__H_celltype_CTL" % trait)
    rc, o, e = sh(NICE + [str(MAGMA), "--gene-results", str(raw), "--gene-covar", str(covar),
                          "--model", "condition-hide=Average", "direction=greater", "--out", str(outp)], timeout=3600)
    gsa = Path(str(outp) + ".gsa.out")
    if rc != 0 or not gsa.exists():
        say("[PROPERTY] %s FAILED rc=%d %s" % (trait, rc, (o + e).strip()[-300:]))
        return None
    rows = [l for l in open(gsa) if not l.startswith("#")]
    hdr = rows[0].split()
    df = pd.DataFrame([dict(zip(hdr, l.split())) for l in rows[1:] if l.strip()])
    for c in ["NGENES", "BETA", "BETA_STD", "SE", "P"]:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["label"] = df["VARIABLE"].map(keymap)
    df = df.dropna(subset=["label"]).copy()
    df["q_bh"] = bh(df["P"].values)
    df["trait"] = trait
    df = df.sort_values("P").reset_index(drop=True)
    say("[PROPERTY_TOP] %s | HECA celltype (CTL)" % trait)
    for _, r in df.head(8).iterrows():
        say("    %-26s beta=%8.4f se=%7.4f p=%.2e q=%.3f" % (r["label"], r["BETA"], r["SE"], r["P"], r["q_bh"]))
    return df


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.14 | HECA SETS + MAGMA GENE-PROPERTY (G4 check 3) | " + TS)
    say("=" * 72)
    pb, gt, genes = build_pseudobulk()
    sym2 = entrez_map()
    say("[GENE_MAP] %d/%d HECA genes mapped to unique autosomal NCBI37.3 genes" % (sum(g in sym2 for g in genes), len(genes)))
    covf, keymap, agg = specificity(pb, gt, genes, sym2)
    present = [t for t in G4_CHECK3 if t in keymap.values()]
    say("[G4_CHECK3_SETS_PRESENT] %s" % (", ".join(present) or "NONE"))

    say("")
    say("-" * 72)
    say("[SECTION] MAGMA GENE-PROPERTY ON HECA TYPES")
    res = []
    for trait in TRAITS:
        df = gene_property(trait, covf, keymap)
        if df is None:
            continue
        res.append(df)
        if trait == "GCST90205183":
            sub = df[df["label"].isin(G4_CHECK3)].copy()
            if len(sub):
                sub["holm"] = holm(sub["P"].values)
                for _, r in sub.sort_values("P").iterrows():
                    rank = int(df.index[df["label"] == r["label"]][0]) + 1
                    say("[G4_CHECK3] %-16s beta=%8.4f se=%7.4f p=%.2e Holm=%.3f rank %d/%d" % (
                        r["label"], r["BETA"], r["SE"], r["P"], r["holm"], rank, len(df)))
                rep = bool((sub["holm"] < 0.05).any())
                verdict = ("HECA dStromal ENRICHMENT REPRODUCED (check 3 passed)" if rep
                           else "HECA dStromal ENRICHMENT NOT REPRODUCED (check 3 failed)")
                say("[G4_CHECK3_VERDICT] " + verdict)
                SUM["results"]["G4_check3"] = {"verdict": verdict,
                                               "rows": sub[["label", "BETA", "SE", "P", "holm"]].to_dict("records")}
        if trait == "HYPERTENSION":
            hit = [x for x in df["label"].head(5) if x in VASCULAR]
            say("[POSITIVE_CONTROL_HTN] top-5 vascular/perivascular: %d (%s)" % (len(hit), ", ".join(hit)))
            SUM["gates"]["PC_HTN_TOP5"] = len(hit)
        if trait == "LEIOMYOMA":
            if (df["label"] == "uSMCs").any():
                r = df[df["label"] == "uSMCs"].iloc[0]
                rank = int(df.index[df["label"] == "uSMCs"][0]) + 1
                say("[POSITIVE_CONTROL_LEIOMYOMA] uSMCs rank %d/%d p=%.2e q=%.3f" % (rank, len(df), r["P"], r["q_bh"]))
                SUM["gates"]["PC_LEIOMYOMA_uSMC_rank"] = rank
            else:
                say("[POSITIVE_CONTROL_LEIOMYOMA] uSMCs not eligible")
    if res:
        pd.concat(res).to_csv(OUT / "heca_magma_property_results.csv", index=False)
    with open(OUT / "p9_v014_heca_magma_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    for k, v in SUM["gates"].items():
        say("[GATE] %-26s %s" % (k, v))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] G4 check 3 rule fixed in CHANGELOG v0.2.6 before this run")
    say("[DONE] PAPER9_V014_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
