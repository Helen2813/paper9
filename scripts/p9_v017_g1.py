#!/usr/bin/env python3
# Paper 9 | v0.0.17 | G1 PRIMARY TEST (rules CHANGELOG v0.2.7 + v0.2.9): cluster-level pseudobulk with lineage
#                    agreement filter, lesion-vs-eutopic logFC per lineage, MAGMA gene-property; secondary contrasts,
#                    matched FinnGen subtypes, specificity controls; HECA all-donor sensitivity for G4.
# (history) v0.0.16 | (0) G4 check 4: do a few genome-wide loci drive the decidual tilt? (rules CHANGELOG v0.2.8)
#                    + v015 G1 data preparation re-run with fixes (exact kNN, per-sample cache, all archive members)
# Paper 9 | v0.0.15 | G1 data preparation (no genetics, no tests):
#   (1) GSE213216 (Tan 2022) per-sample 10x matrices: QC, per-sample clustering, label transfer from
#       GSE266265 Tan-nomenclature centroids, independent canonical-marker lineage check,
#       pseudobulk sample x label with tissue context.
#   (2) HECA GitHub Endometriosis_fGWAS folder: which GWAS / annotation the published result used.
# Runs in ~/Desktop/paper9/.venv, low-heat (2 threads, nice). Rules: CHANGELOG v0.2.7.

import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("NUMBA_NUM_THREADS", "2")
try:
    os.nice(10)
except Exception:
    pass

import datetime as dt
import gzip
import io
import json
import re
import shutil
import subprocess
import sys
import tarfile
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
VENV_PY = REPO / ".venv" / "bin" / "python"
TAR = RAW / "geo" / "GSE213216" / "GSE213216_RAW.tar"
SERIES_MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213216/matrix/GSE213216_series_matrix.txt.gz"
SERIES_MATRIX = RAW / "geo" / "GSE213216" / "GSE213216_series_matrix.txt.gz"
X_DIR = INTERIM / "gse213216_extract"
OUTPB = PROC / "gse213216_pseudobulk_clusters"
REF_PB = PROC / "gse266265_pseudobulk"
REF_SPEC = PROC / "genesets" / "gse266265_v007" / "SA_refined_ALL" / "specificity_by_type.csv.gz"
HECA_GH = "https://api.github.com/repos/ventolab/HECA-Human-Endometrial-Cell-Atlas/contents/Endometriosis_fGWAS"
HECA_GH_DIR = EXT / "heca_fgwas_repo"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "g1" / ("p9_v017_g1_" + TS)
LOG = REPO / "logs" / ("p9_v017_g1_" + TS + ".log")

# Pre-registered processing rules (CHANGELOG v0.2.7)
MIN_GENES = 200
MAX_PCT_MT = 25.0
MIN_CELLS_SAMPLE = 200
N_HVG = 2000
N_PCS = 30
LEIDEN_RES = 1.0
N_MARKERS_PER_TYPE = 30
MIN_REF_DONOR_CELLS = 10
LINEAGE_MARKERS = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
    "stromal": ["PDGFRA", "DCN", "COL1A1", "COL3A1", "LUM"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5"],
    "myeloid": ["CD68", "LYZ", "CD14", "C1QA", "TYROBP"],
    "lymphoid": ["CD3E", "CD2", "NKG7", "GNLY", "MS4A1"],
    "smooth_muscle": ["ACTA2", "MYH11", "DES", "CNN1"],
    "mesothelial": ["MSLN", "UPK3B", "KRT5", "CALB2"],
}

for d in (OUT, LOG.parent, X_DIR, OUTPB, HECA_GH_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V017_G1", "timestamp_utc": TS, "gates": {}}


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def curl(url, dst=None, max_time=120):
    cmd = ["curl", "-sSL", "--fail", "--max-time", str(max_time), "-A", "paper9"]
    if dst:
        cmd += ["-o", str(dst)]
    cmd.append(url)
    p = subprocess.run(cmd, capture_output=True)
    if p.returncode != 0:
        return None
    return True if dst else p.stdout


# ----------------------------------------------------------------------------
# 0. Dependencies (leiden clustering)
# ----------------------------------------------------------------------------
def ensure_leiden():
    try:
        import leidenalg  # noqa
        import igraph  # noqa
        return True
    except Exception:
        p = subprocess.run([str(VENV_PY), "-m", "pip", "install", "-q", "leidenalg", "igraph"], capture_output=True, text=True)
        say("[PIP] leidenalg igraph rc=%d" % p.returncode)
        try:
            import importlib
            importlib.invalidate_caches()
            import leidenalg  # noqa
            import igraph  # noqa
            return True
        except Exception as e:
            say("[LEIDEN] unavailable (%s) -> k-means fallback" % type(e).__name__)
            return False


# ----------------------------------------------------------------------------
# 1. Sample metadata
# ----------------------------------------------------------------------------
def sample_metadata():
    if not SERIES_MATRIX.exists():
        curl(SERIES_MATRIX_URL, SERIES_MATRIX)
    txt = gzip.open(SERIES_MATRIX, "rt", errors="ignore").read().split("!series_matrix_table_begin")[0]
    rows = {}
    for line in txt.splitlines():
        if line.startswith("!Sample_"):
            tag = line.split("\t")[0]
            vals = [v.strip().strip('"') for v in line.split("\t")[1:]]
            rows.setdefault(tag, []).append(vals)
    n = len(rows["!Sample_geo_accession"][0])
    md = pd.DataFrame({"gsm": rows["!Sample_geo_accession"][0], "title": rows["!Sample_title"][0],
                       "source": rows["!Sample_source_name_ch1"][0]})
    for vals in rows.get("!Sample_characteristics_ch1", []):
        for i, v in enumerate(vals):
            if ":" in v:
                k, x = v.split(":", 1)
                md.loc[i, k.strip()] = x.strip()
    for tag in ["!Sample_library_strategy", "!Sample_description", "!Sample_molecule_ch1"]:
        if tag in rows:
            md[tag.replace("!Sample_", "")] = rows[tag][0]
    md.to_csv(OUT / "gse213216_sample_metadata.csv", index=False)
    say("[METADATA] samples=%d columns=%s" % (n, ", ".join(md.columns)))
    return md


def context_of(row):
    t = str(row.get("tissue", "")) + " " + str(row.get("source", ""))
    t = t.lower()
    if "eutopic" in t:
        return "eutopic"
    if "endometrioma" in t:
        return "endometrioma"
    if "no endometriosis" in t:
        return "peritoneum_unaffected"
    if "ovary" in t:
        return "ovary_unaffected"
    if "endometriosis" in t:
        return "peritoneal_lesion"
    return "other"


# ----------------------------------------------------------------------------
# 2. Reference centroids (GSE266265, Tan nomenclature, all donors)
# ----------------------------------------------------------------------------
def reference():
    pb = np.load(REF_PB / "pseudobulk_counts.npy", mmap_mode="r")
    gt = pd.read_csv(REF_PB / "pseudobulk_groups.csv")
    genes = pd.read_csv(REF_PB / "genes.csv")["symbol"].astype(str).tolist()
    S = pd.read_csv(REF_SPEC, index_col=0)
    markers = set()
    for c in S.columns:
        markers.update(S[c].sort_values(ascending=False).index[:N_MARKERS_PER_TYPE])
    gidx = {g: i for i, g in enumerate(genes)}
    markers = sorted(m for m in markers if m in gidx)
    cols = [gidx[m] for m in markers]
    cent = {}
    lineage = {}
    sub = gt[gt["n_cells"] >= MIN_REF_DONOR_CELLS]
    for t in S.columns:
        rows = sub.index[sub["refined"] == t].to_numpy()
        if len(rows) == 0:
            continue
        mats = np.asarray(pb[rows, :], dtype=np.float64)
        lib = mats.sum(axis=1, keepdims=True)
        cpm = mats / np.where(lib > 0, lib, 1) * 1e6
        cent[t] = np.log2(cpm[:, cols] + 1).mean(axis=0)
        lineage[t] = sub.loc[rows, "main"].mode().iloc[0]
    C = pd.DataFrame(cent, index=markers)
    say("[REFERENCE] types=%d marker genes=%d" % (C.shape[1], C.shape[0]))
    return C, lineage


# ----------------------------------------------------------------------------
# 3. Per-sample processing
# ----------------------------------------------------------------------------
def extract_samples():
    done = X_DIR / "extract_v2.done"
    if done.exists():
        return sorted(p for p in X_DIR.iterdir() if p.is_dir())
    with tarfile.open(TAR, "r:") as outer:
        for m in outer.getmembers():
            if not m.name.endswith(".tar.gz"):
                gsm = m.name.split("_")[0]
                d = X_DIR / gsm
                d.mkdir(exist_ok=True)
                target = d / Path(m.name).name
                if not target.exists() and m.isfile():
                    with outer.extractfile(m) as src, open(target, "wb") as dst:
                        shutil.copyfileobj(src, dst)
                say("    [MEMBER] %s (not tar.gz, %d bytes)" % (m.name, m.size))
                continue
            gsm = m.name.split("_")[0]
            d = X_DIR / gsm
            if d.exists() and (any(d.rglob("*matrix.mtx*")) or any(d.rglob("*.h5"))):
                continue
            d.mkdir(exist_ok=True)
            fobj = outer.extractfile(m)
            with tarfile.open(fileobj=io.BytesIO(fobj.read()), mode="r:gz") as inner:
                for im in inner.getmembers():
                    if im.name.startswith("/") or ".." in Path(im.name).parts:
                        continue
                    inner.extract(im, d)
            say("    [EXTRACT] %s %s" % (gsm, m.name))
    done.write_text("ok")
    return sorted(p for p in X_DIR.iterdir() if p.is_dir())


def find_10x(d):
    mtx = sorted(d.rglob("*matrix.mtx*"))
    if mtx:
        name = mtx[0].name
        return ("mtx", mtx[0].parent, name[:name.index("matrix.mtx")])
    h5 = sorted(d.rglob("*.h5"))
    if h5:
        return ("h5", h5[0], "")
    return None


def cluster_labels(adata, use_leiden):
    import scanpy as sc
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat")
    ad = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(ad, max_value=10)
    sc.tl.pca(ad, n_comps=min(N_PCS, ad.n_obs - 1, ad.n_vars - 1))
    if use_leiden:
        sc.pp.neighbors(ad, n_neighbors=15, transformer="sklearn")
        sc.tl.leiden(ad, resolution=LEIDEN_RES, key_added="cl", flavor="igraph", n_iterations=2, directed=False)
        return ad.obs["cl"].astype(str).values
    from sklearn.cluster import KMeans
    k = max(5, min(30, ad.n_obs // 300))
    return KMeans(n_clusters=k, n_init=5, random_state=0).fit_predict(ad.obsm["X_pca"]).astype(str)


def process_sample(gsm, d, meta_row, C, lineage, use_leiden):
    import scanpy as sc
    import scipy.sparse as sp
    tenx = find_10x(d)
    if tenx is None:
        return None, "no_10x_matrix (files: %s)" % ", ".join(sorted(p.name for p in d.rglob("*") if p.is_file())[:6])
    try:
        if tenx[0] == "mtx":
            adata = sc.read_10x_mtx(tenx[1], var_names="gene_symbols", make_unique=True, prefix=tenx[2] or None)
        else:
            adata = sc.read_10x_h5(tenx[1])
            adata.var_names_make_unique()
    except Exception as e:
        return None, "read_failed:%s:%s" % (type(e).__name__, str(e)[:120])
    n0 = adata.n_obs
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    adata = adata[(adata.obs["n_genes_by_counts"] >= MIN_GENES) & (adata.obs["pct_counts_mt"] < MAX_PCT_MT)].copy()
    if adata.n_obs < MIN_CELLS_SAMPLE:
        return None, "too_few_cells(%d)" % adata.n_obs
    counts = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    cl = cluster_labels(adata, use_leiden)
    genes = np.array(adata.var_names)
    gpos = {g: i for i, g in enumerate(genes)}
    common = [g for g in C.index if g in gpos]
    Cc = C.loc[common]
    rows = []
    label_of = {}
    for c in sorted(set(cl)):
        idx = np.where(cl == c)[0]
        tot = np.asarray(counts[idx].sum(axis=0)).ravel()
        cpm = tot / max(tot.sum(), 1) * 1e6
        prof = np.log2(cpm[[gpos[g] for g in common]] + 1)
        r = Cc.apply(lambda col: pd.Series(prof).corr(pd.Series(col.values), method="spearman"))
        r = r.sort_values(ascending=False)
        best, second = r.index[0], r.index[1]
        lm = {}
        for lin, ms in LINEAGE_MARKERS.items():
            ii = [gpos[m] for m in ms if m in gpos]
            lm[lin] = float(np.mean(np.log2(cpm[ii] + 1))) if ii else 0.0
        lin_marker = max(lm, key=lm.get)
        label_of[c] = best
        rows.append({"gsm": gsm, "cluster": c, "n_cells": len(idx), "label": best, "r": float(r.iloc[0]),
                     "margin": float(r.iloc[0] - r.iloc[1]), "second": second,
                     "lineage_ref": lineage.get(best), "lineage_marker": lin_marker})
    pb = {}
    for r in rows:
        idx = np.where(cl == r["cluster"])[0]
        pb["cl" + str(r["cluster"])] = (np.asarray(counts[idx].sum(axis=0)).ravel(), len(idx))
    info = {"gsm": gsm, "cells_raw": n0, "cells_qc": int(adata.n_obs), "clusters": len(set(cl))}
    return (pb, genes, rows, info), "ok"


def lineage_agree(ref, mk):
    if ref is None:
        return False
    ref = str(ref).lower()
    if mk == "smooth_muscle":
        return ref in ("stromal",)
    if mk == "mesothelial":
        return ref in ("stromal", "epithelial")
    if mk == "lymphoid":
        return ref in ("lymphocyte",)
    return ref == mk


# ----------------------------------------------------------------------------
# 4. HECA fGWAS materials
# ----------------------------------------------------------------------------
def heca_fgwas():
    say("")
    say("-" * 72)
    say("[SECTION] HECA Endometriosis_fGWAS materials (GitHub)")
    txt = curl(HECA_GH)
    if not txt:
        say("[HECA_GH] not reachable")
        return
    try:
        items = json.loads(txt)
    except Exception:
        say("[HECA_GH] unexpected response")
        return
    hits = []
    for it in items:
        name, size, url = it.get("name"), it.get("size", 0), it.get("download_url")
        say("    [FILE] %-60s %8d" % (name, size or 0))
        if url and size and size < 3_000_000 and re.search(r"\.(ipynb|py|R|r|md|txt|sh|yaml|yml)$", name):
            dst = HECA_GH_DIR / name
            if curl(url, dst):
                t = dst.read_text(errors="ignore")
                for line in t.splitlines():
                    if re.search(r"gwas|sumstat|summary stat|rahmioglu|sakaue|finngen|ukb|23andme|gcst|\.bed|window|fgwas|marker", line, re.I):
                        s = line.strip().strip('",').strip()
                        if 8 < len(s) < 300:
                            hits.append((name, s))
    seen = set()
    for name, s in hits:
        if s in seen:
            continue
        seen.add(s)
        say("    [HECA_FGWAS] %s :: %s" % (name[:30], s[:240]))
        if len(seen) >= 80:
            break


MAGMA = REPO / "tools" / "magma" / "magma_bin"
MAGMA_OUT = INTERIM / "magma"
GENE_LOC_B37 = REPO / "tools" / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
HECA_COVAR = PROC / "magma_covar" / "H_celltype_CTL.covar.txt"
EXCL_WINDOW = 1_000_000
N_PERM = 10000


def g4_check4():
    say("")
    say("-" * 72)
    say("[SECTION] G4 CHECK 4: genome-wide loci vs polygenic signal in HECA decidual types (GCST90205183)")
    gout = MAGMA_OUT / "GCST90205183_35_10.genes.out"
    cov = pd.read_csv(HECA_COVAR, sep="\t")
    cov["GENE"] = cov["GENE"].astype(str)
    types = [c for c in cov.columns if c.startswith("C_")]
    g = pd.read_csv(gout, sep=r"\s+")
    g["GENE"] = g["GENE"].astype(str)
    bonf = 0.05 / len(g)
    top = g[g["P"] < bonf]
    loc = {}
    for line in open(GENE_LOC_B37):
        x = line.split()
        if len(x) >= 6:
            loc[x[0]] = (x[1], int(x[2]), int(x[3]), x[5])
    say("[TOP_GENES] Bonferroni MAGMA genes = %d: %s" % (len(top), ", ".join(loc.get(e, ("", 0, 0, e))[3] for e in top["GENE"])))
    # (a) specificity of the top genes per type, empirical p vs random gene sets of the same size
    S = cov.set_index("GENE")[types]
    pct = S.rank(pct=True)
    inu = [e for e in top["GENE"] if e in pct.index]
    obs = pct.loc[inu].mean()
    rng = np.random.default_rng(20261106)
    univ = pct.index.to_numpy()
    null = np.empty((N_PERM, len(types)))
    vals = pct.to_numpy()
    for b in range(N_PERM):
        idx = rng.choice(len(univ), size=len(inu), replace=False)
        null[b] = vals[idx].mean(axis=0)
    emp = ((null >= obs.to_numpy()).sum(axis=0) + 1) / (N_PERM + 1)
    a = pd.DataFrame({"type": [t[2:] for t in types], "mean_pct": obs.values, "p_emp": emp}).sort_values("mean_pct", ascending=False)
    a["rank"] = np.arange(1, len(a) + 1)
    say("[CHECK4a] top genes in universe = %d; mean specificity percentile per HECA type (10,000 random sets)" % len(inu))
    for _, r in a.head(6).iterrows():
        say("    %-26s mean_pct=%.3f p_emp=%.4f rank %d" % (r["type"], r["mean_pct"], r["p_emp"], r["rank"]))
    for _, r in a[a["type"].str.startswith("dStromal")].iterrows():
        say("    [dStromal] %-16s mean_pct=%.3f p_emp=%.4f rank %d/%d" % (r["type"], r["mean_pct"], r["p_emp"], r["rank"], len(a)))
    a.to_csv(OUT / "G4_check4a_top_gene_specificity.csv", index=False)
    # (b) MAGMA gene-property after removing genes within +-1 Mb of any Bonferroni gene
    excl = set()
    tops = [loc[e] for e in top["GENE"] if e in loc]
    for e, (c, s, t_, sym) in loc.items():
        for (c2, s2, t2, _) in tops:
            if c == c2 and s <= t2 + EXCL_WINDOW and t_ >= s2 - EXCL_WINDOW:
                excl.add(e)
                break
    sub = cov[~cov["GENE"].isin(excl)]
    subf = OUT / "H_celltype_CTL.minus_top_loci.covar.txt"
    sub.to_csv(subf, sep="\t", index=False)
    say("[CHECK4b] genes removed (+-1 Mb of top genes) = %d; remaining = %d" % (len(cov) - len(sub), len(sub)))
    res = {}
    for name, cf in [("full", HECA_COVAR), ("minus_top_loci", subf)]:
        outp = OUT / ("GCST90205183__HECA__" + name)
        p = subprocess.run(["nice", "-n", "10", str(MAGMA), "--gene-results", str(MAGMA_OUT / "GCST90205183_35_10.genes.raw"),
                            "--gene-covar", str(cf), "--model", "condition-hide=Average", "direction=greater", "--out", str(outp)],
                           capture_output=True, text=True)
        gsa = Path(str(outp) + ".gsa.out")
        if not gsa.exists():
            say("[CHECK4b] MAGMA failed for %s: %s" % (name, (p.stdout + p.stderr)[-300:]))
            return
        rows = [l for l in open(gsa) if not l.startswith("#")]
        hdr = rows[0].split()
        df = pd.DataFrame([dict(zip(hdr, l.split())) for l in rows[1:] if l.strip()])
        for c in ["BETA", "SE", "P"]:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df = df[df["VARIABLE"].str.startswith("C_")].sort_values("P").reset_index(drop=True)
        df["rank"] = np.arange(1, len(df) + 1)
        res[name] = df.set_index("VARIABLE")
    for t in ["C_dStromal_early", "C_dStromal_mid", "C_dStromal_late"]:
        f, m = res["full"].loc[t], res["minus_top_loci"].loc[t]
        say("[CHECK4b] %-16s full: beta=%.3f p=%.3f rank %d | minus top loci: beta=%.3f p=%.3f rank %d" % (
            t[2:], f["BETA"], f["P"], f["rank"], m["BETA"], m["P"], m["rank"]))
    say("[CHECK4b] top-3 after removal: " + ", ".join("%s(p=%.3f)" % (i[2:], r["P"]) for i, r in res["minus_top_loci"].head(3).iterrows()))
    SUM["gates"]["CHECK4_done"] = True


LINEAGES = ["stromal", "epithelial", "endothelial", "myeloid", "lymphocyte"]
MIN_CELLS_PATIENT = 10
MIN_PATIENTS = 3
G1_TRAITS_SECONDARY = ["ENDO_PERITONEUM", "ENDO_OVARY", "ENDO_ALL", "HYPERTENSION", "LEIOMYOMA"]
COVAR_G1 = PROC / "magma_covar" / "g1"


def holm_adj(p):
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


def patient_profiles(M, G, genes, lineage, context, patients=None):
    sub = G[(G["agree"]) & (G["lineage_ref"] == lineage) & (G["context"] == context)]
    if patients is not None:
        sub = sub[sub["patient"].isin(patients)]
    prof = {}
    for pid, g in sub.groupby("patient"):
        n = int(g["n_cells"].sum())
        if n < MIN_CELLS_PATIENT:
            continue
        v = M[g.index.to_numpy()].sum(axis=0)
        prof[pid] = np.log2(v / max(v.sum(), 1) * 1e6 + 1)
    return prof


def magma_property(trait, covf, test_col, cond_col, direction, tag):
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    if not raw.exists():
        return None
    outp = COVAR_G1 / ("%s__%s__%s" % (trait, tag, direction))
    p = subprocess.run(["nice", "-n", "10", str(MAGMA), "--gene-results", str(raw), "--gene-covar", str(covf),
                        "--model", "condition-hide=" + cond_col, "direction=" + direction, "--out", str(outp)],
                       capture_output=True, text=True)
    gsa = Path(str(outp) + ".gsa.out")
    if not gsa.exists():
        say("    [MAGMA_FAIL] %s %s %s" % (trait, tag, (p.stdout + p.stderr)[-200:]))
        return None
    rows = [l for l in open(gsa) if not l.startswith("#")]
    hdr = rows[0].split()
    for l in rows[1:]:
        d = dict(zip(hdr, l.split()))
        if d.get("VARIABLE") == test_col:
            return {"beta": float(d["BETA"]), "se": float(d["SE"]), "p": float(d["P"]), "ngenes": int(d["NGENES"])}
    return None


def contrast_covariates(M, G, genes, sym2, ctx_a, ctx_b, tag, patients=None):
    COVAR_G1.mkdir(parents=True, exist_ok=True)
    mapped = np.array([g in sym2 for g in genes])
    out = {}
    for lin in LINEAGES:
        pa = patient_profiles(M, G, genes, lin, ctx_a, patients)
        pb = patient_profiles(M, G, genes, lin, ctx_b, patients)
        if len(pa) < MIN_PATIENTS or len(pb) < MIN_PATIENTS:
            say("    [%s] %-11s patients %s=%d %s=%d -> SKIP (<%d)" % (tag, lin, ctx_a, len(pa), ctx_b, len(pb), MIN_PATIENTS))
            continue
        A = np.vstack(list(pa.values())).mean(axis=0)
        B = np.vstack(list(pb.values())).mean(axis=0)
        expressed = (np.maximum(2 ** A - 1, 2 ** B - 1) >= 1.0) & mapped
        D = (A - B)[expressed]
        AV = ((A + B) / 2)[expressed]
        g = np.array(genes)[expressed]
        df = pd.DataFrame({"GENE": [sym2[x] for x in g], "D": D, "AVG": AV})
        covf = COVAR_G1 / ("%s__%s.covar.txt" % (tag, lin))
        df.to_csv(covf, sep="\t", index=False, float_format="%.6g")
        out[lin] = (covf, len(pa), len(pb), int(expressed.sum()))
        say("    [%s] %-11s patients %s=%d %s=%d genes=%d" % (tag, lin, ctx_a, len(pa), ctx_b, len(pb), int(expressed.sum())))
    return out


def run_contrast(cov, trait, tag, primary=False):
    res = []
    for lin, (covf, na, nb, ng) in cov.items():
        g = magma_property(trait, covf, "D", "AVG", "greater", tag + "_" + lin)
        l = magma_property(trait, covf, "D", "AVG", "less", tag + "_" + lin)
        if g is None or l is None:
            continue
        res.append({"trait": trait, "contrast": tag, "lineage": lin, "n_a": na, "n_b": nb, "genes": ng,
                    "beta": g["beta"], "se": g["se"], "p_greater": g["p"], "p_less": l["p"]})
    if not res:
        return pd.DataFrame()
    df = pd.DataFrame(res)
    df["holm_greater"] = holm_adj(df["p_greater"].values)
    df["holm_less"] = holm_adj(df["p_less"].values)
    for _, r in df.iterrows():
        say("    [%s] %-15s %-11s beta=%8.4f se=%7.4f p_greater=%.2e Holm=%.3f | p_less=%.2e Holm=%.3f" % (
            "G1" if primary else "SEC", trait[:15], r["lineage"], r["beta"], r["se"], r["p_greater"], r["holm_greater"],
            r["p_less"], r["holm_less"]))
    return df


def g1_tests(M, G, genes):
    say("")
    say("-" * 72)
    say("[SECTION] G1 PRIMARY: peritoneal lesion vs eutopic, per lineage, GCST90205183 (Holm over lineages)")
    sym2 = entrez_unique()
    G = G.reset_index(drop=True)
    say("[G1_FILTER] clusters kept (lineage agreement) %d / %d; cells kept %d / %d" % (
        int(G["agree"].sum()), len(G), int(G.loc[G["agree"], "n_cells"].sum()), int(G["n_cells"].sum())))
    cov = contrast_covariates(M, G, genes, sym2, "peritoneal_lesion", "eutopic", "LES_vs_EU")
    prim = run_contrast(cov, "GCST90205183", "LES_vs_EU", primary=True)
    allres = [prim]
    if len(prim):
        sup = prim[(prim["holm_greater"] < 0.05) & (prim["beta"] > 0)]
        alt = prim[prim["holm_less"] < 0.05]
        if len(sup):
            verdict = "G1 SUPPORTED (lesion-upregulated programs enriched): " + ", ".join(sup["lineage"])
        elif len(alt):
            verdict = "G1 ALTERNATIVE (source-defect: eutopic-upregulated programs enriched): " + ", ".join(alt["lineage"])
        else:
            verdict = "G1 NOT SUPPORTED (no lineage contrast enriched in either direction)"
        say("[G1_VERDICT] " + verdict)
        SUM["results"] = SUM.get("results", {})
        SUM["results"]["G1_verdict"] = verdict
        SUM["results"]["G1_primary"] = prim.to_dict("records")

    say("")
    say("[SECTION] G1 SECONDARY (descriptive): matched subtypes, specificity controls, endometrioma, paired")
    for trait in G1_TRAITS_SECONDARY:
        allres.append(run_contrast(cov, trait, "LES_vs_EU"))
    cov2 = contrast_covariates(M, G, genes, sym2, "endometrioma", "eutopic", "OMA_vs_EU")
    for trait in ["GCST90205183", "ENDO_OVARY", "LEIOMYOMA", "HYPERTENSION"]:
        allres.append(run_contrast(cov2, trait, "OMA_vs_EU"))
    both = set(G[(G["context"] == "peritoneal_lesion")]["patient"]) & set(G[(G["context"] == "eutopic")]["patient"])
    say("    [PAIRED] patients with both lesion and eutopic samples: %d" % len(both))
    if len(both) >= MIN_PATIENTS:
        cov3 = contrast_covariates(M, G, genes, sym2, "peritoneal_lesion", "eutopic", "LES_vs_EU_paired", patients=both)
        allres.append(run_contrast(cov3, "GCST90205183", "LES_vs_EU_paired"))
    out = pd.concat([d for d in allres if len(d)], ignore_index=True) if any(len(d) for d in allres) else pd.DataFrame()
    out.to_csv(OUT / "G1_results.csv", index=False)


def heca_all_donors():
    say("")
    say("-" * 72)
    say("[SECTION] G4 SENSITIVITY: HECA all donors (any pathology, hormones included), MAGMA, GCST90205183")
    import h5py
    f = h5py.File(HECA_FILE, "r")

    def cat(name):
        g = f["obs"][name]
        if isinstance(g, h5py.Group) and "categories" in g:
            cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["categories"][:]], dtype=object)
            codes = g["codes"][:]
        else:
            cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f["obs"]["__categories"][name][:]], dtype=object)
            codes = g[:]
        o = np.empty(len(codes), dtype=object)
        o[codes >= 0] = cats[codes[codes >= 0]]
        return o
    donor, ctype = cat("genotype"), cat("celltype")
    keys = pd.Series(donor).astype(str) + "||" + pd.Series(ctype).astype(str)
    codes, uniq = pd.factorize(keys)
    layer = f["layers"]["counts"]
    var = f["var"]
    genes = [x.decode() if isinstance(x, bytes) else x for x in var[var.attrs.get("_index", "_index")][:]]
    ng, Gn = len(genes), len(uniq)
    indptr = layer["indptr"][:]
    acc = np.zeros(Gn * ng)
    n = len(codes)
    for a in range(0, n, 20000):
        b = min(n, a + 20000)
        s, e = int(indptr[a]), int(indptr[b])
        idx = layer["indices"][s:e].astype(np.int64)
        val = layer["data"][s:e].astype(np.float64)
        rg = np.repeat(codes[a:b].astype(np.int64), np.diff(indptr[a:b + 1]))
        acc += np.bincount(rg * ng + idx, weights=val, minlength=Gn * ng)
    pbm = acc.reshape(Gn, ng)
    gt = pd.DataFrame({"group": list(uniq), "n_cells": pd.Series(codes).value_counts().sort_index().values})
    gt["donor"] = gt["group"].str.split("||", regex=False).str[0]
    gt["celltype"] = gt["group"].str.split("||", regex=False).str[1]
    sub = gt[gt["n_cells"] >= 10]
    agg = sub.groupby("celltype").agg(nd=("donor", "nunique"), nc=("n_cells", "sum"))
    ok = agg[(agg["nd"] >= 3) & (agg["nc"] >= 50)].index.tolist()
    sym2 = entrez_unique()
    mapped = np.array([g in sym2 for g in genes])
    means = {}
    for t in ok:
        rows = sub.index[sub["celltype"] == t].to_numpy()
        m = pbm[rows]
        means[t] = (m / np.maximum(m.sum(axis=1, keepdims=True), 1) * 1e6).mean(axis=0)
    Mx = pd.DataFrame(means, index=genes)[mapped]
    Mx = Mx[(Mx >= 1).any(axis=1)]
    S = Mx.div(Mx.sum(axis=1), axis=0)
    cov = S.copy()
    cov.columns = ["C_" + re.sub(r"[^A-Za-z0-9]+", "_", c) for c in S.columns]
    cov["Average"] = np.log2(Mx.mean(axis=1) + 1)
    cov.insert(0, "GENE", [sym2[g] for g in cov.index])
    covf = COVAR_G1 / "H_celltype_ALL.covar.txt"
    COVAR_G1.mkdir(parents=True, exist_ok=True)
    cov.to_csv(covf, sep="\t", index=False, float_format="%.6g")
    outp = COVAR_G1 / "GCST90205183__HECA_ALL"
    subprocess.run(["nice", "-n", "10", str(MAGMA), "--gene-results", str(MAGMA_OUT / "GCST90205183_35_10.genes.raw"),
                    "--gene-covar", str(covf), "--model", "condition-hide=Average", "direction=greater", "--out", str(outp)],
                   capture_output=True, text=True)
    rows = [l for l in open(str(outp) + ".gsa.out") if not l.startswith("#")]
    hdr = rows[0].split()
    df = pd.DataFrame([dict(zip(hdr, l.split())) for l in rows[1:] if l.strip()])
    df["P"] = pd.to_numeric(df["P"])
    df["BETA"] = pd.to_numeric(df["BETA"])
    df = df[df["VARIABLE"].str.startswith("C_")].sort_values("P").reset_index(drop=True)
    say("[HECA_ALL] donors=%d types=%d" % (gt["donor"].nunique(), len(ok)))
    for i, r in df.head(6).iterrows():
        say("    %-26s beta=%.3f p=%.3e rank %d" % (r["VARIABLE"][2:], r["BETA"], r["P"], i + 1))
    for i, r in df[df["VARIABLE"].str.contains("dStromal")].iterrows():
        say("    [dStromal] %-16s beta=%.3f p=%.3e rank %d/%d" % (r["VARIABLE"][2:], r["BETA"], r["P"], i + 1, len(df)))


HECA_FILE = RAW / "heca" / "endometriumAtlasV2_cells_with_counts.h5ad"


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.17 | GSE213216 PREPARATION (G1 data) + HECA fGWAS check | " + TS)
    say("=" * 72)
    pass
    import scanpy as sc
    sc.settings.verbosity = 0
    use_leiden = ensure_leiden()
    md = sample_metadata()
    md["context"] = md.apply(context_of, axis=1)
    pid_col = next((c for c in md.columns if c.lower().startswith("patient")), None)
    say("[CONTEXTS] " + " | ".join("%s=%d" % (k, v) for k, v in md["context"].value_counts().items()))
    C, lineage = reference()

    say("")
    say("-" * 72)
    say("[SECTION] EXTRACT + PER-SAMPLE QC, CLUSTERING, LABEL TRANSFER")
    dirs = extract_samples()
    say("[SAMPLES_EXTRACTED] %d" % len(dirs))
    all_genes = None
    groups, cl_rows, infos, mats = [], [], [], []
    t0 = time.time()
    for d in dirs:
        gsm = d.name
        row = md[md["gsm"] == gsm]
        if len(row) == 0:
            say("    [SKIP] %s not in metadata" % gsm)
            continue
        row = row.iloc[0]
        desc = (str(row.get("description", "")) + " " + str(row.get("library_strategy", ""))).lower()
        if "visium" in desc or "spatial" in desc or str(row.get("library_strategy", "")).upper() == "OTHER":
            say("    [SKIP] %s spatial / non-scRNA" % gsm)
            continue
        cache = OUTPB / "samples_v2" / (gsm + ".npz")
        cache.parent.mkdir(parents=True, exist_ok=True)
        if cache.exists():
            z = np.load(cache, allow_pickle=True)
            pb = {str(l): (z["counts"][i], int(z["ncells"][i])) for i, l in enumerate(z["labels"])}
            res = (pb, z["genes"], list(z["rows"]), dict(z["info"].item()))
            status = "cached"
        else:
            try:
                res, status = process_sample(gsm, d, row, C, lineage, use_leiden)
            except Exception as e:
                res, status = None, "error:%s:%s" % (type(e).__name__, str(e)[:160])
            if res is not None:
                pb0, genes0, rows0, info0 = res
                labs = sorted(pb0)
                np.savez_compressed(cache, labels=np.array(labs, dtype=object),
                                    counts=np.vstack([pb0[l][0] for l in labs]),
                                    ncells=np.array([pb0[l][1] for l in labs]),
                                    genes=np.array(genes0, dtype=object),
                                    rows=np.array(rows0, dtype=object), info=np.array(info0, dtype=object))
        if res is None:
            say("    [SKIP] %s %s" % (gsm, status))
            continue
        pb, genes, rows, info = res
        info.update({"context": row["context"], "patient": row.get(pid_col) if pid_col else None})
        infos.append(info)
        cl_rows += rows
        if all_genes is None:
            all_genes = list(C.index)  # placeholder, replaced below
        rmap = {"cl" + str(r["cluster"]): r for r in rows}
        for key, (vec, n) in pb.items():
            r = rmap[key]
            groups.append({"gsm": gsm, "patient": info["patient"], "context": row["context"], "cluster": key,
                           "label": r["label"], "lineage_ref": r["lineage_ref"], "lineage_marker": r["lineage_marker"],
                           "agree": lineage_agree(r["lineage_ref"], r["lineage_marker"]), "n_cells": n})
            mats.append((genes, vec))
        say("    [SAMPLE] %s %-22s cells %d->%d clusters=%d labels=%d %.0f min" % (
            gsm, row["context"], info["cells_raw"], info["cells_qc"], info["clusters"], len(pb), (time.time() - t0) / 60))

    # Align all pseudobulk vectors to the union of gene symbols
    union = sorted(set().union(*[set(g) for g, _ in mats])) if mats else []
    gpos = {g: i for i, g in enumerate(union)}
    M = np.zeros((len(mats), len(union)), dtype=np.float64)
    for k, (g, v) in enumerate(mats):
        M[k, [gpos[x] for x in g]] = v
    np.save(OUTPB / "pseudobulk_counts.npy", M)
    pd.DataFrame(groups).to_csv(OUTPB / "pseudobulk_groups.csv", index=False)
    pd.Series(union).to_csv(OUTPB / "genes.csv", index=False, header=["symbol"])
    cl = pd.DataFrame(cl_rows)
    cl.to_csv(OUT / "cluster_label_transfer.csv", index=False)
    pd.DataFrame(infos).to_csv(OUT / "sample_qc.csv", index=False)

    say("")
    say("-" * 72)
    say("[SECTION] SUMMARY")
    if len(cl):
        cl["agree"] = [lineage_agree(a, b) for a, b in zip(cl["lineage_ref"], cl["lineage_marker"])]
        w = cl["n_cells"]
        say("[LABEL_TRANSFER] clusters=%d median r=%.3f median margin=%.3f" % (len(cl), cl["r"].median(), cl["margin"].median()))
        say("[LINEAGE_AGREEMENT] cell-weighted agreement reference vs canonical markers = %.1f%%" % (100 * (w * cl["agree"]).sum() / w.sum()))
        dis = cl[~cl["agree"]].groupby(["lineage_ref", "lineage_marker"])["n_cells"].sum().sort_values(ascending=False).head(10)
        for (a, b), n in dis.items():
            say("    [DISAGREE] ref=%-12s marker=%-14s cells=%d" % (a, b, n))
        SUM["gates"]["LINEAGE_AGREEMENT_PCT"] = round(100 * (w * cl["agree"]).sum() / w.sum(), 1)
    gdf = pd.DataFrame(groups)
    if len(gdf):
        tab = gdf.groupby(["context", "lineage_ref"])["n_cells"].sum().unstack(fill_value=0)
        say("[CELLS_BY_CONTEXT_x_LINEAGE]")
        for line in tab.to_string().splitlines():
            say("    " + line)
        pat = gdf[gdf["n_cells"] >= 10].groupby(["context", "lineage_ref"])["patient"].nunique().unstack(fill_value=0)
        say("[PATIENTS_WITH_>=10_CELLS_BY_CONTEXT_x_LINEAGE]")
        for line in pat.to_string().splitlines():
            say("    " + line)
        top = gdf.groupby(["context", "label"])["n_cells"].sum().reset_index().sort_values("n_cells", ascending=False)
        for ctx in ["peritoneal_lesion", "endometrioma", "eutopic", "peritoneum_unaffected", "ovary_unaffected"]:
            t = top[top["context"] == ctx].head(10)
            say("[TOP_LABELS] %-22s %s" % (ctx, " | ".join("%s=%d" % (a, b) for a, b in zip(t["label"], t["n_cells"]))))
    SUM["gates"]["SAMPLES_PROCESSED"] = len(infos)
    try:
        g1_tests(M, pd.DataFrame(groups), union)
    except Exception as e:
        import traceback
        say("[STEP_ERROR] g1 %s: %s" % (type(e).__name__, str(e)[:300]))
        say(traceback.format_exc()[-1500:])
    try:
        heca_all_donors()
    except Exception as e:
        say("[STEP_ERROR] heca_all_donors %s: %s" % (type(e).__name__, str(e)[:300]))
    with open(OUT / "p9_v017_g1_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    for k, v in SUM["gates"].items():
        say("[GATE] %-24s %s" % (k, v))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] G1 primary rule fixed in CHANGELOG v0.2.7 (+ v0.2.9 filters) before this run")
    say("[DONE] PAPER9_V017_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
