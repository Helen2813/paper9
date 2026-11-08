#!/usr/bin/env python3
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
OUTPB = PROC / "gse213216_pseudobulk"
REF_PB = PROC / "gse266265_pseudobulk"
REF_SPEC = PROC / "genesets" / "gse266265_v007" / "SA_refined_ALL" / "specificity_by_type.csv.gz"
HECA_GH = "https://api.github.com/repos/ventolab/HECA-Human-Endometrial-Cell-Atlas/contents/Endometriosis_fGWAS"
HECA_GH_DIR = EXT / "heca_fgwas_repo"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "g1_prepare" / ("p9_v015_gse213216_prepare_" + TS)
LOG = REPO / "logs" / ("p9_v015_gse213216_prepare_" + TS + ".log")

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
SUM = {"task": "PAPER9_V015_GSE213216_PREPARE", "timestamp_utc": TS, "gates": {}}


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
    done = X_DIR / "extract.done"
    if done.exists():
        return sorted(p for p in X_DIR.iterdir() if p.is_dir())
    with tarfile.open(TAR, "r:") as outer:
        for m in outer.getmembers():
            if not m.name.endswith(".tar.gz"):
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
        sc.pp.neighbors(ad, n_neighbors=15)
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
    labels = np.array([label_of[c] for c in cl])
    pb = {}
    for lab in sorted(set(labels)):
        idx = np.where(labels == lab)[0]
        pb[lab] = (np.asarray(counts[idx].sum(axis=0)).ravel(), len(idx))
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


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.15 | GSE213216 PREPARATION (G1 data) + HECA fGWAS check | " + TS)
    say("=" * 72)
    try:
        heca_fgwas()
    except Exception as e:
        say("[STEP_ERROR] heca_fgwas %s: %s" % (type(e).__name__, str(e)[:200]))
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
        res, status = process_sample(gsm, d, row, C, lineage, use_leiden)
        if res is None:
            say("    [SKIP] %s %s" % (gsm, status))
            continue
        pb, genes, rows, info = res
        info.update({"context": row["context"], "patient": row.get(pid_col) if pid_col else None})
        infos.append(info)
        cl_rows += rows
        if all_genes is None:
            all_genes = list(C.index)  # placeholder, replaced below
        for lab, (vec, n) in pb.items():
            groups.append({"gsm": gsm, "patient": info["patient"], "context": row["context"], "label": lab,
                           "lineage_ref": lineage.get(lab), "n_cells": n})
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
    with open(OUT / "p9_v015_gse213216_prepare_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    for k, v in SUM["gates"].items():
        say("[GATE] %-24s %s" % (k, v))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] Data preparation only; G1 test rule fixed in CHANGELOG v0.2.7 before any genetics")
    say("[DONE] PAPER9_V015_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
