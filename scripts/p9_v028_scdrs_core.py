#!/usr/bin/env python3
# Paper 9 | v0.0.28 | Ochoa-faithful scDRS (rules CHANGELOG v0.2.18): lesion_ochoa context, M2like group, BH, donor-level tests.
# (from v0.0.25) Option A: Ochoa 2025 claim (lesion macrophages, scDRS) re-tested on GSE213216 (Fonseca 2023)
#   with our cell labels, two endometriosis GWAS and unrelated-trait controls. Rules: CHANGELOG v0.2.16.
#   (1) per-sample re-processing as v016 -> per-cell labels (subsampled)  (2) scDRS per trait and context
#   (3) scDRS association test per lineage and fine label.
# Runs in ~/Desktop/paper9/.venv, low-heat (2 threads, nice 10). ~60-90 min.

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
import json
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
VENV_PY = REPO / ".venv" / "bin" / "python"
SERIES_MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213216/matrix/GSE213216_series_matrix.txt.gz"
SERIES_MATRIX = RAW / "geo" / "GSE213216" / "GSE213216_series_matrix.txt.gz"
X_DIR = INTERIM / "gse213216_extract"
REF_PB = PROC / "gse266265_pseudobulk"
REF_SPEC = PROC / "genesets" / "gse266265_v007" / "SA_refined_ALL" / "specificity_by_type.csv.gz"
GENE_LOC_B37 = REPO / "tools" / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
CELLS_DIR = PROC / "gse213216_cells_v028"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "scdrs" / ("p9_v028_scdrs_" + TS)
LOG = REPO / "logs" / ("p9_v028_scdrs_core_" + TS + ".log")

# Processing rules (CHANGELOG v0.2.7, unchanged)
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
# scDRS rules (CHANGELOG v0.2.16)
CAP_PER_SAMPLE_LABEL = 400
SEED = 20261106
TRAITS = ["GCST90205183", "KOLLER_EUR", "HYPERTENSION", "T2D", "LEIOMYOMA"]
N_TOP = 1000
N_CTRL = 1000
MIN_GROUP_CELLS = 100
CONTEXTS = {"lesion_ochoa": ["peritoneal_lesion", "peritoneum_unaffected", "endometrioma"],
            "lesion_ours": ["peritoneal_lesion", "endometrioma"], "eutopic": ["eutopic"]}
M2_MARKERS = ["CD163", "MRC1", "LYVE1", "F13A1", "STAB1", "SELENOP", "FOLR2"]
MIN_DONOR_CELLS = 10

for d in (OUT, LOG.parent, CELLS_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V028_SCDRS", "timestamp_utc": TS, "gates": {}, "results": {}}


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


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


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


# ----------------------------------------------------------------------------
# Functions reused verbatim from v018 (same processing rules)
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
# 1. Per-sample processing -> per-cell labels (subsampled)
# ----------------------------------------------------------------------------
def process_cells(gsm, d, C, lineage, use_leiden, rng):
    import scanpy as sc
    import scipy.sparse as sp
    tenx = find_10x(d)
    if tenx is None:
        return None, "no_10x_matrix"
    if tenx[0] == "mtx":
        adata = sc.read_10x_mtx(tenx[1], var_names="gene_symbols", make_unique=True, prefix=tenx[2] or None)
    else:
        adata = sc.read_10x_h5(tenx[1])
        adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    adata = adata[(adata.obs["n_genes_by_counts"] >= MIN_GENES) & (adata.obs["pct_counts_mt"] < MAX_PCT_MT)].copy()
    if adata.n_obs < MIN_CELLS_SAMPLE:
        return None, "too_few_cells"
    counts = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    cl = cluster_labels(adata, use_leiden)
    genes = np.array(adata.var_names)
    gpos = {g: i for i, g in enumerate(genes)}
    common = [g for g in C.index if g in gpos]
    Cc = C.loc[common]
    lab = np.empty(adata.n_obs, dtype=object)
    lin = np.empty(adata.n_obs, dtype=object)
    agree = np.zeros(adata.n_obs, dtype=bool)
    for c in sorted(set(cl)):
        idx = np.where(cl == c)[0]
        tot = np.asarray(counts[idx].sum(axis=0)).ravel()
        cpm = tot / max(tot.sum(), 1) * 1e6
        prof = np.log2(cpm[[gpos[g] for g in common]] + 1)
        r = Cc.apply(lambda col: pd.Series(prof).corr(pd.Series(col.values), method="spearman")).sort_values(ascending=False)
        lm = {}
        for ln, ms in LINEAGE_MARKERS.items():
            ii = [gpos[m] for m in ms if m in gpos]
            lm[ln] = float(np.mean(np.log2(cpm[ii] + 1))) if ii else 0.0
        best = r.index[0]
        lab[idx] = best
        lin[idx] = lineage.get(best)
        agree[idx] = lineage_agree(lineage.get(best), max(lm, key=lm.get))
    keep = []
    for l in pd.unique(lab[agree]):
        idx = np.where(agree & (lab == l))[0]
        if len(idx) > CAP_PER_SAMPLE_LABEL:
            idx = rng.choice(idx, CAP_PER_SAMPLE_LABEL, replace=False)
        keep.extend(sorted(idx.tolist()))
    keep = np.array(sorted(keep), dtype=int)
    sub = counts[keep].astype(np.float32)
    obs = pd.DataFrame({"gsm": gsm, "label": lab[keep], "lineage": lin[keep],
                        "n_genes": np.asarray((counts[keep] > 0).sum(axis=1)).ravel()})
    return (sub, genes, obs, int(adata.n_obs), int(agree.sum())), "ok"


def build_cells():
    import scipy.sparse as sp
    f = CELLS_DIR / "cells.npz"
    if f.exists():
        z = np.load(f, allow_pickle=True)
        X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
        obs = pd.read_csv(CELLS_DIR / "obs.csv")
        genes = list(pd.read_csv(CELLS_DIR / "genes.csv")["symbol"].astype(str))
        say("[CELLS] cached: cells=%d genes=%d" % X.shape)
        return X, obs, genes
    md = sample_metadata()
    md["context"] = md.apply(context_of, axis=1)
    pid_col = next((c for c in md.columns if c.lower().startswith("patient")), None)
    C, lineage = reference()
    use_leiden = ensure_leiden()
    rng = np.random.default_rng(SEED)
    want = set(sum(CONTEXTS.values(), []))
    mats, obss, union = [], [], {}
    t0 = time.time()
    for d in sorted(p for p in X_DIR.iterdir() if p.is_dir()):
        row = md[md["gsm"] == d.name]
        if not len(row):
            continue
        row = row.iloc[0]
        desc = (str(row.get("description", "")) + " " + str(row.get("library_strategy", ""))).lower()
        if row["context"] not in want or "visium" in desc or "spatial" in desc or str(row.get("library_strategy", "")).upper() == "OTHER":
            continue
        try:
            res, st = process_cells(d.name, d, C, lineage, use_leiden, rng)
        except Exception as e:
            res, st = None, "error:%s" % type(e).__name__
        if res is None:
            say("    [SKIP] %s %s" % (d.name, st))
            continue
        sub, genes, obs, nqc, nagree = res
        obs["context"] = row["context"]
        obs["patient"] = row.get(pid_col) if pid_col else None
        for g in genes:
            if g not in union:
                union[g] = len(union)
        mats.append((sub, genes))
        obss.append(obs)
        say("    [CELLS] %s %-18s qc=%d agree=%d kept=%d %.0f min" % (d.name, row["context"], nqc, nagree, sub.shape[0], (time.time() - t0) / 60))
    glist = list(union.keys())
    blocks = []
    for sub, genes in mats:
        coo = sub.tocoo()
        cols = np.array([union[g] for g in genes])[coo.col]
        blocks.append(sp.csr_matrix((coo.data, (coo.row, cols)), shape=(sub.shape[0], len(glist))))
    X = sp.vstack(blocks).tocsr()
    obs = pd.concat(obss, ignore_index=True)
    np.savez_compressed(f, data=X.data, indices=X.indices, indptr=X.indptr, shape=np.array(X.shape))
    obs.to_csv(CELLS_DIR / "obs.csv", index=False)
    pd.Series(glist).to_csv(CELLS_DIR / "genes.csv", index=False, header=["symbol"])
    say("[CELLS] built: cells=%d genes=%d" % X.shape)
    return X, obs, glist


# ----------------------------------------------------------------------------
# 2-3. scDRS
# ----------------------------------------------------------------------------
def ensure_scdrs():
    try:
        import scdrs  # noqa
        return True
    except Exception:
        p = subprocess.run([str(VENV_PY), "-m", "pip", "install", "-q", "scdrs==1.0.2"], capture_output=True, text=True)
        say("[PIP] scdrs==1.0.2 rc=%d %s" % (p.returncode, p.stderr.strip()[-200:]))
        try:
            import importlib
            importlib.invalidate_caches()
            import scdrs  # noqa
            return True
        except Exception as e:
            say("[SCDRS] unavailable: %s %s" % (type(e).__name__, str(e)[:200]))
            return False


def gene_sets(present):
    sym = {}
    for line in open(GENE_LOC_B37):
        x = line.split()
        if len(x) >= 6:
            sym[x[0]] = x[5]
    out = {}
    for t in TRAITS:
        f = MAGMA_OUT / (t + "_35_10.genes.out")
        if not f.exists():
            say("[GS] %s genes.out missing" % t)
            continue
        g = pd.read_csv(f, sep=r"\s+")
        g["SYM"] = g["GENE"].astype(str).map(sym)
        g = g[g["SYM"].isin(present)].sort_values("ZSTAT", ascending=False).drop_duplicates("SYM").head(N_TOP)
        out[t] = (list(g["SYM"]), list(np.minimum(g["ZSTAT"].values, 10.0)))
        say("[GS] %-13s %d genes (z %.2f .. %.2f)" % (t, len(g), g["ZSTAT"].min(), g["ZSTAT"].max()))
    return out


def assoc(norm, ctrl, idx):
    q_obs = np.quantile(norm[idx], 0.95)
    q_ctrl = np.quantile(ctrl[idx], 0.95, axis=0)
    mcp = (1 + np.sum(q_ctrl >= q_obs)) / (1 + ctrl.shape[1])
    mcz = (q_obs - q_ctrl.mean()) / max(q_ctrl.std(), 1e-12)
    return mcp, mcz


def is_macro(label):
    s = str(label)
    return s.startswith("MPhi") or s.startswith("M$\\Phi$") or s.startswith("MΦ")


def run_scdrs(X, obs, genes):
    import anndata as ad
    import scanpy as sc
    import scdrs
    rows = []
    gs_all = gene_sets(set(genes))
    for ctx, members in CONTEXTS.items():
        m = obs["context"].isin(members).to_numpy()
        a = ad.AnnData(X=X[m].copy(), obs=obs.loc[m].reset_index(drop=True))
        a.var_names = genes
        a.obs_names = ["c%d" % i for i in range(a.n_obs)]
        sc.pp.filter_genes(a, min_cells=10)
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        cov = pd.DataFrame(index=a.obs_names, data={"const": 1.0, "n_genes": a.obs["n_genes"].values.astype(float)})
        scdrs.preprocess(a, cov=cov, n_mean_bin=20, n_var_bin=20, copy=False)
        say("")
        say("[SCDRS] context=%s cells=%d genes=%d" % (ctx, a.n_obs, a.n_vars))
        groups = {}
        for lin in sorted(a.obs["lineage"].dropna().unique()):
            groups["lineage:" + lin] = np.where(a.obs["lineage"].values == lin)[0]
        vc = a.obs["label"].value_counts()
        for l in vc.index[vc >= MIN_GROUP_CELLS]:
            groups["label:" + l] = np.where(a.obs["label"].values == l)[0]
        mac = np.where(a.obs["label"].map(is_macro).values)[0]
        if len(mac) >= MIN_GROUP_CELLS:
            groups["set:MACRO_ALL"] = mac
        mye = np.where(a.obs["lineage"].values == "myeloid")[0]
        mk = [g for g in M2_MARKERS if g in set(a.var_names)]
        if len(mye) >= 2 * MIN_GROUP_CELLS and mk:
            Xm = a[mye, mk].X
            Xm = Xm.toarray() if hasattr(Xm, "toarray") else np.asarray(Xm)
            zs = (Xm - Xm.mean(axis=0)) / np.where(Xm.std(axis=0) > 0, Xm.std(axis=0), 1)
            sc_m2 = zs.mean(axis=1)
            groups["set:M2like"] = mye[sc_m2 > np.median(sc_m2)]
            say("    [M2like] markers=%s myeloid=%d M2like=%d" % (",".join(mk), len(mye), len(groups["set:M2like"])))
        gsm = a.obs["gsm"].values
        for t, (gl, gw) in gs_all.items():
            keep = [i for i, g in enumerate(gl) if g in set(a.var_names)]
            t0 = time.time()
            res = scdrs.score_cell(data=a, gene_list=[gl[i] for i in keep], gene_weight=[gw[i] for i in keep],
                                   ctrl_match_key="mean_var", n_ctrl=N_CTRL, weight_opt="vs", return_ctrl_raw_score=False,
                                   return_ctrl_norm_score=True, verbose=False)
            norm = res["norm_score"].to_numpy(dtype=np.float32)
            ctrl = res[[c for c in res.columns if c.startswith("ctrl_norm_score_")]].to_numpy(dtype=np.float32)
            for gname, idx in groups.items():
                p, z = assoc(norm, ctrl, idx)
                ing = np.zeros(len(norm), dtype=bool)
                ing[idx] = True
                diffs = []
                for s in pd.unique(gsm):
                    sel = np.where((gsm == s) & ing)[0]
                    if len(sel) >= MIN_DONOR_CELLS:
                        # donor-level z: group mean of the disease score vs the same mean for the matched control gene sets
                        cm = ctrl[sel].mean(axis=0)
                        diffs.append(float((norm[sel].mean() - cm.mean()) / max(cm.std(), 1e-12)))
                dp = np.nan
                if len(diffs) >= 3:
                    from scipy.stats import ttest_1samp
                    tt = ttest_1samp(diffs, 0.0)
                    dp = float(tt.pvalue / 2 if tt.statistic > 0 else 1 - tt.pvalue / 2)
                rows.append({"context": ctx, "trait": t, "group": gname, "n_cells": len(idx), "mc_p": p, "mc_z": z,
                             "donor_n": len(diffs), "donor_pos": int(np.sum(np.array(diffs) > 0)) if diffs else 0, "donor_p": dp})
            del res, ctrl
            say("    [SCDRS] %-13s %s scored in %.1f min" % (t, ctx, (time.time() - t0) / 60))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "scdrs_group_tests.csv", index=False)
    return df


def report(df):
    say("")
    say("-" * 72)
    say("[SECTION] scDRS association (MC p, 1,000 controls; BH across labels) and donor-level tests")
    for ctx in CONTEXTS:
        sub = df[df["context"] == ctx]
        if not len(sub):
            continue
        say("[CONTEXT] " + ctx)
        lin = sub[sub["group"].str.startswith("lineage:")]
        piv = lin.pivot_table(index="group", columns="trait", values="mc_p")
        for g, r in piv.iterrows():
            say("    [LIN] %-12s %s" % (g.replace("lineage:", "")[:12], " ".join("%s=%.3f" % (t[:6], r.get(t, np.nan)) for t in TRAITS if t in piv.columns)))
        lab = sub[sub["group"].str.startswith("label:")].copy()
        for t in TRAITS:
            s = lab[lab["trait"] == t].copy()
            if len(s):
                s["q"] = bh(s["mc_p"].values)
                s = s.sort_values("mc_p")
                say("    [TOP] %-13s %s" % (t[:13], "; ".join("%s p=%.3f q=%.2f d=%.2g(%d/%d)" % (
                    g.replace("label:", "")[:18], p, q, dp, dpos, dn) for g, p, q, dp, dpos, dn in
                    zip(s["group"].head(5), s["mc_p"].head(5), s["q"].head(5), s["donor_p"].head(5), s["donor_pos"].head(5), s["donor_n"].head(5)))))
        mac_groups = [g for g in sub["group"].unique() if g in ("set:MACRO_ALL", "set:M2like", "lineage:myeloid") or
                      (g.startswith("label:") and is_macro(g[6:]))]
        for t in TRAITS:
            s = sub[(sub["trait"] == t) & (sub["group"].isin(mac_groups))].copy()
            if not len(s):
                continue
            s["holm"] = holm(s["mc_p"].values)
            say("    [MACRO] %-13s %s" % (t[:13], "; ".join("%s p=%.3f H=%.3f donor_p=%.3g(%d/%d)" % (
                g.split(":", 1)[1][:16], p, h, dp, dpos, dn) for g, p, h, dp, dpos, dn in
                zip(s["group"], s["mc_p"], s["holm"], s["donor_p"], s["donor_pos"], s["donor_n"]))))
            ok = s[(s["holm"] < 0.05) & (s["donor_p"] < 0.05)]
            SUM["results"]["%s_%s_macro_hits" % (ctx, t)] = list(ok["group"])
            SUM["results"]["%s_%s_macro_holm" % (ctx, t)] = dict(zip(s["group"], s["holm"]))
    r = SUM["results"]
    hits = set(r.get("lesion_ochoa_GCST90205183_macro_hits", [])) | set(r.get("lesion_ochoa_KOLLER_EUR_macro_hits", []))
    spec = [g for g in hits if r.get("lesion_ochoa_HYPERTENSION_macro_holm", {}).get(g, 1) >= 0.05
            and r.get("lesion_ochoa_T2D_macro_holm", {}).get(g, 1) >= 0.05]
    say("[OCHOA_FAITHFUL_VERDICT] reproduced: %s (%s) | endometriosis-specific: %s (%s)" % (
        "YES" if hits else "NO", ", ".join(sorted(hits)) or "-", "YES" if spec else "NO", ", ".join(sorted(spec)) or "-"))


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.28 | Ochoa-faithful scDRS + donor-level tests | " + TS)
    say("=" * 72)
    t0 = time.time()
    if not ensure_scdrs():
        say("[FATAL] scdrs not available")
    else:
        try:
            X, obs, genes = build_cells()
            say("[CELLS_BY_CONTEXT] " + " | ".join("%s=%d" % (k, v) for k, v in obs["context"].value_counts().items()))
            df = run_scdrs(X, obs, genes)
            report(df)
        except Exception as e:
            import traceback
            say("[STEP_ERROR] %s: %s" % (type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v028_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.18 before this run")
    say("[DONE] PAPER9_V028_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
