#!/usr/bin/env python3
# Paper 9 | v0.0.52 | real peritoneal mesothelium vs endometrial-type stroma (rules CHANGELOG v0.2.42)
#   M1 GSE213216 raw 10x matrices (eutopic endometrium, peritoneal lesions, endometriomas, unaffected peritoneum):
#      marker-based cell groups (mesothelium by a strict multi-marker rule; stroma split by tissue context;
#      perivascular / smooth muscle; endothelium; epithelium; immune), sample-level pseudobulk, sample-paired DE
#   M2 up-type marker families (UP, UPN, SPEC; plus ABS) and MAGMA on six GWAS; H-MESO; trait contrast
#   M3 robustness of the GSE266265 eF3 / dS3 signal to generic gene properties (mean expression, breadth, LOEUF)
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). ~20-40 min (reads 40+ raw matrices once,
# then cached); no downloads.

import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
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
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RAW = REPO / "data" / "raw"
INTERIM = REPO / "data" / "interim"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
RES = REPO / "results"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
X_DIR = INTERIM / "gse213216_extract"
CELLS213 = PROC / "gse213216_cells_v028"
SERIES_MATRIX = RAW / "geo" / "GSE213216" / "GSE213216_series_matrix.txt.gz"
PB_DIR = PROC / "gse213216_pseudobulk_v052"
DE266 = PROC / "genesets" / "gse266265_v050_de"
GNOMAD_FILE = RAW / "gnomad" / "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = RES / "review" / ("p9_v052_" + TS)
WORK = INTERIM / "magma_v052"
LOG = REPO / "logs" / ("p9_v052_mesothelium_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.42)
MIN_GENES, MAX_PCT_MT = 200, 25.0
MESO_ANY = ["MSLN", "UPK3B", "CALB2", "ITLN1", "PRG4"]
MESO_CORE = ["UPK3B", "ITLN1"]
MESO_EXCLUDE = ["PTPRC", "PECAM1", "EPCAM"]
LINEAGE = {"epithelium": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
           "stroma": ["PDGFRA", "DCN", "COL1A1", "COL3A1", "LUM"],
           "endothelium": ["PECAM1", "VWF", "CDH5", "CLDN5"],
           "immune": ["PTPRC", "CD68", "LYZ", "CD3E", "NKG7"],
           "perivascular / smooth muscle": ["ACTA2", "MYH11", "DES", "CNN1", "RGS5"]}
CONTEXTS = ["eutopic", "peritoneal_lesion", "endometrioma", "peritoneum_unaffected"]
MIN_CELLS_UNIT = 10
MIN_SAMPLES = 4
MIN_CPM = 1.0
LFC, PADJ, MIN_SET, N_MATCH = 1.2, 0.01, 20, 400
FAMS = ["ABS", "UP", "UPN", "SPEC"]
TRAITS = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL", "LEIOMYOMA", "HYPERTENSION", "T2D"]
ENDO_STROMA = ["stroma (eutopic)", "stroma (peritoneal_lesion)", "stroma (endometrioma)"]
SEED = 20261108

for d in (OUT, WORK, LOG.parent, PB_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V052_MESOTHELIUM", "timestamp_utc": TS, "results": {}}
ERR = []


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


def gene_loc():
    sym2ent, dup = {}, set()
    for line in open(GENE_LOC_B37):
        t = line.split()
        if len(t) >= 6 and t[1] in [str(c) for c in range(1, 23)]:
            if t[5] in sym2ent:
                dup.add(t[5])
            sym2ent[t[5]] = t[0]
    for d in dup:
        sym2ent.pop(d, None)
    return sym2ent


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def parse_full(gsa):
    rows = [l for l in open(gsa) if not l.startswith("#") and l.strip()]
    hdr = rows[0].split()
    out = {}
    for l in rows[1:]:
        d = dict(zip(hdr, l.split()))
        k = d.get("FULL_NAME", d.get("VARIABLE"))
        try:
            out[k] = {"p": float(d["P"]), "beta": float(d.get("BETA", "nan")), "beta_std": float(d.get("BETA_STD", "nan")), "se": float(d.get("SE", "nan"))}
        except ValueError:
            continue
    return out


def magma_covar(raw, cf, tag, cond):
    outp = WORK / tag
    gsa = Path(str(outp) + ".gsa.out")
    model = (["condition-hide=" + ",".join(cond)] if cond else []) + ["direction=greater"]
    if gsa.exists():
        gsa.unlink()
    rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--gene-covar", cf, "--model"] + model + ["--out", outp], timeout=7200)
    if gsa.exists():
        return parse_full(gsa)
    say("    [MAGMA_FAIL] %s %s" % (tag, (o + e).strip()[-300:]))
    return None


def tsf(t, df):
    try:
        from scipy.stats import t as T
        return 2.0 * T.sf(np.abs(t), df)
    except Exception:
        return np.array([math.erfc(abs(v) / math.sqrt(2)) for v in t])


# ----------------------------------------------------------------------------
# M1 sample contexts and marker-based groups
# ----------------------------------------------------------------------------
def contexts():
    ctx, pat = {}, {}
    f = CELLS213 / "obs.csv"
    if f.exists():
        O = pd.read_csv(f, usecols=["gsm", "context", "patient"])
        O = O.drop_duplicates("gsm")
        ctx.update(O.set_index("gsm")["context"].to_dict())
        pat.update(O.set_index("gsm")["patient"].to_dict())
    if SERIES_MATRIX.exists():
        txt = gzip.open(SERIES_MATRIX, "rt", errors="ignore").read().split("!series_matrix_table_begin")[0]
        rows = {}
        for line in txt.splitlines():
            if line.startswith("!Sample_"):
                rows.setdefault(line.split("\t")[0], []).append([v.strip().strip('"') for v in line.split("\t")[1:]])
        gsm = rows["!Sample_geo_accession"][0]
        src = rows.get("!Sample_source_name_ch1", [[""] * len(gsm)])[0]
        lib = rows.get("!Sample_library_strategy", [[""] * len(gsm)])[0]
        chars = rows.get("!Sample_characteristics_ch1", [])
        for i, g in enumerate(gsm):
            if g in ctx:
                continue
            t = (src[i] + " " + " ".join(c[i] for c in chars if i < len(c))).lower()
            if "visium" in t or "spatial" in t or str(lib[i]).upper() == "OTHER":
                continue
            c = ("eutopic" if "eutopic" in t else "endometrioma" if "endometrioma" in t else "peritoneum_unaffected" if "no endometriosis" in t
                 else "ovary_unaffected" if "ovary" in t else "peritoneal_lesion" if "endometriosis" in t else "other")
            ctx[g] = c
    return ctx, pat


def find_10x(d):
    mtx = sorted(d.rglob("*matrix.mtx*"))
    if mtx:
        name = mtx[0].name
        return ("mtx", mtx[0].parent, name[:name.index("matrix.mtx")])
    h5 = sorted(d.rglob("*.h5"))
    if h5:
        return ("h5", h5[0], "")
    return None


def build_pseudobulk():
    f_pb, f_meta = PB_DIR / "pseudobulk.csv.gz", PB_DIR / "units.csv"
    if f_pb.exists() and f_meta.exists():
        say("[M1] pseudobulk cached (%s)" % PB_DIR.name)
        return pd.read_csv(f_pb, index_col=0), pd.read_csv(f_meta)
    import scanpy as sc
    import scipy.sparse as sp
    ctx, pat = contexts()
    dirs = sorted(p for p in X_DIR.iterdir() if p.is_dir())
    units, cols = [], {}
    tot = {"cells": 0, "meso": 0}
    t0 = time.time()
    for d in dirs:
        c = ctx.get(d.name)
        if c not in CONTEXTS:
            say("    [SKIP] %s context %s" % (d.name, c))
            continue
        tx = find_10x(d)
        if tx is None:
            say("    [SKIP] %s no 10x matrix" % d.name)
            continue
        try:
            A = sc.read_10x_mtx(tx[1], var_names="gene_symbols", make_unique=True, prefix=tx[2] or None) if tx[0] == "mtx" else sc.read_10x_h5(tx[1])
        except Exception as e:
            say("    [SKIP] %s read error %s" % (d.name, type(e).__name__))
            continue
        A.var_names_make_unique()
        X = A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)
        genes = np.array(A.var_names)
        gi = {g: i for i, g in enumerate(genes)}
        ng = np.asarray((X > 0).sum(axis=1)).ravel()
        lib = np.asarray(X.sum(axis=1)).ravel()
        mt = np.array([g.upper().startswith("MT-") for g in genes])
        pmt = np.asarray(X[:, mt].sum(axis=1)).ravel() / np.maximum(lib, 1) * 100 if mt.any() else np.zeros(len(lib))
        keep = (ng >= MIN_GENES) & (pmt < MAX_PCT_MT)
        X, lib = X[keep], lib[keep]
        if X.shape[0] < 100:
            say("    [SKIP] %s only %d cells after QC" % (d.name, X.shape[0]))
            continue

        mk = [g for g in dict.fromkeys(MESO_ANY + MESO_CORE + MESO_EXCLUDE) if g in gi]
        Mk = X[:, [gi[g] for g in mk]].toarray() > 0 if mk else np.zeros((X.shape[0], 0), dtype=bool)
        mpos = {g: j for j, g in enumerate(mk)}

        def det(g):
            return Mk[:, mpos[g]] if g in mpos else np.zeros(X.shape[0], dtype=bool)
        n_any = sum(det(g).astype(int) for g in MESO_ANY)
        core = np.zeros(X.shape[0], dtype=bool)
        for g in MESO_CORE:
            core |= det(g)
        excl = np.zeros(X.shape[0], dtype=bool)
        for g in MESO_EXCLUDE:
            excl |= det(g)
        meso = (n_any >= 3) & core & ~excl
        norm = sp.diags(1e4 / np.maximum(lib, 1)) @ X
        score = {}
        for ln, ms in LINEAGE.items():
            ii = [gi[m] for m in ms if m in gi]
            score[ln] = np.asarray(np.log1p(norm[:, ii].toarray()).mean(axis=1)).ravel() if ii else np.zeros(X.shape[0])
        S = pd.DataFrame(score)
        best = S.idxmax(axis=1).values
        best[S.max(axis=1).values <= 0] = "unassigned"
        grp = np.where(meso, "mesothelium", best).astype(object)
        grp[grp == "stroma"] = "stroma (%s)" % c
        tot["cells"] += X.shape[0]
        tot["meso"] += int(meso.sum())
        cols["%s|__total__" % d.name] = pd.Series(np.asarray(X.sum(axis=0)).ravel(), index=genes)
        units.append({"unit": "%s|__total__" % d.name, "gsm": d.name, "context": c, "patient": pat.get(d.name), "group": "__total__", "n_cells": X.shape[0]})
        for g in sorted(set(grp)):
            if g == "unassigned":
                continue
            idx = np.where(grp == g)[0]
            if len(idx) < MIN_CELLS_UNIT:
                continue
            key = "%s|%s" % (d.name, g)
            cols[key] = pd.Series(np.asarray(X[idx].sum(axis=0)).ravel(), index=genes)
            units.append({"unit": key, "gsm": d.name, "context": c, "patient": pat.get(d.name), "group": g, "n_cells": len(idx)})
        say("    [M1] %s %-22s cells %5d | mesothelium %4d | %s (%.0f min)" % (
            d.name, c, X.shape[0], int(meso.sum()), ", ".join("%s %d" % (g[:10], (grp == g).sum()) for g in sorted(set(grp)) if g != "unassigned"),
            (time.time() - t0) / 60))
    PB = pd.DataFrame(cols).fillna(0.0)
    PB = PB.groupby(level=0).sum()
    U = pd.DataFrame(units)
    PB.to_csv(PB_DIR / "pseudobulk.csv.gz")
    U.to_csv(PB_DIR / "units.csv", index=False)
    say("[M1] built: %d units, %d genes; %d QC cells, %d mesothelial cells" % (PB.shape[1], PB.shape[0], tot["cells"], tot["meso"]))
    return PB, U


def de(PB, U):
    TU = U[U["group"] == "__total__"].set_index("gsm")["unit"]
    tot = PB[TU.values].copy()
    tot.columns = TU.index
    U = U[U["group"] != "__total__"]
    rows, cpm_t = [], {}
    elig = U.groupby("group")["gsm"].nunique()
    groups = [g for g in elig.index if elig[g] >= MIN_SAMPLES]
    say("[M1] units per group (samples): %s" % ", ".join("%s %d" % (g, elig[g]) for g in elig.index))
    say("[M1] eligible groups (>= %d samples with >= %d cells): %s" % (MIN_SAMPLES, MIN_CELLS_UNIT, ", ".join(groups)))
    genes = PB.index.values
    for g in groups:
        u = U[U["group"] == g]
        X = PB[u["unit"]].values.T
        Rr = tot[u["gsm"]].values.T - X
        Ct = X / np.maximum(X.sum(axis=1, keepdims=True), 1) * 1e6
        Cr = Rr / np.maximum(Rr.sum(axis=1, keepdims=True), 1) * 1e6
        cpm_t[g] = Ct.mean(axis=0)
        ex = (Ct.mean(axis=0) >= MIN_CPM) | (Cr.mean(axis=0) >= MIN_CPM)
        D = np.log2(Ct[:, ex] + 1) - np.log2(Cr[:, ex] + 1)
        n = D.shape[0]
        m, sd = D.mean(axis=0), D.std(axis=0, ddof=1)
        tt = m / np.maximum(sd / math.sqrt(n), 1e-9)
        p = tsf(tt, n - 1)
        rows.append(pd.DataFrame({"group": g, "gene": genes[ex], "lfc": m, "score": tt, "p": p, "padj": bh(p), "n_samples": n}))
    DE = pd.concat(rows, ignore_index=True)
    CPM = pd.DataFrame(cpm_t, index=genes)
    CPM = CPM[(CPM >= MIN_CPM).any(axis=1)]
    DE.to_csv(OUT / "m1_de_all.csv.gz", index=False)
    CPM.to_csv(OUT / "m1_mean_cpm.csv.gz")
    return DE, CPM


def families(DE, CPM, sym2ent):
    uni = sorted({g for g in DE["gene"].astype(str) if g in sym2ent})
    uset = set(uni)
    meta, cols = {}, {}
    for g, d in DE.groupby("group", sort=False):
        d = d.assign(gene=d["gene"].astype(str))
        sig = d["padj"] < PADJ
        sel = {"ABS": d.loc[(d["lfc"].abs() > LFC) & sig, "gene"], "UP": d.loc[(d["lfc"] > LFC) & sig, "gene"],
               "UPN": d.sort_values("score", ascending=False)["gene"].head(N_MATCH)}
        for fam, gl in sel.items():
            gl = [x for x in gl if x in uset]
            if len(gl) >= MIN_SET:
                k = "V%03d" % (len(meta) + 1)
                meta[k], cols[k] = (fam, g), gl
    T = pd.DataFrame({k: np.isin(uni, gl).astype(int) for k, gl in cols.items()}, index=[sym2ent[g] for g in uni])
    T = T[~T.index.duplicated()]
    C = CPM.copy()
    C.index = C.index.astype(str)
    C = C.groupby(level=0).mean()
    C = C.loc[[g for g in C.index if g in sym2ent]]
    S = C.div(C.sum(axis=1), axis=0)
    spec, sd = {}, {}
    for i, c in enumerate(S.columns):
        k = "S%02d" % (i + 1)
        spec[k], sd[k] = c, S[c].values
    sd["Average"] = np.log2(C.mean(axis=1).values + 1.0)
    Ts = pd.DataFrame(sd, index=[sym2ent[g] for g in C.index])
    Ts = Ts[~Ts.index.duplicated()]
    return T, Ts, meta, spec


def ranked(res, meta, spec, fam):
    it = [(spec[k], v) for k, v in res.items() if k in spec] if fam == "SPEC" else [(meta[k][1], v) for k, v in res.items() if k in meta and meta[k][0] == fam]
    it.sort(key=lambda x: x[1]["p"])
    q = bh([v["p"] for _, v in it]) if it else []
    return [(l, v, qq) for (l, v), qq in zip(it, q)]


def m2(sym2ent):
    PB, U = build_pseudobulk()
    if PB is None or not len(U):
        say("[M2] no pseudobulk")
        return
    say("[M1] contexts: %s" % ", ".join("%s %d samples" % (k, v) for k, v in U.groupby("context")["gsm"].nunique().items()))
    DE, CPM = de(PB, U)
    say("")
    say("-" * 72)
    say("[SECTION] M1 group identity check (mean CPM)")
    for g in ["MSLN", "UPK3B", "CALB2", "ITLN1", "PRG4", "KRT19", "EPCAM", "PDGFRA", "DCN", "ESR1", "PGR", "CD10", "MME", "WT1", "HOXA10", "IGF1", "SFRP1", "PTPRC"]:
        if g in CPM.index:
            say("    %-7s %s" % (g, " | ".join("%s %.0f" % (c[:24], CPM.loc[g, c]) for c in CPM.columns)))
    T, Ts, meta, spec = families(DE, CPM, sym2ent)
    cf, cfs = WORK / "sets.txt", WORK / "spec.txt"
    T.to_csv(cf, sep="\t", index_label="GENE")
    Ts.to_csv(cfs, sep="\t", float_format="%.6g", index_label="GENE")
    res = {}
    for t in TRAITS:
        raw = MAGMA_OUT / (t + "_35_10.genes.raw")
        if not raw.exists():
            continue
        r0 = magma_covar(raw, cf, "m_sets_%s" % t, [])
        r1 = magma_covar(raw, cfs, "m_spec_%s" % t, ["Average"])
        if r0 and r1:
            res[t] = {fam: ranked(r1 if fam == "SPEC" else r0, meta, spec, fam) for fam in FAMS}
    rows = []
    for t, F in res.items():
        for fam, it in F.items():
            for i, (l, v, q) in enumerate(it):
                rows.append({"trait": t, "family": fam, "rank": i + 1, "label": l, "p": v["p"], "q": q, "beta_std": v["beta_std"], "beta": v["beta"], "se": v["se"]})
    A = pd.DataFrame(rows)
    A.to_csv(OUT / "m2_rankings.csv", index=False)
    say("")
    say("-" * 72)
    say("[SECTION] M2 rankings of marker-based groups (GSE213216; p / BH q)")
    for t in res:
        say("  %s" % t)
        for fam in FAMS:
            say("    %-4s %s" % (fam, "; ".join("%s %.0e/%.2f" % (l, v["p"], q) for l, v, q in res[t][fam])))
    k = A[(A["trait"] == "KOLLER_EUR") & (A["family"] == "UPN")].set_index("label")
    if len(k):
        es = [g for g in ENDO_STROMA if g in k.index]
        best_es = min(es, key=lambda g: k.loc[g, "rank"]) if es else None
        rm = k.loc["mesothelium", "rank"] if "mesothelium" in k.index else None
        rp = k.loc["stroma (peritoneum_unaffected)", "rank"] if "stroma (peritoneum_unaffected)" in k.index else None
        ok = best_es is not None and rm is not None and k.loc[best_es, "rank"] < rm and (rp is None or k.loc[best_es, "rank"] < rp)
        say("[H_MESO_VERDICT] KOLLER_EUR UPN: best endometrial-type stroma %s rank %s (p %.1e) | mesothelium rank %s (p %s) | peritoneal fibroblasts rank %s -> %s" % (
            best_es, k.loc[best_es, "rank"] if best_es else "-", k.loc[best_es, "p"] if best_es else float("nan"), rm,
            ("%.1e" % k.loc["mesothelium", "p"]) if rm is not None else "-", rp, "SUPPORTED" if ok else "NOT SUPPORTED"))
        SUM["results"]["h_meso"] = {"supported": bool(ok), "best_endometrial_stroma": best_es, "mesothelium_rank": None if rm is None else int(rm)}
        if best_es:
            est = {}
            for t in TRAITS:
                d = A[(A["trait"] == t) & (A["family"] == "UPN") & (A["label"] == best_es)]
                if len(d):
                    r = d.iloc[0]
                    est[t] = (r["beta_std"], abs(r["se"] * r["beta_std"] / r["beta"]) if r["beta"] else np.nan)
            parts = []
            for c in ("T2D", "HYPERTENSION", "LEIOMYOMA"):
                if "KOLLER_EUR" in est and c in est:
                    z = (est["KOLLER_EUR"][0] - est[c][0]) / math.sqrt(est["KOLLER_EUR"][1] ** 2 + est[c][1] ** 2)
                    parts.append("KOLLER>%s z %.2f p %.3f" % (c[:6], z, 0.5 * math.erfc(z / math.sqrt(2))))
            say("  [M2] trait contrast, UPN %s: %s" % (best_es, " | ".join(parts)))


# ----------------------------------------------------------------------------
# M3 GSE266265 eF3 / dS3 vs generic gene properties
# ----------------------------------------------------------------------------
def m3(sym2ent):
    say("")
    say("-" * 72)
    say("[SECTION] M3 GSE266265 UPN eF3 / dS3 conditioned on mean expression, breadth and LOEUF (as v045)")
    if not (DE266 / "de_all.csv.gz").exists():
        say("[M3] v050 DE tables missing")
        return
    DE = pd.read_csv(DE266 / "de_all.csv.gz")
    CPM = pd.read_csv(DE266 / "mean_cpm_by_type.csv.gz", index_col=0)
    CPM.index = CPM.index.astype(str)
    C = CPM.groupby(level=0).mean()
    C = C[(C >= 1.0).any(axis=1)]
    lo = None
    if GNOMAD_FILE.exists():
        G = pd.read_csv(GNOMAD_FILE, sep="\t", compression="gzip")
        sc_ = next(c for c in ("gene", "gene_symbol") if c in G.columns)
        lc = next(c for c in ("oe_lof_upper", "LOEUF") if c in G.columns)
        lo = pd.to_numeric(G.set_index(sc_)[lc], errors="coerce").dropna().groupby(level=0).min()
    uni = [g for g in C.index if g in sym2ent and (lo is None or g in lo.index)]
    cols = {"Average": np.log2(C.loc[uni].mean(axis=1).values + 1), "Breadth": (C.loc[uni] >= 1).mean(axis=1).values}
    cov = ["Average", "Breadth"]
    if lo is not None:
        cols["LOEUF"] = lo.reindex(uni).values
        cov.append("LOEUF")
    tys = ["eF3", "dS3", "Prv-MYH11", "VSMC"]
    for ty in tys:
        d = DE[DE["group"] == ty]
        if len(d):
            cols["T_" + re.sub(r"[^A-Za-z0-9]", "_", ty)] = np.isin(uni, list(d.sort_values("score", ascending=False)["gene"].astype(str).head(N_MATCH))).astype(int)
    T = pd.DataFrame(cols, index=[sym2ent[g] for g in uni])
    T = T[~T.index.duplicated()]
    cf = WORK / "m3.txt"
    T.to_csv(cf, sep="\t", float_format="%.6g", index_label="GENE")
    out = {}
    for t in ("KOLLER_EUR", "GCST90205183", "ENDO_ALL"):
        raw = MAGMA_OUT / (t + "_35_10.genes.raw")
        if not raw.exists():
            continue
        b = magma_covar(raw, cf, "m3_%s_b" % t, [])
        c = magma_covar(raw, cf, "m3_%s_c" % t, cov)
        if b and c:
            parts = []
            for k in [x for x in T.columns if x.startswith("T_")]:
                parts.append("%s %.1e -> %.1e" % (k[2:], b[k]["p"], c[k]["p"]))
                out[(t, k[2:])] = c[k]["p"]
            say("  [M3] %-12s base -> conditioned: %s" % (t, " | ".join(parts)))
    ok = out.get(("KOLLER_EUR", "eF3"), 1) < 0.05
    say("[H_EF3_ROBUST_VERDICT] KOLLER_EUR UPN eF3 keeps p < 0.05 after conditioning on %s -> %s" % (", ".join(cov), "SUPPORTED" if ok else "NOT SUPPORTED"))
    SUM["results"]["h_ef3_robust"] = ok


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.52 | peritoneal mesothelium vs endometrial-type stroma; eF3 robustness | " + TS)
    say("=" * 72)
    t0 = time.time()
    sym2ent = gene_loc()
    for name, fn in (("M3", lambda: m3(sym2ent)), ("M1_M2", lambda: m2(sym2ent))):
        try:
            fn()
        except Exception as e:
            import traceback
            ERR.append(name)
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v052_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.42 before this run")
    if ERR:
        say("[INCOMPLETE] steps with errors: " + ", ".join(ERR))
    else:
        say("[DONE] PAPER9_V052_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
