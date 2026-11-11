#!/usr/bin/env python3
# Paper 9 | v0.0.65 | independent lesion cohort (GSE179640, Tan et al. 2022) + scope of non-European Koller files
#   (rules CHANGELOG v0.2.55, fixed before the run)
#   Z  Koller et al. Zenodo record: list the files and flag non-European / multi-ancestry summary statistics (no download).
#   T0 GSE179640 download (series matrix + RAW tar of per-sample 10x matrices), tissue context per sample from GEO
#      metadata by a fixed rule; organoids and unknown samples excluded.
#   T1 marker-based cell groups exactly as v052 (strict mesothelium rule; stroma split by tissue context), sample pseudobulk.
#   R1 H-LESION-REP: genes of the KOLLER_EUR stromal loci stroma-enriched within peritoneal lesions (and endometriomas if
#      evaluable) beyond expression-matched random genes (v060 method).
#   R2 H-HOST: stromal-locus genes higher in lesion stroma than in stroma of lesion-adjacent peritoneum (host fibroblasts).
#   R3 H-MESO-REP: MAGMA ranking of the marker-based groups (v052 method): best endometrial-type stroma above mesothelium
#      and above adjacent-peritoneum stroma for KOLLER_EUR (size-matched up-only).
# Runs in ~/Desktop/paper9/.venv, low heat (MAGMA single thread, nice 10). Download size printed (likely 2-6 GB);
# reading the matrices ~20-40 min, MAGMA ~10-20 min. Re-running re-uses the download and the pseudobulk cache.

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
import shutil
import subprocess
import tarfile
import time
import urllib.request
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
GEO_DIR = RAW / "geo" / "GSE179640"
X_DIR = INTERIM / "gse179640_extract"
PB_DIR = PROC / "gse179640_pseudobulk_v065"
LOCUS_DIR = RES / "locus_atlas"
GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179640/"
ZENODO_API = "https://zenodo.org/api/records/18983492"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = RES / "review" / ("p9_v065_" + TS)
WORK = INTERIM / "magma_v065"
LOG = REPO / "logs" / ("p9_v065_tan_replication_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.55); cell-group rules identical to v052
MIN_GENES, MAX_PCT_MT = 200, 25.0
MESO_ANY = ["MSLN", "UPK3B", "CALB2", "ITLN1", "PRG4"]
MESO_CORE = ["UPK3B", "ITLN1"]
MESO_EXCLUDE = ["PTPRC", "PECAM1", "EPCAM"]
LINEAGE = {"epithelium": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
           "stroma": ["PDGFRA", "DCN", "COL1A1", "COL3A1", "LUM"],
           "endothelium": ["PECAM1", "VWF", "CDH5", "CLDN5"],
           "immune": ["PTPRC", "CD68", "LYZ", "CD3E", "NKG7"],
           "perivascular / smooth muscle": ["ACTA2", "MYH11", "DES", "CNN1", "RGS5"]}
CONTEXTS = ["eutopic", "control_endometrium", "peritoneal_lesion", "peritoneum_adjacent", "endometrioma"]
CTX_NAME = {"eutopic": "eutopic endometrium", "control_endometrium": "control endometrium", "peritoneal_lesion": "peritoneal lesions",
            "peritoneum_adjacent": "lesion-adjacent peritoneum", "endometrioma": "endometriomas"}
LESIONS = ["peritoneal_lesion", "endometrioma"]
ENDO_STROMA = ["stroma (eutopic)", "stroma (control_endometrium)", "stroma (peritoneal_lesion)", "stroma (endometrioma)"]
LINEAGES = ["stroma", "perivascular / smooth muscle", "endothelium", "epithelium", "immune", "mesothelium"]
COMP_TO_LINEAGE = {"stroma": "stroma", "perivascular/SM": "perivascular / smooth muscle", "endothelium": "endothelium",
                   "epithelium": "epithelium", "immune": "immune"}
COMP_ORDER = ["stroma", "perivascular/SM", "endothelium", "epithelium", "immune"]
TRAIT_LOCI = "KOLLER_EUR"
TRAITS = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL", "LEIOMYOMA", "HYPERTENSION", "T2D"]
MIN_CELLS_UNIT = 10
MIN_KEYS = 4
MIN_GENES_SET = 3
MIN_CPM = 1.0
N_NULL, N_BINS = 10000, 10
LFC, PADJ, MIN_SET, N_MATCH = 1.2, 0.01, 20, 400
FAMS = ["ABS", "UP", "UPN", "SPEC"]
SEED = 20261112

for d in (OUT, WORK, LOG.parent, PB_DIR, GEO_DIR, X_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V065_TAN_REPLICATION", "timestamp_utc": TS, "results": {}}


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


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "-"
    return "%.3f" % p if p >= 0.001 else "%.0e" % p


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
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def tsf2(t, df):
    from scipy.stats import t as T
    return 2.0 * T.sf(np.abs(t), df)


def tsf1(t, df):
    from scipy.stats import t as T
    return float(T.sf(t, df))


def parse_full(gsa):
    rows = [l for l in open(gsa) if not l.startswith("#") and l.strip()]
    hdr = rows[0].split()
    out = {}
    for l in rows[1:]:
        d = dict(zip(hdr, l.split()))
        k = d.get("FULL_NAME", d.get("VARIABLE"))
        try:
            out[k] = {"p": float(d["P"]), "beta": float(d.get("BETA", "nan")), "beta_std": float(d.get("BETA_STD", "nan")),
                      "se": float(d.get("SE", "nan"))}
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


def fetch(url, dst, timeout=6 * 3600):
    if dst.exists() and dst.stat().st_size > 0:
        return True
    part = dst.with_name(dst.name + ".part")
    rc, o, e = sh(["curl", "-L", "-sS", "--fail", "--retry", "3", "--connect-timeout", "30", "-C", "-", "-A", "paper9", "-o", part, url],
                  timeout=timeout)
    if rc == 0 and part.exists() and part.stat().st_size > 0:
        part.rename(dst)
        return True
    say("  [DL] FAILED %s rc=%s %s" % (url.split("/")[-1], rc, (e or o).strip()[:200]))
    return False


def listing(url):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "paper9"}), timeout=60) as r:
            html = r.read().decode("utf-8", "ignore")
        return sorted(set(re.findall(r'href="([^"?/][^"]*)"', html)))
    except Exception as e:
        say("  [LIST] %s failed: %s" % (url, type(e).__name__))
        return []


# ----------------------------------------------------------------------------
# Z  Koller Zenodo record: which ancestries are public
# ----------------------------------------------------------------------------
def part_z():
    say("")
    say("-" * 72)
    say("[SECTION] Z Koller et al. Zenodo record (doi:10.5281/zenodo.18983492): files and ancestries")
    try:
        with urllib.request.urlopen(urllib.request.Request(ZENODO_API, headers={"User-Agent": "paper9"}), timeout=60) as r:
            rec = json.loads(r.read().decode("utf-8"))
    except Exception as e:
        say("[Z] record not reachable: %s" % type(e).__name__)
        return
    files = rec.get("files", []) or rec.get("entries", [])
    flag = re.compile(r"EAS|AFR|AMR|SAS|CSA|MID|ASIAN|AFRICAN|HISPANIC|LATIN|MULTI|TRANS|META|ALL[_.-]?ANC|MAMA|MR-?MEGA", re.I)
    say("[Z] title: %s" % str(rec.get("metadata", {}).get("title", "-"))[:150])
    rows = []
    for f in files:
        k = f.get("key") or f.get("filename") or "?"
        sz = f.get("size") or f.get("filesize") or 0
        rows.append({"file": k, "size_mb": round(sz / 1e6, 1), "non_eur_or_multi": bool(flag.search(k))})
        say("  [Z] %-70s %9.1f MB %s" % (k[:70], sz / 1e6, "<- non-European / multi-ancestry?" if flag.search(k) else ""))
    pd.DataFrame(rows).to_csv(OUT / "z_koller_zenodo_files.csv", index=False)
    n = sum(r["non_eur_or_multi"] for r in rows)
    say("[Z] %d files; %d look non-European or multi-ancestry (by file name; check the record description)" % (len(rows), n))
    SUM["results"]["Z"] = {"n_files": len(rows), "n_flagged": n}


# ----------------------------------------------------------------------------
# T0 GSE179640 download and sample contexts
# ----------------------------------------------------------------------------
def series_samples():
    mats = [h for h in listing(GEO_BASE + "matrix/") if h.endswith("_series_matrix.txt.gz")]
    if not mats:
        mats = ["GSE179640_series_matrix.txt.gz"]
    rows = []
    for m in mats:
        f = GEO_DIR / m
        if not fetch(GEO_BASE + "matrix/" + m, f, timeout=1800):
            continue
        txt = gzip.open(f, "rt", errors="ignore").read().split("!series_matrix_table_begin")[0]
        R = {}
        for line in txt.splitlines():
            if line.startswith("!Sample_"):
                R.setdefault(line.split("\t")[0], []).append([v.strip().strip('"') for v in line.split("\t")[1:]])
        gsm = R.get("!Sample_geo_accession", [[]])[0]
        for i, g in enumerate(gsm):
            get = lambda k: " | ".join(c[i] for c in R.get(k, []) if i < len(c))
            rows.append({"gsm": g, "title": get("!Sample_title"), "source": get("!Sample_source_name_ch1"),
                         "chars": get("!Sample_characteristics_ch1"), "strategy": get("!Sample_library_strategy"),
                         "molecule": get("!Sample_molecule_ch1")})
    return pd.DataFrame(rows).drop_duplicates("gsm") if rows else pd.DataFrame()


def context_of(r):
    # fixed rule (v0.2.55); first match wins
    t = (" ".join([str(r["title"]), str(r["source"]), str(r["chars"])])).lower()
    tt = re.sub(r"[^a-z0-9]+", " ", t)
    if "organoid" in t or re.search(r"\beo[r]?\d", tt) or "spatial" in t or "imaging" in t:
        return "excluded"
    if "adjacent" in t or re.search(r"\becpa\b", tt):
        return "peritoneum_adjacent"
    if re.search(r"\beco\b", tt) or "ovar" in t or "endometrioma" in t:
        return "endometrioma"
    if re.search(r"\becp\b", tt) or "peritone" in t or "ectopic" in t:
        return "peritoneal_lesion"
    if re.search(r"\bctrl\b", tt) or "control" in t or "healthy" in t or "no endometriosis" in t:
        return "control_endometrium"
    if re.search(r"\beue\b", tt) or "eutopic" in t or "endometri" in t:
        return "eutopic"
    return "other"


def patient_of(r):
    for part in str(r["chars"]).split("|"):
        m = re.match(r"\s*(patient|donor|individual|subject|participant|patient id)\s*[:=]\s*(.+)", part.strip(), re.I)
        if m:
            return m.group(2).strip()
    m = re.search(r"(?:patient|pt|donor)[ _-]?([A-Za-z]*\d+)", str(r["title"]), re.I)
    return m.group(1) if m else None


def part_t0():
    say("")
    say("-" * 72)
    say("[SECTION] T0 GSE179640 (Tan et al. 2022) download and tissue contexts")
    S = series_samples()
    if S.empty:
        say("[T0] series matrix not available")
        return None
    S["context"] = [context_of(r) for _, r in S.iterrows()]
    S["patient"] = [patient_of(r) for _, r in S.iterrows()]
    S.to_csv(OUT / "t0_samples.csv", index=False)
    for _, r in S.iterrows():
        say("  [T0] %-11s %-20s patient %-8s | %s | %s | %s" % (r["gsm"], r["context"], r["patient"], r["title"][:28], r["source"][:28], r["chars"][:70]))
    say("[T0] samples per context: %s" % ", ".join("%s %d" % (k, v) for k, v in S["context"].value_counts().items()))
    sup = listing(GEO_BASE + "suppl/")
    say("[T0] supplementary files: %s" % (", ".join(sup) or "listing failed"))
    tar = GEO_DIR / "GSE179640_RAW.tar"
    if not any(X_DIR.iterdir()):
        if "GSE179640_RAW.tar" not in sup and sup:
            say("[T0] no RAW tar in the series supplement")
            return S
        t0 = time.time()
        if not fetch(GEO_BASE + "suppl/GSE179640_RAW.tar", tar):
            return S
        say("[T0] RAW tar %.0f MB in %.0f min" % (tar.stat().st_size / 1e6, (time.time() - t0) / 60))
        with tarfile.open(tar) as tf:
            for m in tf.getmembers():
                if m.isfile() and re.match(r"^GSM\d+_", Path(m.name).name):
                    g = Path(m.name).name.split("_")[0]
                    if S.loc[S["gsm"] == g, "context"].isin(CONTEXTS).any():
                        (X_DIR / g).mkdir(exist_ok=True)
                        with tf.extractfile(m) as src, open(X_DIR / g / Path(m.name).name, "wb") as dst:
                            shutil.copyfileobj(src, dst)
        say("[T0] extracted %d samples to %s" % (sum(1 for p in X_DIR.iterdir() if p.is_dir()), X_DIR.relative_to(REPO)))
    return S


def find_10x(d):
    h5 = sorted(d.rglob("*.h5"))
    h5f = [p for p in h5 if "filtered" in p.name.lower()]
    if h5f or h5:
        return ("h5", (h5f or h5)[0], "")
    mtx = sorted(d.rglob("*matrix.mtx*"))
    if mtx:
        name = mtx[0].name
        return ("mtx", mtx[0].parent, name[:name.index("matrix.mtx")])
    return None


# ----------------------------------------------------------------------------
# T1 marker-based cell groups and pseudobulk (v052 rules)
# ----------------------------------------------------------------------------
def build_pseudobulk(S):
    f_pb, f_meta = PB_DIR / "pseudobulk.csv.gz", PB_DIR / "units.csv"
    if f_pb.exists() and f_meta.exists():
        say("[T1] pseudobulk cached (%s)" % PB_DIR.name)
        return pd.read_csv(f_pb, index_col=0), pd.read_csv(f_meta)
    import scanpy as sc
    import scipy.sparse as sp
    ctx = S.set_index("gsm")["context"].to_dict()
    pat = S.set_index("gsm")["patient"].to_dict()
    units, cols = [], {}
    tot = {"cells": 0, "meso": 0}
    t0 = time.time()
    for d in sorted(p for p in X_DIR.iterdir() if p.is_dir()):
        c = ctx.get(d.name)
        if c not in CONTEXTS:
            continue
        tx = find_10x(d)
        if tx is None:
            say("    [SKIP] %s no 10x matrix" % d.name)
            continue
        try:
            A = sc.read_10x_h5(tx[1]) if tx[0] == "h5" else sc.read_10x_mtx(tx[1], var_names="gene_symbols", make_unique=True, prefix=tx[2] or None)
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
        Sc = pd.DataFrame(score)
        best = Sc.idxmax(axis=1).values
        best[Sc.max(axis=1).values <= 0] = "unassigned"
        grp = np.where(meso, "mesothelium", best).astype(object)
        grp[grp == "stroma"] = "stroma (%s)" % c
        tot["cells"] += X.shape[0]
        tot["meso"] += int(meso.sum())
        cols["%s|__total__" % d.name] = pd.Series(np.asarray(X.sum(axis=0)).ravel(), index=genes)
        units.append({"unit": "%s|__total__" % d.name, "gsm": d.name, "context": c, "patient": pat.get(d.name), "group": "__total__",
                      "n_cells": X.shape[0]})
        for g in sorted(set(grp)):
            if g == "unassigned":
                continue
            idx = np.where(grp == g)[0]
            if len(idx) < MIN_CELLS_UNIT:
                continue
            key = "%s|%s" % (d.name, g)
            cols[key] = pd.Series(np.asarray(X[idx].sum(axis=0)).ravel(), index=genes)
            units.append({"unit": key, "gsm": d.name, "context": c, "patient": pat.get(d.name), "group": g, "n_cells": len(idx)})
        say("    [T1] %s %-20s cells %5d | mesothelium %4d | %s (%.0f min)" % (
            d.name, c, X.shape[0], int(meso.sum()), ", ".join("%s %d" % (g[:10], (grp == g).sum()) for g in sorted(set(grp)) if g != "unassigned"),
            (time.time() - t0) / 60))
    if not units:
        return None, None
    PB = pd.DataFrame(cols).fillna(0.0)
    PB = PB.groupby(level=0).sum()
    U = pd.DataFrame(units)
    PB.to_csv(PB_DIR / "pseudobulk.csv.gz")
    U.to_csv(PB_DIR / "units.csv", index=False)
    say("[T1] built: %d units, %d genes; %d QC cells, %d mesothelial cells" % (PB.shape[1], PB.shape[0], tot["cells"], tot["meso"]))
    return PB, U


def keyed(U):
    U = U.copy()
    U["lineage"] = U["group"].map(lambda g: "stroma" if str(g).startswith("stroma") else g)
    U["key"] = [("P:" + str(p)) if pd.notna(p) and str(p) not in ("", "nan", "None") else ("S:" + str(g)) for p, g in zip(U["patient"], U["gsm"])]
    return U


# ----------------------------------------------------------------------------
# R1 / R2 locus-gene localization (v060 method)
# ----------------------------------------------------------------------------
def context_tables(PB, U, ctx):
    T = U[(U["context"] == ctx) & (U["group"] == "__total__")]
    S = U[(U["context"] == ctx) & (U["group"] != "__total__") & (U["n_cells"] >= MIN_CELLS_UNIT)]
    genes = PB.index.values
    tot_cpm = {}
    for _, r in T.iterrows():
        x = PB[r["unit"]].values.astype(float)
        tot_cpm[r["gsm"]] = x / max(x.sum(), 1.0) * 1e6
    lfc, cpm = {}, {}
    for lin in LINEAGES:
        s = S[S["lineage"] == lin]
        pl, pc = {}, {}
        for _, r in s.iterrows():
            tu = T[T["gsm"] == r["gsm"]]
            if tu.empty:
                continue
            x = PB[r["unit"]].values.astype(float)
            rest = np.maximum(PB[tu["unit"].iloc[0]].values.astype(float) - x, 0.0)
            if rest.sum() <= 0:
                continue
            c_l = x / max(x.sum(), 1.0) * 1e6
            c_r = rest / rest.sum() * 1e6
            pl.setdefault(r["key"], []).append(np.log2(c_l + 1) - np.log2(c_r + 1))
            pc.setdefault(r["key"], []).append(c_l)
        if pl:
            lfc[lin] = pd.DataFrame({k: np.mean(v, axis=0) for k, v in pl.items()}, index=genes)
            cpm[lin] = pd.DataFrame({k: np.mean(v, axis=0) for k, v in pc.items()}, index=genes)
    mean_tot = pd.Series(np.mean(list(tot_cpm.values()), axis=0), index=genes) if tot_cpm else pd.Series(0.0, index=genes)
    return lfc, cpm, mean_tot, len(tot_cpm)


def matched_null(mean_tot, exclude, genes_obs, stat, rng):
    expr = mean_tot[mean_tot >= MIN_CPM]
    bins = pd.qcut(np.log10(expr + 1e-9).rank(method="first"), N_BINS, labels=False)
    vals = stat.reindex(expr.index).values
    idx = {g: i for i, g in enumerate(expr.index)}
    pool_i = {b: np.array([idx[g] for g in bins.index[bins == b] if g not in exclude]) for b in range(N_BINS)}
    want = [int(bins.get(g)) for g in genes_obs if g in bins.index]
    if not want:
        return None
    parts = []
    for b in sorted(set(want)):
        k = want.count(b)
        if len(pool_i[b]) < 10 * k:
            return None
        parts.append(pool_i[b][rng.integers(0, len(pool_i[b]), size=(N_NULL, k))])
    return np.nanmean(vals[np.concatenate(parts, axis=1)], axis=1)


def load_loci(PB):
    fs = sorted(LOCUS_DIR.glob("p9_v042_*/l1_loci_compartments.csv"))
    if not fs:
        return None
    L = pd.read_csv(fs[-1])
    L = L[L["trait"] == TRAIT_LOCI].copy()
    L["gene"] = L["expressed_gene_used"].astype(str).str.strip()
    L = L[L["compartment"].isin(COMP_ORDER) & (L["gene"] != "") & (L["gene"] != "nan")].drop_duplicates("gene")
    return L[L["gene"].isin(PB.index)]


def part_r1(PB, U, L):
    say("")
    say("-" * 72)
    say("[SECTION] R1 locus genes by lineage inside lesions (GSE179640; own lineage vs rest of sample; donor means)")
    rng = np.random.default_rng(SEED)
    exclude = set(L["gene"])
    rows = []
    for ctx in CONTEXTS:
        lfc, cpm, mean_tot, n_s = context_tables(PB, U, ctx)
        if not n_s:
            continue
        say("[R1] %-27s samples %d; donor keys per lineage: %s" % (CTX_NAME[ctx], n_s, ", ".join(
            "%s %d" % (lin.split(" ")[0], lfc[lin].shape[1]) for lin in LINEAGES if lin in lfc)))
        top = pd.DataFrame({lin: cpm[lin].mean(axis=1) for lin in cpm}) if cpm else pd.DataFrame()
        for comp in COMP_ORDER:
            lin = COMP_TO_LINEAGE[comp]
            G = list(L.loc[L["compartment"] == comp, "gene"])
            if lin not in lfc or lfc[lin].shape[1] < MIN_KEYS:
                say("  [R1] %-27s %-16s not evaluable (%d donor keys)" % (CTX_NAME[ctx], comp, lfc[lin].shape[1] if lin in lfc else 0))
                continue
            M = lfc[lin]
            Gx = [g for g in G if mean_tot.get(g, 0) >= MIN_CPM]
            if len(Gx) < MIN_GENES_SET:
                continue
            stat = M.mean(axis=1)
            obs = float(stat.loc[Gx].mean())
            null = matched_null(mean_tot, exclude, Gx, stat, rng)
            if null is None:
                continue
            p = (1 + int((null >= obs).sum())) / (N_NULL + 1)
            n_top = sum(1 for g in Gx if (not top.empty and top.loc[g].idxmax() == lin))
            n = M.shape[1]
            gp = []
            for g in Gx:
                v = M.loc[g].values
                sd = float(np.std(v, ddof=1))
                gp.append(tsf1(float(np.mean(v)) / max(sd / math.sqrt(n), 1e-9), n - 1))
            rows.append({"context": ctx, "compartment": comp, "n_genes": len(Gx), "n_keys": n, "observed": obs,
                         "null_median": float(np.median(null)), "null_p95": float(np.percentile(null, 95)), "p_perm": p,
                         "n_top_lineage": n_top, "n_gene_p05": sum(1 for x in gp if x < 0.05)})
            say("  [R1] %-27s %-16s genes %d, donors %2d: mean log2 ratio %+.2f (null median %+.2f, 95th %+.2f) perm p %s | top lineage %d/%d | gene p<0.05 %d/%d" % (
                CTX_NAME[ctx], comp, len(Gx), n, obs, np.median(null), np.percentile(null, 95), fmt_p(p), n_top, len(Gx),
                sum(1 for x in gp if x < 0.05), len(Gx)))
            if comp == "stroma":
                say("       stromal-locus genes: " + ", ".join("%s %+.1f%s" % (g, stat.loc[g], "" if (top.empty or top.loc[g].idxmax() == lin)
                                                                               else " (top %s)" % top.loc[g].idxmax().split(" ")[0]) for g in Gx))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "r1_set_tests.csv", index=False)
    st = R[(R["compartment"] == "stroma")] if len(R) else R
    ev = {c: st[st["context"] == c] for c in LESIONS}
    have = [c for c in LESIONS if len(ev[c])]
    ok = [bool(ev[c]["p_perm"].iloc[0] < 0.05) for c in have]
    verdict = ("NOT EVALUABLE" if "peritoneal_lesion" not in have else "SUPPORTED" if all(ok) else "NOT SUPPORTED")
    say("[H_LESION_REP_VERDICT] stromal-locus genes stroma-enriched beyond matched random genes in %s -> %s" % (
        ", ".join("%s (%s)" % (CTX_NAME[c], "yes" if o else "no") for c, o in zip(have, ok)) or "no evaluable lesion context", verdict))
    SUM["results"]["R1"] = {"verdict": verdict, "sets": R.to_dict(orient="records")}


def part_r2(PB, U, L):
    say("")
    say("-" * 72)
    say("[SECTION] R2 stromal-locus genes: lesion stroma vs stroma of lesion-adjacent peritoneum (host fibroblasts)")
    rng = np.random.default_rng(SEED + 1)
    G = list(L.loc[L["compartment"] == "stroma", "gene"])
    tabs = {}
    for ctx in CONTEXTS:
        _, cpm, _, _ = context_tables(PB, U, ctx)
        if "stroma" in cpm:
            tabs[ctx] = np.log2(cpm["stroma"] + 1)
    ref = "peritoneum_adjacent"
    if ref not in tabs or tabs[ref].shape[1] < MIN_KEYS:
        say("[H_HOST_VERDICT] adjacent-peritoneum stroma not evaluable (%d donor keys) -> NOT EVALUABLE" % (tabs[ref].shape[1] if ref in tabs else 0))
        SUM["results"]["R2"] = {"verdict": "NOT EVALUABLE"}
        return
    B = tabs[ref]
    rows, verdict = [], "NOT EVALUABLE"
    for ctx in ["eutopic", "control_endometrium"] + LESIONS:
        if ctx not in tabs or tabs[ctx].shape[1] < MIN_KEYS:
            continue
        A = tabs[ctx]
        shared = sorted(set(A.columns) & set(B.columns))
        paired = len(shared) >= MIN_KEYS
        diff = (A[shared] - B[shared]).mean(axis=1) if paired else A.mean(axis=1) - B.mean(axis=1)
        mean_all = ((2 ** A - 1).mean(axis=1) + (2 ** B - 1).mean(axis=1)) / 2
        Gx = [g for g in G if mean_all.get(g, 0) >= MIN_CPM]
        obs = float(diff.loc[Gx].mean())
        null = matched_null(mean_all, set(L["gene"]), Gx, diff, rng)
        p = (1 + int((null >= obs).sum())) / (N_NULL + 1) if null is not None else np.nan
        for g in Gx:
            rows.append({"context": ctx, "gene": g, "diff_log2": float(diff.loc[g]), "paired": paired, "n_ctx": A.shape[1], "n_ref": B.shape[1],
                         "shared": len(shared)})
        say("  [R2] %-27s (%d donors) vs adjacent peritoneum (%d), %s: mean diff %+.2f log2 (matched null median %+.2f) p %s | %s" % (
            CTX_NAME[ctx], A.shape[1], B.shape[1], "paired on %d donors" % len(shared) if paired else "unpaired", obs,
            float(np.median(null)) if null is not None else float("nan"), fmt_p(p), ", ".join("%s %+.1f" % (g, diff.loc[g]) for g in Gx)))
        SUM["results"].setdefault("R2", {})[ctx] = {"mean_diff": obs, "p": p, "paired": paired}
        if ctx == "peritoneal_lesion":
            verdict = "SUPPORTED" if (np.isfinite(p) and p < 0.05) else "NOT SUPPORTED"
    pd.DataFrame(rows).to_csv(OUT / "r2_host_stroma.csv", index=False)
    say("[H_HOST_VERDICT] peritoneal-lesion stroma vs adjacent-peritoneum stroma -> %s" % verdict)
    SUM["results"].setdefault("R2", {})["verdict"] = verdict


# ----------------------------------------------------------------------------
# R3 genetic ranking of marker-based groups (v052 method)
# ----------------------------------------------------------------------------
def de(PB, U):
    TU = U[U["group"] == "__total__"].set_index("gsm")["unit"]
    tot = PB[TU.values].copy()
    tot.columns = TU.index
    U = U[U["group"] != "__total__"]
    rows, cpm_t = [], {}
    elig = U.groupby("group")["gsm"].nunique()
    groups = [g for g in elig.index if elig[g] >= MIN_KEYS]
    say("[R3] samples per group: %s" % ", ".join("%s %d" % (g, elig[g]) for g in elig.index))
    say("[R3] eligible groups (>= %d samples with >= %d cells): %s" % (MIN_KEYS, MIN_CELLS_UNIT, ", ".join(groups)))
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
        p = tsf2(tt, n - 1)
        rows.append(pd.DataFrame({"group": g, "gene": genes[ex], "lfc": m, "score": tt, "p": p, "padj": bh(p), "n_samples": n}))
    if not rows:
        return None, None
    DE = pd.concat(rows, ignore_index=True)
    CPM = pd.DataFrame(cpm_t, index=genes)
    CPM = CPM[(CPM >= MIN_CPM).any(axis=1)]
    DE.to_csv(OUT / "r3_de_all.csv.gz", index=False)
    CPM.to_csv(OUT / "r3_mean_cpm.csv.gz")
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
    Sx = C.div(C.sum(axis=1), axis=0)
    spec, sd = {}, {}
    for i, c in enumerate(Sx.columns):
        k = "S%02d" % (i + 1)
        spec[k], sd[k] = c, Sx[c].values
    sd["Average"] = np.log2(C.mean(axis=1).values + 1.0)
    Ts = pd.DataFrame(sd, index=[sym2ent[g] for g in C.index])
    Ts = Ts[~Ts.index.duplicated()]
    return T, Ts, meta, spec


def ranked(res, meta, spec, fam):
    it = ([(spec[k], v) for k, v in res.items() if k in spec] if fam == "SPEC"
          else [(meta[k][1], v) for k, v in res.items() if k in meta and meta[k][0] == fam])
    it.sort(key=lambda x: x[1]["p"])
    q = bh([v["p"] for _, v in it]) if it else []
    return [(l, v, qq) for (l, v), qq in zip(it, q)]


def part_r3(PB, U, sym2ent):
    say("")
    say("-" * 72)
    say("[SECTION] R3 MAGMA ranking of marker-based groups (GSE179640; p / BH q)")
    DE, CPM = de(PB, U)
    if DE is None:
        say("[H_MESO_REP_VERDICT] no eligible groups -> NOT EVALUABLE")
        return
    for g in ["MSLN", "UPK3B", "ITLN1", "CALB2", "EPCAM", "PDGFRA", "DCN", "ESR1", "PGR", "MME", "WT1", "HOXA10", "IGF1", "SFRP1", "PTPRC"]:
        if g in CPM.index:
            say("    %-7s %s" % (g, " | ".join("%s %.0f" % (c[:26], CPM.loc[g, c]) for c in CPM.columns)))
    T, Ts, meta, spec = families(DE, CPM, sym2ent)
    cf, cfs = WORK / "sets.txt", WORK / "spec.txt"
    T.to_csv(cf, sep="\t", index_label="GENE")
    Ts.to_csv(cfs, sep="\t", float_format="%.6g", index_label="GENE")
    res = {}
    for t in TRAITS:
        raw = MAGMA_OUT / (t + "_35_10.genes.raw")
        if not raw.exists():
            continue
        r0 = magma_covar(raw, cf, "r3_sets_%s" % t, [])
        r1 = magma_covar(raw, cfs, "r3_spec_%s" % t, ["Average"])
        if r0 and r1:
            res[t] = {fam: ranked(r1 if fam == "SPEC" else r0, meta, spec, fam) for fam in FAMS}
    rows = []
    for t, F in res.items():
        say("  %s" % t)
        for fam, it in F.items():
            say("    %-4s %s" % (fam, "; ".join("%s %.0e/%.2f" % (l, v["p"], q) for l, v, q in it)))
            for i, (l, v, q) in enumerate(it):
                rows.append({"trait": t, "family": fam, "rank": i + 1, "label": l, "p": v["p"], "q": q, "beta_std": v["beta_std"]})
    A = pd.DataFrame(rows)
    A.to_csv(OUT / "r3_rankings.csv", index=False)
    k = A[(A["trait"] == "KOLLER_EUR") & (A["family"] == "UPN")].set_index("label") if len(A) else pd.DataFrame()
    es = [g for g in ENDO_STROMA if g in k.index]
    rm = k.loc["mesothelium", "rank"] if "mesothelium" in k.index else None
    rp = k.loc["stroma (peritoneum_adjacent)", "rank"] if "stroma (peritoneum_adjacent)" in k.index else None
    if not es or rm is None:
        v = "NOT EVALUABLE"
        say("[H_MESO_REP_VERDICT] KOLLER_EUR UPN: %s -> %s" % ("mesothelium not an eligible group" if rm is None else "no endometrial-type stroma group", v))
    else:
        best = min(es, key=lambda g: k.loc[g, "rank"])
        ok = k.loc[best, "rank"] < rm and (rp is None or k.loc[best, "rank"] < rp)
        v = "SUPPORTED" if ok else "NOT SUPPORTED"
        say("[H_MESO_REP_VERDICT] KOLLER_EUR UPN: best endometrial-type stroma %s rank %d (p %.1e) | mesothelium rank %d (p %.1e) | adjacent-peritoneum stroma rank %s -> %s" % (
            best, k.loc[best, "rank"], k.loc[best, "p"], rm, k.loc["mesothelium", "p"], "-" if rp is None else int(rp), v))
    SUM["results"]["R3"] = {"verdict": v}


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.65 | independent lesion cohort GSE179640 (Tan et al. 2022) + Koller ancestry scope | " + TS)
    say("=" * 72)
    say("[DISK] free %.1f GB" % (shutil.disk_usage(str(REPO)).free / 1e9))
    part_z()
    S = part_t0()
    if S is None:
        say("[FATAL] no GSE179640 metadata")
        _fh.close()
        return
    PB, U = build_pseudobulk(S)
    if PB is None:
        say("[FATAL] no usable GSE179640 matrices")
        _fh.close()
        return
    U = keyed(U)
    say("[T1] donor keys per context: %s" % ", ".join("%s %d" % (c, U[(U["context"] == c) & (U["group"] == "__total__")]["key"].nunique())
                                                     for c in CONTEXTS))
    L = load_loci(PB)
    if L is None or L.empty:
        say("[FATAL] locus table missing")
    else:
        part_r1(PB, U, L)
        part_r2(PB, U, L)
    part_r3(PB, U, gene_loc())
    (OUT / "summary.json").write_text(json.dumps(SUM, indent=1, default=str), encoding="utf-8")
    say("")
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V065_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
