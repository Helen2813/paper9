#!/usr/bin/env python3
# Paper 9 | v0.0.72 | GSE183771 (Vrljicak et al. 2023, endometrial scATAC) audit and cell-type chromatin annotations
#   Rules: CHANGELOG v0.2.61 (fixed before this run). No GWAS summary statistics are read by this script.
#   A  GEO metadata: samples, characteristics, data processing, supplementary files (series and sample level)
#   B  download of peak-by-cell matrices (no fragment / BAM / bigWig files), reading, common peak universe
#   C  lift-over of the peak universe GRCh38 -> GRCh37 (UCSC chain, >= 95% of bases remapped)          gate G2
#   D  lineage labels: LSI + k-means clusters labelled by marker-gene accessibility (fixed panel);
#      author labels are used instead if a cell-annotation file covers >= 90% of cells;
#      check against the published composition                                                        gate G1
#   E  biopsy-level pseudobulk per lineage (technical replicate libraries merged)                      gate G3
#   F  lineage-specific peak annotations (top 5 / 10 / 20% by biopsy-paired specificity) + all-peaks control
#   G  independent check with endometrial stromal culture ATAC (GSE246585 peaks, all time points pooled)  gate G4
#   H  thin annotations on the 1000 Genomes EUR panel and partitioned LD scores (only if G1-G4 pass)
# Runs in ~/Desktop/paper9/.venv, low heat (2 threads, nice 10). Download size is printed before downloading
# (stops above 8 GB). Reading ~20-40 min, LD scores ~20-40 min. Re-running reuses downloads and caches.

import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
try:
    os.nice(10)
except Exception:
    pass

import datetime as dt
import gzip
import json
import re
import shutil
import subprocess
import tarfile
import time
import urllib.request
from bisect import bisect_right
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

# ----------------------------------------------------------------------------
# Hardcoded locations
# ----------------------------------------------------------------------------
HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RAW = REPO / "data" / "raw"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
LDSC_PY = TOOLS / "ldsc_venv" / "bin" / "python"
LDSC = TOOLS / "ldsc_py3" / "ldsc.py"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
REF = RAW / "ldsc_ref" / "unpacked"
PLINK_DIR = REF / "1000G_EUR_Phase3_plink"
FRQ_DIR = REF / "1000G_Phase3_frq"
HM3_LIST = RAW / "ldsc_ref" / "hm3_no_MHC.list.txt"
GSE = "GSE183771"
GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183771/"
GEO_DIR = RAW / "geo" / GSE
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
CHAIN = RAW / "liftover" / "hg38ToHg19.over.chain.gz"
ENSC_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246585/suppl/"
            "GSE246585_Decidualized_timecourse_ATACseq_counts_matrix.txt.gz")
ENSC_FILE = RAW / "geo" / "GSE246585" / "GSE246585_Decidualized_timecourse_ATACseq_counts_matrix.txt.gz"
CACHE = PROC / "gse183771_v072"
LD_DIR = PROC / "ldscores" / "gse183771_v072"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "chromatin" / ("p9_v072_" + TS)
LOG = REPO / "logs" / ("p9_v072_scatac_audit_" + TS + ".log")
NICE = ["nice", "-n", "10"]
BLAS_THREADS = "2"
L2_ENV = dict(os.environ, OPENBLAS_NUM_THREADS=BLAS_THREADS, OMP_NUM_THREADS=BLAS_THREADS,
              VECLIB_MAXIMUM_THREADS=BLAS_THREADS, MKL_NUM_THREADS=BLAS_THREADS, NUMEXPR_NUM_THREADS=BLAS_THREADS)
LD_WIND_CM = 1.0

# ----------------------------------------------------------------------------
# Fixed rules (CHANGELOG v0.2.61)
# ----------------------------------------------------------------------------
PUBLISHED = {"stroma": 16851, "epithelium": 16075, "ciliated": 430, "immune": 2386, "endothelium": 435}  # 36,177 nuclei
MARKERS = {
    "stroma": ["COL3A1", "COL1A1", "COL1A2", "COL6A3", "PDGFRA", "DCN", "HAND2", "MMP11"],
    "epithelium": ["ELF3", "EPCAM", "KRT8", "KRT18", "CDH1", "CLDN4", "PAX8", "SOX17"],
    "ciliated": ["DNAH5", "FOXJ1", "CDHR3", "RSPH1", "CAPS", "DNAI1"],
    "endothelium": ["PECAM1", "VWF", "CDH5", "CLDN5", "FLT1", "EMCN", "TIE1", "KDR"],
    "immune": ["PTPRC", "KLRD1", "NKG7", "CD2", "CD247", "CD14", "CD68", "CSF1R"],
}
LABELS = ["stroma", "epithelium", "ciliated", "endothelium", "immune"]
LINEAGES = ["stroma", "epithelium", "endothelium", "immune"]       # ciliated pooled into epithelium from E on
GA_UPSTREAM = 2000
MIN_PEAK_COUNTS_CELL = 500
MAX_CELL_FOLD_MEDIAN = 5.0
N_LSI_PEAKS = 25000
N_LSI_SUB = 10000
N_COMP = 30
K_CLUST = 30
K_CLUST_FALLBACK = 60
LABEL_MIN_SCORE = 0.5
LABEL_MIN_MARGIN = 0.25
AUTHOR_MIN_COVER = 0.90
G1_TOL_PP = 10.0
G1_MIN_RARE = 100
G1_MAX_UNASSIGNED = 0.10
MIN_CELLS_UNIT = 20
MIN_UNITS = 6
MIN_PEAK_TOTAL = 20
TOP_Q = [5, 10, 20]
LIFT_MIN_MATCH = 0.95
LIFT_MAX_FOLD = 2.0
G2_MIN_LIFT = 0.95
G4_ALPHA = 0.05
MAX_DOWNLOAD_GB = 8.0
MAX_UNKNOWN_FILE_MB = 300
SEED = 20261111
SKIP_RE = re.compile(r"fragment|\.bam$|\.bai$|possorted|\.bw$|bigwig|\.tbi$|\.cram|\.fastq|\.fq\.gz", re.I)
KEEP_RE = re.compile(r"matrix\.mtx|barcode|peaks?[._]|features|\.h5$|peak_bc|singlecell|meta|annot|cluster|cell_?type", re.I)
ANNOT_RE = re.compile(r"meta|annot|cluster|cell_?type", re.I)
STD_CHR = re.compile(r"^chr([0-9]{1,2}|X|Y)$")

OUT.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)
for d in (GEO_DIR, CHAIN.parent, ENSC_FILE.parent, CACHE, LD_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V072_SCATAC_AUDIT", "timestamp_utc": TS, "gates": {}, "results": {}}


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def sh(cmd, timeout=None, env=None):
    try:
        p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout, env=env)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"
    except FileNotFoundError as e:
        return 127, "", str(e)


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "-"
    return "%.3f" % p if p >= 0.001 else "%.0e" % p


def https(u):
    return re.sub(r"^ftp://", "https://", u.strip())


def fetch(url, dst, timeout=6 * 3600):
    if dst.exists() and dst.stat().st_size > 0:
        return True
    dst.parent.mkdir(parents=True, exist_ok=True)
    part = dst.with_name(dst.name + ".part")
    rc, o, e = sh(["curl", "-L", "-sS", "--fail", "--retry", "3", "--connect-timeout", "30", "-C", "-", "-A", "paper9",
                   "-o", part, url], timeout=timeout)
    if rc == 0 and part.exists() and part.stat().st_size > 0:
        part.rename(dst)
        return True
    say("  [DL] FAILED %s rc=%s %s" % (url.split("/")[-1], rc, (e or o).strip()[:160]))
    return False


def listing(url):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "paper9"}), timeout=60) as r:
            html = r.read().decode("utf-8", "ignore")
    except Exception as e:
        say("  [LIST] %s failed: %s" % (url, type(e).__name__))
        return None
    rows = []
    for m in re.finditer(r'<a href="([^"?/][^"]*)">[^<]*</a>\s*([0-9]{4}-[0-9]{2}-[0-9]{2} [0-9:]+)?\s*([0-9.]+[KMG]?)?', html):
        if not m.group(1).startswith("http"):
            rows.append((m.group(1), m.group(3) or ""))
    return rows


def size_of(s):
    m = re.match(r"^([0-9.]+)([KMG]?)$", s or "")
    if not m:
        return None
    return float(m.group(1)) * {"": 1, "K": 1e3, "M": 1e6, "G": 1e9}[m.group(2)]


def head_bytes(url):
    rc, o, e = sh(["curl", "-sSIL", "--max-time", "60", "-A", "paper9", url], timeout=90)
    if rc != 0:
        return None
    size = re.findall(r"(?i)content-length:\s*(\d+)", o)
    return float(size[-1]) if size else None


def opener(f):
    return gzip.open(f, "rt", errors="ignore") if str(f).endswith(".gz") else open(f, "rt", errors="ignore")


# ----------------------------------------------------------------------------
# A  GEO metadata
# ----------------------------------------------------------------------------
def series_samples():
    names = [n for n, _ in (listing(GEO_BASE + "matrix/") or []) if n.endswith("_series_matrix.txt.gz")]
    if not names:
        names = sorted(p.name for p in GEO_DIR.glob("*_series_matrix.txt.gz")) or [GSE + "_series_matrix.txt.gz"]
    rows, proc = [], []
    for m in names:
        f = GEO_DIR / m
        if not fetch(GEO_BASE + "matrix/" + m, f, timeout=1800):
            continue
        txt = gzip.open(f, "rt", errors="ignore").read().split("!series_matrix_table_begin")[0]
        R = {}
        for line in txt.splitlines():
            if line.startswith("!Sample_") or line.startswith("!Series_"):
                R.setdefault(line.split("\t")[0], []).append([v.strip().strip('"') for v in line.split("\t")[1:]])
        gsm = R.get("!Sample_geo_accession", [[]])[0]
        for i, g in enumerate(gsm):
            def get(k, sep=" | "):
                return sep.join(c[i] for c in R.get(k, []) if i < len(c) and c[i] and c[i] != "NONE")
            sup = [https(c[i]) for k in R if k.startswith("!Sample_supplementary_file") for c in R[k]
                   if i < len(c) and c[i] and c[i] != "NONE"]
            rows.append({"gsm": g, "title": get("!Sample_title"), "source": get("!Sample_source_name_ch1"),
                         "chars": get("!Sample_characteristics_ch1"), "description": get("!Sample_description"),
                         "strategy": get("!Sample_library_strategy"), "supp": sup})
            if not proc:
                proc = [c[i] for c in R.get("!Sample_data_processing", []) if i < len(c) and c[i]]
        SUM["results"]["series_title"] = " ".join(R.get("!Series_title", [[""]])[0])
    S = pd.DataFrame(rows).drop_duplicates("gsm").reset_index(drop=True) if rows else pd.DataFrame()
    return S, proc


def unit_of(r):
    # biopsy key: a donor / individual / biopsy characteristic if present, else the title without replicate tokens
    for part in str(r["chars"]).split("|"):
        m = re.match(r"\s*(donor|individual|patient|subject|participant|biopsy)[ _]?(id)?\s*[:=]\s*(.+)", part.strip(), re.I)
        if m:
            return m.group(3).strip()
    t = str(r["title"]).strip()
    t = re.sub(r"[\s_.-]*(rep|replicate|run|tech|lane|library|lib)[\s_.-]*\d+$", "", t, flags=re.I)
    t = re.sub(r"[\s_.-]+[12]$", "", t)
    return t


def part_a():
    say("")
    say("-" * 72)
    say("[SECTION] A GEO metadata of %s" % GSE)
    S, proc = series_samples()
    if S.empty:
        say("[A] series matrix not available -> stop")
        return None
    say("[A] %s" % SUM["results"].get("series_title", "")[:150])
    S["unit"] = [unit_of(r) for _, r in S.iterrows()]
    for _, r in S.iterrows():
        say("  [A] %-11s %-22s unit %-10s %-10s | %s" % (r["gsm"], r["title"][:22], r["unit"][:10], r["strategy"][:10],
                                                      r["chars"][:120]))
        for u in r["supp"]:
            say("        file: %s" % u.split("/")[-1][:100])
    for line in proc[:15]:
        say("  [A] processing: %s" % line[:170])
    sup = listing(GEO_BASE + "suppl/")
    say("[A] series supplementary: %s" % ("listing failed" if sup is None else
                                         (", ".join("%s (%s)" % (a, b) for a, b in sup) or "none")))
    n_units = S["unit"].nunique()
    say("[A] %d GSM, %d biopsy units (replicate libraries merged by the unit rule)" % (len(S), n_units))
    for u, g in S.groupby("unit"):
        if len(g) > 1:
            say("  [A] unit %s <- %s" % (u, ", ".join(g["gsm"])))
    S.drop(columns=["supp"]).to_csv(OUT / "a_samples.csv", index=False)
    SUM["results"]["A"] = {"n_gsm": int(len(S)), "n_units": int(n_units)}
    return S, (sup or [])


# ----------------------------------------------------------------------------
# B  downloads, reading, peak universe
# ----------------------------------------------------------------------------
def plan_downloads(S, sup):
    plan, skipped = [], []
    for _, r in S.iterrows():
        for u in r["supp"]:
            name = u.split("/")[-1]
            dst = GEO_DIR / r["gsm"] / name
            if SKIP_RE.search(name):
                skipped.append(name)
                continue
            size = dst.stat().st_size if dst.exists() else head_bytes(u)
            if not KEEP_RE.search(name) and (size is None or size > MAX_UNKNOWN_FILE_MB * 1e6):
                skipped.append(name)
                continue
            plan.append((u, dst, size))
    for name, s in sup:
        if name.endswith("_RAW.tar") or SKIP_RE.search(name) or name == "filelist.txt":
            continue
        size = size_of(s)
        if size is not None and size > 2e9:
            skipped.append(name)
            continue
        plan.append((GEO_BASE + "suppl/" + name, GEO_DIR / "series" / name, size))
    if not any(r["supp"] for _, r in S.iterrows()):
        raw = [(n, s) for n, s in sup if n.endswith("_RAW.tar")]
        if raw:
            plan.append((GEO_BASE + "suppl/" + raw[0][0], GEO_DIR / raw[0][0], size_of(raw[0][1])))
    return plan, skipped


def extract_raw(tar):
    done = GEO_DIR / ".raw_extracted_v072"
    if done.exists():
        return
    with tarfile.open(tar) as tf:
        for m in tf.getmembers():
            name = Path(m.name).name
            if not m.isfile() or SKIP_RE.search(name) or not re.match(r"^GSM\d+_", name):
                continue
            g = name.split("_")[0]
            (GEO_DIR / g).mkdir(exist_ok=True)
            with tf.extractfile(m) as src, open(GEO_DIR / g / name, "wb") as dst:
                shutil.copyfileobj(src, dst)
    done.write_text("ok\n")


def parse_coord(s):
    m = re.match(r"^(chr)?([0-9]{1,2}|X|Y|M|MT)[:_-](\d+)[-_:](\d+)$", str(s).strip())
    if not m:
        return None
    return "chr" + m.group(2), int(m.group(3)), int(m.group(4))


def read_peaks_file(f):
    with opener(f) as h:
        first = [next(h, "") for _ in range(3)]
    cols = first[0].rstrip("\n").split("\t")
    if len(cols) >= 3 and re.match(r"^\d+$", cols[1].strip()) and re.match(r"^\d+$", cols[2].strip()):
        D = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"], comment="#",
                        dtype={"chrom": str})
    else:
        T = pd.read_csv(f, sep="\t", header=None, dtype=str)
        col = next((j for j in range(T.shape[1]) if parse_coord(T.iloc[0, j])), None)
        if col is None:
            raise ValueError("no peak coordinates in %s" % f.name)
        P = [parse_coord(x) for x in T.iloc[:, col]]
        D = pd.DataFrame([p if p else ("NA", -1, -1) for p in P], columns=["chrom", "start", "end"])
    D["chrom"] = [c if str(c).startswith("chr") else "chr" + str(c) for c in D["chrom"]]
    return D.reset_index(drop=True)


def read_mtx(f):
    # Matrix Market coordinate file, read in chunks (single thread, int32)
    with opener(f) as h:
        header = None
        skip = 0
        for line in h:
            skip += 1
            if not line.startswith("%"):
                header = [int(x) for x in line.split()[:3]]
                break
    nr, nc, nnz = header
    rows, cols, vals = [], [], []
    for ch in pd.read_csv(f, sep=r"\s+", header=None, skiprows=skip, chunksize=5_000_000, dtype=np.float64,
                          engine="c"):
        rows.append(ch.iloc[:, 0].values.astype(np.int32) - 1)
        cols.append(ch.iloc[:, 1].values.astype(np.int32) - 1)
        vals.append(ch.iloc[:, 2].values.astype(np.float32) if ch.shape[1] > 2 else np.ones(len(ch), np.float32))
    M = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nr, nc)).tocsr()
    return M


def read_h5(f):
    import h5py
    with h5py.File(f, "r") as h:
        g = h["matrix"]
        shape = tuple(g["shape"][()])
        M = sp.csc_matrix((g["data"][()], g["indices"][()], g["indptr"][()]), shape=shape)
        bc = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g["barcodes"][()]])
        fid = g["features"]["id"][()] if "features" in g else g["genes"][()]
        ids = [x.decode() if isinstance(x, bytes) else str(x) for x in fid]
    P = [parse_coord(x) for x in ids]
    peaks = pd.DataFrame([p if p else ("NA", -1, -1) for p in P], columns=["chrom", "start", "end"])
    return M.T.tocsr().astype(np.float32), bc, peaks


def find_lib(d):
    files = [p for p in sorted(d.rglob("*")) if p.is_file() and not p.name.endswith(".part")]
    mtx = [p for p in files if re.search(r"\.mtx(\.gz)?$", p.name)]
    mtx = [p for p in mtx if "raw" not in p.name.lower()] or mtx
    h5 = [p for p in files if p.name.endswith(".h5")]
    h5 = [p for p in h5 if "raw" not in p.name.lower()] or h5
    bc = [p for p in files if "barcode" in p.name.lower() and not p.name.endswith(".h5")]
    pk = [p for p in files if re.search(r"peak|feature|region", p.name.lower()) and re.search(r"\.(bed|tsv|txt)(\.gz)?$", p.name)]
    if mtx and bc and pk:
        return {"kind": "mtx", "mtx": mtx[0], "bc": bc[0], "pk": pk[0]}
    if h5:
        return {"kind": "h5", "h5": h5[0]}
    return None


def read_peaks_only(info):
    if info["kind"] == "h5":
        import h5py
        with h5py.File(info["h5"], "r") as h:
            g = h["matrix"]
            fid = g["features"]["id"][()] if "features" in g else g["genes"][()]
            ids = [x.decode() if isinstance(x, bytes) else str(x) for x in fid]
        P = [parse_coord(x) for x in ids]
        return pd.DataFrame([p if p else ("NA", -1, -1) for p in P], columns=["chrom", "start", "end"])
    return read_peaks_file(info["pk"])


def read_lib(info):
    if info["kind"] == "h5":
        return read_h5(info["h5"])
    M = read_mtx(info["mtx"])
    bc = pd.read_csv(info["bc"], sep="\t", header=None, dtype=str).iloc[:, 0].values
    peaks = read_peaks_file(info["pk"])
    if M.shape == (len(peaks), len(bc)):
        M = M.T.tocsr()
    elif M.shape != (len(bc), len(peaks)):
        raise ValueError("matrix %s vs %d barcodes x %d peaks" % (M.shape, len(bc), len(peaks)))
    return M.astype(np.float32), bc, peaks


def merge_intervals(D):
    out = []
    for c, g in D.groupby("chrom", sort=False):
        a = g[["start", "end"]].to_numpy(np.int64)
        a = a[np.argsort(a[:, 0], kind="stable")]
        cur_s, cur_e = a[0]
        for s, e in a[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                out.append((c, cur_s, cur_e))
                cur_s, cur_e = s, e
        out.append((c, cur_s, cur_e))
    U = pd.DataFrame(out, columns=["chrom", "start", "end"])
    U["k"] = [int(x[3:]) if x[3:].isdigit() else {"X": 23, "Y": 24}.get(x[3:], 25) for x in U["chrom"]]
    return U.sort_values(["k", "start"]).drop(columns="k").reset_index(drop=True)


def map_to_universe(peaks, U):
    idx = np.full(len(peaks), -1, dtype=np.int64)
    for c, g in U.groupby("chrom", sort=False):
        sel = np.where(peaks["chrom"].values == c)[0]
        if len(sel) == 0:
            continue
        st = g["start"].to_numpy()
        en = g["end"].to_numpy()
        ps = peaks["start"].to_numpy()[sel]
        i = np.searchsorted(st, ps, side="right") - 1
        ok = (i >= 0) & (ps < en[np.clip(i, 0, None)])
        idx[sel[ok]] = g.index.to_numpy()[i[ok]]
    return idx


def library_units(S, libs):
    # library id -> (gsm, unit)
    m = {}
    for lib in libs:
        g = lib.split("#")[0]
        row = S[S["gsm"] == g]
        unit = row["unit"].iloc[0] if len(row) else lib
        m[lib] = (g, unit if "#" not in lib else lib)
    return m


def part_b(S, sup):
    say("")
    say("-" * 72)
    say("[SECTION] B downloads, peak matrices and the common peak universe")
    f_u = CACHE / "universe_hg38.csv.gz"
    f_libs = CACHE / "libraries.csv"
    if f_u.exists() and f_libs.exists():
        U = pd.read_csv(f_u)
        L = pd.read_csv(f_libs)
        say("[B] cached: %d libraries, %d cells, %d universe peaks (%s)" % (len(L), int(L["cells_kept"].sum()), len(U), CACHE.name))
        return U, L
    plan, skipped = plan_downloads(S, sup)
    known = [s for _, _, s in plan if s is not None]
    total = sum(known) / 1e9
    say("[B] download plan: %d files, %.2f GB known size (%d of unknown size); skipped %d: %s" % (
        len(plan), total, len(plan) - len(known), len(skipped), ", ".join(sorted(set(skipped)))[:300]))
    if total > MAX_DOWNLOAD_GB:
        say("[B] plan exceeds %.0f GB -> stop (no download)" % MAX_DOWNLOAD_GB)
        return None
    t0 = time.time()
    for u, dst, s in plan:
        if fetch(u, dst) and dst.name.endswith("_RAW.tar"):
            extract_raw(dst)
    say("[B] downloads done in %.0f min" % ((time.time() - t0) / 60))
    libs = {}
    for _, r in S.iterrows():
        d = GEO_DIR / r["gsm"]
        if d.exists():
            info = find_lib(d)
            if info:
                libs[r["gsm"]] = info
            else:
                say("  [B] %s: no peak matrix found (%s)" % (r["gsm"], ", ".join(p.name for p in d.iterdir())[:200]))
    mode = "per-library"
    if not libs and (GEO_DIR / "series").exists():
        info = find_lib(GEO_DIR / "series")
        if info:
            libs["SERIES"] = info
            mode = "combined"
    say("[B] mode: %s (%d matrices)" % (mode, len(libs)))
    if not libs:
        say("[B] no readable peak-by-cell matrix -> stop")
        return None
    # pass 1: peak tables only, then the universe
    peak_tabs, keep_tabs = {}, {}
    for k, info in libs.items():
        peaks = read_peaks_only(info)
        keep_p = np.array([bool(STD_CHR.match(c)) for c in peaks["chrom"]]) & (peaks["start"].values >= 0)
        peak_tabs[k], keep_tabs[k] = peaks, keep_p
    kept = [p[m] for p, m in zip(peak_tabs.values(), keep_tabs.values())]
    first = kept[0]
    same = all(len(p) == len(first) and (p[["chrom", "start", "end"]].values == first[["chrom", "start", "end"]].values).all()
               for p in kept)
    U = merge_intervals(pd.concat(kept, ignore_index=True))
    if same:
        say("[B] all matrices share one peak set: %d peaks on standard chromosomes -> universe %d (sorted, overlaps merged)" % (
            len(first), len(U)))
    else:
        say("[B] per-library peak sets differ -> merged universe of %d peaks (from %d library peaks)" % (
            len(U), sum(len(p) for p in kept)))
    w = (U["end"] - U["start"]).to_numpy()
    say("[B] universe peak width: median %d bp, 99th percentile %d bp, max %d bp" % (
        int(np.median(w)), int(np.percentile(w, 99)), int(w.max())))
    if mode == "combined":
        say("[B] CAUTION combined matrix: library = barcode suffix; GSM mapping and replicate merging are not possible")
    U.to_csv(f_u, index=False)
    # pass 2: one matrix at a time -> universe columns, cell QC, cache
    rows = []
    df_cells = np.zeros(len(U), dtype=np.int64)
    tot_counts = np.zeros(len(U), dtype=np.float64)
    for k, info in libs.items():
        t1 = time.time()
        M, bc, peaks = read_lib(info)
        keep_p = keep_tabs[k]
        if len(peaks) != len(keep_p):
            raise ValueError("%s: peak table changed between passes" % k)
        idx = map_to_universe(peaks, U)
        ok = keep_p & (idx >= 0)
        Map = sp.csr_matrix((np.ones(ok.sum(), np.float32), (np.where(ok)[0], idx[ok])), shape=(len(peaks), len(U)))
        X = (M @ Map).tocsr()
        del M
        tot = np.asarray(X.sum(1)).ravel()
        say("  [B] %-12s %6d barcodes x %7d peaks, nnz %.1fM, %.1f%% of peaks in the universe (%.1f min)" % (
            k, X.shape[0], len(peaks), X.nnz / 1e6, 100 * float(ok.mean()), (time.time() - t1) / 60))
        if mode == "combined":
            suf = np.array([re.search(r"-(\d+)$", b).group(1) if re.search(r"-(\d+)$", b) else "1" for b in bc])
            groups = [(k + "#" + s, np.where(suf == s)[0]) for s in sorted(set(suf), key=lambda x: int(x))]
        else:
            groups = [(k, np.arange(len(bc)))]
        for lib, rows_i in groups:
            t = tot[rows_i]
            base = t[t >= MIN_PEAK_COUNTS_CELL]
            med = float(np.median(base)) if len(base) else 0.0
            keep = rows_i[(t >= MIN_PEAK_COUNTS_CELL) & (t <= MAX_CELL_FOLD_MEDIAN * max(med, 1))]
            Xk = X[keep]
            sp.save_npz(CACHE / ("lib_%s.npz" % lib.replace("#", "_")), Xk)
            pd.Series(bc[keep]).to_csv(CACHE / ("lib_%s_barcodes.txt.gz" % lib.replace("#", "_")), index=False, header=False)
            df_cells += np.asarray((Xk > 0).sum(0)).ravel()
            tot_counts += np.asarray(Xk.sum(0)).ravel()
            rows.append({"lib": lib, "file_key": lib.replace("#", "_"), "barcodes": len(rows_i), "cells_kept": len(keep),
                         "median_counts": med, "frac_in_universe": float(ok.mean())})
            say("  [B] %-14s barcodes %6d -> cells %6d (peak counts >= %d and <= %.0fx median %.0f)" % (
                lib, len(rows_i), len(keep), MIN_PEAK_COUNTS_CELL, MAX_CELL_FOLD_MEDIAN, med))
        del X
    L = pd.DataFrame(rows)
    U["n_cells"] = df_cells
    U["total"] = tot_counts
    U.to_csv(f_u, index=False)
    L["mode"] = mode
    L.to_csv(f_libs, index=False)
    say("[B] %d libraries, %d cells kept (published: 36,177 nuclei)" % (len(L), int(L["cells_kept"].sum())))
    SUM["results"]["B"] = {"mode": mode, "n_libraries": int(len(L)), "cells_kept": int(L["cells_kept"].sum()),
                           "universe_peaks": int(len(U))}
    return U, L


# ----------------------------------------------------------------------------
# C  lift-over GRCh38 -> GRCh37 (chain file; minMatch on remapped bases)
# ----------------------------------------------------------------------------
def load_chain(f):
    blocks = {}
    cur = None
    with gzip.open(f, "rt") as h:
        for line in h:
            t = line.split()
            if not t:
                continue
            if t[0] == "chain":
                tname, tstart = t[2], int(t[5])
                qname, qsize, qstrand, qstart = t[7], int(t[8]), t[9], int(t[10])
                cid = int(t[12]) if len(t) > 12 else len(blocks)
                cur = [tname, tstart, qname, qsize, qstrand, qstart, cid]
                continue
            size = int(t[0])
            tname, tpos, qname, qsize, qstrand, qpos, cid = cur
            blocks.setdefault(tname, []).append((tpos, tpos + size, qname, qpos, qstrand, qsize, cid))
            if len(t) == 3:
                cur[1] = tpos + size + int(t[1])
                cur[5] = qpos + size + int(t[2])
    out = {}
    for c, b in blocks.items():
        b.sort(key=lambda x: x[0])
        ts = np.array([x[0] for x in b], dtype=np.int64)
        te = np.array([x[1] for x in b], dtype=np.int64)
        out[c] = {"ts": ts, "te": te, "maxend": np.maximum.accumulate(te), "b": b}
    return out


def lift_one(ch, c, s, e):
    if c not in ch:
        return None, "chrom"
    B = ch[c]
    i = bisect_right(B["ts"], s) - 1
    hits = []
    j = i
    while j >= 0 and B["maxend"][j] > s:
        if B["te"][j] > s and B["ts"][j] < e:
            hits.append(j)
        j -= 1
    j = i + 1
    while j < len(B["ts"]) and B["ts"][j] < e:
        hits.append(j)
        j += 1
    if not hits:
        return None, "unmapped"
    by = {}
    for j in hits:
        ts, te, qn, qp, qs, qz, cid = B["b"][j]
        os_, oe = max(s, ts), min(e, te)
        if oe <= os_:
            continue
        a, b = qp + (os_ - ts), qp + (oe - ts)
        if qs == "-":
            a, b = qz - b, qz - a
        d = by.setdefault((cid, qn), [0, None, None])
        d[0] += oe - os_
        d[1] = a if d[1] is None else min(d[1], a)
        d[2] = b if d[2] is None else max(d[2], b)
    if not by:
        return None, "unmapped"
    (cid, qn), (nb, a, b) = max(by.items(), key=lambda kv: kv[1][0])
    if nb < LIFT_MIN_MATCH * (e - s):
        return None, "partial"
    if (b - a) > LIFT_MAX_FOLD * (e - s):
        return None, "split"
    return (qn, a, b), "ok"


def part_c(U):
    say("")
    say("-" * 72)
    say("[SECTION] C lift-over of the peak universe GRCh38 -> GRCh37")
    f = CACHE / "universe_lifted.csv.gz"
    if f.exists():
        V = pd.read_csv(f)
    else:
        if not fetch(CHAIN_URL, CHAIN, timeout=1800):
            say("[C] chain file not available -> stop")
            return None
        ch = load_chain(CHAIN)
        res = [lift_one(ch, c, int(s), int(e)) for c, s, e in zip(U["chrom"], U["start"], U["end"])]
        V = U.copy()
        V["chrom37"] = [r[0][0] if r[0] else "" for r in res]
        V["start37"] = [r[0][1] if r[0] else -1 for r in res]
        V["end37"] = [r[0][2] if r[0] else -1 for r in res]
        V["lift"] = [r[1] for r in res]
        V.to_csv(f, index=False)
    auto = V["chrom"].str.match(r"^chr([1-9]|1[0-9]|2[0-2])$")
    same = V["chrom37"] == V["chrom"]
    V["lifted"] = (V["lift"] == "ok") & auto & same
    frac = float(V.loc[auto, "lifted"].mean())
    say("[C] autosomal peaks %d; lifted %d (%.2f%%); reasons: %s; moved to another chromosome %d" % (
        int(auto.sum()), int(V["lifted"].sum()), 100 * frac,
        ", ".join("%s %d" % (k, v) for k, v in V.loc[auto, "lift"].value_counts().items()),
        int(((V["lift"] == "ok") & auto & ~same).sum())))
    ln = (V["end"] - V["start"]).values
    ln37 = (V["end37"] - V["start37"]).values
    ok = V["lifted"].values
    if ok.any():
        say("[C] peak length GRCh38 median %d bp; |length change| > 10 bp in %.2f%% of lifted peaks" % (
            int(np.median(ln[ok])), 100 * float(np.mean(np.abs(ln37[ok] - ln[ok]) > 10))))
    g2 = frac >= G2_MIN_LIFT
    SUM["gates"]["G2_LIFT"] = bool(g2)
    SUM["results"]["C"] = {"autosomal": int(auto.sum()), "lifted": int(V["lifted"].sum()), "frac": frac}
    say("[C] G2 (>= %.0f%% lifted): %s" % (100 * G2_MIN_LIFT, "PASS" if g2 else "FAIL"))
    return V


# ----------------------------------------------------------------------------
# D  lineage labels
# ----------------------------------------------------------------------------
def gene_regions():
    G = {}
    for line in open(GENE_LOC_B37):
        t = line.split()
        if len(t) >= 6:
            G[t[5]] = (t[1], int(t[2]), int(t[3]), t[4])
    return G


def marker_matrix(V):
    G = gene_regions()
    genes, cols = [], []
    chrom37 = V["chrom37"].fillna("").astype(str).values
    s37 = V["start37"].values
    e37 = V["end37"].values
    lifted = V["lifted"].values
    for lin in LABELS:
        for g in MARKERS[lin]:
            if g not in G:
                continue
            c, a, b, strand = G[g]
            lo, hi = (a - GA_UPSTREAM, b) if strand != "-" else (a, b + GA_UPSTREAM)
            sel = np.where(lifted & (chrom37 == "chr" + c) & (s37 < hi) & (e37 > lo - 1))[0]
            if len(sel) == 0:
                continue
            genes.append((lin, g, len(sel)))
            cols.append(sel)
    P = sp.lil_matrix((len(V), len(cols)), dtype=np.float32)
    for j, sel in enumerate(cols):
        P[sel, j] = 1.0
    return P.tocsr(), genes


def load_lib(fk):
    X = sp.load_npz(CACHE / ("lib_%s.npz" % fk)).tocsr()
    bc = pd.read_csv(CACHE / ("lib_%s_barcodes.txt.gz" % fk), header=None, dtype=str).iloc[:, 0].values
    return X, bc


def tfidf(Xs, idf):
    B = (Xs > 0).astype(np.float32)
    rs = np.asarray(B.sum(1)).ravel()
    rs[rs == 0] = 1
    T = sp.diags(1.0 / rs) @ B @ sp.diags(idf.astype(np.float32))
    T.data = np.log1p(T.data * 1e4)
    return T.tocsr()


def author_labels(L):
    files = [p for p in GEO_DIR.rglob("*") if p.is_file() and ANNOT_RE.search(p.name)
             and re.search(r"\.(csv|tsv|txt)(\.gz)?$", p.name) and not p.name.endswith(".part")]
    if not files:
        say("[D] no author cell-annotation file in the GEO supplement")
        return None
    pat = {"ciliated": r"cili", "endothelium": r"endothel|\bec\d?\b|vascular|lymphatic", "immune": r"immun|\bnk\b|unk|lymph|macro|mono|myelo|t.?cell|\bic\d?\b",
           "stroma": r"strom|fibro", "epithelium": r"epith|lumin|glan"}
    best = None
    for f in files:
        try:
            T = pd.read_csv(f, sep=None, engine="python", dtype=str)
        except Exception as e:
            say("  [D] %s unreadable: %s" % (f.name, type(e).__name__))
            continue
        say("  [D] annotation file %s: %d rows, columns %s" % (f.name, len(T), ", ".join(T.columns[:12])))
        for c in T.columns:
            v = T[c].dropna().astype(str).str.lower()
            if v.empty or v.nunique() > 200:
                continue
            hit = sum(bool(re.search("|".join(pat.values()), x)) for x in v.unique())
            if hit >= 2 and (best is None or hit > best[2]):
                best = (f, c, hit, T)
    if best is None:
        return None
    f, c, _, T = best
    bcol = max(T.columns, key=lambda k: T[k].astype(str).str.contains(r"[ACGT]{12,}", regex=True).mean())
    lab = {}
    for b, x in zip(T[bcol].astype(str), T[c].astype(str).str.lower()):
        lin = next((k for k, p in pat.items() if re.search(p, x)), "unassigned")
        lab[re.sub(r"^.*?([ACGT]{12,}-?\d*).*$", r"\1", b)] = lin
    say("[D] author labels: %s column '%s' (barcodes in '%s'); %s" % (f.name, c, bcol, ", ".join(
        "%s %d" % (k, v) for k, v in pd.Series(list(lab.values())).value_counts().items())))
    return lab


def part_d(V, L):
    say("")
    say("-" * 72)
    say("[SECTION] D lineage labels (LSI + k-means, fixed marker panel)")
    f_lab = CACHE / "cell_labels.csv.gz"
    if f_lab.exists():
        C = pd.read_csv(f_lab)
        say("[D] cached labels: %d cells" % len(C))
    else:
        from sklearn.cluster import KMeans
        from sklearn.decomposition import TruncatedSVD
        P, genes = marker_matrix(V)
        say("[D] marker genes with lifted peaks: %d of %d (%s)" % (len(genes), sum(len(v) for v in MARKERS.values()),
             "; ".join("%s %d" % (lin, sum(1 for g in genes if g[0] == lin)) for lin in LABELS)))
        top = np.argsort(-V["n_cells"].values, kind="stable")[:N_LSI_PEAKS]
        n_all = int(L["cells_kept"].sum())
        idf = np.log1p(n_all / (1.0 + V["n_cells"].values[top]))
        rng = np.random.default_rng(SEED)
        frac = min(1.0, N_LSI_SUB / max(n_all, 1))
        subs, meta = [], []
        for _, r in L.iterrows():
            X, bc = load_lib(r["file_key"])
            n = X.shape[0]
            pick = np.sort(rng.choice(n, size=max(1, int(round(frac * n))), replace=False)) if n else np.array([], int)
            subs.append(tfidf(X[pick][:, top], idf))
        svd = TruncatedSVD(n_components=N_COMP, algorithm="randomized", n_iter=5, random_state=SEED)
        svd.fit(sp.vstack(subs).tocsr())
        Z, GA, rows = [], [], []
        for _, r in L.iterrows():
            X, bc = load_lib(r["file_key"])
            tot = np.asarray(X.sum(1)).ravel()
            Z.append(svd.transform(tfidf(X[:, top], idf)))
            ga = np.asarray((X @ P).todense()) if P.shape[1] else np.zeros((X.shape[0], 0))
            GA.append(np.log1p(ga / np.maximum(tot, 1)[:, None] * 1e4))
            rows.append(pd.DataFrame({"lib": r["lib"], "barcode": bc, "total": tot}))
        Z = np.vstack(Z)
        GA = np.vstack(GA)
        C = pd.concat(rows, ignore_index=True)
        r1 = np.corrcoef(Z[:, 0], np.log1p(C["total"].values))[0, 1]
        keepc = list(range(1, N_COMP)) if abs(r1) > 0.5 else list(range(N_COMP))
        say("[D] LSI: %d peaks, %d components fitted on %d cells; component 1 vs depth r = %.2f -> %s" % (
            len(top), N_COMP, sum(s.shape[0] for s in subs), r1, "dropped" if abs(r1) > 0.5 else "kept"))
        Zn = Z[:, keepc] / np.maximum(np.linalg.norm(Z[:, keepc], axis=1, keepdims=True), 1e-12)
        gl = np.array([g[0] for g in genes])
        for k in (K_CLUST, K_CLUST_FALLBACK):
            km = KMeans(n_clusters=min(k, len(C)), n_init=10, random_state=SEED).fit(Zn)
            C["cluster"] = km.labels_
            mean = pd.DataFrame(GA).groupby(C["cluster"].values).mean().values
            z = (mean - mean.mean(0)) / np.maximum(mean.std(0), 1e-9)
            score = np.column_stack([z[:, gl == lin].mean(1) if (gl == lin).any() else np.full(len(z), -9) for lin in LABELS])
            o = np.argsort(-score, axis=1)
            topv = score[np.arange(len(score)), o[:, 0]]
            sec = score[np.arange(len(score)), o[:, 1]]
            lab = np.where((topv >= LABEL_MIN_SCORE) & (topv - sec >= LABEL_MIN_MARGIN), np.array(LABELS)[o[:, 0]], "unassigned")
            if set(LABELS) <= set(lab) or k == K_CLUST_FALLBACK:
                break
            say("[D] k = %d leaves %s without a cluster -> fallback k = %d" % (
                k, ", ".join(sorted(set(LABELS) - set(lab))), K_CLUST_FALLBACK))
        C["label_marker"] = lab[C["cluster"].values]
        for cl in range(len(lab)):
            m = C["cluster"].values == cl
            share = C.loc[m, "lib"].value_counts(normalize=True)
            say("  [D] cluster %2d %6d cells -> %-11s score %5.2f margin %5.2f | largest library %s %.0f%%" % (
                cl, int(m.sum()), lab[cl], topv[cl], topv[cl] - sec[cl], share.index[0] if len(share) else "-",
                100 * (share.iloc[0] if len(share) else 0)))
        au = author_labels(L)
        C["label_author"] = [au.get(re.sub(r"^.*?([ACGT]{12,}-?\d*).*$", r"\1", b), "missing") if au else "missing"
                             for b in C["barcode"]]
        C.to_csv(f_lab, index=False)
    cover = float((C["label_author"] != "missing").mean())
    use = "label_author" if cover >= AUTHOR_MIN_COVER else "label_marker"
    C["label"] = C[use]
    if cover > 0:
        ct = pd.crosstab(C["label_marker"], C["label_author"])
        say("[D] marker labels vs author labels (author coverage %.0f%%):" % (100 * cover))
        for line in ct.to_string().splitlines():
            say("    " + line)
    say("[D] labels used: %s" % ("author" if use == "label_author" else "marker panel"))
    n = len(C)
    comp = C["label"].value_counts()
    pub_n = sum(PUBLISHED.values())
    say("[D] composition (ours vs published): %s" % "; ".join(
        "%s %d (%.1f%% vs %.1f%%)" % (k, int(comp.get(k, 0)), 100 * comp.get(k, 0) / n, 100 * PUBLISHED[k] / pub_n) for k in LABELS)
        + "; unassigned %d (%.1f%%)" % (int(comp.get("unassigned", 0)), 100 * comp.get("unassigned", 0) / n))
    fs = 100 * comp.get("stroma", 0) / n
    fe = 100 * (comp.get("epithelium", 0) + comp.get("ciliated", 0)) / n
    ps = 100 * PUBLISHED["stroma"] / pub_n
    pe = 100 * (PUBLISHED["epithelium"] + PUBLISHED["ciliated"]) / pub_n
    g1 = (abs(fs - ps) <= G1_TOL_PP and abs(fe - pe) <= G1_TOL_PP and comp.get("immune", 0) >= G1_MIN_RARE
          and comp.get("endothelium", 0) >= G1_MIN_RARE and comp.get("unassigned", 0) / n <= G1_MAX_UNASSIGNED)
    SUM["gates"]["G1_COMPOSITION"] = bool(g1)
    SUM["results"]["D"] = {"labels": use, "n_cells": n, "composition": {k: int(v) for k, v in comp.items()}}
    say("[D] G1 (stroma and epithelium within %.0f pp of published, immune and endothelium >= %d cells, unassigned <= %.0f%%): %s" % (
        G1_TOL_PP, G1_MIN_RARE, 100 * G1_MAX_UNASSIGNED, "PASS" if g1 else "FAIL"))
    return C


# ----------------------------------------------------------------------------
# E  biopsy-level pseudobulk
# ----------------------------------------------------------------------------
def part_e(S, V, L, C):
    say("")
    say("-" * 72)
    say("[SECTION] E biopsy-level pseudobulk per lineage")
    lin_of = C["label"].replace({"ciliated": "epithelium"}).values
    um = library_units(S, L["lib"].tolist())
    keys, vecs, ncell = [], [], []
    for _, r in L.iterrows():
        X, bc = load_lib(r["file_key"])
        lab = lin_of[(C["lib"] == r["lib"]).values]
        unit = um[r["lib"]][1]
        for lin in LINEAGES:
            m = lab == lin
            if m.sum() == 0:
                continue
            keys.append((unit, lin))
            vecs.append(np.asarray(X[m].sum(0)).ravel())
            ncell.append(int(m.sum()))
    K = pd.DataFrame(keys, columns=["unit", "lineage"])
    K["cells"] = ncell
    PB = {}
    NC = {}
    for (u, lin), idx in K.groupby(["unit", "lineage"]).groups.items():
        PB[(u, lin)] = np.sum([vecs[i] for i in idx], axis=0)
        NC[(u, lin)] = int(K.loc[idx, "cells"].sum())
    units = sorted(set(k[0] for k in PB))
    tab = pd.DataFrame([[NC.get((u, lin), 0) for lin in LINEAGES] for u in units], index=units, columns=LINEAGES)
    say("[E] cells per biopsy unit and lineage:")
    for line in tab.to_string().splitlines():
        say("    " + line)
    ok = {lin: [u for u in units if NC.get((u, lin), 0) >= MIN_CELLS_UNIT] for lin in LINEAGES}
    use = [lin for lin in LINEAGES if len(ok[lin]) >= MIN_UNITS]
    say("[E] units with >= %d cells: %s" % (MIN_CELLS_UNIT, ", ".join("%s %d" % (lin, len(ok[lin])) for lin in LINEAGES)))
    g3 = "stroma" in use and "epithelium" in use
    SUM["gates"]["G3_UNITS"] = bool(g3)
    SUM["results"]["E"] = {"units": len(units), "units_per_lineage": {k: len(v) for k, v in ok.items()}, "lineages_used": use}
    dropped = [lin for lin in LINEAGES if lin not in use]
    say("[E] G3 (stroma and epithelium in >= %d units): %s%s" % (
        MIN_UNITS, "PASS" if g3 else "FAIL", ("; dropped as comparators: " + ", ".join(dropped)) if dropped else ""))
    tab.to_csv(OUT / "e_cells_per_unit.csv")
    return PB, NC, ok, use


# ----------------------------------------------------------------------------
# F  lineage-specific annotations
# ----------------------------------------------------------------------------
def part_f(V, PB, ok, use):
    say("")
    say("-" * 72)
    say("[SECTION] F lineage-specific peak annotations (biopsy-paired specificity)")
    allm = V["lifted"].values & (V["total"].values >= MIN_PEAK_TOTAL)
    idx = np.where(allm)[0]
    say("[F] ALL (control) = lifted autosomal peaks with >= %d counts: %d peaks, %.1f Mb" % (
        MIN_PEAK_TOTAL, len(idx), float((V["end37"].values[idx] - V["start37"].values[idx]).sum()) / 1e6))
    logc = {}
    for k, v in PB.items():
        if k[1] in use and k[0] in ok[k[1]]:
            logc[k] = np.log2(v[idx] / max(v.sum(), 1) * 1e6 + 1)
    T = pd.DataFrame({"chrom37": V["chrom37"].values[idx], "start37": V["start37"].values[idx], "end37": V["end37"].values[idx],
                      "chrom": V["chrom"].values[idx], "start": V["start"].values[idx], "end": V["end"].values[idx]})
    sets = {"ALL": np.ones(len(idx), bool)}
    for lin in use:
        D = []
        for u in ok[lin]:
            others = [logc[(u, o)] for o in use if o != lin and (u, o) in logc]
            if others:
                D.append(logc[(u, lin)] - np.mean(others, axis=0))
        D = np.array(D)
        n = D.shape[0]
        sd = D.std(0, ddof=1) if n > 1 else np.zeros(D.shape[1])
        s0 = float(np.median(sd))
        t = D.mean(0) / np.sqrt((sd ** 2 + s0 ** 2) / n)
        T["t_" + lin] = t
        T["mean_d_" + lin] = D.mean(0)
        order = np.argsort(-t, kind="stable")
        for q in TOP_Q:
            m = np.zeros(len(idx), bool)
            m[order[:int(round(len(idx) * q / 100.0))]] = True
            sets["%s_S%d" % (lin, q)] = m
        say("[F] %-11s %d units; s0 %.3f; top-10%% t >= %.2f; mean log2 difference in top 10%% %.2f" % (
            lin, n, s0, float(t[order[int(round(len(idx) * 0.1)) - 1]]), float(D.mean(0)[order[:int(round(len(idx) * 0.1))]].mean())))
    for k, m in sets.items():
        T[k] = m.astype(np.int8)
    lin10 = ["%s_S10" % lin for lin in use]
    J = pd.DataFrame([[float((T[a] & T[b]).sum() / max((T[a] | T[b]).sum(), 1)) for b in lin10] for a in lin10],
                     index=lin10, columns=lin10)
    say("[F] Jaccard overlap of the top-10% sets:")
    for line in J.round(3).to_string().splitlines():
        say("    " + line)
    T.to_csv(CACHE / "annotations_v072.csv.gz", index=False)
    T.to_csv(OUT / "f_annotations.csv.gz", index=False)
    SUM["results"]["F"] = {"n_all": int(len(idx)), "sets": [k for k in sets]}
    return T, sets


# ----------------------------------------------------------------------------
# G  stromal culture ATAC check (GSE246585)
# ----------------------------------------------------------------------------
def overlaps(qc, qs, qe, rc, rs, re_):
    # True where query interval overlaps any reference interval (same chromosome)
    out = np.zeros(len(qc), bool)
    R = pd.DataFrame({"c": rc, "s": rs, "e": re_}).sort_values(["c", "s"])
    for c, g in R.groupby("c"):
        s = g["s"].to_numpy()
        e = np.maximum.accumulate(g["e"].to_numpy())
        sel = np.where(qc == c)[0]
        if len(sel) == 0:
            continue
        i = np.searchsorted(s, qe[sel], side="left") - 1
        out[sel] = (i >= 0) & (e[np.clip(i, 0, None)] > qs[sel])
    return out


def read_ensc():
    if not fetch(ENSC_URL, ENSC_FILE, timeout=1800):
        return None
    D = pd.read_csv(ENSC_FILE, sep="\t", comment="#", dtype=str)
    say("[G] GSE246585 matrix %d rows x %d columns; first columns: %s" % (D.shape[0], D.shape[1], ", ".join(D.columns[:6])))
    cl = {c.lower(): c for c in D.columns}
    ck = next((cl[k] for k in ("chr", "chrom", "seqnames", "chromosome") if k in cl), None)
    ek = next((cl[k] for k in ("end", "stop") if k in cl), None)
    if ck and "start" in cl and ek:
        P = pd.DataFrame({"chrom": D[ck], "start": pd.to_numeric(D[cl["start"]], errors="coerce"),
                          "end": pd.to_numeric(D[ek], errors="coerce")})
        if pd.to_numeric(D[cl["start"]], errors="coerce").min() >= 1:
            P["start"] = P["start"] - 1
    else:
        col = next((c for c in D.columns[:3] if parse_coord(D[c].iloc[0])), None)
        if col is None and parse_coord(D.columns[0]):
            D = pd.read_csv(ENSC_FILE, sep="\t", comment="#", dtype=str, header=None)
            col = D.columns[0]
        if col is None:
            say("[G] peak coordinates not recognised -> check skipped")
            return None
        P = pd.DataFrame([parse_coord(x) or ("NA", -1, -1) for x in D[col]], columns=["chrom", "start", "end"])
    P["chrom"] = [c if str(c).startswith("chr") else "chr" + str(c) for c in P["chrom"]]
    P = P.dropna()
    P = P[P["chrom"].str.match(r"^chr([1-9]|1[0-9]|2[0-2])$") & (P["start"] >= 0)].reset_index(drop=True)
    P[["start", "end"]] = P[["start", "end"]].astype(np.int64)
    return P


def part_g(V, T, use):
    say("")
    say("-" * 72)
    say("[SECTION] G check with endometrial stromal culture ATAC (GSE246585, all time points pooled)")
    from scipy.stats import fisher_exact
    P = read_ensc()
    if P is None or P.empty:
        SUM["gates"]["G4_ENSC"] = False
        return
    say("[G] %d stromal culture peaks (autosomes), median width %d bp" % (len(P), int(np.median(P["end"] - P["start"]))))
    a = V[V["lifted"]]
    f38 = overlaps(P["chrom"].values, P["start"].values, P["end"].values, a["chrom"].values, a["start"].values, a["end"].values).mean()
    f37 = overlaps(P["chrom"].values, P["start"].values, P["end"].values, a["chrom37"].values, a["start37"].values, a["end37"].values).mean()
    build = "GRCh38" if f38 >= f37 else "GRCh37"
    say("[G] culture peaks overlapping the scATAC universe: as GRCh38 %.1f%%, as GRCh37 %.1f%% -> treated as %s" % (
        100 * f38, 100 * f37, build))
    c, s, e = (("chrom", "start", "end") if build == "GRCh38" else ("chrom37", "start37", "end37"))
    hit = overlaps(T[c].values, T[s].values, T[e].values, P["chrom"].values, P["start"].values, P["end"].values)
    base = hit.mean()
    rows = []
    for lin in use:
        m = T["%s_S10" % lin].values.astype(bool)
        rows.append({"lineage": lin, "n": int(m.sum()), "frac_culture": float(hit[m].mean()), "vs_all": float(hit[m].mean() / max(base, 1e-12))})
        say("  [G] %-11s top-10%% peaks %6d: %.1f%% in stromal culture peaks (ALL %.1f%%; ratio %.2f)" % (
            lin, m.sum(), 100 * hit[m].mean(), 100 * base, hit[m].mean() / max(base, 1e-12)))
    R = pd.DataFrame(rows).sort_values("frac_culture", ascending=False)
    R.to_csv(OUT / "g_culture_overlap.csv", index=False)
    ms = T["stroma_S10"].values.astype(bool)
    me = T["epithelium_S10"].values.astype(bool)
    tab = [[int(hit[ms].sum()), int((~hit[ms]).sum())], [int(hit[me].sum()), int((~hit[me]).sum())]]
    orr, p = fisher_exact(tab, alternative="greater")
    g4 = (R["lineage"].iloc[0] == "stroma") and p < G4_ALPHA
    say("[G] stroma vs epithelium top-10%%: OR %.2f, one-sided Fisher p %s" % (orr, fmt_p(p)))
    say("[G] G4 (stroma highest overlap and > epithelium at p < %.2f): %s" % (G4_ALPHA, "PASS" if g4 else "FAIL"))
    SUM["gates"]["G4_ENSC"] = bool(g4)
    SUM["results"]["G"] = {"build": build, "frac38": float(f38), "frac37": float(f37), "or_stroma_vs_epi": float(orr), "p": float(p),
                           "overlap": rows}


# ----------------------------------------------------------------------------
# H  thin annotations and partitioned LD scores (v047 pattern)
# ----------------------------------------------------------------------------
def write_annotations(T, keys):
    pd.Series(keys).to_csv(LD_DIR / "columns.txt", index=False, header=False)
    for chrom in range(1, 23):
        f = LD_DIR / ("annot.%d.annot" % chrom)
        if f.exists() and (LD_DIR / ("thin.%d.annot" % chrom)).exists():
            continue
        bim = pd.read_csv(PLINK_DIR / ("1000G.EUR.QC.%d.bim" % chrom), sep=r"\s+", header=None,
                          names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
        bp = bim["BP"].to_numpy()
        df = bim[["CHR", "BP", "SNP", "CM"]].copy()
        sub = T[T["chrom37"] == "chr%d" % chrom]
        for k in keys:
            g = sub[sub[k] == 1].sort_values("start37")
            if g.empty:
                df[k] = 0
                continue
            s = g["start37"].to_numpy()
            e = np.maximum.accumulate(g["end37"].to_numpy())
            i = np.searchsorted(s, bp, side="left") - 1          # peaks with start < BP (BED start 0-based)
            df[k] = ((i >= 0) & (e[np.clip(i, 0, None)] >= bp)).astype(np.int8)
        df.to_csv(f, sep="\t", index=False)
        df[keys].to_csv(LD_DIR / ("thin.%d.annot" % chrom), sep="\t", index=False)


def l2_cmd(c, thin):
    cmd = NICE + [str(LDSC_PY), str(LDSC), "--l2", "--bfile", str(PLINK_DIR / ("1000G.EUR.QC.%d" % c)),
                  "--ld-wind-cm", str(LD_WIND_CM), "--annot", str(LD_DIR / (("thin.%d.annot" if thin else "annot.%d.annot") % c)),
                  "--out", str(LD_DIR / ("multi.%d" % c))]
    return cmd + ["--thin-annot"] if thin else cmd


def frq_maf(c):
    f = next(iter(FRQ_DIR.rglob("1000G.EUR.QC.%d.frq" % c)), None)
    if f is None:
        raise FileNotFoundError("frq for chr%d" % c)
    return pd.read_csv(f, sep=r"\s+")[["SNP", "MAF"]]


def write_M(c):
    mf, m5 = LD_DIR / ("multi.%d.l2.M" % c), LD_DIR / ("multi.%d.l2.M_5_50" % c)
    if mf.exists() and m5.exists() and mf.stat().st_size > 0 and m5.stat().st_size > 0:
        return
    an = pd.read_csv(LD_DIR / ("annot.%d.annot" % c), sep="\t")
    keys = [k for k in an.columns if k not in ("CHR", "BP", "SNP", "CM")]
    mm = an[["SNP"]].merge(frq_maf(c), on="SNP", how="left")
    common = mm["MAF"].fillna(0).values > 0.05
    mf.write_text("\t".join(str(int(x)) for x in an[keys].sum(0)) + "\n")
    m5.write_text("\t".join(str(int(x)) for x in an.loc[common, keys].sum(0)) + "\n")


def run_l2():
    if not LDSC_PY.exists() or not LDSC.exists():
        say("[H] LDSC not found at %s -> LD scores not computed" % LDSC_PY)
        return False
    snps = LD_DIR / "print_snps.txt"
    if not snps.exists():
        with open(HM3_LIST) as f, open(snps, "w") as o:
            for line in f:
                t = line.split()
                if t and t[0].startswith("rs"):
                    o.write(t[0] + "\n")
    t0 = time.time()
    todo = [c for c in range(1, 23) if not (LD_DIR / ("multi.%d.l2.ldscore.gz" % c)).exists()]
    thin = (LD_DIR / "thin.mode").exists()
    if 22 in todo:
        # as v047: full annotation file first, thin annotation if that fails
        msgs = []
        for thin in (False, True):
            rc, o, e = sh(l2_cmd(22, thin), timeout=3 * 3600, env=L2_ENV)
            if (LD_DIR / "multi.22.l2.ldscore.gz").exists():
                break
            msgs += (o + e).strip().splitlines()[-6:]
        if thin:
            (LD_DIR / "thin.mode").write_text("1")
        todo = [c for c in todo if c != 22]
        if not (LD_DIR / "multi.22.l2.ldscore.gz").exists():
            say("    [L2] chr22 FAILED in both annotation modes")
            for line in msgs:
                say("       " + line[:200])
            return False
        write_M(22)
    for c in todo:
        errf = open(LD_DIR / ("multi.%d.stderr.txt" % c), "w")
        p = subprocess.run(l2_cmd(c, thin), stdout=subprocess.DEVNULL, stderr=errf, text=True, env=L2_ENV)
        errf.close()
        if (LD_DIR / ("multi.%d.l2.ldscore.gz" % c)).exists():
            write_M(c)
        else:
            say("    [L2] chr%-2d rc=%d FAILED" % (c, p.returncode))
            logf = LD_DIR / ("multi.%d.log" % c)
            for line in (logf.read_text(errors="ignore") if logf.exists() else "").strip().splitlines()[-8:]:
                say("       " + line[:200])
            return False
    for c in range(1, 23):
        write_M(c)
    say("[H] LD scores: all 22 chromosomes present (%s annotation; %.0f min)" % ("thin" if thin else "full", (time.time() - t0) / 60))
    return True


def split_ldscores(keys):
    marker = LD_DIR / "split.done"
    if marker.exists():
        return
    hm3 = set(open(LD_DIR / "print_snps.txt").read().split())
    for c in range(1, 23):
        ld = pd.read_csv(LD_DIR / ("multi.%d.l2.ldscore.gz" % c), sep="\t")
        ld = ld[ld["SNP"].isin(hm3)].reset_index(drop=True)
        M = open(LD_DIR / ("multi.%d.l2.M" % c)).read().split()
        M550 = open(LD_DIR / ("multi.%d.l2.M_5_50" % c)).read().split()
        cols = [x for x in ld.columns if x not in ("CHR", "SNP", "BP")]
        if len(cols) != len(keys):
            raise RuntimeError("chr%d: %d LD columns vs %d keys" % (c, len(cols), len(keys)))
        for j, k in enumerate(keys):
            d = LD_DIR / "sets" / k
            d.mkdir(parents=True, exist_ok=True)
            ld[["CHR", "SNP", "BP", cols[j]]].rename(columns={cols[j]: k + "L2"}).to_csv(
                d / ("%s.%d.l2.ldscore.gz" % (k, c)), sep="\t", index=False, compression="gzip")
            open(d / ("%s.%d.l2.M" % (k, c)), "w").write(M[j] + "\n")
            open(d / ("%s.%d.l2.M_5_50" % (k, c)), "w").write(M550[j] + "\n")
    marker.write_text("ok")


def part_h(T, sets):
    say("")
    say("-" * 72)
    say("[SECTION] H thin annotations on the 1000 Genomes EUR panel and partitioned LD scores")
    keys = list(sets.keys())
    old = LD_DIR / "columns.txt"
    if old.exists() and old.read_text().split() != keys:
        say("[H] annotation columns changed -> previous LD files removed")
        for p in LD_DIR.iterdir():
            if p.is_file():
                p.unlink()
    write_annotations(T, keys)
    hm3 = set(open(HM3_LIST).read().split()) if HM3_LIST.exists() else set()
    tot = {k: 0 for k in keys}
    hm = {k: 0 for k in keys}
    nsnp = 0
    for c in range(1, 23):
        an = pd.read_csv(LD_DIR / ("annot.%d.annot" % c), sep="\t")
        nsnp += len(an)
        inhm = an["SNP"].isin(hm3).values
        for k in keys:
            tot[k] += int(an[k].sum())
            hm[k] += int(an.loc[inhm, k].sum())
    for k in keys:
        say("  [H] %-16s %8d reference SNPs (%.2f%%), %7d HapMap3 SNPs" % (k, tot[k], 100 * tot[k] / max(nsnp, 1), hm[k]))
    ok = run_l2()
    if ok:
        split_ldscores(keys)
        pd.DataFrame({"set_key": keys}).to_csv(LD_DIR / "manifest.csv", index=False)
    SUM["gates"]["LDSCORES"] = bool(ok)
    SUM["results"]["H"] = {"snps": tot, "hm3": hm}


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def finish(t0):
    with open(OUT / "p9_v072_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[GATES] " + ", ".join("%s %s" % (k, "PASS" if v else "FAIL") for k, v in SUM["gates"].items()))
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.61 before this run; no GWAS data read")
    say("[DONE] PAPER9_V072_COMPLETE")
    _fh.close()


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.72 | GSE183771 scATAC audit and lineage chromatin annotations (no GWAS) | " + TS)
    say("=" * 72)
    t0 = time.time()
    a = part_a()
    if a is None:
        return finish(t0)
    S, sup = a
    b = part_b(S, sup)
    if b is None:
        return finish(t0)
    U, L = b
    V = part_c(U)
    if V is None:
        return finish(t0)
    C = part_d(V, L)
    PB, NC, ok, use = part_e(S, V, L, C)
    if not SUM["gates"].get("G3_UNITS"):
        return finish(t0)
    T, sets = part_f(V, PB, ok, use)
    part_g(V, T, use)
    gates = all(SUM["gates"].get(k) for k in ("G1_COMPOSITION", "G2_LIFT", "G3_UNITS", "G4_ENSC"))
    if gates:
        part_h(T, sets)
    else:
        say("")
        say("[H] not run: a gate failed (annotations kept for inspection; no LD scores)")
    finish(t0)


if __name__ == "__main__":
    main()
