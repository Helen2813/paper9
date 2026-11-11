#!/usr/bin/env python3
# Paper 9 | v0.0.55 | HECA donor structure and LDSC rg (rules CHANGELOG v0.2.45)
#   E1 donor composition of every HECA status x label group (variant P = v030 / v041 cells; variant H = exact Koller cohort)
#   E2 labels pooled over status, donor-balanced cell sample (<= 50 cells per donor and type), Wilcoxon within compartment,
#      UPN / UP sets, MAGMA; leave-one-donor-out over all donors
#   E3 donor-paired pseudobulk DE with labels pooled over status (v054 pseudobulk cache), UPN, MAGMA; leave-one-donor-out
#   A3 LDSC cross-trait rg and intercepts (v054 used the weights folder as reference LD scores, which has no .M files)
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, LDSC 2 threads, nice 10). ~20-40 min; no downloads.

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
import traceback
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
RES = REPO / "results"
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
LDSC_PY = TOOLS / "ldsc_venv" / "bin" / "python"
LDSC = TOOLS / "ldsc_py3" / "ldsc.py"
REF = RAW / "ldsc_ref" / "unpacked"
MUNGED = PROC / "munged"
HECA_FILE = RAW / "heca" / "endometriumAtlasV2_cells_with_counts.h5ad"
HECA_DE = PROC / "genesets" / "heca_koller_v041_de"
GSE_DE = PROC / "genesets" / "gse266265_v050_de"
PB266 = PROC / "gse266265_pseudobulk"
GS_P = PROC / "genesets" / "heca_koller_v031" / "P"
LOCUS_DIR = RES / "locus_atlas"
STATE = RES / "benchmark" / "v041_replicates.csv"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = RES / "review" / ("p9_v055_" + TS)
WORK = INTERIM / "magma_v055"
CACHE = INTERIM / "v055_cache"
CACHE54 = INTERIM / "v054_cache"
LOG = REPO / "logs" / ("p9_v055_donors_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
ENV2 = dict(os.environ, OPENBLAS_NUM_THREADS="2", OMP_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
AUTOS = [str(c) for c in range(1, 23)]

# Fixed (CHANGELOG v0.2.44 / v0.2.45)
SEED = 20261109
N_ROT = 1000
ROT_MIN_FRAC = 0.05
N_PERM = 200
N_BOOT = 10000
MAX_COLS = 1500
LFC, PADJ, MIN_SET, N_MATCH = 1.2, 0.01, 20, 400
TRAITS_C = ["KOLLER_EUR", "GCST90205183", "LEIOMYOMA", "HYPERTENSION", "T2D"]
PAIRS = [("KOLLER_EUR", "T2D"), ("KOLLER_EUR", "HYPERTENSION"), ("KOLLER_EUR", "LEIOMYOMA"),
         ("GCST90205183", "T2D"), ("GCST90205183", "HYPERTENSION"), ("GCST90205183", "LEIOMYOMA")]
ENDO_PAIRS = [("KOLLER_EUR", "GCST90205183"), ("KOLLER_EUR", "ENDO_ALL"), ("GCST90205183", "ENDO_ALL")]
ROT_SETS = [("HECA", "UPN", "dStromal_mid/ENDO"), ("HECA", "UP", "dStromal_mid/ENDO"), ("HECA", "UPN", "dStromal_early/ENDO"),
            ("GSE", "UPN", "eF3"), ("GSE", "UPN", "dS3")]
PRIMARY_C = [("HECA", "UPN", "dStromal_mid/ENDO"), ("GSE", "UPN", "eF3")]
RG_TRAITS = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL", "LEIOMYOMA", "HYPERTENSION", "T2D"]
RG_FIRST = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL"]
TRAITS_LOO = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL"]
LOO_SHARE = 0.90
LOO_HECA_MAX_FAIL = 1
# GSE266265 (as v050)
MIN_CELLS_DONOR, MIN_DONORS, MIN_CPM = 10, 10, 1.0
GROUPS_G = [("endometrial stroma", r"^(eF\d|dS\d|fib C7)"), ("perivascular / MSC-like", r"^(Prv|VSMC)"), ("mesothelial label", r"^mesothelial"),
            ("epithelium", r"^(glandular|lumenal|ciliated|mid-secretory|MUC5B|TP63)"), ("endothelium", r"^(EC-|LEC)")]
COMP_G = [("stroma", r"^(eF\d|dS\d|fib C7)"), ("perivascular/SM", r"^(Prv|VSMC)"),
          ("epithelium", r"^(glandular|lumenal|ciliated|mid-secretory|MUC5B|TP63)"), ("endothelium", r"^(EC-|LEC)"),
          ("immune", r"^(MPhi|T|CD4|CD8|NK|pNK|cDC|pre-cDC|pDC|mDC|DC3|B cell|plasma|ILC|mast|monocytes|granulocytes|CTL)")]
# HECA (as v041 variant P)
DONOR_COL, TYPE_COL, PATHOLOGY_COL, STAGE_COL, EXCLUDE_STAGE = "genotype", "celltype", "Endometrial_pathology", "Binary Stage", "Hormones"
IMMUNE_RE = re.compile(r"Immune|^uM|uNK|^T_|^B_|Mast|DC|Mono|Macro|Neutro|Plasma|ILC|Lymphoid|Myeloid")
NOT_IMMUNE_RE = re.compile(r"Lymphatic")
MIN_CELLS_TYPE, CAP_CELLS_TYPE, SEED_HECA = 30, 400, 20261106
MIN_CELLS_PB, MIN_DONORS_PB = 10, 3
STROMA_H = r"^(dStromal|eStromal|sHormones|dHormones|Fibroblast)"
DECID_H = r"^dStromal"
PERIV_H = r"^(ePV|mPV|uSMC)"
# E (CHANGELOG v0.2.45)
TRAITS_E = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL", "LEIOMYOMA", "HYPERTENSION"]
FAMS_E = ("UPN", "UP")
BAL_CAP, BAL_MIN_CELLS, BAL_MIN_DONORS, SEED_BAL = 50, 10, 3, 20261110
E3_MIN_DONORS = 4
LOO_MAX_FAIL = 1
VASC_H = r"^(Arterial|Venous|Lymphatic|ePV|mPV|uSMC)"
# Locus test (as v042 / v044)
GENE_WIN_BP, DROP_WIN_BP, PCT_ASSIGN = 250_000, 500_000, 0.90

for d in (OUT, WORK, CACHE, LOG.parent):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V055_DONORS", "timestamp_utc": TS, "results": {}}
ERR = []


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def sh(cmd, timeout=None, env=None):
    try:
        p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout, env=env or ENV1)
        return p.returncode, p.stdout, p.stderr
    except FileNotFoundError as e:
        return 127, "", str(e)
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"


def gene_loc():
    ent2sym, sym2ent, pos, dup = {}, {}, {}, set()
    for line in open(GENE_LOC_B37):
        t = line.split()
        if len(t) >= 6 and t[1] in AUTOS:
            ent2sym[t[0]] = t[5]
            pos[t[0]] = (int(t[1]), int(t[2]), int(t[3]))
            if t[5] in sym2ent:
                dup.add(t[5])
            sym2ent[t[5]] = t[0]
    for d in dup:
        sym2ent.pop(d, None)
    return ent2sym, sym2ent, pos


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def phi_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2)) if np.isfinite(z) else np.nan


def tdist_sf(t, df):
    try:
        from scipy.stats import t as T
        return 2.0 * T.sf(np.abs(t), df)
    except Exception:
        return np.array([math.erfc(v / math.sqrt(2)) for v in np.abs(t)])


def wilson(k, n, z=1.959964):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def clopper(k, n, a=0.05):
    from scipy.stats import beta
    lo = 0.0 if k == 0 else float(beta.ppf(a / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a / 2, k + 1, n - k))
    return lo, hi


def mcnemar_exact(b, c):
    n = b + c
    if n == 0:
        return 1.0
    from scipy.stats import binomtest
    return float(binomtest(min(b, c), n, 0.5).pvalue)


def latest(base, pattern):
    fs = sorted(base.glob(pattern))
    return fs[-1] if fs else None


def short(label, n=22):
    return str(label).replace("SOX9_functionalis_", "SOX9f_").replace("dStromal_", "dS_")[:n]


def clean(label):
    return re.sub(r"\$|\\|\{|\}|_\{|\^", "", str(label))


def first_match(label, groups, default):
    for g, rx in groups:
        if re.search(rx, str(label)):
            return g
    return default


def prefix_of(dirname):
    d = REF / dirname
    for p in sorted(d.rglob("*.l2.ldscore.gz")) if d.exists() else []:
        m = re.match(r"^(.*\.)1\.l2\.ldscore\.gz$", p.name)
        if m:
            return str(p.parent / m.group(1))
    return None


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


def se_std(v):
    b, bs, se = v["beta"], v["beta_std"], v["se"]
    return abs(se * bs / b) if (np.isfinite(b) and b != 0 and np.isfinite(bs)) else np.nan


def magma_covar(raw, cf, tag, cond=None):
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


def magma_sets(raw, sf, tag):
    outp = WORK / tag
    gsa = Path(str(outp) + ".gsa.out")
    for model in (["--model", "condition-hide=ALLGENES"], []):
        if gsa.exists():
            gsa.unlink()
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--set-annot", sf] + model + ["--out", outp], timeout=7200)
        if gsa.exists():
            return parse_full(gsa)
    say("    [MAGMA_FAIL] sets %s %s" % (tag, (o + e).strip()[-300:]))
    return None


def covar_runs(raw, cols, uni_ent, tag, chunk=MAX_COLS):
    # cols: dict name -> int8 indicator over uni_ent; written in chunks of <= MAX_COLS columns
    names = list(cols)
    res = {}
    for i in range(0, len(names), chunk):
        part = names[i:i + chunk]
        cf = WORK / ("%s_c%02d.txt" % (tag, i // chunk))
        if not cf.exists():
            M = np.column_stack([cols[k] for k in part]).astype(np.int8)
            pd.DataFrame(M, index=uni_ent, columns=part).to_csv(cf, sep="\t", index_label="GENE")
        r = magma_covar(raw, cf, "%s_%s_c%02d" % (tag, raw.name.split("_35_10")[0], i // chunk))
        if r is None:
            return None
        res.update(r)
    return res


def drop_covar_files(tag):
    for f in WORK.glob(tag + "_c*.txt"):
        try:
            f.unlink()
        except OSError:
            pass


def raw_of(t):
    f = MAGMA_OUT / (t + "_35_10.genes.raw")
    return f if f.exists() else None


# ----------------------------------------------------------------------------
# Marker sets from cached DE tables
# ----------------------------------------------------------------------------
def heca_sets(sym2ent):
    DE = pd.read_csv(HECA_DE / "de_all.csv.gz")
    DE["gene"] = DE["gene"].astype(str)
    uni = sorted({g for g in DE["gene"] if g in sym2ent})
    return DE, uni


def fam_sets(d, uset, fams=("UP", "UPN")):
    out = {}
    sig = d["padj"] < PADJ if "padj" in d else None
    for fam in fams:
        if fam == "UP":
            g = d.loc[(d["lfc"] > LFC) & sig, "gene"]
        elif fam == "UPN":
            g = d.sort_values("score", ascending=False)["gene"].head(N_MATCH)
        else:
            continue
        gl = [x for x in g.astype(str) if x in uset]
        if len(gl) >= MIN_SET:
            out[fam] = gl
    return out


def gse_sets(sym2ent):
    DE = pd.read_csv(GSE_DE / "de_all.csv.gz")
    DE["gene"] = DE["gene"].astype(str)
    uni = sorted({g for g in DE["gene"] if g in sym2ent})
    return DE, uni


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / max(1, len(a | b))



# ----------------------------------------------------------------------------
# Helpers reused from v054
# ----------------------------------------------------------------------------
def parse_rg_log(log):
    txt = log.read_text().splitlines()
    i = next((k for k, l in enumerate(txt) if "Summary of Genetic Correlation Results" in l), None)
    if i is None:
        return []
    hdr = txt[i + 1].split()
    rows = []
    for l in txt[i + 2:]:
        t = l.split()
        if len(t) != len(hdr):
            break
        d = dict(zip(hdr, t))
        rows.append(d)
    return rows


def rank_items(res, keys):
    it = sorted([(lab, res[k]) for k, lab in keys if k in res], key=lambda x: x[1]["p"])
    return it


def dec(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def obs_cat(f, name):
    import h5py
    obs = f["obs"]
    g = obs[name]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = np.array([dec(c) for c in g["categories"][:]], dtype=object)
        codes = g["codes"][:]
    elif "__categories" in obs and name in obs["__categories"]:
        cats = np.array([dec(c) for c in obs["__categories"][name][:]], dtype=object)
        codes = g[:]
    else:
        v = g[:]
        return np.array([dec(x) for x in v], dtype=object)
    out = np.empty(len(codes), dtype=object)
    ok = codes >= 0
    out[ok] = cats[codes[ok]]
    out[~ok] = None
    return out


def find_dataset_col(f, cols):
    for c in cols:
        try:
            v = obs_cat(f, c)
        except Exception:
            continue
        if v.dtype != object:
            continue
        cats = pd.unique(v[pd.notna(v)])
        if len(cats) <= 200:
            hit = [x for x in cats if re.search(r"14039|Marec", str(x))]
            if hit:
                return c, hit
    return None, None


def status_of(path):
    status = np.full(len(path), None, dtype=object)
    status[path == "C"] = "CTL"
    is_case = np.array([bool(re.match(r"^E", str(x))) or "endometriosis" in str(x).lower() for x in path])
    status[is_case] = "ENDO"
    return status


def wilcoxon_de(X, ctype_sub, gnames, s, cp):
    import anndata as ad
    import scanpy as sc
    A = ad.AnnData(X=X.copy())
    A.var_names = list(gnames)
    A.obs["type"] = pd.Categorical(ctype_sub)
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    sc.tl.rank_genes_groups(A, "type", method="wilcoxon", reference="rest", use_raw=False, pts=False)
    parts = []
    for l in A.obs["type"].cat.categories:
        df = sc.get.rank_genes_groups_df(A, group=l)
        parts.append(pd.DataFrame({"group": "%s/%s" % (l, s), "label": l, "status": s, "compartment": cp, "gene": df["names"].astype(str).values,
                                   "score": df["scores"].values, "lfc": df["logfoldchanges"].values, "padj": df["pvals_adj"].values}))
    return pd.concat(parts, ignore_index=True)




# ----------------------------------------------------------------------------
# A3. LDSC cross-trait rg (fixed: LD scores with .M files as reference, HapMap3 weights as regression weights)
# ----------------------------------------------------------------------------
def part_a3():
    say("")
    say("-" * 72)
    say("[SECTION] A3 LDSC cross-trait genetic correlation and intercept (reference LDscore, HapMap3 no-MHC weights)")
    ref, w = prefix_of("LDscore"), prefix_of("1000G_Phase3_weights_hm3_no_MHC")
    if not ref or not w or not LDSC.exists() or not LDSC_PY.exists():
        say("[A3] LDSC, reference LD scores or weights missing (ref %s, w %s); skipped" % (ref, w))
        return
    avail = [t for t in RG_TRAITS if (MUNGED / (t + ".sumstats.gz")).exists()]
    say("[A3] reference %s | weights %s | traits %s" % (Path(ref).parent.name, Path(w).parent.name, ", ".join(avail)))
    name_of = {str(MUNGED / (t + ".sumstats.gz")): t for t in avail}
    rows = []
    t0 = time.time()
    for t in RG_FIRST:
        if t not in avail:
            continue
        others = [x for x in avail if x != t and not (x in RG_FIRST and RG_FIRST.index(x) < RG_FIRST.index(t))]
        if not others:
            continue
        out = CACHE / ("rg_" + t)
        log = Path(str(out) + ".log")
        if not (log.exists() and "Summary of Genetic Correlation Results" in log.read_text(errors="ignore")):
            cmd = NICE + [LDSC_PY, LDSC, "--rg", ",".join(str(MUNGED / (x + ".sumstats.gz")) for x in [t] + others),
                          "--ref-ld-chr", ref, "--w-ld-chr", w, "--out", out]
            rc, o, e = sh(cmd, timeout=6 * 3600, env=ENV2)
            txt = log.read_text(errors="ignore") if log.exists() else (o + e)
            if "Summary of Genetic Correlation Results" not in txt:
                say("  [A3] %s rc=%d; last lines of the LDSC log:" % (t, rc))
                for line in (txt + "\n" + e).strip().splitlines()[-15:]:
                    say("       " + line[:160])
                continue
        for d in parse_rg_log(log):
            a = name_of.get(d.get("p1"), Path(d.get("p1", "")).name.replace(".sumstats.gz", ""))
            b = name_of.get(d.get("p2"), Path(d.get("p2", "")).name.replace(".sumstats.gz", ""))

            def f(k):
                try:
                    return float(d.get(k, "nan"))
                except ValueError:
                    return np.nan
            rows.append({"trait_a": a, "trait_b": b, "rg": f("rg"), "rg_se": f("se"), "rg_p": f("p"), "gcov_int": f("gcov_int"),
                         "gcov_int_se": f("gcov_int_se"), "h2_obs_b": f("h2_obs"), "h2_int_b": f("h2_int")})
    if not rows:
        say("[A3] no LDSC rg results")
        return
    G = pd.DataFrame(rows)
    G.to_csv(OUT / "a3_ldsc_rg.csv", index=False)
    for _, r in G.iterrows():
        say("    %-13s %-13s rg %+.3f (se %.3f, p %.1e) | cross-trait intercept %+.3f (se %.3f)" % (
            r["trait_a"][:13], r["trait_b"][:13], r["rg"], r["rg_se"], r["rg_p"], r["gcov_int"], r["gcov_int_se"]))
    say("[A3] %d pairs (%.0f min); intercept = sample overlap / shared confounding, rg = shared genetics" % (len(G), (time.time() - t0) / 60))
    SUM["results"]["a3"] = G.to_dict("records")


# ----------------------------------------------------------------------------
# E. HECA donor structure
# ----------------------------------------------------------------------------
def heca_open():
    import h5py
    f = h5py.File(HECA_FILE, "r")
    cols = [k for k in f["obs"].keys() if k not in ("_index", "__categories")]
    dataset_col, dataset_vals = find_dataset_col(f, cols)
    O = {"donor": obs_cat(f, DONOR_COL), "ctype": obs_cat(f, TYPE_COL), "stage": obs_cat(f, STAGE_COL)}
    O["status"] = status_of(obs_cat(f, PATHOLOGY_COL))
    base = pd.notna(O["donor"]) & pd.notna(O["ctype"]) & pd.notna(O["status"])
    if dataset_col:
        base &= np.isin(obs_cat(f, dataset_col), dataset_vals)
    O["keep_H"] = base
    O["keep_P"] = base & (O["stage"] != EXCLUDE_STAGE)
    return f, O


def heca_genes(f, sym2ent):
    import h5py
    var = f["var"]
    gobj = var[dec(var.attrs.get("_index", "_index"))]
    if isinstance(gobj, h5py.Group):
        gobj = gobj["values"] if "values" in gobj else gobj[list(gobj.keys())[0]]
    genes = np.array([dec(x) for x in gobj[:]], dtype=object)
    gidx = np.where(np.array([g in sym2ent for g in genes]))[0]
    return genes, gidx


def read_rows(f, want):
    layer = f["layers"]["counts"]
    indptr = layer["indptr"][:]
    rowmat = {}
    step = 10000
    n = len(indptr) - 1
    for a in range(0, n, step):
        b = min(n, a + step)
        lo, hi = np.searchsorted(want, a), np.searchsorted(want, b)
        if hi <= lo:
            continue
        s0, e0 = int(indptr[a]), int(indptr[b])
        idx = np.asarray(layer["indices"][s0:e0])
        val = np.asarray(layer["data"][s0:e0], dtype=np.float32)
        loc = indptr[a:b + 1] - s0
        for i in want[lo:hi]:
            r = i - a
            rowmat[int(i)] = (idx[loc[r]:loc[r + 1]].copy(), val[loc[r]:loc[r + 1]].copy())
        del idx, val
    return rowmat


def comp_of(l):
    return "immune" if IMMUNE_RE.search(l) and not NOT_IMMUNE_RE.search(l) else "endometrial"


def part_e1(O):
    say("")
    say("-" * 72)
    say("[SECTION] E1 donor composition of the HECA status x label groups (all cells; variant P = v030 / v041 cells, H = exact Koller cohort)")
    rows = []
    for v in ("P", "H"):
        k = O["keep_" + v]
        df = pd.DataFrame({"status": O["status"][k], "label": O["ctype"][k], "donor": O["donor"][k].astype(str)})
        say("[E1] variant %s: %d cells; donors CTL %d, ENDO %d" % (v, len(df), df[df["status"] == "CTL"]["donor"].nunique(), df[df["status"] == "ENDO"]["donor"].nunique()))
        for (s, l), g in df.groupby(["status", "label"]):
            if len(g) < MIN_CELLS_TYPE:
                continue
            c = g["donor"].value_counts()
            share = c / c.sum()
            rows.append({"variant": v, "group": "%s/%s" % (l, s), "label": l, "status": s, "cells": int(c.sum()), "donors": int(len(c)),
                         "donors_ge10": int((c >= 10).sum()), "top_donor": c.index[0], "top_share": float(share.iloc[0]),
                         "effective_donors": float(1.0 / (share ** 2).sum())})
    E = pd.DataFrame(rows)
    E.to_csv(OUT / "e1_heca_donor_composition.csv", index=False)
    res = {}
    for v in ("P", "H"):
        e = E[E["variant"] == v]
        if not len(e):
            continue
        say("[E1] variant %s: %d groups with >= %d cells | top donor share > 0.8 in %d (%.0f%%), > 0.5 in %d (%.0f%%) | median effective donors %.1f" % (
            v, len(e), MIN_CELLS_TYPE, int((e["top_share"] > 0.8).sum()), 100 * (e["top_share"] > 0.8).mean(), int((e["top_share"] > 0.5).sum()),
            100 * (e["top_share"] > 0.5).mean(), e["effective_donors"].median()))
        st = e[e["label"].astype(str).str.contains(STROMA_H)].sort_values(["label", "status"])
        for _, r in st.iterrows():
            say("    %-26s %5d cells, %2d donors (%2d with >= 10), top %-8s %3.0f%%, effective donors %.1f" % (
                r["group"][:26], r["cells"], r["donors"], r["donors_ge10"], str(r["top_donor"])[:8], 100 * r["top_share"], r["effective_donors"]))
        res[v] = {"groups": len(e), "top_gt_0.8": int((e["top_share"] > 0.8).sum()), "top_gt_0.5": int((e["top_share"] > 0.5).sum()),
                  "median_eff": float(e["effective_donors"].median())}
    SUM["results"]["e1"] = res


def bal_plan(O):
    rng = np.random.default_rng(SEED_BAL)
    k = O["keep_P"]
    labels = sorted(pd.unique(O["ctype"][k]))
    plan, comp_rows = {}, []
    for cp in ("endometrial", "immune"):
        pick, types = [], []
        for l in labels:
            if comp_of(l) != cp:
                continue
            idx_l = np.where(k & (O["ctype"] == l))[0]
            dn = O["donor"][idx_l].astype(str)
            sel = []
            for d in sorted(set(dn)):
                ii = idx_l[dn == d]
                if len(ii) < BAL_MIN_CELLS:
                    continue
                if len(ii) > BAL_CAP:
                    ii = np.sort(rng.choice(ii, BAL_CAP, replace=False))
                sel.append(ii)
            if len(sel) >= BAL_MIN_DONORS and sum(len(x) for x in sel) >= MIN_CELLS_TYPE:
                pick.append(np.concatenate(sel))
                types.append(l)
                comp_rows.append({"label": l, "compartment": cp, "donors": len(sel), "cells": int(sum(len(x) for x in sel)),
                                  "top_share": max(len(x) for x in sel) / sum(len(x) for x in sel)})
        if len(types) >= 2:
            plan[cp] = (np.sort(np.concatenate(pick)), types)
    return plan, pd.DataFrame(comp_rows)


def eligible(ct, dn, min_donors):
    d = pd.DataFrame({"t": ct, "d": dn})
    nd = d.groupby("t")["d"].nunique()
    nc = d["t"].value_counts()
    return [t for t in nd.index if nd[t] >= min_donors and nc[t] >= MIN_CELLS_TYPE]


def part_e2_de(f, O, sym2ent):
    import scipy.sparse as sp
    f_de, f_cp = CACHE / "e2_balanced_de.csv.gz", CACHE / "e2_balanced_composition.csv"
    if f_de.exists() and f_cp.exists():
        say("[E2] cached donor-balanced DE (%s)" % CACHE.name)
        return pd.read_csv(f_de), pd.read_csv(f_cp)
    t00 = time.time()
    plan, CP = bal_plan(O)
    genes, gidx = heca_genes(f, sym2ent)
    gnames = genes[gidx]
    want = np.unique(np.concatenate([v[0] for v in plan.values()]))
    say("[E2] donor-balanced sample: %d cells, %d types (%s); <= %d cells per donor and type, donors with >= %d cells, types with >= %d donors" % (
        len(want), sum(len(v[1]) for v in plan.values()), ", ".join("%s %d" % (k, len(v[1])) for k, v in plan.items()), BAL_CAP, BAL_MIN_CELLS, BAL_MIN_DONORS))
    rowmat = read_rows(f, want)
    donors_all = sorted({str(x) for v in plan.values() for x in O["donor"][v[0]]})
    parts = []
    for cp, (allidx, types) in plan.items():
        data = [rowmat[int(i)][1] for i in allidx]
        ind = [rowmat[int(i)][0] for i in allidx]
        ptr = np.concatenate([[0], np.cumsum([len(x) for x in ind])])
        X = sp.csr_matrix((np.concatenate(data), np.concatenate(ind), ptr), shape=(len(allidx), len(genes)))[:, gidx]
        ct = np.asarray(O["ctype"][allidx]).astype(str)
        dn = np.asarray(O["donor"][allidx]).astype(str)
        present = set(dn)
        full_de = None
        for itn in ["FULL"] + donors_all:
            if itn != "FULL" and itn not in present and full_de is not None:
                parts.append(full_de.assign(iteration=itn))
                continue
            m = np.ones(len(allidx), dtype=bool) if itn == "FULL" else (dn != itn)
            okt = eligible(ct[m], dn[m], BAL_MIN_DONORS)
            if len(okt) < 2:
                continue
            m &= np.isin(ct, okt)
            t0 = time.time()
            D = wilcoxon_de(X[np.where(m)[0]], ct[m], gnames, "ALL", cp)
            D["group"] = D["label"]
            D["iteration"] = itn
            parts.append(D)
            if itn == "FULL":
                full_de = D
            say("    [E2] %-11s without %-10s %5d cells, %2d types (%.1f min)" % (cp, itn if itn != "FULL" else "(none)", int(m.sum()), len(okt), (time.time() - t0) / 60))
        del X
    del rowmat
    DE = pd.concat(parts, ignore_index=True)
    DE[["iteration", "group", "label", "compartment", "gene", "score", "lfc", "padj"]].to_csv(f_de, index=False)
    CP.to_csv(f_cp, index=False)
    say("[E2] DE cached (%.1f min)" % ((time.time() - t00) / 60))
    return DE, CP


def pb_de(PB, PBG, genes_sym, min_donors):
    # donor-paired pseudobulk DE: type vs the rest of the same donor's cells in the same compartment; moderated paired t
    PBG = PBG.assign(comp=PBG["label"].astype(str).map(comp_of))
    out = {}
    for cp, g in PBG.groupby("comp"):
        tot = {d: PB[gg.index.to_numpy()].astype(np.float64).sum(axis=0) for d, gg in g.groupby("donor")}
        ok = g[g["n_cells"] >= MIN_CELLS_PB]
        for l, r in ok.groupby("label"):
            if r["donor"].nunique() < min_donors:
                continue
            X = PB[r.index.to_numpy()].astype(np.float64)
            Rr = np.vstack([tot[d] for d in r["donor"]]) - X
            Ct = X / np.maximum(X.sum(axis=1, keepdims=True), 1) * 1e6
            Cr = Rr / np.maximum(Rr.sum(axis=1, keepdims=True), 1) * 1e6
            expr = (Ct.mean(axis=0) >= MIN_CPM) | (Cr.mean(axis=0) >= MIN_CPM)
            D = np.log2(Ct[:, expr] + 1) - np.log2(Cr[:, expr] + 1)
            n = D.shape[0]
            m, se = D.mean(axis=0), D.std(axis=0, ddof=1) / math.sqrt(n)
            s0 = float(np.median(se)) if len(se) else 0.0
            out[str(l)] = (pd.DataFrame({"gene": genes_sym[expr], "score": m / (se + s0)}), n)
    return out


def rank_rows(res, meta, itn, fam, t, extra=None):
    keys = [(k, m[2]) for k, m in meta.items() if m[0] == itn and m[1] == fam]
    items = rank_items(res, keys)
    if not items:
        return None
    labs = [l for l, _ in items]
    i_ds = next((i for i, l in enumerate(labs) if l == "dStromal_mid"), None)
    bs = next((i for i, l in enumerate(labs) if re.search(STROMA_H, l)), None)
    bp = next((i for i, l in enumerate(labs) if re.search(PERIV_H, l)), None)
    bv = next((i for i, l in enumerate(labs) if re.search(VASC_H, l)), None)
    return {"trait": t, "family": fam, "left_out": itn, "first": labs[0], "first_p": items[0][1]["p"],
            "first_stromal": bool(re.search(STROMA_H, labs[0])), "first_decidual": bool(re.search(DECID_H, labs[0])),
            "first_vascular": bool(re.search(VASC_H, labs[0])), "first_smooth_muscle": bool(re.search(r"^(uSMC|mPV|ePV)", labs[0])),
            "rank_dStromal_mid": (i_ds + 1) if i_ds is not None else np.nan, "best_stromal_rank": (bs + 1) if bs is not None else np.nan,
            "best_stromal": labs[bs] if bs is not None else "", "p_best_stromal": items[bs][1]["p"] if bs is not None else np.nan,
            "best_perivascular_rank": (bp + 1) if bp is not None else np.nan, "best_vascular_rank": (bv + 1) if bv is not None else np.nan,
            "top3": "; ".join("%s %.0e" % (short(l, 18), x["p"]) for l, x in items[:3]), **(extra or {})}


def report_design(B, design, verdict_tag):
    out = {}
    for t in TRAITS_E:
        for fam in FAMS_E:
            g = B[(B["trait"] == t) & (B["family"] == fam)]
            f = g[g["left_out"] == "FULL"]
            l = g[g["left_out"] != "FULL"]
            if not len(f):
                continue
            f = f.iloc[0]
            say("  %-13s %-3s full: %s | LOO (n=%d): stromal first %d, decidual first %d, perivascular above stroma %d | best stromal p up to %.1e" % (
                t[:13], fam, f["top3"], len(l), int(l["first_stromal"].sum()), int(l["first_decidual"].sum()),
                int((l["best_perivascular_rank"] < l["best_stromal_rank"]).sum()), l["p_best_stromal"].max() if len(l) else np.nan))
            if t == "KOLLER_EUR":
                for _, r in l[~l["first_stromal"]].iterrows():
                    say("       without %-10s %s" % (r["left_out"], r["top3"]))
            out["%s_%s" % (t, fam)] = {"full_first": f["first"], "full_stromal": bool(f["first_stromal"]), "n": len(l),
                                       "stromal_first": int(l["first_stromal"].sum()), "decidual_first": int(l["first_decidual"].sum()),
                                       "perivascular_above": int((l["best_perivascular_rank"] < l["best_stromal_rank"]).sum())}
    v = out.get("KOLLER_EUR_UPN", {})
    ok = bool(v) and v["full_stromal"] and v["stromal_first"] >= v["n"] - LOO_MAX_FAIL
    say("[%s] KOLLER_EUR UPN (%s): full first %s; stromal state first in %s/%s leave-one-donor-out runs (required >= %s) -> %s" % (
        verdict_tag, design, v.get("full_first", "-"), v.get("stromal_first", "-"), v.get("n", "-"), (v["n"] - LOO_MAX_FAIL) if v else "-",
        "SUPPORTED" if ok else "NOT SUPPORTED"))
    for t, rx, name in (("LEIOMYOMA", r"^(uSMC|mPV|ePV)", "smooth muscle / perivascular"), ("HYPERTENSION", VASC_H, "vascular / perivascular")):
        f = B[(B["trait"] == t) & (B["left_out"] == "FULL")]
        if len(f):
            say("  [E_CTRL] %s (%s expected): %s" % (t, name, " | ".join("%s %s%s" % (r["family"], short(r["first"], 18), " (yes)" if re.search(rx, r["first"]) else " (no)")
                                                                        for _, r in f.iterrows())))
    return ok, out


def part_e2(f, O, sym2ent):
    say("")
    say("-" * 72)
    say("[SECTION] E2 HECA labels pooled over status, donor-balanced cell sample, Wilcoxon within compartment; leave-one-donor-out (all donors)")
    DEh, uni = heca_sets(sym2ent)
    uset = set(uni)
    ents = [sym2ent[x] for x in uni]
    pos_of = {g: i for i, g in enumerate(uni)}
    DE, CP = part_e2_de(f, O, sym2ent)
    st = CP[CP["label"].astype(str).str.contains(STROMA_H)]
    say("[E2] stromal types in the balanced sample: %s" % "; ".join("%s %d donors / %d cells (top %.0f%%)" % (r["label"], r["donors"], r["cells"], 100 * r["top_share"])
                                                                  for _, r in st.iterrows()))
    CP.to_csv(OUT / "e2_balanced_composition.csv", index=False)
    DE["gene"] = DE["gene"].astype(str)
    cols, meta = {}, {}
    iters = ["FULL"] + sorted(x for x in DE["iteration"].unique() if x != "FULL")
    for itn in iters:
        for grp, d in DE[DE["iteration"] == itn].groupby("group", sort=False):
            for fam, gl in fam_sets(d, uset, FAMS_E).items():
                k = "E%05d" % (len(cols) + 1)
                v = np.zeros(len(uni), dtype=np.int8)
                v[[pos_of[x] for x in gl]] = 1
                cols[k], meta[k] = v, (itn, fam, str(grp))
    res = {}
    drop_covar_files("e2")
    for t in TRAITS_E:
        raw = raw_of(t)
        if raw is not None:
            r = covar_runs(raw, cols, ents, "e2")
            if r:
                res[t] = r
    drop_covar_files("e2")
    rows = [x for t, r in res.items() for fam in FAMS_E for itn in iters for x in [rank_rows(r, meta, itn, fam, t)] if x]
    B = pd.DataFrame(rows)
    B.to_csv(OUT / "e2_balanced_pooled_loo.csv", index=False)
    if not len(B):
        say("[E2] no MAGMA results")
        return
    ok, out = report_design(B, "donor-balanced, pooled status", "H_HECA_BAL_VERDICT")
    SUM["results"]["e2"] = {"supported": ok, "by": out}


def part_e3(sym2ent):
    say("")
    say("-" * 72)
    say("[SECTION] E3 HECA donor-paired pseudobulk DE with labels pooled over status (v054 pseudobulk cache; >= %d donors with >= %d cells); leave-one-donor-out" % (
        E3_MIN_DONORS, MIN_CELLS_PB))
    f_pb, f_pbg, f_gen = CACHE54 / "heca_pb_counts.npy", CACHE54 / "heca_pb_groups.csv", CACHE54 / "heca_pb_genes.csv"
    if not (f_pb.exists() and f_pbg.exists() and f_gen.exists()):
        say("[E3] v054 pseudobulk cache missing; skipped")
        return
    PB, PBG = np.load(f_pb), pd.read_csv(f_pbg)
    PBG["donor"] = PBG["donor"].astype(str)
    genes_sym = np.array(pd.read_csv(f_gen)["gene"].astype(str))
    DEh, uni = heca_sets(sym2ent)
    uset = set(uni)
    ents = [sym2ent[x] for x in uni]
    pos_of = {g: i for i, g in enumerate(uni)}
    donors = sorted(PBG["donor"].unique())
    cols, meta, nd = {}, {}, {}
    for itn in ["FULL"] + donors:
        sub = PBG if itn == "FULL" else PBG[PBG["donor"] != itn]
        D = pb_de(PB, sub, genes_sym, E3_MIN_DONORS)
        for grp, (d, n) in D.items():
            if itn == "FULL":
                nd[grp] = n
            gl = [x for x in d.sort_values("score", ascending=False)["gene"].head(N_MATCH).astype(str) if x in uset]
            if len(gl) >= MIN_SET:
                k = "P%05d" % (len(cols) + 1)
                v = np.zeros(len(uni), dtype=np.int8)
                v[[pos_of[x] for x in gl]] = 1
                cols[k], meta[k] = v, (itn, "UPN", grp)
    say("[E3] %d donors; eligible types in the full data %d; stromal: %s" % (
        len(donors), len(nd), ", ".join("%s n=%d" % (k, v) for k, v in nd.items() if re.search(STROMA_H, k))))
    res = {}
    drop_covar_files("e3")
    for t in TRAITS_E:
        raw = raw_of(t)
        if raw is not None:
            r = covar_runs(raw, cols, ents, "e3")
            if r:
                res[t] = r
    drop_covar_files("e3")
    rows = [x for t, r in res.items() for itn in ["FULL"] + donors for x in [rank_rows(r, meta, itn, "UPN", t)] if x]
    B = pd.DataFrame(rows)
    B.to_csv(OUT / "e3_pooled_pseudobulk_loo.csv", index=False)
    if not len(B):
        say("[E3] no MAGMA results")
        return
    ok, out = report_design(B[B["family"] == "UPN"], "donor-paired pseudobulk, pooled status", "H_HECA_PB2_VERDICT")
    SUM["results"]["e3"] = {"supported": ok, "by": out}


# ----------------------------------------------------------------------------
def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.55 | HECA donor structure (composition, donor-balanced pooled labels, pooled pseudobulk) and LDSC rg | " + TS)
    say("=" * 72)
    say("[POLICY] rules fixed in CHANGELOG v0.2.45 before this run")
    t0 = time.time()
    ent2sym, sym2ent, pos = gene_loc()
    f, O = None, None
    try:
        if HECA_FILE.exists():
            f, O = heca_open()
    except Exception as e:
        say("[STEP_ERROR] HECA open %s: %s" % (type(e).__name__, str(e)[:300]))
    steps = [("E1", lambda: part_e1(O) if O is not None else say("[E1] HECA missing; skipped")), ("E3", lambda: part_e3(sym2ent)),
             ("E2", lambda: part_e2(f, O, sym2ent) if O is not None else say("[E2] HECA missing; skipped")), ("A3", part_a3)]
    for name, fn in steps:
        try:
            fn()
        except Exception as e:
            ERR.append("%s: %s" % (name, e))
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1500:])
    if f is not None:
        f.close()
    with open(OUT / "p9_v055_summary.json", "w") as fh:
        json.dump(SUM, fh, indent=1, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V055_COMPLETE" + (" with %d step errors" % len(ERR) if ERR else ""))
    _fh.close()


if __name__ == "__main__":
    main()
