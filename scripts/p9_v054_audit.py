#!/usr/bin/env python3
# Paper 9 | v0.0.54 | final statistical audit for the manuscript (rules CHANGELOG v0.2.44)
#   D  benchmark intervals: Wilson CIs per scenario x family, paired bootstrap CIs and exact McNemar tests of family differences
#   A  trait contrasts with GWAS dependence: genome-wide gene-z correlation, a paired genome-rotation null of the BETA_STD
#      difference (1000 shifts, same shifts for every trait), LDSC cross-trait rg and intercepts
#   C  cross-atlas locus test: loci classified by GSE266265 removed from the HECA decidual set, loci classified by HECA
#      removed from the GSE266265 eF3 set (200 random locus subsets each)
#   B  donor stability: leave-one-donor-out for GSE266265 eF3 (donor-paired DE) and for HECA dStromal_mid/ENDO (cell-level
#      Wilcoxon on the fixed v041 cell sample); HECA donor-paired pseudobulk DE as a donor-level alternative
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, LDSC 2 threads, nice 10). ~1-2 h; no downloads; HECA step cached.

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
OUT = RES / "review" / ("p9_v054_" + TS)
WORK = INTERIM / "magma_v054"
CACHE = INTERIM / "v054_cache"
LOG = REPO / "logs" / ("p9_v054_audit_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
ENV2 = dict(os.environ, OPENBLAS_NUM_THREADS="2", OMP_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
AUTOS = [str(c) for c in range(1, 23)]

# Fixed (CHANGELOG v0.2.44)
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
# Locus test (as v042 / v044)
GENE_WIN_BP, DROP_WIN_BP, PCT_ASSIGN = 250_000, 500_000, 0.90

for d in (OUT, WORK, CACHE, LOG.parent):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V054_AUDIT", "timestamp_utc": TS, "results": {}}
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
# D. Benchmark intervals
# ----------------------------------------------------------------------------
def part_d():
    say("")
    say("-" * 72)
    say("[SECTION] D  benchmark intervals (v041 replicates; Wilson 95% CI; paired bootstrap of family differences)")
    if not STATE.exists():
        say("[D] %s missing; skipped" % STATE.relative_to(REPO))
        return
    R = pd.read_csv(STATE)
    for c in ("any_q05", "top1_exact", "top3_exact", "top1_group", "related_miss", "true_q05"):
        if c in R:
            R[c] = R[c].astype(str).str.lower().isin(["true", "1", "1.0"])
    R = R.drop_duplicates(subset=["scenario", "rep", "family"], keep="last")
    nul = R[R["scenario"] == "NULLROT"]
    cal = {}
    say("[D] family-wise false-positive rate (any BH q < 0.05 within a family) under %d genome rotations of the real KOLLER_EUR z:" % nul["rep"].nunique())
    for fam, g in nul.groupby("family", sort=False):
        k, n = int(g["any_q05"].sum()), len(g)
        lo, hi = clopper(k, n)
        cal[fam] = {"k": k, "n": n, "rate": k / n, "cp_lo": lo, "cp_hi": hi}
        say("    %-6s %3d/%d = %.3f (Clopper-Pearson 95%% CI %.3f-%.3f) -> %s" % (
            fam, k, n, k / n, lo, hi, "upper bound below 0.05" if hi < 0.05 else ("CI includes 0.05" if lo <= 0.05 else "above 0.05")))
    T = R[R["scenario"] != "NULLROT"]
    rows = []
    for (sc, fam), g in T.groupby(["scenario", "family"], sort=False):
        n = len(g)
        for m in ("top1_exact", "top1_group", "true_q05"):
            k = int(g[m].sum())
            lo, hi = wilson(k, n)
            rows.append({"scenario": sc, "family": fam, "metric": m, "k": k, "n": n, "rate": k / n, "wilson_lo": lo, "wilson_hi": hi})
    W = pd.DataFrame(rows)
    W.to_csv(OUT / "d_benchmark_rates_ci.csv", index=False)
    rng = np.random.default_rng(SEED)
    comps = [("UPN", "ABSN"), ("UP", "ABS"), ("SPEC", "ABSN"), ("UPN", "SPEC"), ("DOWNN", "UPN"), ("ABSN", "UPN")]
    drows = []
    for sc, g in T.groupby("scenario", sort=False):
        for m in ("top1_exact", "top1_group"):
            P = g.pivot_table(index="rep", columns="family", values=m, aggfunc="first")
            for a, b in comps:
                if a not in P or b not in P:
                    continue
                if sc.startswith("NEG") != (a in ("DOWNN", "ABSN") and b == "UPN"):
                    continue
                x = P[[a, b]].dropna().astype(float)
                if not len(x):
                    continue
                d = (x[a] - x[b]).to_numpy()
                bs = d[rng.integers(0, len(d), size=(N_BOOT, len(d)))].mean(axis=1)
                b10 = int(((x[a] == 1) & (x[b] == 0)).sum())
                b01 = int(((x[a] == 0) & (x[b] == 1)).sum())
                drows.append({"scenario": sc, "metric": m, "family_a": a, "family_b": b, "n": len(d), "rate_a": x[a].mean(), "rate_b": x[b].mean(),
                              "diff": d.mean(), "boot_lo": float(np.percentile(bs, 2.5)), "boot_hi": float(np.percentile(bs, 97.5)),
                              "a_only": b10, "b_only": b01, "mcnemar_p": mcnemar_exact(b10, b01)})
    Dd = pd.DataFrame(drows)
    Dd.to_csv(OUT / "d_benchmark_paired_differences.csv", index=False)
    say("[D] exact top-1 recovery (rate, Wilson 95% CI) and paired differences (bootstrap 95% CI; exact McNemar p):")
    for sc in T["scenario"].unique():
        w = W[(W["scenario"] == sc) & (W["metric"] == "top1_exact")].set_index("family")
        say("  %-18s %s" % (sc, " | ".join("%s %.2f (%.2f-%.2f)" % (f, w.loc[f, "rate"], w.loc[f, "wilson_lo"], w.loc[f, "wilson_hi"])
                                           for f in ("ABS", "UP", "DOWN", "ABSN", "UPN", "DOWNN", "SPEC") if f in w.index)))
        dd = Dd[(Dd["scenario"] == sc) & (Dd["metric"] == "top1_exact")]
        if len(dd):
            say("  %-18s %s" % ("", " | ".join("%s-%s %+.2f (%+.2f to %+.2f; p %.3f)" % (r["family_a"], r["family_b"], r["diff"], r["boot_lo"], r["boot_hi"], r["mcnemar_p"])
                                              for _, r in dd.iterrows())))
    SUM["results"]["d"] = {"fwer": cal, "n_truth_reps": int(T.groupby("scenario")["rep"].nunique().min()) if len(T) else 0}


# ----------------------------------------------------------------------------
# A. Trait contrasts that respect GWAS dependence
# ----------------------------------------------------------------------------
def read_gene_z(t):
    f = MAGMA_OUT / (t + "_35_10.genes.out")
    if not f.exists():
        return None
    g = pd.read_csv(f, sep=r"\s+", dtype={"GENE": str})
    return g.set_index("GENE")["ZSTAT"].astype(float)


def part_a_genez():
    say("")
    say("-" * 72)
    say("[SECTION] A1 genome-wide correlation of MAGMA gene z between GWAS (all genes; includes shared genetics and sample overlap)")
    Z = {t: read_gene_z(t) for t in set(RG_TRAITS)}
    rows = []
    for a, b in ENDO_PAIRS + PAIRS:
        if Z.get(a) is None or Z.get(b) is None:
            continue
        j = pd.concat([Z[a], Z[b]], axis=1, join="inner").dropna()
        if len(j) < 100:
            continue
        r = float(np.corrcoef(j.iloc[:, 0], j.iloc[:, 1])[0, 1])
        rs = float(j.rank().corr().iloc[0, 1])
        rows.append({"trait_a": a, "trait_b": b, "n_genes": len(j), "pearson": r, "spearman": rs})
        say("    %-13s %-13s genes %5d  r %.3f  rho %.3f" % (a[:13], b[:13], len(j), r, rs))
    pd.DataFrame(rows).to_csv(OUT / "a1_gene_z_correlation.csv", index=False)
    SUM["results"]["a1"] = rows


def rot_columns(x_sorted, n_rot, rng):
    n = len(x_sorted)
    lo, hi = max(1, int(ROT_MIN_FRAC * n)), max(2, int((1 - ROT_MIN_FRAC) * n))
    offs = rng.integers(lo, hi, size=n_rot)
    cols = {"OBS": x_sorted.astype(np.int8)}
    for i, o in enumerate(offs):
        cols["R%04d" % (i + 1)] = np.roll(x_sorted, int(o)).astype(np.int8)
    return cols, offs


def part_a_rotation(sym2ent, pos):
    say("")
    say("-" * 72)
    say("[SECTION] A2 trait contrasts with a paired rotation null (%d genome shifts per set, identical for every trait)" % N_ROT)
    DEh, uni_h = heca_sets(sym2ent)
    DEg, uni_g = gse_sets(sym2ent)
    rng = np.random.default_rng(SEED)
    est = {}
    rows = []
    for atlas, fam, lab in ROT_SETS:
        DE, uni = (DEh, uni_h) if atlas == "HECA" else (DEg, uni_g)
        d = DE[DE["group"] == lab]
        if not len(d):
            say("  [A2] %s %s %s: group not found; skipped" % (atlas, fam, lab))
            continue
        fs = fam_sets(d, set(uni), (fam,))
        if fam not in fs:
            say("  [A2] %s %s %s: set < %d genes; skipped" % (atlas, fam, lab, MIN_SET))
            continue
        ents = [sym2ent[g] for g in uni]
        order = sorted(range(len(ents)), key=lambda i: (pos[ents[i]][0], pos[ents[i]][1]))
        ents_sorted = [ents[i] for i in order]
        member = set(sym2ent[g] for g in fs[fam])
        x = np.array([1 if e in member else 0 for e in ents_sorted], dtype=np.int8)
        cols, offs = rot_columns(x, N_ROT, rng)
        tag = "rot_%s_%s_%s" % (atlas, fam, re.sub(r"[^A-Za-z0-9]", "", lab))
        drop_covar_files(tag)
        t0 = time.time()
        for t in TRAITS_C:
            raw = raw_of(t)
            if raw is None:
                continue
            r = covar_runs(raw, cols, ents_sorted, tag, chunk=N_ROT + 1)
            if not r or "OBS" not in r:
                continue
            rb = np.array([r.get("R%04d" % (i + 1), {}).get("beta_std", np.nan) for i in range(N_ROT)])
            est[(atlas, fam, lab, t)] = (r["OBS"], rb)
        drop_covar_files(tag)
        say("  [A2] %s %-4s %-20s %d genes in set, %d-gene universe; %d traits (%.1f min)" % (
            atlas, fam, short(lab, 20), int(x.sum()), len(x), sum(1 for k in est if k[:3] == (atlas, fam, lab)), (time.time() - t0) / 60))
        for t in TRAITS_C:
            if (atlas, fam, lab, t) in est:
                o, rb = est[(atlas, fam, lab, t)]
                sd = float(np.nanstd(rb, ddof=1))
                say("       %-13s BETA_STD %.4f (MAGMA se %.4f, p %.1e) | rotation mean %.4f sd %.4f (sd / MAGMA se %.2f; empirical p %.3f)" % (
                    t[:13], o["beta_std"], se_std(o), o["p"], float(np.nanmean(rb)), sd, sd / se_std(o) if se_std(o) > 0 else np.nan,
                    (1 + int(np.sum(rb >= o["beta_std"]))) / (1 + int(np.isfinite(rb).sum()))))
        for a, b in PAIRS:
            if (atlas, fam, lab, a) not in est or (atlas, fam, lab, b) not in est:
                continue
            (oa, ra), (ob, rb_) = est[(atlas, fam, lab, a)], est[(atlas, fam, lab, b)]
            ok = np.isfinite(ra) & np.isfinite(rb_)
            r_ab = float(np.corrcoef(ra[ok], rb_[ok])[0, 1]) if ok.sum() > 10 else np.nan
            sa, sb = se_std(oa), se_std(ob)
            diff = oa["beta_std"] - ob["beta_std"]
            se_ind = math.sqrt(sa ** 2 + sb ** 2)
            var_cov = sa ** 2 + sb ** 2 - 2 * r_ab * sa * sb if np.isfinite(r_ab) else np.nan
            se_cov = math.sqrt(var_cov) if np.isfinite(var_cov) and var_cov > 0 else np.nan
            dr = (ra - rb_)[ok]
            z_rot = (diff - float(dr.mean())) / float(dr.std(ddof=1)) if len(dr) > 10 else np.nan
            p_emp = (1 + int(np.sum(dr >= diff))) / (1 + len(dr))
            p_ind, p_cov = phi_sf(diff / se_ind), phi_sf(diff / se_cov) if np.isfinite(se_cov) else np.nan
            rows.append({"atlas": atlas, "family": fam, "label": lab, "trait_a": a, "trait_b": b, "beta_std_a": oa["beta_std"], "beta_std_b": ob["beta_std"],
                         "se_a": sa, "se_b": sb, "diff": diff, "r_rot": r_ab, "se_independent": se_ind, "p_independent": p_ind,
                         "se_covariance": se_cov, "p_covariance": p_cov, "sd_rot_diff": float(dr.std(ddof=1)), "z_rot": z_rot, "p_rot_empirical": p_emp,
                         "p_conservative": np.nanmax([p_cov, p_emp])})
    C = pd.DataFrame(rows)
    C.to_csv(OUT / "a2_trait_contrasts_rotation.csv", index=False)
    if not len(C):
        say("[A2] no contrasts computed")
        return
    say("[A2] one-sided p for 'endometriosis > control': independence (v049 / v050) | covariance from the rotations | empirical rotation:")
    for (atlas, fam, lab), g in C.groupby(["atlas", "family", "label"], sort=False):
        say("  %-4s %-4s %-20s" % (atlas, fam, short(lab, 20)))
        for _, r in g.iterrows():
            say("       %-6s>%-6s diff %+.4f | r_rot %+.2f | p_ind %.3f | p_cov %.3f | p_rot %.3f" % (
                r["trait_a"][:6], r["trait_b"][:6], r["diff"], r["r_rot"], r["p_independent"], r["p_covariance"], r["p_rot_empirical"]))
    verdict = {}
    for atlas, fam, lab in PRIMARY_C:
        g = C[(C["atlas"] == atlas) & (C["family"] == fam) & (C["label"] == lab) & (C["trait_a"] == "KOLLER_EUR") & C["trait_b"].isin(["T2D", "HYPERTENSION"])]
        ok = len(g) == 2 and bool((g["p_conservative"] < 0.05).all())
        verdict["%s_%s_%s" % (atlas, fam, lab)] = {"supported": ok, "p_conservative": dict(zip(g["trait_b"], g["p_conservative"]))}
        say("[H_COV_CONTRAST_VERDICT] %s %s %s: KOLLER_EUR > T2D and > hypertension with covariance-aware AND empirical rotation one-sided p < 0.05 -> %s (%s)" % (
            atlas, fam, lab, "SUPPORTED" if ok else "NOT SUPPORTED",
            ", ".join("%s p %.3f" % (r["trait_b"][:6], r["p_conservative"]) for _, r in g.iterrows())))
    SUM["results"]["a2"] = verdict


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


def part_a_rg():
    say("")
    say("-" * 72)
    say("[SECTION] A3 LDSC cross-trait genetic correlation and intercept (HapMap3 weights as reference and regression LD scores)")
    w = prefix_of("1000G_Phase3_weights_hm3_no_MHC")
    if not w or not LDSC.exists() or not LDSC_PY.exists():
        say("[A3] LDSC or weights missing; skipped")
        return
    avail = [t for t in RG_TRAITS if (MUNGED / (t + ".sumstats.gz")).exists()]
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
        if not (log.exists() and "Summary of Genetic Correlation Results" in log.read_text()):
            cmd = NICE + [LDSC_PY, LDSC, "--rg", ",".join(str(MUNGED / (x + ".sumstats.gz")) for x in [t] + others),
                          "--ref-ld-chr", w, "--w-ld-chr", w, "--out", out]
            rc, o, e = sh(cmd, timeout=6 * 3600, env=ENV2)
            if not log.exists():
                say("  [A3] %s rc=%d %s" % (t, rc, (o + e).strip()[-300:]))
                continue
        for d in parse_rg_log(log):
            a = name_of.get(d.get("p1"), Path(d.get("p1", "")).name.replace(".sumstats.gz", ""))
            b = name_of.get(d.get("p2"), Path(d.get("p2", "")).name.replace(".sumstats.gz", ""))
            def f(k):
                try:
                    return float(d.get(k, "nan"))
                except ValueError:
                    return np.nan
            rows.append({"trait_a": a, "trait_b": b, "rg": f("rg"), "rg_se": f("se"), "rg_p": f("p"), "gcov_int": f("gcov_int"), "gcov_int_se": f("gcov_int_se"),
                         "h2_obs_b": f("h2_obs"), "h2_int_b": f("h2_int")})
    if not rows:
        say("[A3] no LDSC rg results")
        return
    G = pd.DataFrame(rows)
    G.to_csv(OUT / "a3_ldsc_rg.csv", index=False)
    for _, r in G.iterrows():
        say("    %-13s %-13s rg %+.3f (se %.3f, p %.1e) | cross-trait intercept %+.3f (se %.3f)" % (
            r["trait_a"][:13], r["trait_b"][:13], r["rg"], r["rg_se"], r["rg_p"], r["gcov_int"], r["gcov_int_se"]))
    say("[A3] %d pairs (%.0f min); the intercept reflects sample overlap / shared confounding, rg the shared genetics" % (len(G), (time.time() - t0) / 60))
    SUM["results"]["a3"] = G.to_dict("records")


# ----------------------------------------------------------------------------
# C. Cross-atlas locus test
# ----------------------------------------------------------------------------
def gse_ratio():
    C = pd.read_csv(GSE_DE / "mean_cpm_by_type.csv.gz", index_col=0)
    C.index = C.index.astype(str)
    C = C.groupby(level=0).mean()
    C = C[(C >= MIN_CPM).any(axis=1)]
    S = C.div(C.sum(axis=1), axis=0)
    comp = {c: first_match(c, COMP_G, "other") for c in S.columns}
    names = [n for n, _ in COMP_G]
    share = pd.DataFrame({n: S[[c for c in S.columns if comp[c] == n]].sum(axis=1) for n in names if any(comp[c] == n for c in S.columns)}, index=S.index)
    return comp, share.rank(pct=True), S


def assign(row):
    r = row.dropna()
    if not len(r):
        return "no expression", np.nan
    b = r.idxmax()
    return (b if r[b] >= PCT_ASSIGN else "broad"), float(r[b])


def part_c(ent2sym, sym2ent, pos):
    say("")
    say("-" * 72)
    say("[SECTION] C  cross-atlas locus test (classification and enrichment from different atlases; %d random subsets each)" % N_PERM)
    lf = latest(LOCUS_DIR, "p9_v042_*/l1_loci_compartments.csv")
    go = MAGMA_OUT / "KOLLER_EUR_35_10.genes.out"
    mf = GS_P / "manifest.csv"
    raw = raw_of("KOLLER_EUR")
    if lf is None or not go.exists() or not mf.exists() or raw is None:
        say("[C] inputs missing (v042 loci / genes.out / v031 sets); skipped")
        return
    L = pd.read_csv(lf)
    K = L[L["trait"] == "KOLLER_EUR"].reset_index(drop=True)
    g = pd.read_csv(go, sep=r"\s+", dtype={"GENE": str, "CHR": str})
    for c in ("CHR", "START", "STOP", "P"):
        g[c] = pd.to_numeric(g[c], errors="coerce")
    g = g.dropna(subset=["CHR", "START", "STOP", "P"])
    g["CHR"] = g["CHR"].astype(int)
    comp, ratio, S = gse_ratio()
    say("[C] GSE266265 compartments (34 eligible types, donor-mean CPM): " + "; ".join(
        "%s %s" % (n, ",".join(short(c, 10) for c in S.columns if comp[c] == n)) for n, _ in COMP_G) +
        "; not assigned: " + ",".join(c for c in S.columns if comp[c] == "other"))
    gse_c, gse_gene = [], []
    for _, r in K.iterrows():
        m = g[(g["CHR"] == int(r["chr"])) & (g["STOP"] >= r["bp"] - GENE_WIN_BP) & (g["START"] <= r["bp"] + GENE_WIN_BP)].sort_values("P")
        syms = [ent2sym.get(str(e), str(e)) for e in m["GENE"]]
        expr = [s for s in syms if s in ratio.index]
        c_, _ = assign(ratio.loc[expr[0]]) if expr else ("no expression", np.nan)
        gse_c.append(c_)
        gse_gene.append(expr[0] if expr else "")
    K["compartment_gse"] = gse_c
    K["gene_gse"] = gse_gene
    K.to_csv(OUT / "c_loci_two_atlases.csv", index=False)
    X = pd.crosstab(K["compartment"], K["compartment_gse"])
    X.to_csv(OUT / "c_crosstab.csv")
    say("[C] KOLLER_EUR loci %d; HECA (rows) x GSE266265 (columns):" % len(K))
    for line in X.to_string().splitlines():
        say("    " + line)
    po = float((K["compartment"] == K["compartment_gse"]).mean())
    pe = float(sum(K["compartment"].value_counts(normalize=True).get(c, 0) * K["compartment_gse"].value_counts(normalize=True).get(c, 0)
                   for c in set(K["compartment"]) | set(K["compartment_gse"])))
    kappa = (po - pe) / (1 - pe) if pe < 1 else np.nan
    sh_ = K[K["compartment"] == "stroma"]
    say("[C] agreement %.2f, Cohen's kappa %.2f; HECA stromal loci also stromal in GSE266265: %d/%d (%s)" % (
        po, kappa, int((sh_["compartment_gse"] == "stroma").sum()), len(sh_),
        ", ".join("%s->%s" % (r["prioritized_gene"], r["compartment_gse"]) for _, r in sh_.iterrows())))
    say("[C] GSE266265 stromal loci: %s" % ", ".join("%s (%s)" % (r["gene_gse"], r["compartment"]) for _, r in K[K["compartment_gse"] == "stroma"].iterrows()))

    def genes_of(ids):
        d = set()
        for _, r in K[K["lead_snp"].isin(ids)].iterrows():
            d |= {e for e, (c, s, en) in pos.items() if c == int(r["chr"]) and en >= r["bp"] - DROP_WIN_BP and s <= r["bp"] + DROP_WIN_BP}
        return d

    def test(focal, allgenes, stroma_ids, other_stroma_ids, tag, seed):
        rng = np.random.default_rng(seed)
        others = [x for x in K["lead_snp"] if x not in stroma_ids]
        if len(others) < len(stroma_ids):
            say("  [C] %s: %d stromal loci but only %d other loci; not evaluable" % (tag, len(stroma_ids), len(others)))
            return None
        versions = {"BASE": focal, "DROPA": [x for x in focal if x not in genes_of(stroma_ids)],
                    "DROPB": [x for x in focal if x not in genes_of(other_stroma_ids)]}
        for i in range(N_PERM):
            pick = list(rng.choice(others, size=min(len(stroma_ids), len(others)), replace=False))
            dg = genes_of(pick)
            versions["R%03d" % (i + 1)] = [x for x in focal if x not in dg]
        sf = WORK / ("c_%s.txt" % tag)
        with open(sf, "w") as o:
            o.write("ALLGENES " + " ".join(allgenes) + "\n")
            for k, gl in versions.items():
                o.write(k + " " + " ".join(gl) + "\n")
        r = magma_sets(raw, sf, "c_" + tag)
        if not r or "BASE" not in r:
            return None
        base = r["BASE"]["p"]
        obs = math.log10(r["DROPA"]["p"]) - math.log10(base)
        null = np.array([math.log10(r[k]["p"]) - math.log10(base) for k in versions if k.startswith("R") and k in r])
        pp = (1 + int((null >= obs).sum())) / (1 + len(null))
        return {"p_base": base, "p_drop": r["DROPA"]["p"], "p_drop_other_atlas": r["DROPB"]["p"], "delta": obs, "null_median": float(np.median(null)),
                "null_p95": float(np.percentile(null, 95)), "n_null": len(null), "p_perm": pp, "n_loci": len(stroma_ids)}

    out = {}
    gse_str = list(K.loc[K["compartment_gse"] == "stroma", "lead_snp"])
    heca_str = list(K.loc[K["compartment"] == "stroma", "lead_snp"])
    man = pd.read_csv(mf)
    key = next((k for k, l, s in zip(man["set_key"], man["label"], man["status"]) if k.startswith("UP_") and l == "dStromal_mid" and s == "ENDO"), None)
    if key and len(gse_str) >= 3:
        focal = [sym2ent[x] for x in open(GS_P / (key + ".txt")).read().split() if x in sym2ent]
        allg = [sym2ent[x] for x in open(GS_P / "ALLGENES.txt").read().split() if x in sym2ent] if (GS_P / "ALLGENES.txt").exists() else \
            sorted({e for k2 in man["set_key"] for e in [sym2ent[x] for x in open(GS_P / (k2 + ".txt")).read().split() if x in sym2ent]})
        res = test(focal, allg, gse_str, heca_str, "heca_dsmid", SEED)
        if res:
            out["a"] = res
            say("  [C-a] HECA UP dStromal_mid/ENDO (v031 sets) without the %d GSE266265-stromal loci: p %.1e -> %.1e (%+.2f log10); random %d-locus subsets median %+.2f, 95th %+.2f -> perm p %.3f | without the 9 HECA-stromal loci (v044 reference) %.1e" % (
                res["n_loci"], res["p_base"], res["p_drop"], res["delta"], res["n_loci"], res["null_median"], res["null_p95"], res["p_perm"], res["p_drop_other_atlas"]))
    else:
        say("  [C-a] not evaluable (set key %s; GSE266265 stromal loci %d)" % (key, len(gse_str)))
    DEg, uni_g = gse_sets(sym2ent)
    d = DEg[DEg["group"] == "eF3"]
    fs = fam_sets(d, set(uni_g), ("UPN",)) if len(d) else {}
    if "UPN" in fs and len(heca_str) >= 3:
        res = test([sym2ent[x] for x in fs["UPN"]], [sym2ent[x] for x in uni_g], heca_str, gse_str, "gse_ef3", SEED + 1)
        if res:
            out["b"] = res
            say("  [C-b] GSE266265 UPN eF3 without the %d HECA-stromal loci: p %.1e -> %.1e (%+.2f log10); random %d-locus subsets median %+.2f, 95th %+.2f -> perm p %.3f | without the GSE266265-stromal loci %.1e" % (
                res["n_loci"], res["p_base"], res["p_drop"], res["delta"], res["n_loci"], res["null_median"], res["null_p95"], res["p_perm"], res["p_drop_other_atlas"]))
    else:
        say("  [C-b] not evaluable")
    pd.DataFrame([{"test": k, **v} for k, v in out.items()]).to_csv(OUT / "c_cross_atlas_removal.csv", index=False)
    na = sum(1 for v in out.values() if v["p_perm"] < 0.05)
    verdict = "SUPPORTED" if na == 2 else ("PARTIAL" if na == 1 else "NOT SUPPORTED")
    say("[H_LOCUS_X_VERDICT] stromal loci from one atlas carry the stromal set of the other more than random loci (both directions, perm p < 0.05): %d/%d -> %s" % (
        na, len(out), verdict if len(out) == 2 else verdict + " (only %d direction evaluable)" % len(out)))
    SUM["results"]["c"] = {"kappa": kappa, "agreement": po, "n_gse_stroma": len(gse_str), "tests": out, "verdict": verdict}


# ----------------------------------------------------------------------------
# B1. GSE266265 leave-one-donor-out (donor-paired DE as v050)
# ----------------------------------------------------------------------------
def gse_prep(G, pb, tot):
    # per eligible type: sums over donors of D = log2(CPM_type + 1) - log2(CPM_rest + 1), D^2, CPM_type, CPM_rest (all genes)
    ok = G[G["n_cells"] >= MIN_CELLS_DONOR]
    stats = {}
    for ty, rows in ok.groupby("refined"):
        X = np.asarray(pb[rows.index.to_numpy(), :], dtype=np.float64)
        Rr = np.vstack([tot[s] for s in rows["sample"]]) - X
        Ct = X / np.maximum(X.sum(axis=1, keepdims=True), 1) * 1e6
        Cr = Rr / np.maximum(Rr.sum(axis=1, keepdims=True), 1) * 1e6
        D = np.log2(Ct + 1) - np.log2(Cr + 1)
        stats[ty] = {"n": len(rows), "s1": D.sum(axis=0), "s2": (D * D).sum(axis=0), "sct": Ct.sum(axis=0), "scr": Cr.sum(axis=0),
                     "row": dict(zip(rows["sample"], rows.index.to_numpy()))}
    return stats


def gse_de_scores(stats, genes, pb, tot, drop=None):
    out = {}
    for ty, st in stats.items():
        n, s1, s2, sct, scr = st["n"], st["s1"], st["s2"], st["sct"], st["scr"]
        if drop is not None and drop in st["row"]:
            x = np.asarray(pb[st["row"][drop], :], dtype=np.float64)
            r = tot[drop] - x
            ct = x / max(x.sum(), 1) * 1e6
            cr = r / max(r.sum(), 1) * 1e6
            dd = np.log2(ct + 1) - np.log2(cr + 1)
            n, s1, s2, sct, scr = n - 1, s1 - dd, s2 - dd * dd, sct - ct, scr - cr
        if n < MIN_DONORS:
            continue
        expr = (sct / n >= MIN_CPM) | (scr / n >= MIN_CPM)
        m = s1 / n
        sd = np.sqrt(np.maximum((s2 - s1 * s1 / n) / (n - 1), 0.0))
        tt = m / np.maximum(sd / math.sqrt(n), 1e-9)
        out[clean(ty)] = pd.DataFrame({"gene": genes[expr], "score": tt[expr]})
    return out


def rank_items(res, keys):
    it = sorted([(lab, res[k]) for k, lab in keys if k in res], key=lambda x: x[1]["p"])
    return it


def part_b1(sym2ent):
    say("")
    say("-" * 72)
    say("[SECTION] B1 GSE266265 leave-one-donor-out (donor-paired DE recomputed without each donor; UPN, N = %d; MAGMA)" % N_MATCH)
    if not (PB266 / "pseudobulk_counts.npy").exists():
        say("[B1] GSE266265 pseudobulk missing; skipped")
        return
    G = pd.read_csv(PB266 / "pseudobulk_groups.csv")
    genes = np.array(pd.read_csv(PB266 / "genes.csv")["symbol"].astype(str))
    pb = np.load(PB266 / "pseudobulk_counts.npy", mmap_mode="r")
    tot = {s: np.asarray(pb[np.sort(np.asarray(idx)), :], dtype=np.float64).sum(axis=0) for s, idx in G.groupby("sample").groups.items()}
    DEc, uni = gse_sets(sym2ent)
    uset = set(uni)
    ents = [sym2ent[x] for x in uni]
    pos_of = {g: i for i, g in enumerate(uni)}
    donors = sorted(G["sample"].unique(), key=str)
    cols, meta, sets_ef3 = {}, {}, {}
    t0 = time.time()
    check = None
    stats = gse_prep(G, pb, tot)
    for it, drop in enumerate([None] + donors):
        S = gse_de_scores(stats, genes, pb, tot, drop)
        for j, (ty, d) in enumerate(sorted(S.items())):
            gl = [x for x in d.sort_values("score", ascending=False)["gene"].head(N_MATCH).astype(str) if x in uset]
            if len(gl) < MIN_SET:
                continue
            k = "L%03dT%02d" % (it, j)
            v = np.zeros(len(uni), dtype=np.int8)
            v[[pos_of[x] for x in gl]] = 1
            cols[k], meta[k] = v, ("FULL" if drop is None else str(drop), ty)
            if ty == "eF3":
                sets_ef3["FULL" if drop is None else str(drop)] = gl
        if drop is None and "FULL" in sets_ef3:
            c = DEc[DEc["group"] == "eF3"]
            ref = [x for x in c.sort_values("score", ascending=False)["gene"].head(N_MATCH).astype(str) if x in uset]
            check = jaccard(ref, sets_ef3["FULL"])
    say("[B1] %d donors; %d leave-one-out DE runs + full (%.1f min); full-data eF3 set vs v050 cache: Jaccard %.3f" % (
        len(donors), len(donors), (time.time() - t0) / 60, check if check is not None else float("nan")))
    res = {}
    drop_covar_files("b1")
    for t in TRAITS_LOO:
        raw = raw_of(t)
        if raw is None:
            continue
        r = covar_runs(raw, cols, ents, "b1")
        if r:
            res[t] = r
    drop_covar_files("b1")
    rows = []
    for t, r in res.items():
        for itn in ["FULL"] + [str(x) for x in donors]:
            keys = [(k, m[1]) for k, m in meta.items() if m[0] == itn]
            items = rank_items(r, keys)
            if not items:
                continue
            labs = [l for l, _ in items]
            grp = [first_match(l, GROUPS_G, "immune / other") for l in labs]
            ef3 = next((i for i, l in enumerate(labs) if l == "eF3"), None)
            bs = next((i for i, g_ in enumerate(grp) if g_ == "endometrial stroma"), None)
            bp = next((i for i, g_ in enumerate(grp) if g_ == "perivascular / MSC-like"), None)
            v = items[ef3][1] if ef3 is not None else {}
            rows.append({"trait": t, "left_out": itn, "first": labs[0], "first_group": grp[0], "first_p": items[0][1]["p"],
                         "rank_eF3": (ef3 + 1) if ef3 is not None else np.nan, "p_eF3": v.get("p", np.nan), "beta_std_eF3": v.get("beta_std", np.nan),
                         "se_eF3": se_std(v) if v else np.nan, "best_stroma_rank": (bs + 1) if bs is not None else np.nan,
                         "best_perivascular_rank": (bp + 1) if bp is not None else np.nan,
                         "jaccard_eF3_vs_full": jaccard(sets_ef3.get(itn, []), sets_ef3.get("FULL", []))})
    B = pd.DataFrame(rows)
    B.to_csv(OUT / "b1_gse266265_loo.csv", index=False)
    verdict = {}
    for t in TRAITS_LOO:
        g = B[(B["trait"] == t) & (B["left_out"] != "FULL")]
        f = B[(B["trait"] == t) & (B["left_out"] == "FULL")]
        if not len(g) or not len(f):
            continue
        f = f.iloc[0]
        n = len(g)
        n_str = int((g["first_group"] == "endometrial stroma").sum())
        n_ef3 = int((g["rank_eF3"] == 1).sum())
        n_pv = int((g["best_perivascular_rank"] < g["best_stroma_rank"]).sum())
        worst = g.sort_values("p_eF3", ascending=False).iloc[0]
        say("  %-13s full: first %s (p %.1e), eF3 rank %s | LOO (n=%d): stroma first %d, eF3 first %d, perivascular above stroma %d | eF3 BETA_STD %.4f [%.4f-%.4f], p up to %.1e (without %s; rank %s) | min Jaccard %.2f" % (
            t[:13], f["first"], f["first_p"], int(f["rank_eF3"]) if np.isfinite(f["rank_eF3"]) else "-", n, n_str, n_ef3, n_pv, f["beta_std_eF3"],
            g["beta_std_eF3"].min(), g["beta_std_eF3"].max(), worst["p_eF3"], worst["left_out"], int(worst["rank_eF3"]) if np.isfinite(worst["rank_eF3"]) else "-",
            g["jaccard_eF3_vs_full"].min()))
        others = g[g["first_group"] != "endometrial stroma"]
        if len(others):
            say("       first not stroma without: %s" % "; ".join("%s (%s)" % (r["left_out"], r["first"]) for _, r in others.head(10).iterrows()))
        verdict[t] = {"n": n, "stroma_first": n_str, "eF3_first": n_ef3, "perivascular_above": n_pv}
    ok = all(t in verdict and verdict[t]["stroma_first"] >= LOO_SHARE * verdict[t]["n"] for t in ("KOLLER_EUR", "GCST90205183"))
    say("[H_LOO_GSE_VERDICT] endometrial stroma first in >= %.0f%% of leave-one-donor-out runs for KOLLER_EUR and GCST90205183 (UPN) -> %s (%s)" % (
        100 * LOO_SHARE, "SUPPORTED" if ok else "NOT SUPPORTED",
        ", ".join("%s %d/%d" % (t[:6], v["stroma_first"], v["n"]) for t, v in verdict.items())))
    SUM["results"]["b1"] = {"supported": ok, "by_trait": verdict, "jaccard_full_vs_cache": check}


# ----------------------------------------------------------------------------
# B2 / B3. HECA leave-one-ENDO-donor-out (cell-level, fixed v041 sample) and donor-paired pseudobulk DE
# ----------------------------------------------------------------------------
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


def heca_compute(sym2ent):
    f_loo, f_pb, f_pbg, f_comp = CACHE / "heca_loo_de.csv.gz", CACHE / "heca_pb_counts.npy", CACHE / "heca_pb_groups.csv", CACHE / "heca_donor_composition.csv"
    f_gen = CACHE / "heca_pb_genes.csv"
    if f_loo.exists() and f_pb.exists() and f_pbg.exists() and f_gen.exists():
        say("[B2] cached HECA leave-one-out DE and pseudobulk (%s)" % CACHE.name)
        return (pd.read_csv(f_loo), np.load(f_pb), pd.read_csv(f_pbg), (pd.read_csv(f_comp) if f_comp.exists() else None),
                np.array(pd.read_csv(f_gen)["gene"].astype(str)))
    import h5py
    import scipy.sparse as sp
    t00 = time.time()
    f = h5py.File(HECA_FILE, "r")
    cols = [k for k in f["obs"].keys() if k not in ("_index", "__categories")]
    dataset_col, dataset_vals = find_dataset_col(f, cols)
    donor = obs_cat(f, DONOR_COL)
    ctype = obs_cat(f, TYPE_COL)
    path = obs_cat(f, PATHOLOGY_COL)
    stage = obs_cat(f, STAGE_COL)
    keep = pd.notna(donor) & pd.notna(ctype) & (stage != EXCLUDE_STAGE)
    if dataset_col:
        keep &= np.isin(obs_cat(f, dataset_col), dataset_vals)
    status = status_of(path)
    keep &= pd.notna(status)
    labels = sorted(pd.unique(ctype[keep]))
    comp = {l: ("immune" if IMMUNE_RE.search(l) and not NOT_IMMUNE_RE.search(l) else "endometrial") for l in labels}
    layer = f["layers"]["counts"]
    var = f["var"]
    gobj = var[dec(var.attrs.get("_index", "_index"))]
    if isinstance(gobj, h5py.Group):
        gobj = gobj["values"] if "values" in gobj else gobj[list(gobj.keys())[0]]
    genes = np.array([dec(x) for x in gobj[:]], dtype=object)
    gidx = np.where(np.array([g in sym2ent for g in genes]))[0]
    gnames = genes[gidx]
    indptr = layer["indptr"][:]
    # identical cell sample to v041 (same seed, same draw order over CTL and ENDO)
    rng = np.random.default_rng(SEED_HECA)
    plan = {}
    for s in ("CTL", "ENDO"):
        for cp in ("endometrial", "immune"):
            pick, ok_types = [], []
            for l in labels:
                if comp[l] != cp:
                    continue
                idx = np.where(keep & (status == s) & (ctype == l))[0]
                if len(idx) < MIN_CELLS_TYPE:
                    continue
                if len(idx) > CAP_CELLS_TYPE:
                    idx = np.sort(rng.choice(idx, CAP_CELLS_TYPE, replace=False))
                pick.append(idx)
                ok_types.append(l)
            if len(ok_types) >= 2:
                plan[(s, cp)] = (np.sort(np.concatenate(pick)), ok_types)
    endo_plans = {k: v for k, v in plan.items() if k[0] == "ENDO"}
    want = np.unique(np.concatenate([v[0] for v in endo_plans.values()])) if endo_plans else np.array([], dtype=int)
    # pseudobulk groups (status, donor, label) over all kept cells
    key = pd.Series(list(zip(status, donor, ctype)))
    gid = np.full(len(donor), -1, dtype=np.int64)
    kk = key[keep]
    codes, uniq = pd.factorize(kk)
    gid[np.where(keep)[0]] = codes
    PBG = pd.DataFrame(list(uniq), columns=["status", "donor", "label"])
    PBG["n_cells"] = np.bincount(codes, minlength=len(PBG))
    PB = np.zeros((len(PBG), len(gidx)), dtype=np.float32)
    say("[B2] HECA variant P: %d cells kept (CTL donors %d, ENDO donors %d); %d pseudobulk units; ENDO cell sample %d cells" % (
        int(keep.sum()), pd.Series(donor[keep & (status == "CTL")]).nunique(), pd.Series(donor[keep & (status == "ENDO")]).nunique(), len(PBG), len(want)))
    rowmat = {}
    step = 10000
    n = len(indptr) - 1
    for a in range(0, n, step):
        b = min(n, a + step)
        gch = gid[a:b]
        lo, hi = np.searchsorted(want, a), np.searchsorted(want, b)
        if hi <= lo and not (gch >= 0).any():
            continue
        s0, e0 = int(indptr[a]), int(indptr[b])
        idx = np.asarray(layer["indices"][s0:e0])
        val = np.asarray(layer["data"][s0:e0], dtype=np.float32)
        loc = indptr[a:b + 1] - s0
        for i in want[lo:hi]:
            r = i - a
            rowmat[int(i)] = (idx[loc[r]:loc[r + 1]].copy(), val[loc[r]:loc[r + 1]].copy())
        rows_ok = np.where(gch >= 0)[0]
        if len(rows_ok):
            M = sp.csr_matrix((val, idx, loc), shape=(b - a, len(genes)))
            ug, inv = np.unique(gch[rows_ok], return_inverse=True)
            Ag = sp.csr_matrix((np.ones(len(rows_ok), dtype=np.float32), (inv, rows_ok)), shape=(len(ug), b - a))
            PB[ug] += (Ag @ M)[:, gidx].toarray()
        del idx, val
    f.close()
    say("[B2] counts read (%.1f min)" % ((time.time() - t00) / 60))
    comp_rows = []
    de_parts = []
    endo_donors = sorted({str(x) for k, (ix, _) in endo_plans.items() for x in donor[ix]})
    for (s, cp), (allidx, ok_types) in endo_plans.items():
        data = [rowmat[int(i)][1] for i in allidx]
        ind = [rowmat[int(i)][0] for i in allidx]
        ptr = np.concatenate([[0], np.cumsum([len(x) for x in ind])])
        X = sp.csr_matrix((np.concatenate(data), np.concatenate(ind), ptr), shape=(len(allidx), len(genes)))[:, gidx]
        ct = np.asarray(ctype[allidx]).astype(str)
        dn = np.asarray(donor[allidx]).astype(str)
        for l in ok_types:
            for dd, c in pd.Series(dn[ct == l]).value_counts().items():
                comp_rows.append({"group": "%s/%s" % (l, s), "donor": dd, "cells": int(c)})
        present = set(dn)
        full_de = None
        for itn in ["FULL"] + endo_donors:
            if itn != "FULL" and itn not in present and full_de is not None:
                de_parts.append(full_de.assign(iteration=itn))
                continue
            m = np.ones(len(allidx), dtype=bool) if itn == "FULL" else (dn != itn)
            cnt = pd.Series(ct[m]).value_counts()
            okt = set(cnt[cnt >= MIN_CELLS_TYPE].index)
            m &= np.isin(ct, list(okt))
            if len(okt) < 2:
                continue
            t0 = time.time()
            D = wilcoxon_de(X[np.where(m)[0]], ct[m], gnames, s, cp)
            D["iteration"] = itn
            de_parts.append(D)
            if itn == "FULL":
                full_de = D
            say("    [B2] ENDO %-11s without %-10s %5d cells, %2d types (%.1f min)" % (cp, itn if itn != "FULL" else "(none)", int(m.sum()), len(okt), (time.time() - t0) / 60))
        del X
    del rowmat
    LOO = pd.concat(de_parts, ignore_index=True) if de_parts else pd.DataFrame()
    LOO[["iteration", "group", "label", "status", "compartment", "gene", "score", "lfc", "padj"]].to_csv(f_loo, index=False)
    np.save(f_pb, PB)
    PBG.to_csv(f_pbg, index=False)
    pd.DataFrame({"gene": [str(x) for x in gnames]}).to_csv(f_gen, index=False)
    CMP = pd.DataFrame(comp_rows)
    CMP.to_csv(f_comp, index=False)
    say("[B2] leave-one-out DE and pseudobulk cached (%.1f min total)" % ((time.time() - t00) / 60))
    return LOO, PB, PBG, CMP, np.array([str(x) for x in gnames])


def pb_de(PB, PBG, genes_sym):
    # donor-paired pseudobulk DE: type vs the rest of the same donor's cells in the same compartment; moderated paired t
    comp = {l: ("immune" if IMMUNE_RE.search(l) and not NOT_IMMUNE_RE.search(l) else "endometrial") for l in PBG["label"].unique()}
    PBG = PBG.assign(comp=PBG["label"].map(comp))
    out = {}
    for (s, cp), g in PBG.groupby(["status", "comp"]):
        tot = {d: PB[gg.index.to_numpy()].astype(np.float64).sum(axis=0) for d, gg in g.groupby("donor")}
        ok = g[g["n_cells"] >= MIN_CELLS_PB]
        for l, r in ok.groupby("label"):
            if r["donor"].nunique() < MIN_DONORS_PB:
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
            out["%s/%s" % (l, s)] = (pd.DataFrame({"gene": genes_sym[expr], "score": m / (se + s0), "lfc": m}), n)
    return out


def heca_rank_rows(res, meta, it_name, fam, t, extra):
    keys = [(k, m[2]) for k, m in meta.items() if m[1] == fam and m[0] in ("CTL", it_name)]
    items = rank_items(res, keys)
    if not items:
        return None
    labs = [l for l, _ in items]
    i_ds = next((i for i, l in enumerate(labs) if l == "dStromal_mid/ENDO"), None)
    bs = next((i for i, l in enumerate(labs) if re.search(STROMA_H, l)), None)
    bp = next((i for i, l in enumerate(labs) if re.search(PERIV_H, l)), None)
    v = items[i_ds][1] if i_ds is not None else {}
    return {"trait": t, "family": fam, "left_out": it_name, "first": labs[0], "first_p": items[0][1]["p"],
            "first_decidual": bool(re.search(DECID_H, labs[0])), "first_stromal": bool(re.search(STROMA_H, labs[0])),
            "rank_dS_mid_ENDO": (i_ds + 1) if i_ds is not None else np.nan, "p_dS_mid_ENDO": v.get("p", np.nan),
            "beta_std_dS_mid_ENDO": v.get("beta_std", np.nan), "se_dS_mid_ENDO": se_std(v) if v else np.nan,
            "best_stromal_rank": (bs + 1) if bs is not None else np.nan, "best_perivascular_rank": (bp + 1) if bp is not None else np.nan,
            "top3": "; ".join("%s %.0e" % (short(l, 18), x["p"]) for l, x in items[:3]), **extra}


def part_b2(sym2ent):
    say("")
    say("-" * 72)
    say("[SECTION] B2 HECA leave-one-ENDO-donor-out (fixed v041 cell sample; Wilcoxon within ENDO compartments; CTL groups unchanged) and B3 donor-paired pseudobulk DE")
    if not HECA_FILE.exists() and not (CACHE / "heca_loo_de.csv.gz").exists():
        say("[B2] HECA file missing; skipped")
        return
    DEh, uni = heca_sets(sym2ent)
    uset = set(uni)
    ents = [sym2ent[x] for x in uni]
    pos_of = {g: i for i, g in enumerate(uni)}
    LOO, PB, PBG, CMP, genes_sym = heca_compute(sym2ent)
    if CMP is not None and len(CMP):
        for grp in ("dStromal_mid/ENDO", "dStromal_early/ENDO"):
            c = CMP[CMP["group"] == grp].sort_values("cells", ascending=False)
            if len(c):
                say("[B2] donor composition of %s: %s (largest share %.0f%%)" % (grp, ", ".join("%s %d" % (r["donor"], r["cells"]) for _, r in c.iterrows()),
                                                                                 100 * c["cells"].max() / c["cells"].sum()))
        CMP.to_csv(OUT / "b2_donor_composition.csv", index=False)
    LOO["gene"] = LOO["gene"].astype(str)
    full = LOO[LOO["iteration"] == "FULL"]
    chk = []
    for grp in ("dStromal_mid/ENDO", "dStromal_early/ENDO"):
        a = full[full["group"] == grp].set_index("gene")["score"]
        b = DEh[DEh["group"] == grp].set_index("gene")["score"]
        j = pd.concat([a, b], axis=1, join="inner")
        if len(j):
            chk.append("%s r %.4f, max |diff| %.2g" % (grp, float(np.corrcoef(j.iloc[:, 0], j.iloc[:, 1])[0, 1]), float((j.iloc[:, 0] - j.iloc[:, 1]).abs().max())))
    say("[B2] full re-run vs v041 cache: %s" % ("; ".join(chk) or "no overlap"))
    cols, meta, sets_ds = {}, {}, {}

    def add(it_name, fam, grp, gl):
        k = "H%05d" % (len(cols) + 1)
        v = np.zeros(len(uni), dtype=np.int8)
        v[[pos_of[x] for x in gl]] = 1
        cols[k], meta[k] = v, (it_name, fam, grp)

    for grp, d in DEh[DEh["status"] == "CTL"].groupby("group", sort=False):
        for fam, gl in fam_sets(d, uset).items():
            add("CTL", fam, grp, gl)
    iters = ["FULL"] + sorted(x for x in LOO["iteration"].unique() if x != "FULL")
    for itn in iters:
        for grp, d in LOO[LOO["iteration"] == itn].groupby("group", sort=False):
            for fam, gl in fam_sets(d, uset).items():
                add(itn, fam, grp, gl)
                if grp == "dStromal_mid/ENDO":
                    sets_ds[(itn, fam)] = gl
    PBD = {}
    if len(genes_sym) == PB.shape[1]:
        PBD = pb_de(PB, PBG, genes_sym)
        for grp, (d, n) in PBD.items():
            gl = [x for x in d.sort_values("score", ascending=False)["gene"].head(N_MATCH).astype(str) if x in uset]
            if len(gl) >= MIN_SET:
                add("PB", "UPN", grp, gl)
                if grp == "dStromal_mid/ENDO":
                    sets_ds[("PB", "UPN")] = gl
        say("[B3] donor-paired pseudobulk DE: %d groups with >= %d donors of >= %d cells (%s)" % (
            len(PBD), MIN_DONORS_PB, MIN_CELLS_PB, ", ".join("%s n=%d" % (short(k, 16), n) for k, (_, n) in list(PBD.items())[:60] if re.search(STROMA_H, k))))
    else:
        say("[B3] gene order of the pseudobulk cache could not be matched; skipped")
    res = {}
    drop_covar_files("b2")
    for t in TRAITS_LOO:
        raw = raw_of(t)
        if raw is None:
            continue
        r = covar_runs(raw, cols, ents, "b2")
        if r:
            res[t] = r
    drop_covar_files("b2")
    rows = []
    for t, r in res.items():
        for fam in ("UPN", "UP"):
            for itn in iters:
                x = heca_rank_rows(r, meta, itn, fam, t, {"jaccard_vs_full": jaccard(sets_ds.get((itn, fam), []), sets_ds.get(("FULL", fam), []))})
                if x:
                    rows.append(x)
        keys = [(k, m[2]) for k, m in meta.items() if m[0] == "PB"]
        items = rank_items(r, keys)
        if items:
            labs = [l for l, _ in items]
            i_ds = next((i for i, l in enumerate(labs) if l == "dStromal_mid/ENDO"), None)
            bs = next((i for i, l in enumerate(labs) if re.search(STROMA_H, l)), None)
            bp = next((i for i, l in enumerate(labs) if re.search(PERIV_H, l)), None)
            rows.append({"trait": t, "family": "UPN_pseudobulk", "left_out": "-", "first": labs[0], "first_p": items[0][1]["p"],
                         "first_decidual": bool(re.search(DECID_H, labs[0])), "first_stromal": bool(re.search(STROMA_H, labs[0])),
                         "rank_dS_mid_ENDO": (i_ds + 1) if i_ds is not None else np.nan, "p_dS_mid_ENDO": items[i_ds][1]["p"] if i_ds is not None else np.nan,
                         "beta_std_dS_mid_ENDO": items[i_ds][1]["beta_std"] if i_ds is not None else np.nan, "se_dS_mid_ENDO": np.nan,
                         "best_stromal_rank": (bs + 1) if bs is not None else np.nan, "best_perivascular_rank": (bp + 1) if bp is not None else np.nan,
                         "top3": "; ".join("%s %.0e" % (short(l, 18), x["p"]) for l, x in items[:5]),
                         "jaccard_vs_full": jaccard(sets_ds.get(("PB", "UPN"), []), sets_ds.get(("FULL", "UPN"), []))})
    B = pd.DataFrame(rows)
    B.to_csv(OUT / "b2_b3_heca_loo_pseudobulk.csv", index=False)
    if not len(B):
        say("[B2] no MAGMA results")
        return
    verdict = {}
    for t in TRAITS_LOO:
        for fam in ("UPN", "UP"):
            g = B[(B["trait"] == t) & (B["family"] == fam)]
            f = g[g["left_out"] == "FULL"]
            l = g[g["left_out"] != "FULL"]
            if not len(f) or not len(l):
                continue
            f = f.iloc[0]
            say("  %-13s %-3s full: %s | LOO (n=%d): decidual first %d, dS_mid/ENDO first %d, stromal first %d, perivascular above stroma %d | dS_mid/ENDO BETA_STD %.4f [%.4f-%.4f], worst p %.1e (rank %s) | min Jaccard %.2f" % (
                t[:13], fam, f["top3"], len(l), int(l["first_decidual"].sum()), int((l["rank_dS_mid_ENDO"] == 1).sum()), int(l["first_stromal"].sum()),
                int((l["best_perivascular_rank"] < l["best_stromal_rank"]).sum()), f["beta_std_dS_mid_ENDO"], l["beta_std_dS_mid_ENDO"].min(),
                l["beta_std_dS_mid_ENDO"].max(), l["p_dS_mid_ENDO"].max(), int(l.loc[l["p_dS_mid_ENDO"].idxmax(), "rank_dS_mid_ENDO"]) if l["p_dS_mid_ENDO"].notna().any() else "-",
                l["jaccard_vs_full"].min()))
            for _, r in l.iterrows():
                say("       without %-10s %s" % (r["left_out"], r["top3"]))
            verdict["%s_%s" % (t, fam)] = {"n": len(l), "decidual_first": int(l["first_decidual"].sum()), "dsmid_first": int((l["rank_dS_mid_ENDO"] == 1).sum()),
                                          "stromal_first": int(l["first_stromal"].sum()), "perivascular_above": int((l["best_perivascular_rank"] < l["best_stromal_rank"]).sum())}
        p = B[(B["trait"] == t) & (B["family"] == "UPN_pseudobulk")]
        if len(p):
            p = p.iloc[0]
            say("  %-13s UPN donor-paired pseudobulk: %s | dS_mid/ENDO rank %s | best stromal %s, best perivascular %s | set Jaccard vs cell-level %.2f" % (
                t[:13], p["top3"], int(p["rank_dS_mid_ENDO"]) if np.isfinite(p["rank_dS_mid_ENDO"]) else "-",
                int(p["best_stromal_rank"]) if np.isfinite(p["best_stromal_rank"]) else "-", int(p["best_perivascular_rank"]) if np.isfinite(p["best_perivascular_rank"]) else "-",
                p["jaccard_vs_full"]))
            verdict["%s_PB" % t] = {"first": p["first"], "first_stromal": bool(p["first_stromal"]), "first_decidual": bool(p["first_decidual"])}
    v = verdict.get("KOLLER_EUR_UPN", {})
    need = (v["n"] - LOO_HECA_MAX_FAIL) if v else 0
    ok = bool(v) and v["decidual_first"] >= need
    say("[H_LOO_HECA_VERDICT] KOLLER_EUR UPN: decidual stromal state first in %s/%s leave-one-ENDO-donor-out runs (required >= %d); perivascular above stroma in %s -> %s" % (
        v.get("decidual_first", "-"), v.get("n", "-"), need, v.get("perivascular_above", "-"), "SUPPORTED" if ok else "NOT SUPPORTED"))
    pbv = verdict.get("KOLLER_EUR_PB", {})
    say("[H_HECA_PB_VERDICT] KOLLER_EUR UPN with donor-paired pseudobulk DE: first %s (stromal %s, decidual %s) -> %s" % (
        pbv.get("first", "-"), pbv.get("first_stromal", "-"), pbv.get("first_decidual", "-"), "SUPPORTED" if pbv.get("first_stromal") else "NOT SUPPORTED"))
    SUM["results"]["b2"] = {"supported": ok, "pb_supported": bool(pbv.get("first_stromal")), "by": verdict}


# ----------------------------------------------------------------------------
def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.54 | final statistical audit: benchmark CIs, dependent-GWAS contrasts, cross-atlas loci, donor stability | " + TS)
    say("=" * 72)
    say("[POLICY] rules fixed in CHANGELOG v0.2.44 before this run; this run answers review #3 (no new discovery analyses)")
    t0 = time.time()
    ent2sym, sym2ent, pos = gene_loc()
    steps = [("D", part_d), ("A1", part_a_genez), ("C", lambda: part_c(ent2sym, sym2ent, pos)), ("A2", lambda: part_a_rotation(sym2ent, pos)),
             ("B1", lambda: part_b1(sym2ent)), ("B2", lambda: part_b2(sym2ent)), ("A3", part_a_rg)]
    for name, fn in steps:
        try:
            fn()
        except Exception as e:
            ERR.append("%s: %s" % (name, e))
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v054_summary.json", "w") as f:
        json.dump(SUM, f, indent=1, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V054_COMPLETE" + (" with %d step errors" % len(ERR) if ERR else ""))
    _fh.close()


if __name__ == "__main__":
    main()
