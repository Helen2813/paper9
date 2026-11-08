#!/usr/bin/env python3
# Paper 9 | v0.0.41 | Matched directional benchmark (rules CHANGELOG v0.2.30): marker rules |log2FC| / up / down,
#   their size-matched versions and continuous specificity, on the same HECA cells (v031 P), under truths that do
#   and do not favour any rule (specificity, depletion, expression level), with 300 calibrated null rotations.
#   Real KOLLER_EUR gene z and correlations; annotation layer rotated along the genome. Resumable.
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). ~3-4 h (overnight); re-run continues.

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
import json
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
LDSC_PY = TOOLS / "ldsc_venv" / "bin" / "python"
LDSC = TOOLS / "ldsc_py3" / "ldsc.py"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
REF = RAW / "ldsc_ref" / "unpacked"
PLINK_DIR = REF / "1000G_EUR_Phase3_plink"
FRQ_DIR = REF / "1000G_Phase3_frq"
HM3_LIST = RAW / "ldsc_ref" / "hm3_no_MHC.list.txt"
MUNGED = PROC / "munged"
HECA_FILE = RAW / "heca" / "endometriumAtlasV2_cells_with_counts.h5ad"
V030_GS = PROC / "genesets" / "heca_koller_v030"
V030_LD = PROC / "ldscores" / "heca_koller_v030"
GS_BASE = PROC / "genesets" / "heca_koller_v031"
LD_BASE = PROC / "ldscores" / "heca_koller_v031"
LD_DIR = V030_LD  # reassigned per family (helpers below read this global)

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "benchmark" / ("p9_v041_" + TS)
LOG = REPO / "logs" / ("p9_v041_benchmark_" + TS + ".log")
NICE = ["nice", "-n", "10"]
N_WORKERS = 1
BLAS_THREADS = "2"
L2_ENV = dict(os.environ, OPENBLAS_NUM_THREADS=BLAS_THREADS, OMP_NUM_THREADS=BLAS_THREADS,
              VECLIB_MAXIMUM_THREADS=BLAS_THREADS, MKL_NUM_THREADS=BLAS_THREADS, NUMEXPR_NUM_THREADS=BLAS_THREADS)
WINDOW_BP = 100_000
LD_WIND_CM = 1.0

# Fixed (CHANGELOG v0.2.20)
DONOR_COL = "genotype"
TYPE_COL = "celltype"
PATHOLOGY_COL = "Endometrial_pathology"
STAGE_COL = "Binary Stage"
EXCLUDE_STAGE = "Hormones"
IMMUNE_RE = re.compile(r"Immune|^uM|uNK|^T_|^B_|Mast|DC|Mono|Macro|Neutro|Plasma|ILC|Lymphoid|Myeloid")
NOT_IMMUNE_RE = re.compile(r"Lymphatic")
MIN_CELLS_TYPE = 30
CAP_CELLS_TYPE = 400
SEED = 20261106
LFC = 1.2
PADJ = 0.01
MIN_SET = 20
TRAITS = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL", "HYPERTENSION", "T2D", "LEIOMYOMA"]
TRAITS_SENS = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL"]
TRAITS_COND = ["KOLLER_EUR", "GCST90205183"]
KOLLER_N = {"ENDO": 22, "CTL": 13}
FOCUS_RE = r"dStromal|SOX9_functionalis"
ENDOTHELIAL_RE = r"^(Venous|Arterial)$"
VARIANTS = {"P": {"hormones": False, "cap": CAP_CELLS_TYPE, "desc": "v030 cells (non-hormone, <= 400 per type)"},
            "H": {"hormones": True, "cap": CAP_CELLS_TYPE, "desc": "hormone-treated samples included"},
            "N": {"hormones": False, "cap": None, "desc": "all cells (no cap)"}}

for d in (OUT, LOG.parent, GS_BASE, LD_BASE):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V041_BENCHMARK", "timestamp_utc": TS, "gates": {}, "results": {}}


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def sh(cmd, cwd=None, timeout=None, env=None):
    try:
        p = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True, timeout=timeout, env=env)
        return p.returncode, p.stdout, p.stderr
    except FileNotFoundError as e:
        return 127, "", str(e)
    except OSError as e:
        return 126, "", str(e)
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"


def safe(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def dec(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def obs_cat(f, name):
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


def gene_map():
    sym2, dup = {}, set()
    for line in open(GENE_LOC_B37):
        t = line.split()
        if len(t) >= 6 and t[1] in [str(c) for c in range(1, 23)]:
            if t[5] in sym2:
                dup.add(t[5])
            sym2[t[5]] = (t[0], int(t[1]), int(t[2]), int(t[3]))
    for d in dup:
        sym2.pop(d, None)
    return sym2


# ----------------------------------------------------------------------------
# LDSC-SEG helpers reused verbatim from v012
# ----------------------------------------------------------------------------
def merged_intervals(gene_rows, chrom):
    iv = sorted((max(0, s - WINDOW_BP), e + WINDOW_BP) for c, s, e in gene_rows if c == chrom)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return np.array(out, dtype=np.int64).reshape(-1, 2)


def write_annotations(all_sets, sym2):
    keys = list(all_sets.keys())
    pd.Series(keys).to_csv(LD_DIR / "columns.txt", index=False, header=False)
    for chrom in range(1, 23):
        f = LD_DIR / ("annot.%d.annot" % chrom)
        old = LD_DIR / ("annot.%d.annot.gz" % chrom)
        if old.exists():
            old.unlink()
        if f.exists():
            continue
        bim = pd.read_csv(PLINK_DIR / ("1000G.EUR.QC.%d.bim" % chrom), sep=r"\s+", header=None,
                          names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
        bp = bim["BP"].to_numpy()
        df = bim[["CHR", "BP", "SNP", "CM"]].copy()
        for k in keys:
            rows = [(sym2[g][1], sym2[g][2], sym2[g][3]) for g in all_sets[k]]
            iv = merged_intervals(rows, chrom)
            if len(iv) == 0:
                df[k] = 0
                continue
            idx = np.searchsorted(iv[:, 0], bp, side="right") - 1
            inside = (idx >= 0) & (bp <= iv[np.clip(idx, 0, None), 1])
            df[k] = inside.astype(np.int8)
        df.to_csv(f, sep="\t", index=False)
        df[keys].to_csv(LD_DIR / ("thin.%d.annot" % chrom), sep="\t", index=False)
    pass


def l2_cmd(c, snps, thin):
    pref = ("thin.%d.annot" if thin else "annot.%d.annot") % c
    cmd = NICE + [str(LDSC_PY), str(LDSC), "--l2", "--bfile", str(PLINK_DIR / ("1000G.EUR.QC.%d" % c)),
           "--ld-wind-cm", str(LD_WIND_CM), "--annot", str(LD_DIR / pref),
           "--out", str(LD_DIR / ("multi.%d" % c))]
    if thin:
        cmd.append("--thin-annot")
    return cmd


def frq_maf(c):
    f = next(iter(FRQ_DIR.rglob("1000G.EUR.QC.%d.frq" % c)), None)
    if f is None:
        raise FileNotFoundError("frq for chr%d" % c)
    fr = pd.read_csv(f, sep=r"\s+")
    return fr[["SNP", "MAF"]]


def write_M(c):
    mf = LD_DIR / ("multi.%d.l2.M" % c)
    m5 = LD_DIR / ("multi.%d.l2.M_5_50" % c)
    if mf.exists() and mf.stat().st_size > 0 and m5.exists() and m5.stat().st_size > 0:
        return
    an = pd.read_csv(LD_DIR / ("annot.%d.annot" % c), sep="\t")
    keys = [k for k in an.columns if k not in ("CHR", "BP", "SNP", "CM")]
    maf = frq_maf(c)
    if len(maf) == len(an) and (maf["SNP"].values == an["SNP"].values).all():
        common = maf["MAF"].values > 0.05
    else:
        mm = an[["SNP"]].merge(maf, on="SNP", how="left")
        common = mm["MAF"].fillna(0).values > 0.05
    M = an[keys].sum(axis=0).astype(int).tolist()
    M550 = an.loc[common, keys].sum(axis=0).astype(int).tolist()
    mf.write_text("\t".join(str(x) for x in M) + "\n")
    m5.write_text("\t".join(str(x) for x in M550) + "\n")


def check_M_against_baseline(c=22):
    base = prefix_of("baseline_v1.2")
    f = Path(str(base) + "%d.l2.M_5_50" % c) if base else None
    if not f or not f.exists():
        say("[M_CHECK] baseline M_5_50 not found")
        return
    ref_first = float(open(f).read().split()[0])
    n_common = int((frq_maf(c)["MAF"].values > 0.05).sum())
    if abs(n_common - ref_first) >= 0.5:
        say("[M_CHECK] chr%d SNPs with MAF>0.05: ours=%d baseline=%.0f DIFFERENT" % (c, n_common, ref_first))
    SUM["gates"]["M_CHECK"] = abs(n_common - ref_first) < 0.5


def run_l2():
    snps = LD_DIR / "print_snps.txt"
    if not snps.exists():
        with open(HM3_LIST) as f, open(snps, "w") as o:
            for line in f:
                t = line.split()
                if t and t[0].startswith("rs"):
                    o.write(t[0] + "\n")
    thin = False
    if not (LD_DIR / "multi.22.l2.ldscore.gz").exists():
        for thin in (False, True):
            t0 = time.time()
            rc, o, e = sh(l2_cmd(22, snps, thin), timeout=3 * 3600, env=L2_ENV)
            ok = (LD_DIR / "multi.22.l2.ldscore.gz").exists()
            if ok:
                break
            logf = LD_DIR / "multi.22.log"
            txt = (logf.read_text(errors="ignore") if logf.exists() else "") + "\n" + e
            for line in txt.strip().splitlines()[-25:]:
                say("       " + line[:220])
        if not (LD_DIR / "multi.22.l2.ldscore.gz").exists():
            SUM["gates"]["LDSCORES"] = False
            return False
    else:
        thin = (LD_DIR / "thin.mode").exists()
    if thin:
        (LD_DIR / "thin.mode").write_text("1")
    if (LD_DIR / "multi.22.l2.ldscore.gz").exists():
        write_M(22)
    check_M_against_baseline(22)
    todo = [c for c in range(1, 23) if not (LD_DIR / ("multi.%d.l2.ldscore.gz" % c)).exists()]
    n_todo = len(todo)
    running = []
    t0 = time.time()
    while todo or running:
        while todo and len(running) < N_WORKERS:
            c = todo.pop(0)
            errf = open(LD_DIR / ("multi.%d.stderr.txt" % c), "w")
            p = subprocess.Popen(l2_cmd(c, snps, thin), stdout=subprocess.DEVNULL, stderr=errf, text=True, env=L2_ENV)
            running.append((c, p, errf))
        time.sleep(5)
        still = []
        for c, p, errf in running:
            if p.poll() is None:
                still.append((c, p, errf))
                continue
            errf.close()
            ok = (LD_DIR / ("multi.%d.l2.ldscore.gz" % c)).exists()
            if ok:
                write_M(c)
            if not ok:
                say("    [L2] chr%-2d rc=%d FAILED %.0f min" % (c, p.returncode, (time.time() - t0) / 60))
                logf = LD_DIR / ("multi.%d.log" % c)
                txt = logf.read_text(errors="ignore") if logf.exists() else ""
                for line in txt.strip().splitlines()[-8:]:
                    say("       " + line[:200])
        running = still
    ok_all = all((LD_DIR / ("multi.%d.l2.ldscore.gz" % c)).exists() for c in range(1, 23))
    if ok_all:
        for c in range(1, 23):
            write_M(c)
    SUM["gates"]["LDSCORES_" + LD_DIR.name] = ok_all
    say("[L2] %s: %d chromosomes computed, all 22 present: %s (%.0f min)" % (LD_DIR.name, n_todo, ok_all, (time.time() - t0) / 60))
    return ok_all


def split_ldscores(keys):
    marker = LD_DIR / "split.done"
    if marker.exists():
        return
    hm3 = set(open(LD_DIR / "print_snps.txt").read().split())
    n_all = n_kept = 0
    for c in range(1, 23):
        ld = pd.read_csv(LD_DIR / ("multi.%d.l2.ldscore.gz" % c), sep="\t")
        n_all += len(ld)
        ld = ld[ld["SNP"].isin(hm3)].reset_index(drop=True)
        n_kept += len(ld)
        M = open(LD_DIR / ("multi.%d.l2.M" % c)).read().split()
        M550 = open(LD_DIR / ("multi.%d.l2.M_5_50" % c)).read().split()
        cols = [x for x in ld.columns if x not in ("CHR", "SNP", "BP")]
        if len(cols) != len(keys):
            raise RuntimeError("chr%d: %d LD columns vs %d keys" % (c, len(cols), len(keys)))
        for j, k in enumerate(keys):
            d = LD_DIR / "sets" / k
            d.mkdir(parents=True, exist_ok=True)
            out = ld[["CHR", "SNP", "BP", cols[j]]].rename(columns={cols[j]: k + "L2"})
            out.to_csv(d / ("%s.%d.l2.ldscore.gz" % (k, c)), sep="\t", index=False, compression="gzip")
            open(d / ("%s.%d.l2.M" % (k, c)), "w").write(M[j] + "\n")
            open(d / ("%s.%d.l2.M_5_50" % (k, c)), "w").write(M550[j] + "\n")
    marker.write_text("ok")



def prefix_of(dirname):
    d = REF / dirname
    for p in sorted(d.rglob("*.l2.ldscore.gz")) if d.exists() else []:
        m = re.match(r"^(.*\.)1\.l2\.ldscore\.gz$", p.name)
        if m:
            return str(p.parent / m.group(1))
    return None


def ldsc_seg(trait, fam):
    out = OUT.parent / "cache" / ("%s__%s" % (trait, fam["family"]))
    out.parent.mkdir(parents=True, exist_ok=True)
    res = Path(str(out) + ".cell_type_results.txt")
    if not res.exists():
        ldcts = LD_DIR / ("%s.ldcts" % fam["family"])
        ctrl = str(LD_DIR / "sets" / fam["control"] / (fam["control"] + "."))
        with open(ldcts, "w") as f:
            for k in fam["sets"]:
                if k == fam["control"]:
                    continue
                f.write("%s\t%s,%s\n" % (k, str(LD_DIR / "sets" / k / (k + ".")), ctrl))
        cmd = NICE + [str(LDSC_PY), str(LDSC), "--h2-cts", str(MUNGED / (trait + ".sumstats.gz")),
               "--ref-ld-chr", BASELINE, "--ref-ld-chr-cts", str(ldcts), "--w-ld-chr", WEIGHTS, "--out", str(out)]
        t0 = time.time()
        rc, o, e = sh(cmd, timeout=6 * 3600, env=L2_ENV)
        if rc != 0 or not res.exists():
            say("    [SEG] %-16s %-22s rc=%d FAILED %.0f min" % (trait, fam["family"], rc, (time.time() - t0) / 60))
            say("       " + (o + e).strip()[-400:])
            return None
    df = pd.read_csv(res, sep="\t")
    lab = fam["manifest"].set_index("set_key")["label"]
    df["label"] = df["Name"].map(lab)
    df["family"] = fam["family"]
    df["trait"] = trait
    df["q_bh"] = bh(df["Coefficient_P_value"].values)
    return df.sort_values("Coefficient_P_value")


import math
# ----------------------------------------------------------------------------
# Fixed (CHANGELOG v0.2.30)
# ----------------------------------------------------------------------------
HECA_GS = PROC / "genesets" / "heca_v014"
DE_DIR = PROC / "genesets" / "heca_koller_v041_de"
WORK = INTERIM / "magma_v041"
WORK.mkdir(parents=True, exist_ok=True)
STATE = REPO / "results" / "benchmark" / "v041_replicates.csv"   # shared across re-runs (resume)
OBS_BETA_DSMID = 0.947
SIM_SEED = 20261108
N_NULL = 300
N_TRUTH = 40
STRENGTH = 2.0
N_MATCH = 400
TARGETS = ["dStromal_mid", "Venous", "ePV_2"]
KINDS = ["POS", "NEG", "EXPR"]      # specificity, depletion (negative specificity), expression level
FAMILIES = ["ABS", "UP", "DOWN", "ABSN", "UPN", "DOWNN", "SPEC"]
STATE_COLS = ["scenario", "kind", "target", "rep", "offset", "family", "top", "top_group", "any_q05", "rank_true", "top1_exact",
              "top3_exact", "top1_group", "related_miss", "true_q05"]
GROUPS = {"stromal": r"^(dStromal|eStromal|Fibroblast)", "endothelial": r"^(Arterial|Venous|Lymphatic)",
          "perivascular/SM": r"^(ePV|mPV|uSMC)"}


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


def group_of(label):
    b = str(label).split("/")[0]
    for g, rx in GROUPS.items():
        if re.search(rx, b):
            return g
    return "other"


def clopper(k, n, a=0.05):
    from scipy.stats import beta
    lo = 0.0 if k == 0 else float(beta.ppf(a / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a / 2, k + 1, n - k))
    return lo, hi


def markers_de(sym2):
    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp
    vkey = "P"
    V = VARIANTS[vkey]
    DE_DIR.mkdir(parents=True, exist_ok=True)
    done = DE_DIR / "de_all.csv.gz"
    if done.exists():
        DE = pd.read_csv(done)
        CPM = pd.read_csv(DE_DIR / "mean_cpm_by_group.csv.gz", index_col=0)
        say("[MARKERS P] cached DE tables: %d groups" % DE["group"].nunique())
        return DE, CPM
    t00 = time.time()
    f = h5py.File(HECA_FILE, "r")
    cols = [k for k in f["obs"].keys() if k not in ("_index", "__categories")]
    dataset_col, dataset_vals = find_dataset_col(f, cols)
    donor = obs_cat(f, DONOR_COL)
    ctype = obs_cat(f, TYPE_COL)
    path = obs_cat(f, PATHOLOGY_COL)
    stage = obs_cat(f, STAGE_COL)
    if V["hormones"]:
        keep = pd.notna(donor) & pd.notna(ctype)
    else:
        keep = pd.notna(donor) & pd.notna(ctype) & (stage != EXCLUDE_STAGE)
    if dataset_col:
        dv = obs_cat(f, dataset_col)
        keep &= np.isin(dv, dataset_vals)
    status = status_of(path)
    keep &= pd.notna(status)
    say("[MARKERS %s] %s: CTL cells %d donors %d | ENDO cells %d donors %d" % (
        vkey, V["desc"], int((keep & (status == "CTL")).sum()), pd.Series(donor[keep & (status == "CTL")]).nunique(),
        int((keep & (status == "ENDO")).sum()), pd.Series(donor[keep & (status == "ENDO")]).nunique()))
    labels = sorted(pd.unique(ctype[keep]))
    comp = {l: ("immune" if IMMUNE_RE.search(l) and not NOT_IMMUNE_RE.search(l) else "endometrial") for l in labels}
    layer = f["layers"]["counts"]
    var = f["var"]
    gobj = var[dec(var.attrs.get("_index", "_index"))]
    if isinstance(gobj, h5py.Group):
        gobj = gobj["values"] if "values" in gobj else gobj[list(gobj.keys())[0]]
    genes = np.array([dec(x) for x in gobj[:]], dtype=object)
    gmask = np.array([g in sym2 for g in genes])
    gidx = np.where(gmask)[0]
    indptr = layer["indptr"][:]
    rng = np.random.default_rng(SEED)
    cap = V["cap"]
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
                if cap is not None and len(idx) > cap:
                    idx = np.sort(rng.choice(idx, cap, replace=False))
                pick.append(idx)
                ok_types.append(l)
            if len(ok_types) < 2:
                say("[MARKERS %s] %s %s: < 2 types, skipped" % (vkey, s, cp))
                continue
            plan[(s, cp)] = (np.sort(np.concatenate(pick)), ok_types)
    want = np.unique(np.concatenate([v[0] for v in plan.values()]))
    rowmat = {}
    step = 20000
    n = len(indptr) - 1
    for a in range(0, n, step):
        b = min(n, a + step)
        lo = np.searchsorted(want, a)
        hi = np.searchsorted(want, b)
        if hi <= lo:
            continue
        s0, e0 = int(indptr[a]), int(indptr[b])
        idx = np.asarray(layer["indices"][s0:e0])
        val = np.asarray(layer["data"][s0:e0], dtype=np.float32)
        loc = indptr[a:b + 1] - s0
        for i in want[lo:hi]:
            r = i - a
            rowmat[int(i)] = (idx[loc[r]:loc[r + 1]], val[loc[r]:loc[r + 1]])
    f.close()
    rows, sets = [], {}
    de_parts, cpm_cols = [], {}
    for (s, cp), (allidx, ok_types) in plan.items():
        data = [rowmat[int(i)][1] for i in allidx]
        ind = [rowmat[int(i)][0] for i in allidx]
        ptr = np.concatenate([[0], np.cumsum([len(x) for x in ind])])
        X = sp.csr_matrix((np.concatenate(data), np.concatenate(ind), ptr), shape=(len(allidx), len(genes)))
        X = X[:, gidx]
        lib = np.asarray(X.sum(axis=1)).ravel()
        lib[lib == 0] = 1
        Xc = sp.diags(1e6 / lib) @ X
        for l in ok_types:
            mk = np.asarray(ctype[allidx] == l)
            cpm_cols["%s/%s" % (l, s)] = np.asarray(Xc[mk].mean(axis=0)).ravel()
        A = ad.AnnData(X=X)
        A.var_names = list(genes[gidx])
        A.obs["type"] = pd.Categorical(ctype[allidx])
        sc.pp.normalize_total(A, target_sum=1e4)
        sc.pp.log1p(A)
        sc.tl.rank_genes_groups(A, "type", method="wilcoxon", reference="rest", use_raw=False, pts=False)
        for l in ok_types:
            df = sc.get.rank_genes_groups_df(A, group=l)
            de_parts.append(pd.DataFrame({"group": "%s/%s" % (l, s), "label": l, "status": s, "compartment": cp, "gene": df["names"].astype(str).values,
                                          "score": df["scores"].values, "lfc": df["logfoldchanges"].values, "padj": df["pvals_adj"].values}))
            sig = df["pvals_adj"] < PADJ
            for mode, sel in (("ABS", (df["logfoldchanges"].abs() > LFC) & sig),
                              ("UP", (df["logfoldchanges"] > LFC) & sig),
                              ("DOWN", (df["logfoldchanges"] < -LFC) & sig)):
                g = list(df.loc[sel, "names"].astype(str))
                key = "%s_%s_%s_%s" % (mode, s, cp[:3], safe(l))[:60]
                rows.append({"set_key": key, "mode": mode, "status": s, "compartment": cp, "label": l, "n_genes": len(g),
                             "n_up": int((sel & (df["logfoldchanges"] > 0)).sum()), "n_down": int((sel & (df["logfoldchanges"] < 0)).sum()),
                             "n_cells_type": int((ctype[allidx] == l).sum()), "kept": len(g) >= MIN_SET})
                if len(g) >= MIN_SET:
                    sets[key] = g
        sets.setdefault("ALLGENES", list(A.var_names))
        del A, X
    del rowmat
    DE = pd.concat(de_parts, ignore_index=True)
    DE.to_csv(DE_DIR / "de_all.csv.gz", index=False)
    CPM = pd.DataFrame(cpm_cols, index=list(genes[gidx]))
    CPM.to_csv(DE_DIR / "mean_cpm_by_group.csv.gz")
    say("[MARKERS P] DE tables %d groups, %d genes (%.1f min)" % (DE["group"].nunique(), len(genes[gidx]), (time.time() - t00) / 60))
    return DE, CPM



# ----------------------------------------------------------------------------
# Families and continuous specificity from the same cells
# ----------------------------------------------------------------------------
def build_families(DE, sym2):
    fams = {}   # key -> list of symbols
    meta = {}   # key -> (family, group label)
    for g, d in DE.groupby("group", sort=False):
        d = d.copy()
        sig = d["padj"] < PADJ
        lab = g
        st, cp, l = d["status"].iloc[0], d["compartment"].iloc[0], d["label"].iloc[0]
        sel = {"ABS": d.loc[(d["lfc"].abs() > LFC) & sig, "gene"], "UP": d.loc[(d["lfc"] > LFC) & sig, "gene"],
               "DOWN": d.loc[(d["lfc"] < -LFC) & sig, "gene"],
               "ABSN": d.reindex(d["score"].abs().sort_values(ascending=False).index)["gene"].head(N_MATCH),
               "UPN": d.sort_values("score", ascending=False)["gene"].head(N_MATCH),
               "DOWNN": d.sort_values("score", ascending=True)["gene"].head(N_MATCH)}
        for fam, genes in sel.items():
            gl = [x for x in genes.astype(str) if x in sym2]
            if len(gl) >= MIN_SET:
                key = ("%s_%s_%s_%s" % (fam, st, cp[:3], safe(l)))[:60]
                fams[key] = gl
                meta[key] = (fam, lab)
    allg = sorted(set(DE["gene"].astype(str)) & set(sym2))
    fams["ALLGENES"] = allg
    return fams, meta


def spec_covariates(CPM, sym2):
    C = CPM[[c for c in CPM.columns]].copy()
    C = C[(C >= 1.0).any(axis=1)]
    C = C.loc[[g for g in C.index if g in sym2]]
    S = C.div(C.sum(axis=1), axis=0)
    names = {"C%03d" % (i + 1): c for i, c in enumerate(S.columns)}
    S.columns = list(names)
    S["Average"] = np.log2(C.mean(axis=1) + 1.0)
    return S, names


# ----------------------------------------------------------------------------
# MAGMA runs
# ----------------------------------------------------------------------------
def read_raw(trait):
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    lines = raw.read_text().splitlines()
    data_idx = [i for i, l in enumerate(lines) if l.strip() and not l.startswith("#")]
    out = pd.read_csv(MAGMA_OUT / (trait + "_35_10.genes.out"), sep=r"\s+", dtype={"GENE": str})
    zmap = dict(zip(out["GENE"].astype(str), out["ZSTAT"]))
    votes = {}
    for i in data_idx[:300]:
        tok = lines[i].split()
        z = zmap.get(tok[0])
        if z is None:
            continue
        for j in range(1, min(len(tok), 12)):
            try:
                if abs(float(tok[j]) - z) < 1e-3:
                    votes[j] = votes.get(j, 0) + 1
            except ValueError:
                pass
    zi = max(votes, key=votes.get)
    genes = [lines[i].split()[0] for i in data_idx]
    z = np.array([float(lines[i].split()[zi]) for i in data_idx])
    return lines, data_idx, zi, genes, z


def write_raw(lines, data_idx, zi, znew, path):
    out = list(lines)
    for k, i in enumerate(data_idx):
        t2 = out[i].split()
        t2[zi] = "%.6g" % znew[k]
        out[i] = " ".join(t2)
    Path(path).write_text("\n".join(out) + "\n")


def parse_gsa(gsa):
    rows = [l for l in open(gsa) if not l.startswith("#") and l.strip()]
    hdr = rows[0].split()
    res = {}
    for l in rows[1:]:
        d = dict(zip(hdr, l.split()))
        try:
            res[d.get("FULL_NAME", d.get("VARIABLE"))] = float(d["P"])
        except (KeyError, ValueError):
            continue
    return res


def magma_sets(raw, sf, tag):
    outp = WORK / tag
    gsa = Path(str(outp) + ".gsa.out")
    for model in (["--model", "condition-hide=ALLGENES"], []):
        if gsa.exists():
            gsa.unlink()
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--set-annot", sf] + model + ["--out", outp], timeout=3600, env=L2_ENV)
        if gsa.exists():
            return parse_gsa(gsa)
    say("    [MAGMA_FAIL] sets %s" % (o + e).strip()[-200:])
    return None


def magma_covar(raw, cf, tag):
    outp = WORK / tag
    gsa = Path(str(outp) + ".gsa.out")
    for key in ("condition-hide=", "condition="):
        if gsa.exists():
            gsa.unlink()
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--gene-covar", cf, "--model", key + "Average", "direction=greater", "--out", outp],
                      timeout=3600, env=L2_ENV)
        if gsa.exists():
            return parse_gsa(gsa)
    say("    [MAGMA_FAIL] covar %s" % (o + e).strip()[-200:])
    return None


# ----------------------------------------------------------------------------
# Benchmark
# ----------------------------------------------------------------------------
def summarize(res_sets, res_cov, meta, cov_names, target):
    out = {}
    for fam in FAMILIES:
        if fam == "SPEC":
            if not res_cov:
                continue
            items = [(cov_names[k], p) for k, p in res_cov.items() if k in cov_names]
        else:
            if not res_sets:
                continue
            items = [(meta[k][1], p) for k, p in res_sets.items() if k in meta and meta[k][0] == fam]
        if not items:
            continue
        items.sort(key=lambda x: x[1])
        labs = [x[0] for x in items]
        p = np.array([x[1] for x in items])
        q = bh(p)
        base = [l.split("/")[0] for l in labs]
        r = {"family": fam, "top": labs[0], "top_group": group_of(labs[0]), "any_q05": bool((q < 0.05).any())}
        if target is not None:
            hit = [i for i, b in enumerate(base) if b == target]
            r["rank_true"] = (hit[0] + 1) if hit else np.nan
            r["top1_exact"] = bool(hit and hit[0] == 0)
            r["top3_exact"] = bool(hit and hit[0] < 3)
            tg = group_of(target)
            r["top1_group"] = group_of(labs[0]) == tg
            r["related_miss"] = (not r["top1_exact"]) and r["top1_group"]
            r["true_q05"] = bool(hit and q[hit[0]] < 0.05)
        out[fam] = r
    return out


def run_benchmark(fams, meta, S_cov, cov_names, sym2):
    say("")
    say("-" * 72)
    say("[SECTION] Benchmark (rotation design; real %s z; families %s)" % ("KOLLER_EUR", ", ".join(FAMILIES)))
    lines, data_idx, zi, genes, z0 = read_raw("KOLLER_EUR")
    n = len(genes)
    idx_of = {g: i for i, g in enumerate(genes)}
    keys = list(fams)
    set_idx = {k: np.array(sorted({idx_of[sym2[x][0]] for x in fams[k] if sym2[x][0] in idx_of}), dtype=int) for k in keys}
    cov_idx = np.array([idx_of.get(sym2[g][0], -1) for g in S_cov.index])
    okc = cov_idx >= 0
    covvals = S_cov.values[okc]
    cov_idx = cov_idx[okc]
    Hs = pd.read_csv(HECA_GS / "specificity_by_type.csv.gz", index_col=0)
    Hm = pd.read_csv(HECA_GS / "mean_cpm_by_type.csv.gz", index_col=0)
    ent = {sym2[g][0]: g for g in Hs.index if g in sym2}
    sym_of = np.array([ent.get(g) for g in genes], dtype=object)
    has = np.array([s is not None for s in sym_of])
    b_std = OBS_BETA_DSMID * float(Hs["dStromal_mid"].std())
    say("[BENCH] genes.raw %d genes; sets %d (%s); continuous covariates %d; effect unit %.4f" % (
        n, len(keys) - 1, ", ".join("%s %d" % (f, sum(1 for k in meta if meta[k][0] == f)) for f in FAMILIES if f != "SPEC"), len(cov_names), b_std))

    def profile(kind, target):
        p_ = np.zeros(n)
        if kind == "EXPR":
            v = np.log2(Hm.loc[list(sym_of[has]), target].to_numpy(dtype=float) + 1.0)
        else:
            v = Hs.loc[list(sym_of[has]), target].to_numpy(dtype=float)
        v = (v - v.mean()) / v.std()
        p_[has] = -v if kind == "NEG" else v
        return p_

    sf, cf, simraw = WORK / "sets_rot.txt", WORK / "covar_rot.txt", WORK / "SIM_35_10.genes.raw"
    realraw = MAGMA_OUT / "KOLLER_EUR_35_10.genes.raw"

    def write_layer(off):
        with open(sf, "w") as o:
            for k in keys:
                o.write(k + " " + " ".join(genes[i] for i in (set_idx[k] + off) % n) + "\n")
        with open(cf, "w") as o:
            o.write("GENE\t" + "\t".join(S_cov.columns) + "\n")
            for gi, row in zip((cov_idx + off) % n, covvals):
                o.write(genes[gi] + "\t" + "\t".join("%.6g" % x for x in row) + "\n")

    done = set()
    if STATE.exists():
        st = pd.read_csv(STATE)
        done = set(zip(st["scenario"], st["rep"]))
    scen = [("NULLROT", None, None)] + [("%s_%s" % (k, t), k, t) for t in TARGETS for k in KINDS if t in Hs.columns]
    rng = np.random.default_rng(SIM_SEED)
    offs = {s[0]: [int(rng.integers(500, n - 500)) for _ in range(N_NULL if s[0] == "NULLROT" else N_TRUTH)] for s in scen}
    t0 = time.time()
    for name, kind, target in scen:
        prof = None if kind is None else STRENGTH * b_std * profile(kind, target)
        todo = [r for r in range(len(offs[name])) if (name, r) not in done]
        for rep in todo:
            off = offs[name][rep]
            write_layer(off)
            if prof is None:
                raw = realraw
            else:
                write_raw(lines, data_idx, zi, z0 + np.roll(prof, off), simraw)
                raw = simraw
            rs = magma_sets(raw, sf, "sim_sets")
            rc = magma_covar(raw, cf, "sim_cov")
            rows = [{"scenario": name, "kind": kind or "NULLROT", "target": target or "-", "rep": rep, "offset": off, **v}
                    for v in summarize(rs, rc, meta, cov_names, target).values()]
            if rows:
                pd.DataFrame(rows).reindex(columns=STATE_COLS).to_csv(STATE, mode="a", header=not STATE.exists(), index=False)
        say("    [BENCH] %-18s %d/%d replicates (%.0f min)" % (name, len(offs[name]), len(offs[name]), (time.time() - t0) / 60))
    return pd.read_csv(STATE)


def report(R):
    say("")
    say("-" * 72)
    say("[SECTION] Results")
    nul = R[R["scenario"] == "NULLROT"]
    say("[CALIB] FWER (any BH q < 0.05 within family) under %d null rotations, 95%% Clopper-Pearson:" % nul["rep"].nunique())
    cal = {}
    for fam in FAMILIES:
        g = nul[nul["family"] == fam]
        if not len(g):
            continue
        k, m = int(g["any_q05"].sum()), len(g)
        lo, hi = clopper(k, m)
        cal[fam] = {"k": k, "n": m, "fwer": k / m, "lo": lo, "hi": hi, "ok": lo <= 0.05}
        se = math.sqrt(max(k / m * (1 - k / m), 1e-12) / m)
        say("    %-6s %3d/%d = %.3f (se %.3f; 95%% CI %.3f-%.3f) -> %s" % (fam, k, m, k / m, se, lo, hi, "consistent with <= 5%" if lo <= 0.05 else "ANTI-CONSERVATIVE"))
    say("[RECOVERY] per truth: top-1 exact type / top-3 exact / top-1 correct group / related-state miss / true type q<0.05")
    T = R[R["scenario"] != "NULLROT"]
    tab = T.groupby(["kind", "target", "family"])[["top1_exact", "top3_exact", "top1_group", "related_miss", "true_q05"]].mean().reset_index()
    tab.to_csv(OUT / "recovery_table.csv", index=False)
    for (kind, target), g in tab.groupby(["kind", "target"], sort=False):
        say("  %s %s" % (kind, target))
        for _, r in g.set_index("family").reindex(FAMILIES).dropna(how="all").reset_index().iterrows():
            say("    %-6s %.2f / %.2f / %.2f / %.2f / %.2f" % (r["family"], r["top1_exact"], r["top3_exact"], r["top1_group"], r["related_miss"], r["true_q05"]))
    # Fixed readings (v0.2.30)
    def t1(kind, target, fam):
        v = tab[(tab["kind"] == kind) & (tab["target"] == target) & (tab["family"] == fam)]["top1_exact"]
        return float(v.iloc[0]) if len(v) else np.nan
    tg = [t for t in TARGETS if t in set(tab["target"])]
    p1 = all(t1("POS", t, "UPN") >= t1("POS", t, "ABSN") for t in tg)
    p2 = all(max(t1("NEG", t, "DOWNN"), t1("NEG", t, "ABSN")) > t1("NEG", t, "UPN") for t in tg)
    p3 = sum(t1("POS", t, "UPN") - t1("POS", t, "ABSN") >= 0.2 for t in tg) >= 2
    p4 = all(max(t1("EXPR", t, f) for f in FAMILIES) <= 0.5 for t in tg)
    say("[BENCH_VERDICT] P1 size-matched up >= size-matched |log2FC| under specificity truths (all targets): %s" % ("YES" if p1 else "NO"))
    say("[BENCH_VERDICT] P2 under depletion truths down or |log2FC| beats up (all targets): %s" % ("YES" if p2 else "NO"))
    say("[BENCH_VERDICT] P3 up advantage survives size matching (>= 0.2 for >= 2 of 3 targets): %s" % ("YES" if p3 else "NO"))
    say("[EXPLORATORY] P4 (reclassified before the run) highest top-1 exact under expression-level truths <= 0.5 for all targets: %s" % ("YES" if p4 else "NO"))
    SUM["results"] = {"calibration": cal, "P1": p1, "P2": p2, "P3": p3, "P4": p4}


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.41 | matched directional benchmark of marker rules | " + TS)
    say("=" * 72)
    t0 = time.time()
    sym2 = gene_map()
    try:
        DE, CPM = markers_de(sym2)
        fams, meta = build_families(DE, sym2)
        # identity check against the v031 thresholded sets
        mf = GS_BASE / "P" / "manifest.csv"
        if mf.exists():
            man = pd.read_csv(mf)
            same = sum(1 for k in man["set_key"] if k in fams and set(fams[k]) == set(open(GS_BASE / "P" / (k + ".txt")).read().split()) & set(sym2))
            say("[CHECK] thresholded sets identical to v031: %d of %d" % (same, len(man)))
        S_cov, cov_names = spec_covariates(CPM, sym2)
        R = run_benchmark(fams, meta, S_cov, cov_names, sym2)
        report(R)
    except Exception as e:
        import traceback
        say("[STEP_ERROR] %s: %s" % (type(e).__name__, str(e)[:300]))
        say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v041_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] prospectively specified follow-up analysis (CHANGELOG v0.2.30-v0.2.31), designed after exploratory results v001-v040")
    say("[DONE] PAPER9_V041_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
