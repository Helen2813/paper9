#!/usr/bin/env python3
# Paper 9 | v0.0.31 | A4/A5 follow-up of the Koller 2026 pipeline (rules CHANGELOG v0.2.21):
#   K1 cohort audit (donors / samples vs Koller 22 + 13); K2 composition of |log2FC| sets (share of down genes);
#   K3 S-LDSC families UP and DOWN; K4 pairwise joint S-LDSC (focus | top set and reverse); K5 sensitivity
#   (hormone-treated included; all cells); K6 FinnGen endometriosis as cohort-matched comparator, positive controls.
# Runs in ~/Desktop/paper9/.venv, low-heat (1 worker, 2 BLAS threads, nice 10). ~2 h; resumable (re-run continues).

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
OUT = REPO / "results" / "koller_repro" / ("p9_v031_" + TS)
LOG = REPO / "logs" / ("p9_v031_koller_resolution_" + TS + ".log")
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
SUM = {"task": "PAPER9_V031_KOLLER_RESOLUTION", "timestamp_utc": TS, "gates": {}, "results": {}}


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


# ----------------------------------------------------------------------------
# Shared helpers
# ----------------------------------------------------------------------------
SAMPLE_RE = re.compile(r"^(sample|sample_id|sampleid|sample_name|library|library_id|libraryid)$", re.I)
VASC_RE = r"^(Arterial|Venous|mPV|ePV_1a|ePV_1b|ePV_2|uSMCs)/"
LEIO_RE = r"^(uSMCs|mPV)/"


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


def is_focus(label):
    return bool(re.search(FOCUS_RE, str(label)))


def short(label, n=22):
    s = str(label)
    s = s.replace("SOX9_functionalis_", "SOX9f_").replace("dStromal_", "dS_").replace("Glandular_secretory", "Gland_sec")
    return s[:n]


# ----------------------------------------------------------------------------
# K1. Cohort audit (donors and samples vs Koller 22 endometriosis / 13 control samples)
# ----------------------------------------------------------------------------
def cohort_audit():
    say("")
    say("-" * 72)
    say("[SECTION] K1 cohort audit")
    f = h5py.File(HECA_FILE, "r")
    cols = [k for k in f["obs"].keys() if k not in ("_index", "__categories")]
    say("[OBS] columns: " + ", ".join(cols))
    dcol, dvals = find_dataset_col(f, cols)
    scol = next((c for c in cols if SAMPLE_RE.match(c)), None)
    if scol is None:
        scol = next((c for c in cols if "sample" in c.lower()), None)
    say("[COHORT] dataset column: %s %s | sample column: %s" % (dcol, dvals, scol or "none (donor used)"))
    donor = obs_cat(f, DONOR_COL)
    path = obs_cat(f, PATHOLOGY_COL)
    stage = obs_cat(f, STAGE_COL)
    smp = obs_cat(f, scol) if scol else donor
    dv = obs_cat(f, dcol) if dcol else np.full(len(donor), "all", dtype=object)
    f.close()
    status = status_of(path)
    mar = np.isin(dv, dvals) if dcol else np.ones(len(donor), bool)
    say("[COHORT] Mareckova '%s' values: %s" % (STAGE_COL, ", ".join("%s=%d" % (k, n) for k, n in pd.Series(stage[mar]).fillna("NA").value_counts().items())))
    rows = []
    for scope, m0 in (("Mareckova", mar), ("all HECA", np.ones(len(donor), bool))):
        for horm in (False, True):
            m = m0 & pd.notna(donor) & pd.notna(status)
            if not horm:
                m &= (stage != EXCLUDE_STAGE)
            r = {"scope": scope, "hormones_included": horm}
            for s in ("ENDO", "CTL"):
                mm = m & (status == s)
                r[s + "_donors"] = int(pd.Series(donor[mm]).nunique())
                r[s + "_samples"] = int(pd.Series(smp[mm]).nunique())
                r[s + "_cells"] = int(mm.sum())
            r["diff_to_koller"] = abs(r["ENDO_samples"] - KOLLER_N["ENDO"]) + abs(r["CTL_samples"] - KOLLER_N["CTL"])
            rows.append(r)
            say("[COHORT] %-9s hormones %-8s ENDO donors %3d samples %3d | CTL donors %3d samples %3d | diff to 22/13 = %d" % (
                scope, "included" if horm else "excluded", r["ENDO_donors"], r["ENDO_samples"], r["CTL_donors"], r["CTL_samples"], r["diff_to_koller"]))
    D = pd.DataFrame(rows)
    D.to_csv(OUT / "k1_cohort.csv", index=False)
    b = D.sort_values("diff_to_koller").iloc[0]
    say("[COHORT] closest to Koller's 22/13 samples: %s, hormones %s" % (b["scope"], "included" if b["hormones_included"] else "excluded"))
    SUM["results"]["k1_cohort"] = rows


# ----------------------------------------------------------------------------
# Markers per variant (v030 code; modes ABS, UP, DOWN)
# ----------------------------------------------------------------------------
def markers(sym2, vkey):
    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp
    V = VARIANTS[vkey]
    gs = GS_BASE / vkey
    gs.mkdir(parents=True, exist_ok=True)
    man_f = gs / "manifest.csv"
    if man_f.exists():
        man = pd.read_csv(man_f)
        sets = {r["set_key"]: open(gs / (r["set_key"] + ".txt")).read().split() for _, r in man.iterrows()}
        say("[MARKERS %s] cached: %d sets" % (vkey, len(man)))
        return man, sets
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
    for (s, cp), (allidx, ok_types) in plan.items():
        data = [rowmat[int(i)][1] for i in allidx]
        ind = [rowmat[int(i)][0] for i in allidx]
        ptr = np.concatenate([[0], np.cumsum([len(x) for x in ind])])
        X = sp.csr_matrix((np.concatenate(data), np.concatenate(ind), ptr), shape=(len(allidx), len(genes)))
        X = X[:, gidx]
        A = ad.AnnData(X=X)
        A.var_names = list(genes[gidx])
        A.obs["type"] = pd.Categorical(ctype[allidx])
        sc.pp.normalize_total(A, target_sum=1e4)
        sc.pp.log1p(A)
        sc.tl.rank_genes_groups(A, "type", method="wilcoxon", reference="rest", use_raw=False, pts=False)
        for l in ok_types:
            df = sc.get.rank_genes_groups_df(A, group=l)
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
    man = pd.DataFrame(rows)
    man = pd.concat([man, pd.DataFrame([{"set_key": "ALLGENES", "mode": "CTRL", "status": "-", "compartment": "-", "label": "all tested genes",
                                         "n_genes": len(sets["ALLGENES"]), "n_up": 0, "n_down": 0, "n_cells_type": 0, "kept": True}])],
                    ignore_index=True)
    for k, g in sets.items():
        (gs / (k + ".txt")).write_text("\n".join(g) + "\n")
    man.to_csv(gs / "manifest_all.csv", index=False)
    man = man[man["kept"]].reset_index(drop=True)
    man.to_csv(man_f, index=False)
    med = {m: int(man.loc[man["mode"] == m, "n_genes"].median()) if (man["mode"] == m).any() else 0 for m in ("ABS", "UP", "DOWN")}
    say("[MARKERS %s] sets kept: ABS %d UP %d DOWN %d; median size ABS %d UP %d DOWN %d (%.1f min)" % (
        vkey, (man["mode"] == "ABS").sum(), (man["mode"] == "UP").sum(), (man["mode"] == "DOWN").sum(), med["ABS"], med["UP"], med["DOWN"],
        (time.time() - t00) / 60))
    return man, sets


def same_as_v030(man, sets):
    m30 = V030_GS / "manifest.csv"
    if not m30.exists() or not (V030_LD / "split.done").exists():
        return False
    man30 = pd.read_csv(m30)
    k30 = list(man30.loc[man30["mode"] == "ABS", "set_key"]) + ["ALLGENES"]
    k31 = list(man.loc[man["mode"] == "ABS", "set_key"]) + ["ALLGENES"]
    if k30 != k31:
        return False
    if (V030_LD / "columns.txt").read_text().split() != k30:
        return False
    for k in k30:
        if open(V030_GS / (k + ".txt")).read().split() != sets[k]:
            return False
    return True


# ----------------------------------------------------------------------------
# S-LDSC family runner (v030 layout)
# ----------------------------------------------------------------------------
def family(vkey, mode, man, sets, sym2, traits, reuse=False):
    global LD_DIR
    keys = list(man.loc[man["mode"] == mode, "set_key"]) + ["ALLGENES"]
    if reuse:
        LD_DIR, fname = V030_LD, "KOLLER_ABS"
    else:
        LD_DIR, fname = LD_BASE / ("%s_%s" % (vkey, mode)), "K31_%s_%s" % (vkey, mode)
    LD_DIR.mkdir(parents=True, exist_ok=True)
    all_sets = {k: [g for g in sets[k] if g in sym2] for k in keys}
    if not (LD_DIR / "split.done").exists():
        kf = LD_DIR / "columns.txt"
        if kf.exists() and kf.read_text().split() != keys:
            for p in LD_DIR.glob("*"):
                if p.is_file():
                    p.unlink()
        write_annotations(all_sets, sym2)
        if not run_l2():
            say("[STOP] %s: LD scores incomplete; re-run to resume" % fname)
            return None, None
        split_ldscores(keys)
        if not reuse:
            for pat in ("annot.*.annot", "thin.*.annot", "multi.*.l2.ldscore.gz"):
                for p in LD_DIR.glob(pat):
                    p.unlink()
    sub = man[man["set_key"].isin(keys)][["set_key", "label"]]
    sub = sub.assign(label=sub["label"] + "/" + sub["set_key"].str.split("_").str[1])
    fam = {"family": fname, "sets": keys, "control": "ALLGENES", "manifest": sub, "dir": LD_DIR, "tag": "%s_%s" % (vkey, mode)}
    res = []
    for t in traits:
        if not (MUNGED / (t + ".sumstats.gz")).exists():
            say("    [SEG] %s sumstats missing" % t)
            continue
        df = ldsc_seg(t, fam)
        if df is not None:
            res.append(df)
    if not res:
        return fam, None
    R = pd.concat(res)
    R["variant"], R["mode"] = vkey, mode
    R.to_csv(OUT / ("sldsc_%s_%s.csv" % (vkey, mode)), index=False)
    say("[FAMILY %s_%s] %d annotations" % (vkey, mode, len(keys) - 1))
    for t in traits:
        d = R[R["trait"] == t].sort_values("Coefficient_P_value")
        if not len(d):
            continue
        top = "; ".join("%s %.0e/%.2f" % (short(r["label"]), r["Coefficient_P_value"], r["q_bh"]) for _, r in d.head(3).iterrows())
        foc = d[d["label"].map(is_focus)].head(2)
        fs = "; ".join("%s %.2g/%.2f" % (short(r["label"]), r["Coefficient_P_value"], r["q_bh"]) for _, r in foc.iterrows())
        say("    %-12s top(p/q): %s | focus: %s" % (t[:12], top, fs))
    return fam, R


def rank_of(R, trait, pattern, k=6):
    d = R[R["trait"] == trait].sort_values("Coefficient_P_value").head(k)
    return [l for l in d["label"] if re.search(pattern, str(l).split("/")[0])]


# ----------------------------------------------------------------------------
# K2 composition and K3 H-ABS reading
# ----------------------------------------------------------------------------
def composition(man, sets, R_abs):
    say("")
    say("-" * 72)
    say("[SECTION] K2 composition of |log2FC| sets (P cells)")
    A = man[man["mode"] == "ABS"].copy()
    A["down_share"] = A["n_down"] / A["n_genes"]
    A["lab"] = A["label"] + "/" + A["status"]
    A.to_csv(OUT / "k2_abs_composition.csv", index=False)
    k = R_abs[R_abs["trait"] == "KOLLER_EUR"].sort_values("Coefficient_P_value")
    top = [n for n in k["Name"] if not is_focus(k.set_index("Name").loc[n, "label"])][:6]
    foc = [n for n in k["Name"] if is_focus(k.set_index("Name").loc[n, "label"])]
    di = A.set_index("set_key")
    say("[COMPOSITION] top non-focus (KOLLER_EUR): " + "; ".join("%s n=%d down %.0f%%" % (short(di.loc[n, "lab"]), di.loc[n, "n_genes"], 100 * di.loc[n, "down_share"]) for n in top if n in di.index))
    say("[COMPOSITION] focus: " + "; ".join("%s n=%d down %.0f%%" % (short(di.loc[n, "lab"]), di.loc[n, "n_genes"], 100 * di.loc[n, "down_share"]) for n in foc if n in di.index))
    endo = A[A["label"].astype(str).str.match(ENDOTHELIAL_RE)]
    oth = A[~A["label"].astype(str).str.match(ENDOTHELIAL_RE)]
    say("[COMPOSITION] down share: endothelial sets %s | other sets median %.0f%% (n=%d)" % (
        ", ".join("%s %.0f%%" % (short(r["lab"]), 100 * r["down_share"]) for _, r in endo.iterrows()), 100 * oth["down_share"].median(), len(oth)))
    jac = []
    for a in foc:
        for b in top:
            sa, sb = set(sets[a]), set(sets[b])
            jac.append({"focus": a, "other": b, "jaccard": len(sa & sb) / max(1, len(sa | sb))})
    J = pd.DataFrame(jac)
    if len(J):
        J.to_csv(OUT / "k2_jaccard.csv", index=False)
        say("[COMPOSITION] Jaccard focus x top non-focus: median %.2f, max %.2f" % (J["jaccard"].median(), J["jaccard"].max()))
    c1 = bool(len(endo)) and bool((endo["down_share"] > 0.5).all())
    SUM["results"]["k2"] = {"endothelial_down_share": dict(zip(endo["lab"], endo["down_share"].round(3))), "all_endothelial_down_gt_half": c1}
    return c1


def h_abs_verdict(c1, R_up, R_down):
    c2 = c3 = None
    if R_up is not None:
        c2 = (len(rank_of(R_up, "KOLLER_EUR", ENDOTHELIAL_RE)) == 0) and (len(rank_of(R_up, "KOLLER_EUR", r"^dStromal")) > 0)
    if R_down is not None:
        c3 = len(rank_of(R_down, "KOLLER_EUR", ENDOTHELIAL_RE)) > 0
    ok = bool(c1) and bool(c2) and bool(c3)
    yn = lambda x: "NA" if x is None else ("YES" if x else "NO")
    say("[H_ABS_VERDICT] endothelial ABS sets down-dominated: %s | UP top 6: no endothelial + dStromal present: %s | DOWN top 6 has endothelial: %s -> %s" % (
        yn(c1), yn(c2), yn(c3), "SUPPORTED" if ok else "NOT SUPPORTED"))
    SUM["results"]["h_abs"] = {"c1": c1, "c2": c2, "c3": c3, "supported": ok}


# ----------------------------------------------------------------------------
# K4 conditional resolution (joint S-LDSC: A given B and B given A)
# ----------------------------------------------------------------------------
def conditional(fam, R):
    if fam is None or R is None:
        return
    say("")
    say("-" * 72)
    say("[SECTION] K4 joint S-LDSC, family %s (cell: p(focus | other) / p(other | focus))" % fam["tag"])
    k = R[R["trait"] == "KOLLER_EUR"].sort_values("Coefficient_P_value")
    lab = dict(zip(k["Name"], k["label"]))
    focus = [n for n in k["Name"] if is_focus(lab[n])]
    cols = [n for n in k["Name"] if not is_focus(lab[n])][:6]
    if not focus or not cols:
        say("[K4] no focus or no column sets")
        return
    pairs = [(a, b) for a in focus for b in cols] + [(b, a) for a in focus for b in cols]
    pref = lambda n: str(fam["dir"] / "sets" / n / (n + "."))
    say("[K4] columns: " + ", ".join("%d=%s" % (i + 1, short(lab[b], 18)) for i, b in enumerate(cols)))
    out_rows = []
    for t in TRAITS_COND:
        if t not in set(R["trait"]):
            continue
        out = OUT.parent / "cache" / ("%s__COND_%s" % (t, fam["family"]))
        res = Path(str(out) + ".cell_type_results.txt")
        if not res.exists():
            ldcts = fam["dir"] / ("COND_%s.ldcts" % fam["family"])
            with open(ldcts, "w") as f:
                for a, b in pairs:
                    f.write("%s__GIVEN__%s\t%s,%s,%s\n" % (a, b, pref(a), pref("ALLGENES"), pref(b)))
            cmd = NICE + [str(LDSC_PY), str(LDSC), "--h2-cts", str(MUNGED / (t + ".sumstats.gz")), "--ref-ld-chr", BASELINE,
                          "--ref-ld-chr-cts", str(ldcts), "--w-ld-chr", WEIGHTS, "--out", str(out)]
            t0 = time.time()
            rc, o, e = sh(cmd, timeout=8 * 3600, env=L2_ENV)
            if rc != 0 or not res.exists():
                say("    [K4_FAIL] %s rc=%d %s" % (t, rc, (o + e).strip()[-300:]))
                continue
            say("    [K4] %s: %d joint fits (%.0f min)" % (t, len(pairs), (time.time() - t0) / 60))
        C = pd.read_csv(res, sep="\t")
        ab = C["Name"].str.split("__GIVEN__", expand=True)
        cp = {(a, b): p for a, b, p in zip(ab[0], ab[1], C["Coefficient_P_value"])}
        d = R[R["trait"] == t].set_index("Name")["Coefficient_P_value"].to_dict()
        say("  %s" % t)
        for a in focus:
            cells = []
            for b in cols:
                pab, pba = cp.get((a, b), np.nan), cp.get((b, a), np.nan)
                pa, pb = d.get(a, np.nan), d.get(b, np.nan)
                if pa < 0.05 and pb < 0.05:
                    oc = {(True, True): "both", (True, False): "focus only", (False, True): "other only", (False, False): "neither"}[(pab < 0.05, pba < 0.05)]
                else:
                    oc = "not eligible"
                out_rows.append({"family": fam["tag"], "trait": t, "focus": lab[a], "other": lab[b], "p_focus": pa, "p_other": pb,
                                 "p_focus_given_other": pab, "p_other_given_focus": pba, "outcome": oc})
                cells.append("%.2f/%.2f" % (pab, pba))
            say("    %-22s p=%.3f | %s" % (short(lab[a]), d.get(a, np.nan), "  ".join(cells)))
        O = pd.DataFrame([r for r in out_rows if r["trait"] == t])
        E = O[O["outcome"] != "not eligible"]
        if not len(E):
            say("    [RESOLUTION_VERDICT] %s %s: no pair with both marginal p < 0.05" % (fam["tag"], t))
            continue
        fs = [short(f) for f, g in E.groupby("focus") if (g["p_focus_given_other"] < 0.05).all()]
        os_ = [short(o) for o, g in E.groupby("other") if (g["p_other_given_focus"] < 0.05).all()]
        cnt = E["outcome"].value_counts().to_dict()
        say("    [RESOLUTION_VERDICT] %s %s: eligible pairs %d (%s) | focus surviving every other: %s | others surviving every focus: %s" % (
            fam["tag"], t, len(E), ", ".join("%s %d" % (k2, v) for k2, v in cnt.items()), ", ".join(fs) or "none", ", ".join(os_) or "none"))
        SUM["results"]["k4_%s_%s" % (fam["tag"], t)] = {"eligible": len(E), "outcomes": cnt, "focus_surviving": fs, "others_surviving": os_}
    if out_rows:
        pd.DataFrame(out_rows).to_csv(OUT / ("k4_joint_%s.csv" % fam["tag"]), index=False)


# ----------------------------------------------------------------------------
# K6 comparator informativeness and Koller reading per variant
# ----------------------------------------------------------------------------
def koller_reading(R, tag):
    k = R[R["trait"] == "KOLLER_EUR"]
    ds = k[k["label"].astype(str).str.contains("dStromal")]
    sx = k[k["label"].astype(str).str.contains("SOX9_functionalis")]
    rep = bool((ds["Coefficient_P_value"] < 0.05).any() and (sx["Coefficient_P_value"] < 0.05).any())
    corr = bool((pd.concat([ds, sx])["q_bh"] < 0.05).any())
    hits = pd.concat([ds, sx])
    hits = list(hits.loc[hits["Coefficient_P_value"] < 0.05, "Name"])
    return rep, corr, hits


def comparators(R_abs, R_up, R_down):
    say("")
    say("-" * 72)
    say("[SECTION] K6 comparator informativeness")
    rep, corr, hits = koller_reading(R_abs, "P_ABS")
    e = R_abs[(R_abs["trait"] == "ENDO_ALL") & (R_abs["Name"].isin(hits))]
    if len(e):
        say("[K6] Koller-reproduced sets in FinnGen ENDO_ALL: " + "; ".join("%s p=%.3f" % (short(r["label"]), r["Coefficient_P_value"]) for _, r in e.iterrows()))
    endo_null = (bool((e["Coefficient_P_value"] >= 0.05).all()) if len(e) else None) if hits else None
    pc = {}
    for tag, R in (("ABS", R_abs), ("UP", R_up), ("DOWN", R_down)):
        if R is None:
            continue
        for t, pat in (("HYPERTENSION", VASC_RE), ("LEIOMYOMA", LEIO_RE)):
            d = R[(R["trait"] == t) & (R["label"].astype(str).str.match(pat))].sort_values("Coefficient_P_value")
            if len(d):
                r = d.iloc[0]
                pc["%s_%s" % (t, tag)] = bool(r["q_bh"] < 0.05)
                say("[K6] positive control %-12s %-4s best %s p=%.1e q=%.2f -> %s" % (t, tag, short(r["label"]), r["Coefficient_P_value"], r["q_bh"],
                                                                                   "DETECTED" if r["q_bh"] < 0.05 else "not detected"))
        e2 = R[(R["trait"] == "ENDO_ALL")].sort_values("Coefficient_P_value")
        if len(e2):
            say("[K6] ENDO_ALL %-4s top: %s" % (tag, "; ".join("%s %.0e/%.2f" % (short(r["label"]), r["Coefficient_P_value"], r["q_bh"]) for _, r in e2.head(3).iterrows())))
    inform = (endo_null is False) and any(pc.values())
    say("[SPECIFICITY_VERDICT] FinnGen endometriosis also null on the reproduced sets: %s | any S-LDSC positive control detected: %s -> comparators informative for specificity: %s" % (
        "NA (no reproduced sets)" if endo_null is None else ("YES" if endo_null else "NO"), "YES" if any(pc.values()) else "NO", "YES" if inform else "NO"))
    SUM["results"]["k6"] = {"endo_all_null": endo_null, "positive_controls": pc, "informative": inform}


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def step(name, fn, *a):
    try:
        return fn(*a)
    except Exception as e:
        import traceback
        say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
        say(traceback.format_exc()[-1500:])
        return None


def main():
    global BASELINE, WEIGHTS
    say("=" * 72)
    say("PAPER 9 | v0.0.31 | Koller pipeline: selection rule, joint resolution, sensitivity | " + TS)
    say("=" * 72)
    t0 = time.time()
    BASELINE = prefix_of("baseline_v1.2")
    WEIGHTS = prefix_of("1000G_Phase3_weights_hm3_no_MHC")
    sym2 = gene_map()
    step("K1", cohort_audit)

    say("")
    say("-" * 72)
    say("[SECTION] K3 S-LDSC families on the v030 cells (cell entries: p / BH q)")
    r = step("markers P", markers, sym2, "P")
    if r is None:
        return finish(t0)
    man, sets = r
    reuse = same_as_v030(man, sets)
    say("[REUSE] ABS sets identical to v030: %s" % ("YES (v030 LD scores and results reused)" if reuse else "NO (recomputed)"))
    fa = step("P_ABS", family, "P", "ABS", man, sets, sym2, TRAITS, reuse) or (None, None)
    fu = step("P_UP", family, "P", "UP", man, sets, sym2, TRAITS) or (None, None)
    fd = step("P_DOWN", family, "P", "DOWN", man, sets, sym2, TRAITS) or (None, None)
    if fa[1] is not None:
        c1 = step("K2", composition, man, sets, fa[1])
        step("H_ABS", h_abs_verdict, c1, fu[1], fd[1])
        step("K4_ABS", conditional, fa[0], fa[1])
    if fu[1] is not None:
        step("K4_UP", conditional, fu[0], fu[1])
    if fa[1] is not None:
        step("K6", comparators, fa[1], fu[1], fd[1])

    say("")
    say("-" * 72)
    say("[SECTION] K5 sensitivity (ABS family)")
    readings = {}
    if fa[1] is not None:
        readings["P"] = koller_reading(fa[1], "P_ABS")[:2]
    for vk in ("H", "N"):
        r = step("markers " + vk, markers, sym2, vk)
        if r is None:
            continue
        fv = step(vk + "_ABS", family, vk, "ABS", r[0], r[1], sym2, TRAITS_SENS) or (None, None)
        if fv[1] is not None:
            readings[vk] = koller_reading(fv[1], vk + "_ABS")[:2]
    say("[KOLLER_SENS_VERDICT] " + " | ".join("%s: nominal %s, BH %s" % (vk, "YES" if a else "NO", "YES" if b else "NO") for vk, (a, b) in readings.items()))
    SUM["results"]["k5"] = readings
    finish(t0)


def finish(t0):
    with open(OUT / "p9_v031_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.21 before this run")
    say("[DONE] PAPER9_V031_COMPLETE")
    _fh.close()


BASELINE = None
WEIGHTS = None

if __name__ == "__main__":
    main()
