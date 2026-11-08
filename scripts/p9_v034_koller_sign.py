#!/usr/bin/env python3
# Paper 9 | v0.0.34 | A5 follow-up on the v031 LD scores (rules CHANGELOG v0.2.23); no new LD scores:
#   K7 sign decomposition: for every cell type, joint S-LDSC of its UP and DOWN sets (UP | DOWN and DOWN | UP);
#   K8 decidual stroma (UP sets) vs every vascular / perivascular / smooth-muscle UP set, both directions.
# Runs in ~/Desktop/paper9/.venv, low-heat (nice 10, 2 BLAS threads). ~30-45 min; cached fits are reused.

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
OUT = REPO / "results" / "koller_repro" / ("p9_v034_" + TS)
LOG = REPO / "logs" / ("p9_v034_koller_sign_" + TS + ".log")
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
SUM = {"task": "PAPER9_V034_KOLLER_SIGN", "timestamp_utc": TS, "gates": {}, "results": {}}


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
# Helpers
# ----------------------------------------------------------------------------
VASC_LABELS = ["Arterial", "Venous", "ePV_1a", "ePV_1b", "ePV_2", "mPV", "uSMCs"]
TRAITS_K = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL"]


def short(label, n=22):
    s = str(label)
    s = s.replace("SOX9_functionalis_", "SOX9f_").replace("dStromal_", "dS_").replace("Glandular_secretory", "Gland_sec")
    return s[:n]


def load_fam(mode, man):
    keys = list(man.loc[man["mode"] == mode, "set_key"]) + ["ALLGENES"]
    d, fname = (V030_LD, "KOLLER_ABS") if mode == "ABS" else (LD_BASE / ("P_" + mode), "K31_P_" + mode)
    if not (d / "split.done").exists():
        say("[STOP] LD scores for %s not found (%s); run v031 first" % (mode, d))
        return None
    sub = man[man["set_key"].isin(keys)][["set_key", "label"]]
    sub = sub.assign(label=sub["label"] + "/" + sub["set_key"].str.split("_").str[1])
    return {"family": fname, "sets": keys, "control": "ALLGENES", "manifest": sub, "dir": d, "tag": "P_" + mode}


def pref(fam, n):
    return str(fam["dir"] / "sets" / n / (n + "."))


def marginal(trait, fam):
    global LD_DIR
    LD_DIR = fam["dir"]
    df = ldsc_seg(trait, fam)
    return {} if df is None else dict(zip(df["Name"], df["Coefficient_P_value"]))


def run_joint(trait, name, lines):
    out = OUT.parent / "cache" / ("%s__%s" % (trait, name))
    out.parent.mkdir(parents=True, exist_ok=True)
    res = Path(str(out) + ".cell_type_results.txt")
    if not res.exists():
        ldcts = LD_BASE / ("%s.ldcts" % name)
        ldcts.write_text("".join(l + "\n" for l in lines))
        cmd = NICE + [str(LDSC_PY), str(LDSC), "--h2-cts", str(MUNGED / (trait + ".sumstats.gz")), "--ref-ld-chr", BASELINE,
                      "--ref-ld-chr-cts", str(ldcts), "--w-ld-chr", WEIGHTS, "--out", str(out)]
        t0 = time.time()
        rc, o, e = sh(cmd, timeout=8 * 3600, env=L2_ENV)
        if rc != 0 or not res.exists():
            say("    [JOINT_FAIL] %s %s rc=%d %s" % (trait, name, rc, (o + e).strip()[-300:]))
            return None
        say("    [JOINT] %s %s: %d fits (%.0f min)" % (trait, name, len(lines), (time.time() - t0) / 60))
    C = pd.read_csv(res, sep="\t")
    return dict(zip(C["Name"], C["Coefficient_P_value"]))


def fp(x):
    return "NA" if x is None or not np.isfinite(x) else ("%.0e" % x if x < 0.001 else "%.3f" % x)


# ----------------------------------------------------------------------------
# K7 sign decomposition (UP | DOWN and DOWN | UP for each cell type and status)
# ----------------------------------------------------------------------------
def k7(man, fA, fU, fD):
    say("")
    say("-" * 72)
    say("[SECTION] K7 sign decomposition of the |log2FC| annotations (joint S-LDSC of a type's UP and DOWN sets)")
    key = {(r["mode"], r["status"], r["compartment"], r["label"]): r["set_key"] for _, r in man.iterrows()}
    types = []
    for _, r in man[man["mode"] == "ABS"].iterrows():
        ku, kd = key.get(("UP", r["status"], r["compartment"], r["label"])), key.get(("DOWN", r["status"], r["compartment"], r["label"]))
        if ku and kd:
            types.append((r["set_key"], ku, kd, "%s/%s" % (r["label"], r["status"])))
    lines = []
    for ka, ku, kd, lab in types:
        lines.append("%s__GIVEN__%s\t%s,%s,%s" % (ku, kd, pref(fU, ku), pref(fU, "ALLGENES"), pref(fD, kd)))
        lines.append("%s__GIVEN__%s\t%s,%s,%s" % (kd, ku, pref(fD, kd), pref(fU, "ALLGENES"), pref(fU, ku)))
    say("[K7] cell types with UP and DOWN sets: %d" % len(types))
    rows = []
    for t in TRAITS_K:
        cp = run_joint(t, "K7_SIGN", lines)
        if cp is None:
            continue
        mA, mU, mD = marginal(t, fA), marginal(t, fU), marginal(t, fD)
        for ka, ku, kd, lab in types:
            pa, pu, pd_ = mA.get(ka, np.nan), mU.get(ku, np.nan), mD.get(kd, np.nan)
            pud, pdu = cp.get("%s__GIVEN__%s" % (ku, kd), np.nan), cp.get("%s__GIVEN__%s" % (kd, ku), np.nan)
            if pa < 0.05:
                cls = {(True, False): "UP-carried", (False, True): "DOWN-carried", (True, True): "both", (False, False): "neither"}[(pud < 0.05, pdu < 0.05)]
            else:
                cls = "ABS not nominal"
            rows.append({"trait": t, "type": lab, "p_abs": pa, "p_up": pu, "p_down": pd_, "p_up_given_down": pud, "p_down_given_up": pdu, "class": cls})
        R = pd.DataFrame([r for r in rows if r["trait"] == t])
        E = R[R["class"] != "ABS not nominal"]
        cnt = E["class"].value_counts().to_dict()
        if t == "KOLLER_EUR":
            panel = R[R["type"].str.split("/").str[0].isin(VASC_LABELS) | R["type"].str.contains(FOCUS_RE)].sort_values("p_abs")
            for _, r in panel.iterrows():
                say("    %-22s ABS %-6s | UP %-6s (|DOWN %-6s) | DOWN %-6s (|UP %-6s) -> %s" % (
                    short(r["type"]), fp(r["p_abs"]), fp(r["p_up"]), fp(r["p_up_given_down"]), fp(r["p_down"]), fp(r["p_down_given_up"]), r["class"]))
        say("  [K7] %-12s ABS-nominal types %d: %s" % (t[:12], len(E), ", ".join("%s %d" % (k, v) for k, v in sorted(cnt.items())) or "-"))
        if t == "KOLLER_EUR":
            endo = E[E["type"].str.match(r"^(Venous|Arterial)/")]
            ds = E[E["type"].str.match(r"^dStromal")]
            c1 = bool(len(endo)) and bool((endo["class"] == "DOWN-carried").all())
            c2 = bool(len(ds)) and bool((ds["class"] == "UP-carried").all())
            say("[H_SIGN_VERDICT] KOLLER_EUR: endothelial ABS signal DOWN-carried (%d/%d): %s | decidual-stroma ABS signal UP-carried (%d/%d): %s -> %s" % (
                int((endo["class"] == "DOWN-carried").sum()), len(endo), "YES" if c1 else "NO", int((ds["class"] == "UP-carried").sum()), len(ds),
                "YES" if c2 else "NO", "SUPPORTED" if (c1 and c2) else "NOT SUPPORTED"))
            SUM["results"]["h_sign"] = {"endothelial_down": c1, "dstromal_up": c2}
        SUM["results"]["k7_" + t] = cnt
    if rows:
        pd.DataFrame(rows).to_csv(OUT / "k7_sign_decomposition.csv", index=False)


# ----------------------------------------------------------------------------
# K8 decidual stroma (UP) vs vascular / perivascular / smooth muscle (UP)
# ----------------------------------------------------------------------------
def k8(man, fU):
    say("")
    say("-" * 72)
    say("[SECTION] K8 decidual stroma vs vascular, perivascular and smooth-muscle programs (UP sets, joint S-LDSC)")
    up = man[man["mode"] == "UP"]
    lab = dict(zip(up["set_key"], up["label"] + "/" + up["status"]))
    rows_ = [k for k in up["set_key"] if str(lab[k]).startswith("dStromal")]
    cols = [k for k in up["set_key"] if str(lab[k]).split("/")[0] in VASC_LABELS]
    if not rows_ or not cols:
        say("[K8] no decidual-stroma or vascular UP sets")
        return
    lines = []
    for a in rows_:
        for b in cols:
            lines.append("%s__GIVEN__%s\t%s,%s,%s" % (a, b, pref(fU, a), pref(fU, "ALLGENES"), pref(fU, b)))
            lines.append("%s__GIVEN__%s\t%s,%s,%s" % (b, a, pref(fU, b), pref(fU, "ALLGENES"), pref(fU, a)))
    say("[K8] rows %d decidual-stroma sets x columns %d vascular sets: %s" % (len(rows_), len(cols), ", ".join(short(lab[b], 14) for b in cols)))
    out = []
    for t in TRAITS_K:
        cp = run_joint(t, "K8_UPVASC", lines)
        if cp is None:
            continue
        m = marginal(t, fU)
        for a in rows_:
            for b in cols:
                out.append({"trait": t, "dstromal": lab[a], "other": lab[b], "p_ds": m.get(a, np.nan), "p_other": m.get(b, np.nan),
                            "p_ds_given_other": cp.get("%s__GIVEN__%s" % (a, b), np.nan), "p_other_given_ds": cp.get("%s__GIVEN__%s" % (b, a), np.nan)})
        O = pd.DataFrame([r for r in out if r["trait"] == t])
        say("  %s" % t)
        for a in rows_:
            g = O[O["dstromal"] == lab[a]]
            lost = g[g["p_ds_given_other"] >= 0.05]
            surv = g[(g["p_other"] < 0.05) & (g["p_other_given_ds"] < 0.05)]
            say("    %-16s p=%-6s survives %2d/%d | lost given: %s | vascular surviving it: %s" % (
                short(lab[a], 16), fp(m.get(a, np.nan)), int((g["p_ds_given_other"] < 0.05).sum()), len(g),
                ", ".join("%s %.2f" % (short(r["other"], 12), r["p_ds_given_other"]) for _, r in lost.iterrows()) or "none",
                ", ".join(short(r["other"], 12) for _, r in surv.iterrows()) or "none"))
        pv = O[O["dstromal"].str.startswith("dStromal_mid") & O["other"].str.startswith("ePV_2")]
        if len(pv):
            sep = bool((pv["p_ds_given_other"] < 0.05).all())
            say("  [K8_VERDICT] %s dStromal_mid given ePV_2 (UP): %s -> G4 prediction (not separable) %s" % (
                t, "; ".join("%s|%s %.3f" % (short(r["dstromal"], 12), short(r["other"], 10), r["p_ds_given_other"]) for _, r in pv.iterrows()),
                "NOT CONFIRMED (separable in S-LDSC UP)" if sep else "CONFIRMED"))
            SUM["results"]["k8_" + t] = {"dsmid_separable_from_epv2": sep}
        else:
            say("  [K8_VERDICT] %s: no dStromal_mid x ePV_2 pair" % t)
    if out:
        pd.DataFrame(out).to_csv(OUT / "k8_dstromal_vs_vascular_up.csv", index=False)


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
    say("PAPER 9 | v0.0.34 | Koller pipeline: sign decomposition and decidual stroma vs vascular | " + TS)
    say("=" * 72)
    t0 = time.time()
    BASELINE = prefix_of("baseline_v1.2")
    WEIGHTS = prefix_of("1000G_Phase3_weights_hm3_no_MHC")
    mf = GS_BASE / "P" / "manifest.csv"
    if not mf.exists():
        say("[STOP] v031 marker sets not found; run v031 first")
        return finish(t0)
    man = pd.read_csv(mf)
    fA, fU, fD = load_fam("ABS", man), load_fam("UP", man), load_fam("DOWN", man)
    if fA and fU and fD:
        step("K7", k7, man, fA, fU, fD)
    if fU:
        step("K8", k8, man, fU)
    finish(t0)


def finish(t0):
    with open(OUT / "p9_v034_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.23 before this run")
    say("[DONE] PAPER9_V034_COMPLETE")
    _fh.close()


BASELINE = None
WEIGHTS = None

if __name__ == "__main__":
    main()
