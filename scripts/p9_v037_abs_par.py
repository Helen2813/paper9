#!/usr/bin/env python3
# Paper 9 | v0.0.37 | A5 simulation, rotation design (rules CHANGELOG v0.2.26): does a lineage-agnostic endometrial
#   parenchymal source (stroma + epithelium, no cell state) surface as "endothelium" under the |log2FC| marker rule?
#   Same design as v036 (real KOLLER_EUR z and correlations; marker sets and planted profile rotated together).
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). ~30-40 min.

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
OUT = REPO / "results" / "koller_repro" / ("p9_v037_" + TS)
LOG = REPO / "logs" / ("p9_v037_abs_par_" + TS + ".log")
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
SUM = {"task": "PAPER9_V037_ABS_PAR", "timestamp_utc": TS, "gates": {}, "results": {}}


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
# Fixed (CHANGELOG v0.2.26)
# ----------------------------------------------------------------------------
HECA_GS = PROC / "genesets" / "heca_v014"
WORK = INTERIM / "magma_v037"
WORK.mkdir(parents=True, exist_ok=True)
OBS_BETA_DSMID = 0.947      # KOLLER_EUR dStromal_mid beta (v024 base model), as in R4 (v027)
SIM_SEED = 20261107
N_SIM = 60
# scenario: (name, strength label, [(profile, strength), ...]); PAR = summed specificity over stromal + epithelial types
SCENARIOS = [("null", 0.0, []), ("PAR", 2.0, [("PAR", 2.0)]), ("PAR", 4.0, [("PAR", 4.0)]),
             ("PAR+dS_mid", 2.0, [("PAR", 2.0), ("dStromal_mid", 2.0)])]
PAR_RE = r"Stromal|Fibroblast|SOX9|Lumin|Gland|MUC5B|KRT5|iliated|HOXA13"
STROMAL_RE = r"^(dStromal|eStromal|Fibroblast)"
VASC = ["Arterial", "Venous", "ePV_1a", "ePV_1b", "ePV_2", "mPV", "uSMCs"]
ENDOTH = ["Arterial", "Venous"]


def group_of(label):
    b = str(label).split("/")[0]
    if b.startswith("dStromal"):
        return "dStromal"
    if b in ENDOTH:
        return "endothelial"
    if b in VASC:
        return "perivascular/SM"
    return "other"


def read_raw(trait):
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    lines = raw.read_text().splitlines()
    out = pd.read_csv(MAGMA_OUT / (trait + "_35_10.genes.out"), sep=r"\s+")
    zmap = dict(zip(out["GENE"].astype(str), out["ZSTAT"]))
    data_idx = [i for i, l in enumerate(lines) if l.strip() and not l.startswith("#")]
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


def magma_sets_run(raw, sf, tag):
    outp = WORK / tag
    gsa = Path(str(outp) + ".gsa.out")
    for model in (["--model", "condition-hide=ALLGENES"], []):
        if gsa.exists():
            gsa.unlink()
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--set-annot", sf] + model + ["--out", outp], timeout=3600, env=L2_ENV)
        if gsa.exists():
            break
    if not gsa.exists():
        say("    [MAGMA_FAIL] %s %s" % (tag, (o + e).strip()[-200:]))
        return None
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


def summarize(res, man):
    out = {}
    for mode in ("ABS", "UP", "DOWN"):
        m = man[man["mode"] == mode]
        keys = [k for k in m["set_key"] if k in res]
        if not keys:
            continue
        lab = dict(zip(m["set_key"], m["label"] + "/" + m["status"]))
        p = np.array([res[k] for k in keys])
        q = bh(p)
        o = np.argsort(p)
        top = lab[keys[o[0]]]
        g = np.array([group_of(lab[k]) for k in keys])
        out[mode] = {"top": top, "top_group": group_of(top), "any_q05": bool((q < 0.05).any()),
                     "endo_q05": bool(((q < 0.05) & (g == "endothelial")).any()), "ds_q05": bool(((q < 0.05) & (g == "dStromal")).any()),
                     "best_ds_rank": int(np.where(g[o] == "dStromal")[0][0]) + 1 if (g == "dStromal").any() else -1,
                     "best_endo_rank": int(np.where(g[o] == "endothelial")[0][0]) + 1 if (g == "endothelial").any() else -1}
    return out


def simulate(man, sets, sym2):
    say("")
    say("-" * 72)
    say("[SECTION] S-ABS rotation design (MAGMA gene-set analysis, condition-hide=ALLGENES; real KOLLER_EUR z, annotation layer rotated)")
    S = pd.read_csv(HECA_GS / "specificity_by_type.csv.gz", index_col=0)
    lines, data_idx, zi, genes, z0 = read_raw("KOLLER_EUR")
    n = len(genes)
    idx_of = {g: i for i, g in enumerate(genes)}
    ent = {sym2[g][0]: g for g in S.index if g in sym2}
    sym_of = np.array([ent.get(g) for g in genes], dtype=object)
    has = np.array([s is not None for s in sym_of])
    b_std = OBS_BETA_DSMID * float(S["dStromal_mid"].std())
    keys = list(man["set_key"])
    set_idx = {k: np.array(sorted({idx_of[sym2[x][0]] for x in sets[k] if x in sym2 and sym2[x][0] in idx_of}), dtype=int) for k in keys}
    say("[SIM] genes.raw %d genes (Z column %d); HECA profiles for %d genes; effect per SD of profile x1 = %.4f" % (n, zi, int(has.sum()), b_std))
    sf = WORK / "sets_rot.txt"
    simraw = WORK / "SIM_35_10.genes.raw"
    realraw = MAGMA_OUT / "KOLLER_EUR_35_10.genes.raw"

    def write_sets(off):
        with open(sf, "w") as o:
            for k in keys:
                ii = (set_idx[k] + off) % n
                o.write(k + " " + " ".join(genes[i] for i in ii) + "\n")

    write_sets(0)
    obs = magma_sets_run(realraw, sf, "observed")
    if obs:
        so = summarize(obs, man)
        say("[SIM] observed KOLLER_EUR (no rotation): " + " | ".join("%s top %s" % (m, short(v["top"])) for m, v in so.items()))
    rng = np.random.default_rng(SIM_SEED)
    rows = []
    t0 = time.time()
    par_cols = [c for c in S.columns if re.search(PAR_RE, c) and c not in ("Arterial", "Venous")]
    say("[SIM] PAR profile = sum of specificity over %d types: %s" % (len(par_cols), ", ".join(par_cols)))

    def std_profile(target):
        prof = np.zeros(n)
        if target == "PAR":
            v = S.loc[list(sym_of[has]), par_cols].sum(axis=1).to_numpy(dtype=float)
        else:
            v = S.loc[list(sym_of[has]), target].to_numpy(dtype=float)
        prof[has] = (v - v.mean()) / v.std()
        return prof

    for name, strength, comps in SCENARIOS:
        if any(c != "PAR" and c not in S.columns for c, _ in comps) or (any(c == "PAR" for c, _ in comps) and not par_cols):
            say("[SIM] %s: profile not available; skipped" % name)
            continue
        prof = None
        for c, st in comps:
            p_ = st * b_std * std_profile(c)
            prof = p_ if prof is None else prof + p_
        for rep in range(N_SIM):
            off = int(rng.integers(500, n - 500))
            write_sets(off)
            if prof is None:
                raw = realraw
            else:
                write_raw(lines, data_idx, zi, z0 + np.roll(prof, off), simraw)
                raw = simraw
            res = magma_sets_run(raw, sf, "sim")
            if not res:
                continue
            for mode, v in summarize(res, man).items():
                rows.append({"truth": name, "strength": strength, "rep": rep, "offset": off, "mode": mode, **v})
        say("    [SIM] truth %-12s x%.0f done (%.0f min)" % (name, strength, (time.time() - t0) / 60))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "s_abs_parenchymal.csv", index=False)
    say("")
    say("[SIM_SUM] share of replicates: top set group (dStromal / endothelial / perivascular+SM / other) | any q<0.05 | endothelial q<0.05 | dStromal q<0.05")
    for (name, strength), g0 in R.groupby(["truth", "strength"], sort=False):
        for mode in ("ABS", "UP", "DOWN"):
            g = g0[g0["mode"] == mode]
            if not len(g):
                continue
            sh_ = g["top_group"].value_counts(normalize=True)
            fq = g["top"].value_counts(normalize=True).head(3)
            say("    truth %-12s x%.0f %-4s n=%-3d top: dS %.2f endo %.2f pv/SM %.2f other %.2f | anyq %.2f endoq %.2f dSq %.2f | rank dS %d endo %d | most often: %s" % (
                name, strength, mode, len(g), sh_.get("dStromal", 0), sh_.get("endothelial", 0), sh_.get("perivascular/SM", 0), sh_.get("other", 0),
                g["any_q05"].mean(), g["endo_q05"].mean(), g["ds_q05"].mean(), int(g["best_ds_rank"].median()), int(g["best_endo_rank"].median()),
                ", ".join("%s %.2f" % (short(k, 16), v) for k, v in fq.items())))

    def share(name, strength, mode, col, val=None):
        g = R[(R["truth"] == name) & (R["strength"] == strength) & (R["mode"] == mode)]
        if not len(g):
            return np.nan
        return float((g[col] == val).mean()) if val is not None else float(g[col].mean())
    calib = {m: share("null", 0.0, m, "any_q05") for m in ("ABS", "UP", "DOWN")}
    e2, e4, e0 = share("PAR", 2.0, "ABS", "top_group", "endothelial"), share("PAR", 4.0, "ABS", "top_group", "endothelial"), share("null", 0.0, "ABS", "top_group", "endothelial")
    say("[H_PAR_VERDICT] null any q<0.05: ABS %.2f, UP %.2f, DOWN %.2f | ABS top endothelial: null %.2f, PAR x2 %.2f, PAR x4 %.2f -> H-PAR (x2 >= null + 0.10): %s" % (
        calib["ABS"], calib["UP"], calib["DOWN"], e0, e2, e4, "SUPPORTED" if (e2 - e0 >= 0.10) else "NOT SUPPORTED"))
    # signature of the observed data: ABS top endothelial and UP top stromal
    W = R.pivot_table(index=["truth", "strength", "rep"], columns="mode", values="top", aggfunc="first").reset_index()
    W["sig"] = W["ABS"].map(lambda x: group_of(x) == "endothelial") & W["UP"].map(lambda x: bool(re.match(STROMAL_RE, str(x))))
    say("[SIGNATURE] share of replicates with ABS top endothelial AND UP top stromal (observed KOLLER pattern): " + "; ".join(
        "%s x%.0f %.2f" % (a, b, g["sig"].mean()) for (a, b), g in W.groupby(["truth", "strength"], sort=False)))
    SUM["results"]["h_par"] = {"null_anyq": calib, "abs_endo_top": {"null": e0, "x2": e2, "x4": e4},
                               "signature": {"%s_x%.0f" % k: float(g["sig"].mean()) for k, g in W.groupby(["truth", "strength"])}}


def short(label, n=22):
    s = str(label)
    s = s.replace("SOX9_functionalis_", "SOX9f_").replace("dStromal_", "dS_").replace("Glandular_secretory", "Gland_sec")
    return s[:n]


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.37 | simulation, rotation design: parenchymal source vs |log2FC| marker sets | " + TS)
    say("=" * 72)
    t0 = time.time()
    mf = GS_BASE / "P" / "manifest.csv"
    if not mf.exists():
        say("[STOP] v031 marker sets not found; run v031 first")
    else:
        man = pd.read_csv(mf)
        sets = {k: open(GS_BASE / "P" / (k + ".txt")).read().split() for k in man["set_key"]}
        try:
            simulate(man, sets, gene_map())
        except Exception as e:
            import traceback
            say("[STEP_ERROR] %s: %s" % (type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v037_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.26 before this run")
    say("[DONE] PAPER9_V037_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
