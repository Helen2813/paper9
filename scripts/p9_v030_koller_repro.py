#!/usr/bin/env python3
# Paper 9 | v0.0.30 | A4: Koller et al. 2026 HECA S-LDSC, method-matched reproduction (rules CHANGELOG v0.2.20):
#   HECA cases / controls separately, each cell type vs all others within the endometrial or immune compartment,
#   |log2FC| > 1.2 and adjusted p < 0.01 (scanpy Wilcoxon), +-100 kb annotations, S-LDSC (LDSC-SEG layout, v012 code),
#   plus MAGMA gene-set analysis of the same sets. Traits KOLLER_EUR (primary), GCST90205183, HYPERTENSION, T2D, LEIOMYOMA.
# Runs in ~/Desktop/paper9/.venv, low-heat (1 worker, 2 BLAS threads, nice 10). Long: LD scores ~4-8 h (resumable).

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
GS_DIR = PROC / "genesets" / "heca_koller_v030"
LD_DIR = PROC / "ldscores" / "heca_koller_v030"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "koller_repro" / ("p9_v030_" + TS)
LOG = REPO / "logs" / ("p9_v030_koller_repro_" + TS + ".log")
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
TRAITS = ["KOLLER_EUR", "GCST90205183", "HYPERTENSION", "T2D", "LEIOMYOMA"]

for d in (OUT, LOG.parent, GS_DIR, LD_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V030_KOLLER_REPRO", "timestamp_utc": TS, "gates": {}, "results": {}}


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
    say("[ANNOT] %d columns x 22 chromosomes written" % len(keys))


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
    say("[M_CHECK] chr%d SNPs with MAF>0.05: ours=%d baseline(base annotation)=%.0f %s" % (
        c, n_common, ref_first, "MATCH" if abs(n_common - ref_first) < 0.5 else "DIFFERENT"))
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
            say("[L2_PROBE] chr22 format=%s rc=%d %s %.1f min (rc=1 expected: port crashes after LD scores, M rebuilt)" % ("thin" if thin else "full", rc, "OK" if ok else "FAILED",
                                                                 (time.time() - t0) / 60))
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
    say("[L2] BLAS threads per worker=%s" % BLAS_THREADS)
    say("[L2] chromosomes to compute: %s (workers=%d, format=%s)" % (todo or "none", N_WORKERS, "thin" if thin else "full"))
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
            say("    [L2] chr%-2d rc=%d %s %.0f min" % (c, p.returncode, "OK (M rebuilt)" if ok else "FAILED", (time.time() - t0) / 60))
            if not ok:
                logf = LD_DIR / ("multi.%d.log" % c)
                txt = logf.read_text(errors="ignore") if logf.exists() else ""
                for line in txt.strip().splitlines()[-8:]:
                    say("       " + line[:200])
        running = still
    ok_all = all((LD_DIR / ("multi.%d.l2.ldscore.gz" % c)).exists() for c in range(1, 23))
    if ok_all:
        for c in range(1, 23):
            write_M(c)
    SUM["gates"]["LDSCORES"] = ok_all
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
    say("[SPLIT] per-set LD score files written for %d sets; rows kept (HapMap3) %d of %d" % (len(keys), n_kept, n_all))


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
        say("    [SEG] %-16s %-22s rc=%d %.0f min" % (trait, fam["family"], rc, (time.time() - t0) / 60))
        if rc != 0 or not res.exists():
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
# 1. HECA markers (cases and controls separately; endometrial and immune separately)
# ----------------------------------------------------------------------------
def markers(sym2):
    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp
    man_f = GS_DIR / "manifest.csv"
    if man_f.exists():
        man = pd.read_csv(man_f)
        sets = {r["set_key"]: open(GS_DIR / (r["set_key"] + ".txt")).read().split() for _, r in man.iterrows()}
        say("[MARKERS] cached: %d sets" % len(man))
        return man, sets
    f = h5py.File(HECA_FILE, "r")
    cols = [k for k in f["obs"].keys() if k not in ("_index", "__categories")]
    dataset_col, dataset_vals = None, None
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
                dataset_col, dataset_vals = c, hit
                break
    donor = obs_cat(f, DONOR_COL)
    ctype = obs_cat(f, TYPE_COL)
    path = obs_cat(f, PATHOLOGY_COL)
    stage = obs_cat(f, STAGE_COL)
    say("[HECA] pathology categories: " + " | ".join("%s=%d" % (k, n) for k, n in pd.Series(path).value_counts().items()))
    keep = pd.notna(donor) & pd.notna(ctype) & (stage != EXCLUDE_STAGE)
    if dataset_col:
        dv = obs_cat(f, dataset_col)
        keep &= np.isin(dv, dataset_vals)
        say("[HECA] dataset column %s values %s -> cells %d" % (dataset_col, dataset_vals, int(keep.sum())))
    else:
        say("[HECA] no dataset column naming E-MTAB-14039 -> all HECA donors (non-hormone)")
    status = np.full(len(path), None, dtype=object)
    status[path == "C"] = "CTL"
    is_case = np.array([bool(re.match(r"^E", str(x))) or "endometriosis" in str(x).lower() for x in path])
    status[is_case] = "ENDO"
    keep &= pd.notna(status)
    for s in ("CTL", "ENDO"):
        m = keep & (status == s)
        say("[HECA] status %s: cells %d donors %d" % (s, int(m.sum()), pd.Series(donor[m]).nunique()))
    labels = sorted(pd.unique(ctype[keep]))
    comp = {l: ("immune" if IMMUNE_RE.search(l) and not NOT_IMMUNE_RE.search(l) else "endometrial") for l in labels}
    say("[COMPARTMENT] immune: " + ", ".join(l for l in labels if comp[l] == "immune"))
    say("[COMPARTMENT] endometrial: " + ", ".join(l for l in labels if comp[l] == "endometrial"))
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
            if len(ok_types) < 2:
                say("[MARKERS] %s %s: < 2 types, skipped" % (s, cp))
                continue
            plan[(s, cp)] = (np.sort(np.concatenate(pick)), ok_types)
    want = np.unique(np.concatenate([v[0] for v in plan.values()]))
    t0 = time.time()
    rowmat = {}
    step = 20000
    n = len(indptr) - 1
    wpos = 0
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
    say("[MARKERS] read %d sampled cells in one pass (%.1f min)" % (len(rowmat), (time.time() - t0) / 60))
    rows, sets = [], {}
    for (s, cp), (allidx, ok_types) in plan.items():
            t0 = time.time()
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
            say("[MARKERS] %s %-11s types=%d cells=%d genes=%d (%.1f min)" % (s, cp, len(ok_types), len(allidx), A.n_vars, (time.time() - t0) / 60))
            for l in ok_types:
                df = sc.get.rank_genes_groups_df(A, group=l)
                for mode, sel in (("ABS", (df["logfoldchanges"].abs() > LFC) & (df["pvals_adj"] < PADJ)),
                                  ("UP", (df["logfoldchanges"] > LFC) & (df["pvals_adj"] < PADJ))):
                    g = list(df.loc[sel, "names"].astype(str))
                    key = "%s_%s_%s_%s" % (mode, s, cp[:3], safe(l))[:60]
                    rows.append({"set_key": key, "mode": mode, "status": s, "compartment": cp, "label": l, "n_genes": len(g),
                                 "kept": len(g) >= MIN_SET})
                    if len(g) >= MIN_SET:
                        sets[key] = g
            allg = list(A.var_names)
            sets.setdefault("ALLGENES", allg)
    man = pd.DataFrame(rows)
    man = pd.concat([man, pd.DataFrame([{"set_key": "ALLGENES", "mode": "CTRL", "status": "-", "compartment": "-", "label": "all tested genes",
                                         "n_genes": len(sets["ALLGENES"]), "kept": True}])], ignore_index=True)
    for k, g in sets.items():
        (GS_DIR / (k + ".txt")).write_text("\n".join(g) + "\n")
    man[man["kept"]].to_csv(man_f, index=False)
    man.to_csv(GS_DIR / "manifest_all.csv", index=False)
    say("[MARKERS] sets kept %d of %d (>= %d genes); ABS median size %d" % (
        int(man["kept"].sum()), len(man), MIN_SET, int(man.loc[man["mode"] == "ABS", "n_genes"].median())))
    return man[man["kept"]].reset_index(drop=True), sets


# ----------------------------------------------------------------------------
# 2. MAGMA gene-set analysis (cheap complement)
# ----------------------------------------------------------------------------
def magma_sets(man, sets, sym2):
    say("")
    say("-" * 72)
    say("[SECTION] MAGMA gene-set analysis (competitive, condition-hide=ALLGENES)")
    sf = GS_DIR / "magma_sets.txt"
    with open(sf, "w") as o:
        for k, g in sets.items():
            ids = [sym2[x][0] for x in g if x in sym2]
            o.write(k + " " + " ".join(ids) + "\n")
    lab = man.set_index("set_key")["label"].to_dict()
    rows = []
    for t in TRAITS:
        raw = MAGMA_OUT / (t + "_35_10.genes.raw")
        if not raw.exists():
            continue
        outp = OUT / ("magma_sets_" + t)
        gsa = Path(str(outp) + ".gsa.out")
        for model in (["--model", "condition-hide=ALLGENES"], []):
            if gsa.exists():
                gsa.unlink()
            rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--set-annot", sf] + model + ["--out", outp], timeout=3600, env=L2_ENV)
            if gsa.exists():
                if not model:
                    say("    [MAGMA_SET] %s: condition on ALLGENES not accepted; unconditioned run used" % t)
                break
        if not gsa.exists():
            say("    [MAGMA_SET_FAIL] %s %s" % (t, (o + e).strip()[-200:]))
            continue
        lines = [l for l in open(gsa) if not l.startswith("#") and l.strip()]
        hdr = lines[0].split()
        for l in lines[1:]:
            d = dict(zip(hdr, l.split()))
            name = d.get("FULL_NAME", d.get("VARIABLE"))
            try:
                rows.append({"trait": t, "set_key": name, "label": lab.get(name, name), "p": float(d["P"]), "beta": float(d["BETA"])})
            except (KeyError, ValueError):
                pass
    R = pd.DataFrame(rows)
    if not len(R):
        return
    R["mode"] = R["set_key"].str.split("_").str[0]
    R["status"] = R["set_key"].str.split("_").str[1]
    out = []
    for (t, mode), g in R.groupby(["trait", "mode"]):
        g = g.copy()
        g["q"] = bh(g["p"].values)
        out.append(g)
        g = g.sort_values("p")
        say("    [MSET] %-13s %-3s top: %s" % (t[:13], mode, "; ".join("%s/%s p=%.1e q=%.2f" % (r["label"][:18], r["status"], r["p"], r["q"])
                                                                   for _, r in g.head(4).iterrows())))
        foc = g[g["label"].str.contains("dStromal|SOX9_functionalis", regex=True)]
        say("    [MSET] %-13s %-3s focus: %s" % (t[:13], mode, "; ".join("%s/%s p=%.3f" % (r["label"][:18], r["status"], r["p"]) for _, r in foc.iterrows())))
    pd.concat(out).to_csv(OUT / "magma_genesets.csv", index=False)


# ----------------------------------------------------------------------------
# 3. S-LDSC (ABS family)
# ----------------------------------------------------------------------------
def sldsc(man, sets, sym2):
    global BASELINE, WEIGHTS
    say("")
    say("-" * 72)
    say("[SECTION] S-LDSC (+-100 kb, 1 cM, baseline v1.2, --h2-cts with all-genes control); family ABS")
    fam_keys = [k for k in man.loc[man["mode"] == "ABS", "set_key"]] + ["ALLGENES"]
    all_sets = {k: [g for g in sets[k] if g in sym2] for k in fam_keys}
    keys_file = LD_DIR / "columns.txt"
    if keys_file.exists() and keys_file.read_text().split() != list(all_sets.keys()):
        say("[ANNOT] set list changed -> rebuilding")
        for p in LD_DIR.glob("*"):
            if p.is_file():
                p.unlink()
    write_annotations(all_sets, sym2)
    if not run_l2():
        say("[STOP] LD scores incomplete; re-run to resume")
        return
    split_ldscores(list(all_sets.keys()))
    BASELINE = prefix_of("baseline_v1.2")
    WEIGHTS = prefix_of("1000G_Phase3_weights_hm3_no_MHC")
    fam = {"family": "KOLLER_ABS", "sets": list(all_sets.keys()), "control": "ALLGENES",
           "manifest": man[man["set_key"].isin(fam_keys)][["set_key", "label"]].assign(
               label=lambda d: d["label"] + "/" + d["set_key"].str.split("_").str[1])}
    allres = []
    for t in TRAITS:
        if not (MUNGED / (t + ".sumstats.gz")).exists():
            say("[SEG] %s sumstats missing" % t)
            continue
        df = ldsc_seg(t, fam)
        if df is None:
            continue
        allres.append(df)
        say("[SEG_TOP] %s" % t)
        for _, r in df.head(6).iterrows():
            say("    %-34s coef=%9.2e p=%.2e q=%.3f" % (str(r["label"])[:34], r["Coefficient"], r["Coefficient_P_value"], r["q_bh"]))
        foc = df[df["label"].astype(str).str.contains("dStromal|SOX9_functionalis", regex=True)]
        say("    [FOCUS] " + "; ".join("%s p=%.3f q=%.2f" % (str(r["label"])[:26], r["Coefficient_P_value"], r["q_bh"]) for _, r in foc.iterrows()))
    if not allres:
        return
    A = pd.concat(allres)
    A.to_csv(OUT / "sldsc_koller_abs.csv", index=False)
    k = A[A["trait"] == "KOLLER_EUR"]
    ds = k[k["label"].astype(str).str.contains("dStromal")]
    sx = k[k["label"].astype(str).str.contains("SOX9_functionalis")]
    rep = (ds["Coefficient_P_value"] < 0.05).any() and (sx["Coefficient_P_value"] < 0.05).any()
    corr = (pd.concat([ds, sx])["q_bh"] < 0.05).any()
    hits = pd.concat([ds, sx])
    hits = list(hits.loc[hits["Coefficient_P_value"] < 0.05, "Name"])
    ctrl = A[(A["trait"].isin(["HYPERTENSION", "T2D"])) & (A["Name"].isin(hits)) & (A["Coefficient_P_value"] < 0.05)]
    spec = rep and len(ctrl) == 0
    say("[KOLLER_REPRO_VERDICT] nominal dStromal + SOX9 functionalis reproduced: %s | survives BH: %s | endometriosis-specific: %s" % (
        "YES" if rep else "NO", "YES" if corr else "NO", "YES" if spec else "NO"))
    SUM["results"]["koller_repro"] = {"nominal": rep, "bh": corr, "specific": spec}


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.30 | Koller 2026 HECA S-LDSC, method-matched reproduction | " + TS)
    say("=" * 72)
    t0 = time.time()
    sym2 = gene_map()
    try:
        man, sets = markers(sym2)
        magma_sets(man, sets, {g: v for g, v in sym2.items()})
        sldsc(man, sets, sym2)
    except Exception as e:
        import traceback
        say("[STEP_ERROR] %s: %s" % (type(e).__name__, str(e)[:300]))
        say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v030_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.20 before this run")
    say("[DONE] PAPER9_V030_COMPLETE")
    _fh.close()


BASELINE = None
WEIGHTS = None

if __name__ == "__main__":
    main()
