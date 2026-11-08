#!/usr/bin/env python3
# Paper 9 | v0.0.7 | (1) MAGMA relink with an older macOS SDK
#                   (2) Gene sets from the v006 GSE266265 pseudobulk (CTL donors; rules CHANGELOG v0.2.3-v0.2.4)
#                   (3) Gene-set SNP annotations (+-100 kb) and partitioned LD scores (1000G EUR, hg19)
#                   (4) LDSC-SEG (--h2-cts) per trait; G4 primary test (decidual stroma, Holm)
# Runs in ~/Desktop/paper9/.venv. Re-run safe: finished chromosomes / traits are reused.

import datetime as dt
import gzip
import json
import math
import os
import re
import stat
import subprocess
import time
from pathlib import Path

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
MAGMA_SRC = TOOLS / "magma" / "source"
SDK_DIR = Path("/Library/Developer/CommandLineTools/SDKs")
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
REF = RAW / "ldsc_ref" / "unpacked"
PLINK_DIR = REF / "1000G_EUR_Phase3_plink"
HM3_LIST = RAW / "ldsc_ref" / "hm3_no_MHC.list.txt"
PB_DIR = PROC / "gse266265_pseudobulk"
GS_DIR = PROC / "genesets" / "gse266265_v007"
LD_DIR = PROC / "ldscores" / "gse266265_v007"
MUNGED = PROC / "munged"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "ldsc_seg" / ("p9_v007_ldsc_seg_" + TS)
LOG = REPO / "logs" / ("p9_v007_ldsc_seg_" + TS + ".log")

# Pre-registered rules (CHANGELOG v0.2.3 / v0.2.4)
CTL_LABEL = "CTL"
MIN_CELLS_PER_DONOR_GROUP = 10
MIN_DONORS_PER_TYPE = 3
MIN_CELLS_PER_TYPE = 50
MIN_CPM_ANY_TYPE = 1.0
TOP_FRACTION = 0.10
WINDOW_BP = 100_000
LD_WIND_CM = 1.0
N_WORKERS = 4
PHASE3 = {"proliferative": "proliferative", "periovulatory": "periovulatory",
          "secretory_early": "secretory", "secretory_mid": "secretory", "secretory_late": "secretory"}
G4_PRIMARY = ["dS1-myofibroblast", "dS2", "dS3", "dS4"]
TRAITS = ["GCST90205183", "ENDO_ALL", "LEIOMYOMA", "ADENOMYOSIS", "HYPERTENSION", "T2D",
          "ENDO_OVARY", "ENDO_PERITONEUM", "ENDO_DEEP"]

for d in (OUT, LOG.parent, GS_DIR, LD_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V007_LDSC_SEG", "timestamp_utc": TS, "gates": {}}


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def sh(cmd, cwd=None, timeout=None, env=None):
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout, env=env)
        return p.returncode, p.stdout, p.stderr
    except FileNotFoundError as e:
        return 127, "", str(e)
    except OSError as e:
        return 126, "", str(e)
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"


def safe(name):
    s = name.replace("$\\Phi$", "Phi").replace("$", "").replace("{", "").replace("}", "").replace("\\", "")
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def holm(p):
    p = np.asarray(p, dtype=float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m)
    run = 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (m - r) * p[i]))
        adj[i] = run
    return adj


# ----------------------------------------------------------------------------
# 1. MAGMA relink with older SDKs
# ----------------------------------------------------------------------------
def magma_relink():
    say("")
    say("-" * 72)
    say("[SECTION] MAGMA RELINK")
    built = MAGMA_SRC / "magma"
    sdks = sorted([p for p in SDK_DIR.glob("*.sdk")], key=lambda p: p.name) if SDK_DIR.exists() else []
    say("[SDKS] " + ", ".join(p.name for p in sdks))
    rc, out, err = sh(["xcrun", "--show-sdk-path"])
    say("[XCRUN_SDK] " + (out or err).strip())
    rc, out, err = sh(["/usr/bin/clang++", "--version"])
    say("[CLANG] " + (out or err).strip().splitlines()[0] if (out or err).strip() else "[CLANG] ?")
    if built.exists():
        rc, o, e = sh([str(built)], timeout=30)
        if re.search(r"v1\.\d+", o + e):
            say("[MAGMA] already linked")
            SUM["gates"]["MAGMA"] = "OK"
            return built
    shim = TOOLS / "cxx_shim"
    shim.mkdir(exist_ok=True)
    candidates = [p for p in sdks if not p.name.startswith("MacOSX27")] + [None]
    for sdk in candidates:
        isys = ("-isysroot %s " % sdk) if sdk else ""
        for name in ["g++", "c++", "clang++"]:
            f = shim / name
            f.write_text("#!/bin/sh\nexec /usr/bin/clang++ %s-DEIGEN_DONT_VECTORIZE -Wno-everything \"$@\"\n" % isys)
            f.chmod(f.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
        env = dict(os.environ, PATH=str(shim) + ":" + os.environ.get("PATH", ""))
        if sdk:
            env["SDKROOT"] = str(sdk)
        if built.exists():
            built.unlink()
        rc, out, err = sh(["make", "CXX=" + str(shim / "g++")], cwd=str(MAGMA_SRC), timeout=7200, env=env)
        ok = built.exists()
        ver = None
        if ok:
            r2, o2, e2 = sh([str(built)], timeout=30)
            m = re.search(r"v1\.\d+\w*", o2 + e2)
            ver = m.group(0) if m else None
        say("[MAGMA_LINK] sdk=%s rc=%d binary=%s version=%s" % (sdk.name if sdk else "default", rc, ok, ver))
        if ver:
            link = TOOLS / "magma" / "magma_bin"
            if link.is_symlink() or link.exists():
                link.unlink()
            link.symlink_to(built)
            SUM["gates"]["MAGMA"] = "OK " + ver
            return built
        tail = (out + err).strip().splitlines()[-3:]
        for l in tail:
            say("       " + l[:200])
    SUM["gates"]["MAGMA"] = "NOT_LINKED"
    say("[MAGMA] not linked; LDSC-SEG is the primary method, MAGMA remains secondary")
    return None


# ----------------------------------------------------------------------------
# 2. Gene sets
# ----------------------------------------------------------------------------
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


def build_family(pb, gt, genes, sym2, family, type_cols, donors):
    sub = gt[gt["status"].isin(donors) & (gt["n_cells"] >= MIN_CELLS_PER_DONOR_GROUP)].copy()
    if len(sub) == 0:
        say("[%s] no groups" % family)
        return None
    sub["type"] = sub[type_cols].astype(str).apply(lambda r: "@".join(r.values), axis=1)
    agg = sub.groupby("type").agg(n_donors=("sample", "nunique"), n_cells=("n_cells", "sum"))
    ok = agg[(agg["n_donors"] >= MIN_DONORS_PER_TYPE) & (agg["n_cells"] >= MIN_CELLS_PER_TYPE)].index.tolist()
    say("")
    say("[%s] donors=%s eligible types %d / %d" % (family, "+".join(donors), len(ok), len(agg)))
    for t, r in agg.sort_values("n_cells", ascending=False).iterrows():
        say("    [TYPE] %-44s donors=%-3d cells=%-7d %s" % (t[:44], r["n_donors"], r["n_cells"], "OK" if t in ok else "DROP"))
    mapped = np.array([g in sym2 for g in genes])
    means = {}
    for t in ok:
        rows = sub.index[sub["type"] == t].to_numpy()
        mats = np.asarray(pb[rows, :], dtype=np.float64)
        lib = mats.sum(axis=1, keepdims=True)
        means[t] = (mats / np.where(lib > 0, lib, 1) * 1e6).mean(axis=0)
    M = pd.DataFrame(means, index=genes)[mapped]
    M = M[(M >= MIN_CPM_ANY_TYPE).any(axis=1)]
    S = M.div(M.sum(axis=1), axis=0)
    k = int(round(TOP_FRACTION * len(S)))
    fam_dir = GS_DIR / family
    fam_dir.mkdir(parents=True, exist_ok=True)
    sets = {}
    rows = []
    for t in S.columns:
        top = list(S[t].sort_values(ascending=False).index[:k])
        key = family[:2] + "_" + safe(t)
        n = 2
        while key in sets:
            key = family[:2] + "_" + safe(t) + "_%d" % n
            n += 1
        sets[key] = top
        rows.append({"family": family, "set_key": key, "label": t, "n_genes": len(top),
                     "n_donors": int(agg.loc[t, "n_donors"]), "n_cells": int(agg.loc[t, "n_cells"])})
        pd.DataFrame([{"ENTREZ": sym2[g][0], "SYMBOL": g, "CHR": sym2[g][1], "START": sym2[g][2], "STOP": sym2[g][3]}
                      for g in top]).to_csv(fam_dir / (key + ".genes.tsv"), sep="\t", index=False)
    ctrl_key = family[:2] + "_CONTROL_ALL_EXPRESSED"
    sets[ctrl_key] = list(S.index)
    rows.append({"family": family, "set_key": ctrl_key, "label": "ALL_EXPRESSED", "n_genes": len(S.index),
                 "n_donors": None, "n_cells": None})
    M.to_csv(fam_dir / "mean_cpm_by_type.csv.gz")
    S.to_csv(fam_dir / "specificity_by_type.csv.gz")
    man = pd.DataFrame(rows)
    man.to_csv(fam_dir / "manifest.csv", index=False)
    say("[%s] sets=%d universe=%d genes, top10%%=%d genes per set" % (family, len(rows) - 1, len(S), k))
    return {"family": family, "sets": sets, "control": ctrl_key, "manifest": man}


# ----------------------------------------------------------------------------
# 3. Annotations and LD scores
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
        f = LD_DIR / ("annot.%d.annot.gz" % chrom)
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
        df.to_csv(f, sep="\t", index=False, compression="gzip")
    say("[ANNOT] %d columns x 22 chromosomes written" % len(keys))


def run_l2():
    snps = LD_DIR / "print_snps.txt"
    if not snps.exists():
        with open(HM3_LIST) as f, open(snps, "w") as o:
            for line in f:
                t = line.split()
                if t and t[0].startswith("rs"):
                    o.write(t[0] + "\n")
    todo = [c for c in range(1, 23) if not (LD_DIR / ("multi.%d.l2.ldscore.gz" % c)).exists()]
    say("[L2] chromosomes to compute: %s (workers=%d)" % (todo or "none", N_WORKERS))
    running = []
    t0 = time.time()
    while todo or running:
        while todo and len(running) < N_WORKERS:
            c = todo.pop(0)
            cmd = [str(LDSC_PY), str(LDSC), "--l2", "--bfile", str(PLINK_DIR / ("1000G.EUR.QC.%d" % c)),
                   "--ld-wind-cm", str(LD_WIND_CM), "--annot", str(LD_DIR / ("annot.%d" % c)),
                   "--out", str(LD_DIR / ("multi.%d" % c)), "--print-snps", str(snps)]
            errf = open(LD_DIR / ("multi.%d.stderr.txt" % c), "w")
            p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=errf, text=True)
            running.append((c, p, errf))
        time.sleep(5)
        still = []
        for c, p, errf in running:
            if p.poll() is None:
                still.append((c, p, errf))
            else:
                errf.close()
                err = (LD_DIR / ("multi.%d.stderr.txt" % c)).read_text(errors="ignore")
                logf = LD_DIR / ("multi.%d.log" % c)
                if logf.exists():
                    err = err + logf.read_text(errors="ignore")[-400:]
                ok = (LD_DIR / ("multi.%d.l2.ldscore.gz" % c)).exists()
                say("    [L2] chr%-2d rc=%d %s %.0f min" % (c, p.returncode, "OK" if ok else "FAILED", (time.time() - t0) / 60))
                if not ok:
                    say("       " + err.strip()[-300:])
        running = still
    ok_all = all((LD_DIR / ("multi.%d.l2.ldscore.gz" % c)).exists() for c in range(1, 23))
    SUM["gates"]["LDSCORES"] = ok_all
    return ok_all


def split_ldscores(keys):
    marker = LD_DIR / "split.done"
    if marker.exists():
        return
    for c in range(1, 23):
        ld = pd.read_csv(LD_DIR / ("multi.%d.l2.ldscore.gz" % c), sep="\t")
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
    say("[SPLIT] per-set LD score files written for %d sets" % len(keys))


# ----------------------------------------------------------------------------
# 4. LDSC-SEG
# ----------------------------------------------------------------------------
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
        cmd = [str(LDSC_PY), str(LDSC), "--h2-cts", str(MUNGED / (trait + ".sumstats.gz")),
               "--ref-ld-chr", BASELINE, "--ref-ld-chr-cts", str(ldcts), "--w-ld-chr", WEIGHTS, "--out", str(out)]
        t0 = time.time()
        rc, o, e = sh(cmd, timeout=6 * 3600)
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


def main():
    global BASELINE, WEIGHTS
    say("=" * 72)
    say("PAPER 9 | v0.0.7 | GENE SETS + LD SCORES + LDSC-SEG | " + TS)
    say("=" * 72)
    try:
        magma_relink()
    except Exception as e:
        say("[STEP_ERROR] magma %s: %s" % (type(e).__name__, str(e)[:200]))

    say("")
    say("-" * 72)
    say("[SECTION] GENE SETS (GSE266265 pseudobulk from v006)")
    pb = np.load(PB_DIR / "pseudobulk_counts.npy", mmap_mode="r")
    gt = pd.read_csv(PB_DIR / "pseudobulk_groups.csv")
    genes = pd.read_csv(PB_DIR / "genes.csv")["symbol"].astype(str).tolist()
    gt["phase3"] = gt["phase5"].map(PHASE3)
    say("[PB] groups=%d genes=%d statuses=%s" % (len(gt), len(genes), sorted(gt["status"].dropna().unique())))
    sym2 = gene_map()
    fams = []
    for family, cols, donors in [("A_refined_CTL", ["refined"], [CTL_LABEL]),
                                 ("B_main_phase3_CTL", ["main", "phase3"], [CTL_LABEL])]:
        f = build_family(pb, gt, genes, sym2, family, cols, donors)
        if f:
            fams.append(f)
    # Sensitivity families: gene sets only (LD scores later if needed)
    for family, cols, donors in [("SA_refined_ALL", ["refined"], [CTL_LABEL, "Endo"]),
                                 ("SB_main_phase5_ALL", ["main", "phase5"], [CTL_LABEL, "Endo"])]:
        build_family(pb, gt, genes, sym2, family, cols, donors)
    g4_present = [t for t in G4_PRIMARY if any(f["family"] == "A_refined_CTL" and t in set(f["manifest"]["label"]) for f in fams)]
    say("[G4_PRIMARY_SETS_PRESENT] %s" % (", ".join(g4_present) or "NONE"))

    say("")
    say("-" * 72)
    say("[SECTION] ANNOTATIONS + PARTITIONED LD SCORES (+-100 kb, 1 cM, HapMap3 print SNPs)")
    all_sets = {}
    for f in fams:
        all_sets.update(f["sets"])
    keys_file = LD_DIR / "columns.txt"
    if keys_file.exists() and keys_file.read_text().split() != list(all_sets.keys()):
        say("[ANNOT] set list changed since last run -> rebuilding annotations and LD scores")
        for p in LD_DIR.glob("*"):
            if p.is_file():
                p.unlink()
    write_annotations(all_sets, sym2)
    if not run_l2():
        say("[STOP] LD scores incomplete; re-run to resume")
        finish()
        return
    split_ldscores(list(all_sets.keys()))

    BASELINE = prefix_of("baseline_v1.2")
    WEIGHTS = prefix_of("1000G_Phase3_weights_hm3_no_MHC")
    say("[BASELINE] %s" % BASELINE)
    say("[WEIGHTS] %s" % WEIGHTS)

    say("")
    say("-" * 72)
    say("[SECTION] LDSC-SEG (--h2-cts), one-sided coefficient test, BH within family")
    allres = []
    for trait in TRAITS:
        if not (MUNGED / (trait + ".sumstats.gz")).exists():
            say("[SEG] %s sumstats missing" % trait)
            continue
        for f in fams:
            df = ldsc_seg(trait, f)
            if df is None:
                continue
            allres.append(df)
            say("[SEG_TOP] %s | %s" % (trait, f["family"]))
            for _, r in df.head(8).iterrows():
                say("    %-34s coef=%9.2e se=%8.2e p=%.2e q=%.3f" % (str(r["label"])[:34], r["Coefficient"],
                                                                   r["Coefficient_std_error"], r["Coefficient_P_value"], r["q_bh"]))
        if allres:
            pd.concat(allres).to_csv(OUT / "ldsc_seg_all_results.csv", index=False)
        if trait == "GCST90205183" and allres:
            g4_test(pd.concat(allres))
    finish()


def g4_test(res):
    say("")
    say("[SECTION] G4 PRIMARY TEST (GCST90205183, family A_refined_CTL, Holm over dS types)")
    sub = res[(res["trait"] == "GCST90205183") & (res["family"] == "A_refined_CTL") & (res["label"].isin(G4_PRIMARY))].copy()
    if len(sub) == 0:
        say("[G4] no primary sets available")
        return
    sub["holm"] = holm(sub["Coefficient_P_value"].values)
    for _, r in sub.iterrows():
        say("[G4] %-20s coef=%9.2e se=%8.2e p=%.2e Holm=%.3f" % (r["label"], r["Coefficient"], r["Coefficient_std_error"],
                                                               r["Coefficient_P_value"], r["holm"]))
    rep = bool((sub["holm"] < 0.05).any())
    verdict = ("DECIDUAL ENRICHMENT REPLICATED (check 1 passed; G4 bold hypothesis not supported at this check)"
               if rep else "DECIDUAL ENRICHMENT NOT REPLICATED (check 1 failed; supports G4)")
    say("[G4_VERDICT] " + verdict)
    sec = res[(res["trait"] == "GCST90205183") & (res["family"] == "B_main_phase3_CTL") & (res["label"] == "stromal@secretory")]
    for _, r in sec.iterrows():
        say("[G4_SECONDARY] stromal@secretory coef=%9.2e p=%.2e q=%.3f" % (r["Coefficient"], r["Coefficient_P_value"], r["q_bh"]))
    sub.to_csv(OUT / "G4_primary.csv", index=False)
    SUM["G4_primary"] = sub[["label", "Coefficient", "Coefficient_std_error", "Coefficient_P_value", "holm"]].to_dict("records")
    SUM["G4_verdict"] = verdict


def finish():
    with open(OUT / "p9_v007_ldsc_seg_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    for k, v in SUM["gates"].items():
        say("[GATE] %-10s %s" % (k, v))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] Gene sets from CTL donors only; G4 primary rule fixed in CHANGELOG v0.2.4 before this run")
    say("[DONE] PAPER9_V007_COMPLETE")
    _fh.close()


BASELINE = None
WEIGHTS = None

if __name__ == "__main__":
    main()
