#!/usr/bin/env python3
# Paper 9 | v0.0.6 | (1) MAGMA build for arm64 (Eigen vectorization off) + sanity gene analysis
#                   (2) GSE266265 donor-level pseudobulk (sample x AnnotationRefined) from the Paper 8 export
#                   (3) Eutopic specificity gene sets on non-endometriosis donors (rules: CHANGELOG v0.2.3)
#                   (4) HECA download with resume
# Runs in ~/Desktop/paper9/.venv (numpy, pandas). No disease effects are computed.

import datetime as dt
import gzip
import json
import os
import re
import stat
import subprocess
import sys
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
MAGMA_DIR = TOOLS / "magma"
MAGMA_OUT = INTERIM / "magma"
GENE_LOC_B37 = MAGMA_DIR / "genes_b37" / "NCBI37.3.gene.loc"
G1000 = MAGMA_DIR / "ref_eur" / "g1000_eur"
GWAS_H = RAW / "gwas" / "GCST90205183" / "harmonised" / "GCST90205183.h.tsv.gz"
GWAS_N = 474160
EXPORT = (REPO / "data" / "external" / "p8_GSE266265_v232c_export").resolve() / "rna_counts_and_metadata"
PB_DIR = PROC / "gse266265_pseudobulk"
GS_DIR = PROC / "genesets" / "gse266265"
HECA_URL = "https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/endometriumAtlasV2_cells_with_counts.h5ad"
HECA_FILE = RAW / "heca" / "endometriumAtlasV2_cells_with_counts.h5ad"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "genesets" / ("p9_v006_genesets_magma_" + TS)
LOG = REPO / "logs" / ("p9_v006_genesets_magma_" + TS + ".log")

# Pre-registered gene-set rules (CHANGELOG v0.2.3)
MIN_CELLS_PER_DONOR_GROUP = 10
MIN_DONORS_PER_TYPE = 3
MIN_CELLS_PER_TYPE = 50
MIN_CPM_ANY_TYPE = 1.0
TOP_FRACTION = 0.10

KNOWN_ENDO_GENES = ["WNT4", "CDC42", "GREB1", "ESR1", "FN1", "VEZT", "ID4", "KDR", "IL1A", "SYNE1",
                    "CDKN2B", "FSHB", "ETAA1", "RSPO3", "SKAP1", "BSN", "MLLT10", "HOXC10"]

for d in (OUT, LOG.parent, MAGMA_OUT, PB_DIR, GS_DIR, HECA_FILE.parent):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V006_GENESETS_MAGMA", "timestamp_utc": TS, "gates": {}}


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


def human(n):
    n = float(n)
    for u in ["B", "KB", "MB", "GB", "TB"]:
        if n < 1024:
            return "%.1f %s" % (n, u)
        n /= 1024.0
    return "%.1f PB" % n


# ----------------------------------------------------------------------------
# 1. MAGMA build on arm64
# ----------------------------------------------------------------------------
def magma_version(binpath):
    for args in (["--version"], []):
        rc, out, err = sh([str(binpath)] + args, timeout=60)
        m = re.search(r"v1\.\d+\w*", out + err)
        if m:
            return m.group(0)
    return None


def build_magma():
    say("")
    say("-" * 72)
    say("[SECTION] MAGMA BUILD (arm64, -DEIGEN_DONT_VECTORIZE)")
    src = MAGMA_DIR / "source"
    built = src / "magma"
    if built.exists() and magma_version(built):
        say("[MAGMA] already built: %s %s" % (built, magma_version(built)))
        return built
    if not (src / "makefile").exists() and not (src / "Makefile").exists():
        say("[MAGMA] source missing")
        return None
    for o in src.rglob("*.o"):
        o.unlink()
    shim = TOOLS / "cxx_shim"
    shim.mkdir(exist_ok=True)
    flags = "-DEIGEN_DONT_VECTORIZE -Wno-everything"
    for name in ["g++", "c++", "clang++"]:
        f = shim / name
        f.write_text("#!/bin/sh\nexec /usr/bin/clang++ %s \"$@\"\n" % flags)
        f.chmod(f.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    env = dict(os.environ, PATH=str(shim) + ":" + os.environ.get("PATH", ""))
    t0 = time.time()
    rc, out, err = sh(["make", "CXX=" + str(shim / "g++")], cwd=str(src), timeout=7200, env=env)
    say("[MAGMA_COMPILE] rc=%d %.0fs" % (rc, time.time() - t0))
    if rc != 0:
        for l in (out + err).strip().splitlines()[-15:]:
            say("       " + l[:220])
    ver = magma_version(built) if built.exists() else None
    say("[MAGMA] %s %s" % ("COMPILED_OK" if ver else "NOT_INSTALLED", ver or ""))
    SUM["gates"]["MAGMA"] = "COMPILED_OK" if ver else "NOT_INSTALLED"
    if ver:
        link = MAGMA_DIR / "magma_bin"
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(built)
        return built
    return None


def magma_sanity(binpath):
    say("")
    say("-" * 72)
    say("[SECTION] MAGMA SANITY: GCST90205183 gene analysis (35kb/10kb)")
    annot = MAGMA_OUT / "annot_b37_35_10"
    if not Path(str(annot) + ".genes.annot").exists():
        rc, out, err = sh([str(binpath), "--annotate", "window=35,10", "--snp-loc", str(G1000) + ".bim",
                           "--gene-loc", str(GENE_LOC_B37), "--out", str(annot)], timeout=3600)
        say("[MAGMA_ANNOTATE] rc=%d" % rc)
        if rc != 0:
            say("       " + (out + err).strip()[-300:])
            return
    pfile = MAGMA_OUT / "GCST90205183.pval.txt"
    if not pfile.exists():
        n = 0
        with gzip.open(GWAS_H, "rt") as f, open(pfile, "w") as o:
            hdr = f.readline().rstrip("\n").split("\t")
            ix = {h: i for i, h in enumerate(hdr)}
            o.write("SNP\tP\n")
            for line in f:
                t = line.rstrip("\n").split("\t")
                rs, p = t[ix["rsid"]], t[ix["p_value"]]
                if rs.startswith("rs") and p not in ("", "NA"):
                    o.write(rs + "\t" + p + "\n")
                    n += 1
        say("[MAGMA_PFILE] %d SNPs" % n)
    outp = MAGMA_OUT / "GCST90205183_35_10"
    gout = Path(str(outp) + ".genes.out")
    if not gout.exists():
        say("[MAGMA_GENE] running (30-90 min, console is quiet)")
        t0 = time.time()
        rc, out, err = sh([str(binpath), "--bfile", str(G1000), "--pval", str(pfile), "N=%d" % GWAS_N,
                           "--gene-annot", str(annot) + ".genes.annot", "--out", str(outp)], timeout=6 * 3600)
        say("[MAGMA_GENE] rc=%d %.0f min" % (rc, (time.time() - t0) / 60))
        if rc != 0:
            say("       " + (out + err).strip()[-400:])
            return
    sym = {}
    for line in open(GENE_LOC_B37):
        t = line.split()
        if len(t) >= 6:
            sym[t[0]] = t[5]
    df = pd.read_csv(gout, sep=r"\s+")
    df["SYMBOL"] = df["GENE"].astype(str).map(sym)
    df = df.sort_values("P")
    bonf = 0.05 / len(df)
    say("[MAGMA_GENES] tested=%d bonferroni(p<%.1e)=%d" % (len(df), bonf, int((df["P"] < bonf).sum())))
    for _, r in df.head(15).iterrows():
        say("    [TOP] %-12s z=%6.2f p=%.2e" % (r["SYMBOL"], r["ZSTAT"], r["P"]))
    ranks = {s: i + 1 for i, s in enumerate(df["SYMBOL"])}
    hits = [(g, ranks[g], float(df.loc[df["SYMBOL"] == g, "P"].iloc[0])) for g in KNOWN_ENDO_GENES if g in ranks]
    say("[MAGMA_KNOWN_LOCI] " + "; ".join("%s rank=%d p=%.1e" % h for h in hits))
    k = sum(1 for h in hits if h[2] < bonf)
    SUM["gates"]["MAGMA_SANITY"] = "PASS" if k >= 5 else "CHECK"
    say("[GATE] MAGMA_SANITY %s (%d known genes Bonferroni-significant)" % (SUM["gates"]["MAGMA_SANITY"], k))


# ----------------------------------------------------------------------------
# 2. GSE266265 export readers
# ----------------------------------------------------------------------------
def find_col_file(idx, name):
    for p in EXPORT.iterdir():
        m = re.match(r"^\d+_meta\.data_(\d+)__(.+)\.(json|npy)$", p.name)
        if m and int(m.group(1)) == idx and m.group(2) == name:
            return p
    return None


def levels_file(idx):
    for p in EXPORT.iterdir():
        if re.match(r"^\d+_meta\.data_%d__levels\.json$" % idx, p.name):
            return p
    return None


def load_col(names, name):
    idx = names.index(name)
    f = find_col_file(idx, name)
    if f is None:
        raise FileNotFoundError("column file for " + name)
    if f.suffix == ".json":
        return np.array(json.loads(f.read_text()), dtype=object)
    arr = np.load(f)
    lv = levels_file(idx)
    if lv is not None:
        levels = json.loads(lv.read_text())
        out = np.empty(arr.shape[0], dtype=object)
        ok = (arr >= 1) & (arr <= len(levels))
        out[ok] = np.array(levels, dtype=object)[arr[ok] - 1]
        out[~ok] = None
        return out
    return arr


def file_by_prefix(tag):
    for p in EXPORT.iterdir():
        if tag in p.name:
            return p
    raise FileNotFoundError(tag)


def build_pseudobulk():
    say("")
    say("-" * 72)
    say("[SECTION] GSE266265 PSEUDOBULK (sample x AnnotationRefined, raw UMI sums)")
    done = PB_DIR / "pseudobulk_counts.npy"
    names = json.loads(file_by_prefix("meta.data_names.json").read_text())
    genes = json.loads(file_by_prefix("counts_Dimnames_0.json").read_text())
    cells = json.loads(file_by_prefix("counts_Dimnames_1.json").read_text())
    rown = json.loads(file_by_prefix("meta.data_row.names.json").read_text())
    say("[EXPORT] genes=%d cells=%d metadata_rows=%d order_match=%s" % (len(genes), len(cells), len(rown), cells == rown))
    if cells != rown:
        raise RuntimeError("cell order mismatch between matrix and metadata")
    sample = load_col(names, "sample")
    status = load_col(names, "EndometriosisStatus")
    phase5 = load_col(names, "MenstrualCyclePhase")
    main = load_col(names, "AnnotationMain")
    refined = load_col(names, "AnnotationRefined")
    dbl = load_col(names, "scDblFinder.class")
    for nm, arr in [("EndometriosisStatus", status), ("MenstrualCyclePhase", phase5), ("scDblFinder.class", dbl),
                    ("AnnotationMain", main)]:
        vc = pd.Series(arr).value_counts(dropna=False)
        say("[VALUES] %-22s %s" % (nm, " | ".join("%s=%d" % (k, v) for k, v in vc.items())[:700]))
    vc = pd.Series(refined).value_counts(dropna=False)
    say("[VALUES] AnnotationRefined n=%d" % len(vc))
    for k, v in vc.items():
        say("    [REFINED] %-40s %d" % (k, v))

    keep = np.array([d is None or str(d).lower() != "doublet" for d in dbl])
    keep &= np.array([r is not None and str(r) not in ("", "NA", "nan") for r in refined])
    say("[CELLS_KEPT] %d of %d (doublets and unannotated removed)" % (int(keep.sum()), len(keep)))

    meta = pd.DataFrame({"sample": sample, "status": status, "phase5": phase5, "main": main, "refined": refined})
    # Sample-level attributes must be constant within sample
    gb = meta.groupby("sample")
    n_status = gb["status"].nunique()
    n_phase = gb["phase5"].nunique()
    s_status = gb["status"].first()
    s_phase = gb["phase5"].first()
    bad = int(((n_status != 1) | (n_phase != 1)).sum())
    say("[SAMPLES] %d; inconsistent status/phase within sample: %d" % (len(n_status), bad))
    say("[SAMPLES_BY_STATUS_PHASE] " + " | ".join("%s/%s=%d" % (a, b, c) for (a, b), c in
        pd.DataFrame({"s": s_status, "p": s_phase}).value_counts().sort_index().items()))
    grp_keys = meta.loc[keep, ["sample", "refined"]].astype(str).agg("||".join, axis=1)
    codes, uniq = pd.factorize(grp_keys)
    group_of_cell = np.full(len(cells), -1, dtype=np.int64)
    group_of_cell[np.where(keep)[0]] = codes
    G = len(uniq)
    say("[GROUPS] sample x refined groups=%d" % G)

    if done.exists():
        say("[PSEUDOBULK] present -> reuse")
    else:
        p = np.load(file_by_prefix("counts_p.npy"))
        i = np.load(file_by_prefix("counts_i.npy"), mmap_mode="r")
        x = np.load(file_by_prefix("counts_x.npy"), mmap_mode="r")
        sample_x = np.asarray(x[:2_000_000])
        say("[COUNTS_CHECK] integer-valued=%s max=%.0f" % (bool(np.all(np.mod(sample_x, 1) == 0)), sample_x.max()))
        ng = len(genes)
        acc = np.zeros(G * ng, dtype=np.float64)
        t0 = time.time()
        ncell = len(cells)
        step = 20000
        for a in range(0, ncell, step):
            b = min(ncell, a + step)
            s, e = int(p[a]), int(p[b])
            ii = np.asarray(i[s:e], dtype=np.int64)
            xx = np.asarray(x[s:e])
            g = np.repeat(group_of_cell[a:b], np.diff(p[a:b + 1]))
            m = g >= 0
            acc += np.bincount(g[m] * ng + ii[m], weights=xx[m], minlength=G * ng)
            if (a // step) % 5 == 0:
                say("    [PB] cells %d/%d %.0fs" % (b, ncell, time.time() - t0))
        np.save(done, acc.reshape(G, ng))
        say("[PSEUDOBULK] built %.0fs" % (time.time() - t0))
    sizes = pd.Series(codes).value_counts().sort_index()
    gt = pd.DataFrame({"group": list(uniq)})
    gt["sample"] = gt["group"].str.split("||", regex=False).str[0]
    gt["refined"] = gt["group"].str.split("||", regex=False).str[1]
    gt["n_cells"] = sizes.values
    gt["status"] = gt["sample"].map(s_status)
    gt["phase5"] = gt["sample"].map(s_phase)
    main_of = meta.loc[keep].groupby("refined")["main"].agg(lambda s: s.value_counts().index[0])
    gt["main"] = gt["refined"].map(main_of)
    gt.to_csv(PB_DIR / "pseudobulk_groups.csv", index=False)
    pd.Series(genes).to_csv(PB_DIR / "genes.csv", index=False, header=["symbol"])
    return np.load(done, mmap_mode="r"), gt, genes


# ----------------------------------------------------------------------------
# 3. Specificity gene sets (non-endometriosis donors)
# ----------------------------------------------------------------------------
def gene_map():
    sym2 = {}
    dup = set()
    for line in open(GENE_LOC_B37):
        t = line.split()
        if len(t) >= 6 and t[1] in [str(c) for c in range(1, 23)]:
            if t[5] in sym2:
                dup.add(t[5])
            sym2[t[5]] = (t[0], t[1], int(t[2]), int(t[3]))
    for d in dup:
        sym2.pop(d, None)
    return sym2


def control_label(values):
    vals = [v for v in pd.Series(values).dropna().unique()]
    ctrl = [v for v in vals if re.search(r"\bno\b|non|control|without", str(v), re.I)]
    return ctrl


def specificity_sets(pb, gt, genes, family, type_col, extra_cols, sym2):
    ctrl = control_label(gt["status"])
    say("[%s] control status values: %s" % (family, ctrl))
    sub = gt[gt["status"].isin(ctrl) & (gt["n_cells"] >= MIN_CELLS_PER_DONOR_GROUP)].copy()
    sub["type"] = sub[[type_col] + extra_cols].astype(str).agg("@".join, axis=1) if extra_cols else sub[type_col].astype(str)
    agg = sub.groupby("type").agg(n_donors=("sample", "nunique"), n_cells=("n_cells", "sum"))
    ok_types = agg[(agg["n_donors"] >= MIN_DONORS_PER_TYPE) & (agg["n_cells"] >= MIN_CELLS_PER_TYPE)].index.tolist()
    say("[%s] types eligible %d / %d" % (family, len(ok_types), len(agg)))
    for t, r in agg.iterrows():
        say("    [TYPE] %-48s donors=%-3d cells=%-7d %s" % (t[:48], r["n_donors"], r["n_cells"], "OK" if t in ok_types else "DROP"))
    mapped = np.array([g in sym2 for g in genes])
    cpm_means = {}
    for t in ok_types:
        rows = sub.index[sub["type"] == t].to_numpy()
        mats = np.asarray(pb[rows, :], dtype=np.float64)
        lib = mats.sum(axis=1, keepdims=True)
        cpm = mats / np.where(lib > 0, lib, 1) * 1e6
        cpm_means[t] = cpm.mean(axis=0)
    M = pd.DataFrame(cpm_means, index=genes)
    M = M[mapped]
    M = M[(M >= MIN_CPM_ANY_TYPE).any(axis=1)]
    S = M.div(M.sum(axis=1), axis=0)
    fam_dir = GS_DIR / family
    fam_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    k = int(round(TOP_FRACTION * len(S)))
    for t in S.columns:
        top = S[t].sort_values(ascending=False).index[:k]
        fname = re.sub(r"[^A-Za-z0-9_.@+-]+", "_", t) + ".genes.txt"
        with open(fam_dir / fname, "w") as f:
            f.write("ENTREZ\tSYMBOL\tCHR\tSTART\tSTOP\n")
            for g in top:
                e, c, s0, s1 = sym2[g]
                f.write("%s\t%s\t%s\t%d\t%d\n" % (e, g, c, s0, s1))
        rows.append({"family": family, "set": t, "file": fname, "n_genes": len(top),
                     "n_donors": int(agg.loc[t, "n_donors"]), "n_cells": int(agg.loc[t, "n_cells"])})
    with open(fam_dir / "ALL_EXPRESSED.genes.txt", "w") as f:
        f.write("ENTREZ\tSYMBOL\tCHR\tSTART\tSTOP\n")
        for g in S.index:
            e, c, s0, s1 = sym2[g]
            f.write("%s\t%s\t%s\t%d\t%d\n" % (e, g, c, s0, s1))
    M.to_csv(fam_dir / "mean_cpm_by_type.csv.gz")
    S.to_csv(fam_dir / "specificity_by_type.csv.gz")
    pd.DataFrame(rows).to_csv(fam_dir / "manifest.csv", index=False)
    say("[%s] sets=%d genes_universe=%d top10%%=%d genes/set -> %s" % (family, len(rows), len(S), k, fam_dir))
    return rows


# ----------------------------------------------------------------------------
# 4. HECA download with resume loop
# ----------------------------------------------------------------------------
def heca_download():
    say("")
    say("-" * 72)
    say("[SECTION] HECA download (14.9 GB, resumable)")
    part = HECA_FILE.with_name(HECA_FILE.name + ".part")
    if HECA_FILE.exists():
        say("[HECA] present %s" % human(HECA_FILE.stat().st_size))
        SUM["gates"]["HECA_DOWNLOAD"] = True
        return
    for attempt in range(1, 16):
        rc, out, err = sh(["curl", "-L", "-sS", "--fail", "-C", "-", "--retry", "5", "--retry-all-errors",
                           "--retry-delay", "10", "--connect-timeout", "30", "-o", str(part), HECA_URL])
        size = part.stat().st_size if part.exists() else 0
        say("[HECA] attempt %d rc=%d have=%s" % (attempt, rc, human(size)))
        if rc == 0:
            part.rename(HECA_FILE)
            break
        time.sleep(20)
    ok = HECA_FILE.exists()
    SUM["gates"]["HECA_DOWNLOAD"] = ok
    say("[HECA] %s" % ("OK " + human(HECA_FILE.stat().st_size) if ok else "INCOMPLETE - re-run to resume"))


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.6 | GENE SETS + MAGMA | " + TS)
    say("=" * 72)
    magma = None
    try:
        magma = build_magma()
    except Exception as e:
        say("[STEP_ERROR] magma_build %s: %s" % (type(e).__name__, str(e)[:200]))
    try:
        pb, gt, genes = build_pseudobulk()
        sym2 = gene_map()
        mapped = sum(1 for g in genes if g in sym2)
        say("[GENE_MAP] %d/%d symbols mapped to unique autosomal NCBI37.3 genes" % (mapped, len(genes)))
        rows = []
        rows += specificity_sets(pb, gt, genes, "A_refined_types", "refined", [], sym2)
        rows += specificity_sets(pb, gt, genes, "B_main_by_phase5", "main", ["phase5"], sym2)
        pd.DataFrame(rows).to_csv(OUT / "geneset_manifest.csv", index=False)
        SUM["gates"]["GENESETS"] = len(rows)
    except Exception as e:
        say("[STEP_ERROR] pseudobulk/genesets %s: %s" % (type(e).__name__, str(e)[:300]))
    try:
        heca_download()
    except Exception as e:
        say("[STEP_ERROR] heca %s: %s" % (type(e).__name__, str(e)[:200]))
    if magma:
        try:
            magma_sanity(magma)
        except Exception as e:
            say("[STEP_ERROR] magma_sanity %s: %s" % (type(e).__name__, str(e)[:200]))
    with open(OUT / "p9_v006_genesets_magma_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    for k, v in SUM["gates"].items():
        say("[GATE] %-14s %s" % (k, v))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] Gene sets from non-endometriosis donors only; no disease effects computed")
    say("[DONE] PAPER9_V006_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
