#!/usr/bin/env python3
# Paper 9 | v0.0.23 | Second (larger, partially overlapping) GWAS: Koller et al. 2026 EUR without 23andMe
#   (Zenodo 10.5281/zenodo.18983492, CC-BY-NC-4.0). Rules: CHANGELOG v0.2.14 (fixed before this run).
#   1 download + format detection + rsID mapping   2 MAGMA gene analysis   3 munge, h2, rg (overlap)
#   4 re-run unchanged: G6, X4a/X4b, G4 HECA M1/M3, G1 primary (+ lymphocyte permutation), next to GCST90205183
# Runs in ~/Desktop/paper9/.venv, low-heat (single thread, nice 10). ~60-90 min (MAGMA gene analysis ~35 min).

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
try:
    os.nice(10)
except Exception:
    pass

import datetime as dt
import gzip
import json
import math
import re
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
EXT = REPO / "data" / "external"
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
G1000 = TOOLS / "magma" / "ref_eur" / "g1000_eur"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
ANNOT = MAGMA_OUT / "annot_b37_35_10.genes.annot"
LDSC_PY = TOOLS / "ldsc_venv" / "bin" / "python"
LDSC_DIR = TOOLS / "ldsc_py3"
REF_UNPACKED = RAW / "ldsc_ref" / "unpacked"
HM3_LIST = RAW / "ldsc_ref" / "hm3_no_MHC.list.txt"
MERGE_ALLELES = INTERIM / "ldsc" / "hm3_merge_alleles.tsv"
MUNGED = PROC / "munged"
HM3_DIR = INTERIM / "sumstats_hm3"
GCST_H = RAW / "gwas" / "GCST90205183" / "harmonised" / "GCST90205183.h.tsv.gz"
GTEX_FILE = EXT / "gtex" / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
HECA_GS = PROC / "genesets" / "heca_v014"
COV_G6 = INTERIM / "magma_v019" / "G6_gtex_v8_52.covar.txt"
COV_X4 = INTERIM / "magma_v022"
COV_C5 = INTERIM / "magma_v020"
COVAR_G1 = PROC / "magma_covar" / "g1"
PB213 = PROC / "gse213216_pseudobulk_clusters"
KDIR = RAW / "gwas" / "Koller2026"
ZENODO_API = "https://zenodo.org/api/records/18983492"
KFILE_NAME = "endometriosis_EUR_wo_23andMe.txt.gz"
README_NAME = "README_GWAS.txt"
TRAIT = "KOLLER_EUR"
REF = "GCST90205183"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "koller" / ("p9_v023_koller_" + TS)
WORK = INTERIM / "magma_v023"
LOG = REPO / "logs" / ("p9_v023_koller_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

REPRO = ["Uterus", "Ovary", "Vagina", "Fallopian Tube", "Cervix - Ectocervix", "Cervix - Endocervix"]
DSTROMAL = ["dStromal_early", "dStromal_mid", "dStromal_late"]
LINEAGES = ["stromal", "epithelial", "endothelial", "myeloid", "lymphocyte"]
N_PERM = 100
PERM_SEED = 20261106

for d in (OUT, WORK, LOG.parent, KDIR, MUNGED, HM3_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V023_KOLLER", "timestamp_utc": TS, "gates": {}, "results": {}}


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def sh(cmd, timeout=None):
    try:
        p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout, env=ENV1)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"


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


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def entrez_unique():
    sym2, dup = {}, set()
    for line in open(GENE_LOC_B37):
        x = line.split()
        if len(x) >= 6 and x[1] in [str(c) for c in range(1, 23)]:
            if x[5] in sym2:
                dup.add(x[5])
            sym2[x[5]] = x[0]
    for d in dup:
        sym2.pop(d, None)
    return sym2


def magma_gp(trait, covf, cond, tag):
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    if not raw.exists() or not Path(covf).exists():
        return None
    outp = WORK / ("%s__%s" % (trait, tag))
    gsa = Path(str(outp) + ".gsa.out")
    keys = ["condition-hide=", "condition="] if cond else [None]
    o = e = ""
    for key in keys:
        if gsa.exists():
            gsa.unlink()
        model = ["direction=greater"] if key is None else [key + cond, "direction=greater"]
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--gene-covar", covf, "--model"] + model + ["--out", outp],
                      timeout=3 * 3600)
        if rc == 0 and gsa.exists():
            rows = [l for l in open(gsa) if not l.startswith("#") and l.strip()]
            hdr = rows[0].split()
            res = {}
            for l in rows[1:]:
                d = dict(zip(hdr, l.split()))
                name = d.get("FULL_NAME", d.get("VARIABLE"))
                try:
                    res[name] = {"beta": float(d["BETA"]), "beta_std": float(d.get("BETA_STD", "nan")),
                                 "se": float(d["SE"]), "p": float(d["P"])}
                except (KeyError, ValueError):
                    continue
            return res
    say("    [MAGMA_FAIL] %s %s %s" % (trait, tag, (o + e).strip()[-200:]))
    return None


# ----------------------------------------------------------------------------
# 1. Download, format detection, rsID mapping
# ----------------------------------------------------------------------------
def zenodo_links():
    p = subprocess.run(["curl", "-sSL", "--fail", "--max-time", "60", "-A", "paper9", ZENODO_API], capture_output=True)
    links = {}
    if p.returncode == 0:
        rec = json.loads(p.stdout.decode("utf-8", "replace"))
        for f in rec.get("files", []):
            links[f.get("key")] = (f.get("links", {}).get("self"), f.get("size", 0))
    return links


def download(name, links):
    dst = KDIR / name
    url, size = links.get(name, (None, 0))
    if dst.exists() and (not size or dst.stat().st_size == size):
        say("[DL] %-40s SKIP_PRESENT %.1f MB" % (name, dst.stat().st_size / 1e6))
        return dst
    urls = [u for u in [url, "https://zenodo.org/records/18983492/files/%s?download=1" % name] if u]
    for u in urls:
        part = dst.with_name(dst.name + ".part")
        rc, o, e = sh(["curl", "-L", "-sS", "--fail", "--retry", "3", "-C", "-", "-A", "paper9", "-o", part, u], timeout=7200)
        if rc == 0 and part.exists() and (not size or part.stat().st_size == size):
            part.rename(dst)
            say("[DL] %-40s OK %.1f MB" % (name, dst.stat().st_size / 1e6))
            return dst
        say("[DL] %-40s FAILED rc=%s %s" % (name, rc, e.strip()[:120]))
    return None


def pick(cols, cands):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c in low:
            return low[c]
    return None


def detect(kfile):
    head = pd.read_csv(kfile, sep=None, engine="python", nrows=20000, compression="infer")
    cols = list(head.columns)
    say("[FORMAT] columns: " + ", ".join(cols))
    for line in head.head(3).to_string(index=False).splitlines():
        say("    " + line[:220])
    c = {
        "id": pick(cols, ["rsid", "snp", "markername", "variant_id", "snpid", "id", "rs", "rs_id"]),
        "chr": pick(cols, ["chr", "chrom", "chromosome", "#chrom", "#chr", "hg38_chr", "chr_b38"]),
        "pos": pick(cols, ["pos", "bp", "position", "base_pair_location", "pos_b38", "hg38_pos", "bp_b38"]),
        "a1": pick(cols, ["allele1", "a1", "effect_allele", "ea", "alt", "tested_allele"]),
        "a2": pick(cols, ["allele2", "a2", "other_allele", "nea", "ref", "non_effect_allele"]),
        "p": pick(cols, ["p-value", "p", "pval", "p_value", "pvalue", "p.value"]),
        "z": pick(cols, ["zscore", "z", "z_score", "stat"]),
        "beta": pick(cols, ["beta", "effect", "b", "log_or", "logor"]),
        "se": pick(cols, ["se", "stderr", "standard_error", "sebeta"]),
        "n": pick(cols, ["n", "weight", "n_total", "neff", "totaln", "n_eff", "samplesize"]),
        "frq": pick(cols, ["freq1", "eaf", "effect_allele_frequency", "frq", "af", "maf"]),
    }
    ids = head[c["id"]].astype(str) if c["id"] else pd.Series([], dtype=str)
    rs_share = float(ids.str.startswith("rs").mean()) if len(ids) else 0.0
    say("[FORMAT] mapped columns: " + ", ".join("%s=%s" % (k, v) for k, v in c.items()) + " | rsID share=%.3f" % rs_share)
    return c, rs_share


def key_of(chrom, pos):
    ch = pd.Series(chrom).astype(str).str.replace("chr", "", regex=False).str.upper().replace({"X": "23", "Y": "24"})
    ch = pd.to_numeric(ch, errors="coerce")
    ps = pd.to_numeric(pd.Series(pos), errors="coerce")
    k = ch * 10_000_000_000 + ps
    return k.to_numpy(dtype="float64")


def build_map(build):
    # Arrays of (chr*1e10 + pos) -> rs number, sorted; b38 from the GCST harmonised file, b37 from the 1000G bim
    keys, rs = [], []
    if build == "b38":
        for ch in pd.read_csv(GCST_H, sep="\t", compression="gzip", chunksize=1_000_000,
                              usecols=lambda x: x in ("hm_rsid", "rsid", "hm_chrom", "hm_pos", "chromosome", "base_pair_location")):
            idc = "hm_rsid" if "hm_rsid" in ch.columns else "rsid"
            cc = "hm_chrom" if "hm_chrom" in ch.columns else "chromosome"
            pc = "hm_pos" if "hm_pos" in ch.columns else "base_pair_location"
            r = ch[idc].astype(str)
            ok = r.str.startswith("rs")
            keys.append(key_of(ch.loc[ok, cc], ch.loc[ok, pc]))
            rs.append(pd.to_numeric(r[ok].str[2:], errors="coerce").to_numpy(dtype="float64"))
    else:
        for ch in pd.read_csv(str(G1000) + ".bim", sep=r"\s+", header=None, usecols=[0, 1, 3], chunksize=2_000_000):
            r = ch[1].astype(str)
            ok = r.str.startswith("rs")
            keys.append(key_of(ch.loc[ok, 0], ch.loc[ok, 3]))
            rs.append(pd.to_numeric(r[ok].str[2:], errors="coerce").to_numpy(dtype="float64"))
    k = np.concatenate(keys)
    r = np.concatenate(rs)
    m = np.isfinite(k) & np.isfinite(r)
    k, r = k[m], r[m]
    o = np.argsort(k, kind="mergesort")
    return k[o], r[o]


def lookup(kmap, rmap, k):
    i = np.searchsorted(kmap, k)
    i = np.clip(i, 0, len(kmap) - 1)
    hit = kmap[i] == k
    out = np.full(len(k), np.nan)
    out[hit] = rmap[i[hit]]
    return out


def prepare(kfile):
    say("")
    say("-" * 72)
    say("[SECTION] 1. Koller 2026 EUR without 23andMe: format and rsID mapping")
    pfile = MAGMA_OUT / (TRAIT + ".pval.txt")
    hfile = HM3_DIR / (TRAIT + ".tsv")
    c, rs_share = detect(kfile)
    if pfile.exists() and hfile.exists():
        say("[PREPARE] present -> reuse")
        return c
    if c["p"] is None or (c["id"] is None and (c["chr"] is None or c["pos"] is None)):
        raise RuntimeError("cannot find p-value / SNP identifier columns")
    hm3 = set()
    for line in open(HM3_LIST, errors="ignore"):
        t = line.split()
        if t and t[0].startswith("rs"):
            hm3.add(t[0])
    use_pos = rs_share < 0.5
    kmap = rmap = None
    if use_pos:
        sample = pd.read_csv(kfile, sep=None, engine="python", nrows=200000, compression="infer")
        if c["chr"] and c["pos"]:
            ks = key_of(sample[c["chr"]], sample[c["pos"]])
        else:
            parts = sample[c["id"]].astype(str).str.replace("chr", "", regex=False).str.split(r"[:_]", expand=True)
            ks = key_of(parts[0], parts[1])
        best = None
        for b in ("b38", "b37"):
            km, rm = build_map(b)
            share = float(np.isfinite(lookup(km, rm, ks)).mean())
            say("[BUILD] %s match share in first 200k rows: %.3f" % (b, share))
            if best is None or share > best[0]:
                best = (share, b, km, rm)
            else:
                del km, rm
        _, b, kmap, rmap = best
        say("[BUILD] using %s positions" % b)
        SUM["gates"]["BUILD"] = b
    n_in = n_out = n_hm3 = 0
    with gzip.open(kfile, "rt") as f0:
        first = f0.readline()
    sep = "\t" if "\t" in first else r"\s+"
    reader = pd.read_csv(kfile, sep=sep, chunksize=500000, compression="infer", low_memory=False)
    with open(pfile, "w") as fp, open(hfile, "w") as fh:
        fp.write("SNP\tP\tN\n" if c["n"] else "SNP\tP\n")
        fh.write("SNP\tA1\tA2\tZ\tP\tN\n" if c["n"] else "SNP\tA1\tA2\tZ\tP\n")
        for ch in reader:
            n_in += len(ch)
            if use_pos:
                if c["chr"] and c["pos"]:
                    k = key_of(ch[c["chr"]], ch[c["pos"]])
                else:
                    parts = ch[c["id"]].astype(str).str.replace("chr", "", regex=False).str.split(r"[:_]", expand=True)
                    k = key_of(parts[0], parts[1])
                rn = lookup(kmap, rmap, k)
                ok = np.isfinite(rn)
                snp = pd.Series(["rs%d" % v for v in rn[ok]], index=ch.index[ok])
                ch = ch.loc[ok]
            else:
                snp = ch[c["id"]].astype(str).str.split(r"[;,]").str[0]
                ok = snp.str.startswith("rs").to_numpy()
                ch, snp = ch.loc[ok], snp[ok]
            p = pd.to_numeric(ch[c["p"]], errors="coerce")
            if c["z"]:
                z = pd.to_numeric(ch[c["z"]], errors="coerce")
            elif c["beta"] and c["se"]:
                z = pd.to_numeric(ch[c["beta"]], errors="coerce") / pd.to_numeric(ch[c["se"]], errors="coerce")
            elif c["beta"]:
                from scipy.stats import norm
                z = np.sign(pd.to_numeric(ch[c["beta"]], errors="coerce")) * norm.isf(p / 2)
            else:
                z = pd.Series(np.nan, index=ch.index)
            n = pd.to_numeric(ch[c["n"]], errors="coerce") if c["n"] else None
            good = p.notna() & (p > 0) & (p <= 1)
            df = pd.DataFrame({"SNP": snp, "P": p})
            if n is not None:
                df["N"] = n
            df = df[good.to_numpy()].drop_duplicates("SNP")
            df.to_csv(fp, sep="\t", header=False, index=False)
            n_out += len(df)
            if c["a1"] and c["a2"]:
                h = pd.DataFrame({"SNP": snp, "A1": ch[c["a1"]].astype(str).str.upper(), "A2": ch[c["a2"]].astype(str).str.upper(),
                                  "Z": z, "P": p})
                if n is not None:
                    h["N"] = n
                h = h[good.to_numpy() & h["SNP"].isin(hm3).to_numpy()].drop_duplicates("SNP")
                h.to_csv(fh, sep="\t", header=False, index=False)
                n_hm3 += len(h)
    say("[PREPARE] rows in=%d | MAGMA SNPs=%d | HapMap3 rows=%d" % (n_in, n_out, n_hm3))
    SUM["gates"]["KOLLER_SNPS"] = n_out
    return c


# ----------------------------------------------------------------------------
# 2-3. MAGMA gene analysis, munge, h2, rg
# ----------------------------------------------------------------------------
def gene_analysis(c):
    raw = MAGMA_OUT / (TRAIT + "_35_10.genes.raw")
    if raw.exists():
        say("[MAGMA_GENE] present")
        return True
    pfile = MAGMA_OUT / (TRAIT + ".pval.txt")
    if c["n"]:
        pv = [str(pfile), "ncol=N"]
        nn = []
    else:
        pv = [str(pfile)]
        nn = ["N=1000000"]
        say("[MAGMA_GENE] WARNING no per-SNP N column; N=1,000,000 used")
    t0 = time.time()
    say("[MAGMA_GENE] %s running (~35 min, quiet)" % TRAIT)
    rc, o, e = sh(NICE + [MAGMA, "--bfile", G1000, "--pval"] + pv + nn + ["--gene-annot", ANNOT,
                  "--out", MAGMA_OUT / (TRAIT + "_35_10")], timeout=6 * 3600)
    say("[MAGMA_GENE] rc=%d %.0f min %s" % (rc, (time.time() - t0) / 60, "" if rc == 0 else (o + e).strip()[-300:]))
    return raw.exists()


def top_genes():
    sym = {}
    for line in open(GENE_LOC_B37):
        x = line.split()
        if len(x) >= 6:
            sym[x[0]] = x[5]
    res = {}
    for t in (REF, TRAIT):
        f = MAGMA_OUT / (t + "_35_10.genes.out")
        if f.exists():
            g = pd.read_csv(f, sep=r"\s+")
            g["GENE"] = g["GENE"].astype(str)
            res[t] = g
    if TRAIT in res:
        g = res[TRAIT]
        bonf = 0.05 / len(g)
        top = g.sort_values("P").head(15)
        say("[MAGMA_TOP] %s Bonferroni genes=%d; top: %s" % (TRAIT, int((g["P"] < bonf).sum()),
            ", ".join("%s(%.0e)" % (sym.get(e, e), p) for e, p in zip(top["GENE"], top["P"]))))
        if REF in res:
            m = res[REF].merge(g, on="GENE", suffixes=("_ref", "_k"))
            r = np.corrcoef(m["ZSTAT_ref"], m["ZSTAT_k"])[0, 1]
            both = int(((m["P_ref"] < bonf) & (m["P_k"] < bonf)).sum())
            say("[MAGMA_OVERLAP] gene Z correlation with %s r=%.3f; Bonferroni in both=%d (ref %d, Koller %d)" % (
                REF, r, both, int((m["P_ref"] < bonf).sum()), int((m["P_k"] < bonf).sum())))


def ld_prefixes():
    def find(dirname):
        d = REF_UNPACKED / dirname
        for p in sorted(d.rglob("*.l2.ldscore.gz")) if d.exists() else []:
            m = re.match(r"^(.*\.)1\.l2\.ldscore\.gz$", p.name)
            if m:
                return str(p.parent / m.group(1))
        return None
    return find("LDscore"), find("1000G_Phase3_weights_hm3_no_MHC")


def ldsc_block():
    say("")
    say("-" * 72)
    say("[SECTION] 3. munge, h2, rg with GCST90205183 and ENDO_ALL (overlap)")
    ldsc_py = LDSC_DIR / "ldsc.py"
    munge_py = LDSC_DIR / "munge_sumstats.py"
    ss = MUNGED / (TRAIT + ".sumstats.gz")
    if not ss.exists():
        src = HM3_DIR / (TRAIT + ".tsv")
        hdr = open(src).readline().rstrip("\n").split("\t")
        n_med = pd.read_csv(src, sep="\t", usecols=["N"])["N"].median() if "N" in hdr else float("nan")
        cmd = [LDSC_PY, munge_py, "--sumstats", src, "--snp", "SNP", "--a1", "A1", "--a2", "A2", "--p", "P",
               "--signed-sumstats", "Z,0", "--merge-alleles", MERGE_ALLELES, "--chunksize", "500000", "--out", MUNGED / TRAIT]
        cmd += ["--N-col", "N"] if np.isfinite(n_med) else ["--N", "1000000"]
        rc, o, e = sh(cmd, timeout=3600)
        say("[MUNGE] rc=%d %s" % (rc, "" if rc == 0 else (o + e).strip()[-300:]))
    if not ss.exists():
        return
    ref, w = ld_prefixes()
    rc, o, e = sh([LDSC_PY, ldsc_py, "--h2", ss, "--ref-ld-chr", ref, "--w-ld-chr", w, "--out", WORK / "h2_koller"], timeout=3600)
    t = Path(str(WORK / "h2_koller") + ".log").read_text(errors="ignore") if Path(str(WORK / "h2_koller") + ".log").exists() else o + e
    for pat in [r"Total Observed scale h2:.*", r"Lambda GC:.*", r"Mean Chi\^2:.*", r"Intercept:.*", r"Ratio.*"]:
        m = re.search(pat, t)
        if m:
            say("    [H2] " + m.group(0).strip())
    others = [x for x in [REF, "ENDO_ALL"] if (MUNGED / (x + ".sumstats.gz")).exists()]
    rc, o, e = sh([LDSC_PY, ldsc_py, "--rg", ",".join(str(MUNGED / (x + ".sumstats.gz")) for x in [TRAIT] + others),
                   "--ref-ld-chr", ref, "--w-ld-chr", w, "--out", WORK / "rg_koller"], timeout=7200)
    t = Path(str(WORK / "rg_koller") + ".log").read_text(errors="ignore") if Path(str(WORK / "rg_koller") + ".log").exists() else o + e
    if "Summary of Genetic Correlation Results" in t:
        block = t.split("Summary of Genetic Correlation Results")[1].strip().splitlines()
        hdr = block[0].split()
        for line in block[1:]:
            x = line.split()
            if len(x) != len(hdr):
                break
            d = dict(zip(hdr, x))
            say("    [RG] %s ~ %-14s rg=%s se=%s gcov_int=%s (%s)" % (TRAIT, Path(d["p2"]).name.replace(".sumstats.gz", ""),
                d.get("rg"), d.get("se"), d.get("gcov_int"), d.get("gcov_int_se")))
    else:
        say("    [RG] failed rc=%d %s" % (rc, (o + e).strip()[-200:]))


# ----------------------------------------------------------------------------
# 4. Re-runs
# ----------------------------------------------------------------------------
def gtex_tissues():
    hdr = gzip.open(GTEX_FILE, "rt").read(200000).splitlines()[2].split("\t")
    return [c for c in hdr[2:] if not c.startswith("Cells")]


def rerun_g6():
    say("")
    say("-" * 72)
    say("[SECTION] 4a. G6 tissues (GTEx v8) and X4, Koller vs GCST90205183")
    tissues = gtex_tissues()
    lab = {"T%02d" % (i + 1): t for i, t in enumerate(tissues)}
    for trait in (REF, TRAIT):
        r = magma_gp(trait, COV_G6, "Average", "G6")
        if not r:
            continue
        df = pd.DataFrame([{"tissue": lab[k], **r[k]} for k in lab if k in r]).sort_values("p").reset_index(drop=True)
        df["q"] = bh(df["p"].values)
        rep = df[df["tissue"].isin(REPRO)].copy()
        rep["holm"] = holm(rep["p"].values)
        say("    [G6] %-13s top: %s" % (trait[:13], "; ".join("%s p=%.1e" % (t, p) for t, p in zip(df["tissue"].head(5), df["p"].head(5)))))
        say("    [G6] %-13s reproductive: %s" % (trait[:13], " | ".join("%s r%d H=%.2g" % (t, i + 1, h) for i, (t, h) in
                                                                      zip(rep.index, zip(rep["tissue"], rep["holm"])))))
        df.to_csv(OUT / ("G6_%s.csv" % trait), index=False)
    fa = COV_X4 / "X4a.covar.txt"
    for trait in (REF, TRAIT):
        r = magma_gp(trait, fa, "Average,SMOOTH", "X4a")
        if r and "MULL" in r:
            say("    [X4a] %-13s MULL|Avg+SMOOTH beta=%.4f p=%.2e | UTERUS beta=%.4f p=%.2e" % (
                trait[:13], r["MULL"]["beta"], r["MULL"]["p"], r.get("UTERUS", {}).get("beta", np.nan), r.get("UTERUS", {}).get("p", np.nan)))
        parts = []
        for k in ("T1", "T2", "T3"):
            rr = magma_gp(trait, COV_X4 / ("X4b_%s.covar.txt" % k), "Average,SMOOTH", "X4b_" + k)
            if rr and "MULL" in rr:
                parts.append("%s BETA_STD=%.3f p=%.1e" % (k, rr["MULL"]["beta_std"], rr["MULL"]["p"]))
        say("    [X4b] %-13s %s" % (trait[:13], " | ".join(parts)))


def rerun_g4():
    say("")
    say("-" * 72)
    say("[SECTION] 4b. G4 HECA (control donors) M1 specificity|Average and M3 log-mean (no conditioning)")
    types = list(pd.read_csv(HECA_GS / "specificity_by_type.csv.gz", index_col=0, nrows=1).columns)
    lab = {"T%02d" % (i + 1): t for i, t in enumerate(types)}
    for mname, covf, cond in (("M1", COV_C5 / "C5_M1_specificity.covar.txt", "Average"),
                              ("M3", COV_C5 / "C5_M23_logmean.covar.txt", None)):
        for trait in (REF, TRAIT):
            r = magma_gp(trait, covf, cond, "C5_" + mname)
            if not r:
                continue
            df = pd.DataFrame([{"type": lab[k], **r[k]} for k in lab if k in r]).sort_values("p").reset_index(drop=True)
            ds = df[df["type"].isin(DSTROMAL)].copy()
            ds["holm"] = holm(ds["p"].values)
            say("    [G4_%s] %-13s top3: %s | dStromal: %s" % (mname, trait[:13],
                "; ".join("%s(p=%.1e)" % (t, p) for t, p in zip(df["type"].head(3), df["p"].head(3))),
                "; ".join("%s r%d p=%.2g H=%.2g" % (t.replace("dStromal_", "dS_"), i + 1, p, h)
                          for i, t, p, h in zip(ds.index, ds["type"], ds["p"], ds["holm"]))))
            df.to_csv(OUT / ("G4_%s_%s.csv" % (mname, trait)), index=False)
            if trait == TRAIT and mname == "M1":
                SUM["results"]["G4_M1_koller_dS_holm_min"] = float(ds["holm"].min()) if len(ds) else None


def unit_profiles(M, G, lineage, context):
    sub = G[(G["lineage_ref"] == lineage) & (G["context"] == context) & (G["agree"])]
    prof = {}
    for pid, g in sub.groupby("patient"):
        if int(g["n_cells"].sum()) < 10:
            continue
        v = M[g.index.to_numpy()].sum(axis=0)
        prof[(pid, context)] = np.log2(v / max(v.sum(), 1) * 1e6 + 1)
    return prof


def write_contrast(A_list, B_list, genes, mapped, sym2, tag):
    A = np.vstack(A_list).mean(axis=0)
    B = np.vstack(B_list).mean(axis=0)
    ex = (np.maximum(2 ** A - 1, 2 ** B - 1) >= 1.0) & mapped
    g = np.asarray(genes)[ex]
    covf = WORK / ("G1_%s.covar.txt" % tag)
    pd.DataFrame({"GENE": [sym2[x] for x in g], "D": (A - B)[ex], "DN": (B - A)[ex], "AVG": ((A + B) / 2)[ex]}).to_csv(
        covf, sep="\t", index=False, float_format="%.6g")
    return covf


def rerun_g1(sym2):
    say("")
    say("-" * 72)
    say("[SECTION] 4c. G1 primary (peritoneal lesion vs eutopic), Koller vs GCST90205183, Holm over lineages")
    rows = []
    for lin in LINEAGES:
        f = COVAR_G1 / ("LES_vs_EU__%s.covar.txt" % lin)
        for trait in (REF, TRAIT):
            r = magma_gp(trait, f, "AVG", "G1_" + lin)
            if r and "D" in r and "DN" in r:
                rows.append({"trait": trait, "lineage": lin, "beta": r["D"]["beta"], "se": r["D"]["se"],
                             "p_greater": r["D"]["p"], "p_less": r["DN"]["p"]})
    df = pd.DataFrame(rows)
    for trait, g in df.groupby("trait"):
        hl, hg = holm(g["p_less"].values), holm(g["p_greater"].values)
        say("    [G1] %-13s " % trait[:13] + " ".join("%s:%+.3f(pL=%.3f,HL=%.2f,HG=%.2f)" % (l[:5], b, pl, a, c)
                                                     for l, b, pl, a, c in zip(g["lineage"], g["beta"], g["p_less"], hl, hg)))
    df.to_csv(OUT / "G1_koller.csv", index=False)
    M = np.load(PB213 / "pseudobulk_counts.npy", mmap_mode="r")
    G = pd.read_csv(PB213 / "pseudobulk_groups.csv").reset_index(drop=True)
    genes = list(pd.read_csv(PB213 / "genes.csv")["symbol"].astype(str))
    mapped = np.array([x in sym2 for x in genes])
    pa = unit_profiles(M, G, "lymphocyte", "peritoneal_lesion")
    pb = unit_profiles(M, G, "lymphocyte", "eutopic")
    base = df[(df["trait"] == TRAIT) & (df["lineage"] == "lymphocyte")]
    if len(pa) >= 3 and len(pb) >= 3 and len(base):
        units = list(pa.values()) + list(pb.values())
        lab = np.array([1] * len(pa) + [0] * len(pb))
        rng = np.random.default_rng(PERM_SEED)
        ts = []
        t0 = time.time()
        for b in range(N_PERM):
            pl = rng.permutation(lab)
            covf = write_contrast([units[j] for j in range(len(units)) if pl[j] == 1],
                                  [units[j] for j in range(len(units)) if pl[j] == 0], genes, mapped, sym2, "perm")
            r = magma_gp(TRAIT, covf, "AVG", "G1_perm")
            if r and "D" in r:
                ts.append(r["D"]["beta"] / r["D"]["se"])
        tobs = base["beta"].iloc[0] / base["se"].iloc[0]
        ts = np.array(ts)
        pe = (1 + np.sum(ts <= tobs)) / (1 + len(ts)) if len(ts) else float("nan")
        say("    [G1_PERM] %s lymphocyte t_obs=%.2f permutation p_less=%.3f (B=%d) vs parametric %.4f | %.0f min" % (
            TRAIT, tobs, pe, len(ts), base["p_less"].iloc[0], (time.time() - t0) / 60))


# ----------------------------------------------------------------------------
def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.23 | Koller 2026 EUR (without 23andMe) as second GWAS | " + TS)
    say("=" * 72)
    t0 = time.time()
    links = zenodo_links()
    say("[ZENODO] files listed: %d" % len(links))
    rd = download(README_NAME, links)
    if rd:
        for line in rd.read_text(errors="ignore").splitlines()[:60]:
            say("    [README] " + line[:200])
    kfile = download(KFILE_NAME, links)
    if kfile is None:
        say("[FATAL] GWAS file not available")
    else:
        try:
            c = prepare(kfile)
            ok = gene_analysis(c)
            SUM["gates"]["KOLLER_GENES_RAW"] = ok
            top_genes()
        except Exception as e:
            import traceback
            say("[STEP_ERROR] prepare/gene %s: %s" % (type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1200:])
        for name, fn in (("ldsc", ldsc_block), ("g6", rerun_g6), ("g4", rerun_g4), ("g1", lambda: rerun_g1(entrez_unique()))):
            try:
                fn()
            except Exception as e:
                import traceback
                say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
                say(traceback.format_exc()[-1000:])
    with open(OUT / "p9_v023_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    for k, v in SUM["gates"].items():
        say("[GATE] %-24s %s" % (k, v))
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.14 before this run (larger, partially overlapping GWAS)")
    say("[DONE] PAPER9_V023_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
