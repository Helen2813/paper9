#!/usr/bin/env python3
# Paper 9 | v0.0.67 | non-European replication of the stromal signal (rules CHANGELOG v0.2.57, fixed before the run)
#   N0 Koller et al. Zenodo: endometriosis_EAS / AFR / AMR summary statistics + README (case numbers printed)
#   N1 MAGMA 1000 Genomes reference panels for EAS / AFR / AMR (links from the MAGMA page), gene annotation per panel,
#      rsID mapping (as v023), gene analysis per ancestry (KOLLER_EAS, KOLLER_AFR, KOLLER_AMR)
#   N2 GSE266265 (60 donors) seven marker families from the cached v050 DE; MAGMA per ancestry; H-NE-EAS (eF3 size-matched
#      up-only), H-NE-ORIGIN (best endometrial-stroma type above perivascular, mesothelial-like and epithelial types);
#      AFR / AMR and a Stouffer combination of the three non-European estimates reported
# Runs in ~/Desktop/paper9/.venv, low heat (MAGMA single thread, nice 10). Downloads ~0.4 GB of summary statistics plus
# the reference panels; ~2-3 h in total. Re-running skips finished steps.

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
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RAW = REPO / "data" / "raw"
INTERIM = REPO / "data" / "interim"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
MDIR = TOOLS / "magma"
GENE_LOC_B37 = MDIR / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
GCST_H = RAW / "gwas" / "GCST90205183" / "harmonised" / "GCST90205183.h.tsv.gz"
KDIR = RAW / "gwas" / "Koller2026"
DE_DIR = PROC / "genesets" / "gse266265_v050_de"
ZENODO_API = "https://zenodo.org/api/records/18983492"
MAGMA_PAGE = "https://ctg.thebluebus.nl/software/magma"
MAGMA_FALLBACK = "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_%s.zip"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "review" / ("p9_v067_" + TS)
WORK = INTERIM / "magma_v067"
LOG = REPO / "logs" / ("p9_v067_non_european_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.57)
ANC = ["EAS", "AFR", "AMR"]
LFC, PADJ, MIN_SET, N_MATCH = 1.2, 0.01, 20, 400
FAMILIES = ["ABS", "UP", "DOWN", "ABSN", "UPN", "DOWNN", "SPEC"]
GROUPS = [("endometrial stroma", r"^(eF\d|dS\d|fib C7)"), ("perivascular / MSC-like", r"^(Prv|VSMC)"), ("mesothelial-like", r"^mesothelial"),
          ("epithelium", r"^(glandular|lumenal|ciliated|mid-secretory|MUC5B|TP63)"), ("endothelium", r"^(EC-|LEC)")]
KEY_TYPE = "eF3"
REF_TRAIT = "KOLLER_EUR"

for d in (OUT, WORK, LOG.parent, KDIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V067_NON_EUROPEAN", "timestamp_utc": TS, "results": {}}


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


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def gene_loc():
    sym2ent, dup = {}, set()
    for line in open(GENE_LOC_B37):
        t = line.split()
        if len(t) >= 6 and t[1] in [str(c) for c in range(1, 23)]:
            if t[5] in sym2ent:
                dup.add(t[5])
            sym2ent[t[5]] = t[0]
    for d in dup:
        sym2ent.pop(d, None)
    return sym2ent


def curl(url, dst, timeout=4 * 3600):
    part = dst.with_name(dst.name + ".part")
    rc, o, e = sh(["curl", "-L", "-sS", "--fail", "--retry", "3", "--connect-timeout", "30", "-C", "-", "-A", "paper9", "-o", part, url],
                  timeout=timeout)
    if rc == 0 and part.exists() and part.stat().st_size > 0:
        part.rename(dst)
        return True
    say("  [DL] FAILED %s rc=%s %s" % (url[:90], rc, (e or o).strip()[:160]))
    return False


# ----------------------------------------------------------------------------
# N0 Koller non-European files
# ----------------------------------------------------------------------------
def part_n0():
    say("")
    say("-" * 72)
    say("[SECTION] N0 Koller et al. non-European summary statistics (Zenodo 10.5281/zenodo.18983492, CC BY-NC 4.0)")
    rc, o, e = sh(["curl", "-sSL", "--fail", "--max-time", "60", "-A", "paper9", ZENODO_API], timeout=120)
    links = {}
    if rc == 0:
        for f in json.loads(o).get("files", []):
            links[f.get("key")] = (f.get("links", {}).get("self"), f.get("size", 0))
    files = {}
    for name in ["README_GWAS.txt"] + ["endometriosis_%s.txt.gz" % a for a in ANC]:
        dst = KDIR / name
        url, size = links.get(name, (None, 0))
        if not (dst.exists() and (not size or dst.stat().st_size == size)):
            for u in [x for x in [url, "https://zenodo.org/records/18983492/files/%s?download=1" % name] if x]:
                if curl(u, dst):
                    break
        if dst.exists():
            files[name] = dst
            say("  [N0] %-34s %.1f MB" % (name, dst.stat().st_size / 1e6))
    rd = files.get("README_GWAS.txt")
    if rd:
        for line in rd.read_text(errors="ignore").splitlines():
            if re.search(r"EAS|AFR|AMR|East|African|Latin|Hispanic|case|control|N_|sample size|column|build|GRCh", line, re.I):
                say("  [README] " + line.strip()[:200])
    return files


# ----------------------------------------------------------------------------
# N1 reference panels, annotation, rsID mapping, gene analysis
# ----------------------------------------------------------------------------
def magma_links():
    rc, o, e = sh(["curl", "-sSL", "--fail", "--max-time", "60", "-A", "paper9", MAGMA_PAGE], timeout=120)
    out = {}
    if rc != 0:
        return out
    cands = []
    for m in re.finditer(r'<a[^>]*href="([^"]+)"[^>]*>(.*?)</a>', o, re.I | re.S):
        href, text = m.group(1), re.sub(r"<[^>]+>", " ", m.group(2)).lower()
        if not re.search(r"\.zip|download|/s/", href, re.I):
            continue
        before = re.sub(r"<[^>]+>", "|", o[max(0, m.start() - 120):m.start()]).lower().split("|")
        label = next((t.strip() for t in reversed(before) if t.strip()), "")
        url = href if href.startswith("http") else "https://ctg.thebluebus.nl" + (href if href.startswith("/") else "/" + href)
        cands.append((url, href.lower() + " " + text, label))
    for a in ("eas", "afr", "amr"):
        hit = next((u for u, ht, lb in cands if "g1000_" + a in ht), None)
        if hit is None:
            rx = {"eas": r"east[ -]asian", "afr": r"^african|[^n] african", "amr": r"american|admixed"}[a]
            hit = next((u for u, ht, lb in cands if re.search(rx, lb) and not (a == "afr" and "american" in lb)), None)
        if hit:
            out[a] = hit
    return out


def reference(anc, links):
    a = anc.lower()
    base = MDIR / ("ref_" + a)
    bed = list(base.rglob("g1000_%s.bed" % a)) if base.exists() else []
    if bed:
        return bed[0].with_suffix("")
    z = MDIR / ("ref_%s.zip" % a)
    if z.exists() and not zipfile.is_zipfile(z):
        z.unlink()
    if not z.exists():
        urls = []
        u0 = links.get(a)
        if u0:
            urls.append(u0)
            if "/s/" in u0 and not u0.rstrip("/").endswith("download"):
                urls.append(u0.rstrip("/") + "/download")
        urls.append(MAGMA_FALLBACK % a)
        ok = False
        for u in urls:
            if curl(u, z) and zipfile.is_zipfile(z):
                ok = True
                break
            if z.exists():
                z.unlink()
        if not ok:
            return None
        say("  [REF] %s panel %.0f MB" % (anc, z.stat().st_size / 1e6))
    try:
        with zipfile.ZipFile(z) as zf:
            zf.extractall(base)
    except Exception as e:
        say("  [REF] %s unzip failed %s" % (anc, type(e).__name__))
        return None
    bed = list(base.rglob("g1000_%s.bed" % a))
    return bed[0].with_suffix("") if bed else None


def pick(cols, cands):
    low = {c.lower(): c for c in cols}
    for c in cands:
        if c in low:
            return low[c]
    return None


def detect(kfile):
    head = pd.read_csv(kfile, sep=None, engine="python", nrows=20000, compression="infer")
    cols = list(head.columns)
    c = {"id": pick(cols, ["rsid", "snp", "markername", "variant_id", "snpid", "id", "rs", "rs_id"]),
         "chr": pick(cols, ["chr", "chrom", "chromosome", "#chrom", "#chr", "hg38_chr", "chr_b38"]),
         "pos": pick(cols, ["pos", "bp", "position", "base_pair_location", "pos_b38", "hg38_pos", "bp_b38"]),
         "p": pick(cols, ["p-value", "p", "pval", "p_value", "pvalue", "p.value"]),
         "n": pick(cols, ["n", "weight", "n_total", "neff", "totaln", "n_eff", "samplesize"])}
    ids = head[c["id"]].astype(str) if c["id"] else pd.Series([], dtype=str)
    rs_share = float(ids.str.startswith("rs").mean()) if len(ids) else 0.0
    say("  [FORMAT] %s: %s | rsID share %.3f" % (kfile.name, ", ".join("%s=%s" % (k, v) for k, v in c.items()), rs_share))
    return c, rs_share


def key_of(chrom, pos):
    ch = pd.Series(chrom).astype(str).str.replace("chr", "", regex=False).str.upper().replace({"X": "23", "Y": "24"})
    ch = pd.to_numeric(ch, errors="coerce")
    ps = pd.to_numeric(pd.Series(pos), errors="coerce")
    return (ch * 10_000_000_000 + ps).to_numpy(dtype="float64")


def build_map(build, bim):
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
        for ch in pd.read_csv(str(bim) + ".bim", sep=r"\s+", header=None, usecols=[0, 1, 3], chunksize=2_000_000):
            r = ch[1].astype(str)
            ok = r.str.startswith("rs")
            keys.append(key_of(ch.loc[ok, 0], ch.loc[ok, 3]))
            rs.append(pd.to_numeric(r[ok].str[2:], errors="coerce").to_numpy(dtype="float64"))
    k, r = np.concatenate(keys), np.concatenate(rs)
    m = np.isfinite(k) & np.isfinite(r)
    o = np.argsort(k[m], kind="mergesort")
    return k[m][o], r[m][o]


def lookup(kmap, rmap, k):
    i = np.clip(np.searchsorted(kmap, k), 0, len(kmap) - 1)
    hit = kmap[i] == k
    out = np.full(len(k), np.nan)
    out[hit] = rmap[i[hit]]
    return out


def prepare(kfile, trait, bfile):
    pfile = MAGMA_OUT / (trait + ".pval.txt")
    c, rs_share = detect(kfile)
    if pfile.exists():
        return c, pfile
    if c["p"] is None or (c["id"] is None and (c["chr"] is None or c["pos"] is None)):
        say("  [PREPARE] %s: no p-value / SNP columns" % trait)
        return None, None
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
            km, rm = build_map(b, bfile)
            share = float(np.isfinite(lookup(km, rm, ks)).mean())
            say("  [BUILD] %s %s match share %.3f" % (trait, b, share))
            if best is None or share > best[0]:
                best = (share, b, km, rm)
        _, b, kmap, rmap = best
    with gzip.open(kfile, "rt") as f0:
        first = f0.readline()
    sep = "\t" if "\t" in first else r"\s+"
    n_in = n_out = 0
    tmp = pfile.with_name(pfile.name + ".part")
    with open(tmp, "w") as fp:
        fp.write("SNP\tP\tN\n" if c["n"] else "SNP\tP\n")
        for ch in pd.read_csv(kfile, sep=sep, chunksize=500000, compression="infer", low_memory=False):
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
            df = pd.DataFrame({"SNP": snp, "P": p.clip(lower=1e-300)})
            if c["n"]:
                df["N"] = pd.to_numeric(ch[c["n"]], errors="coerce")
            df = df[(p.notna() & (p > 0) & (p <= 1)).to_numpy()].drop_duplicates("SNP")
            if c["n"]:
                df = df[df["N"].notna()]
            df.to_csv(fp, sep="\t", header=False, index=False)
            n_out += len(df)
    tmp.rename(pfile)
    say("  [PREPARE] %s rows in %d | MAGMA SNPs %d" % (trait, n_in, n_out))
    return c, pfile


def gene_analysis(anc, files, links):
    trait = "KOLLER_" + anc
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    if raw.exists():
        say("  [MAGMA_GENE] %s present" % trait)
        return raw
    kf = files.get("endometriosis_%s.txt.gz" % anc)
    if kf is None:
        say("  [N1] %s summary statistics missing" % anc)
        return None
    bfile = reference(anc, links)
    if bfile is None:
        say("  [N1] %s reference panel missing (MAGMA page links: %s)" % (anc, links.get(anc.lower(), "none")))
        return None
    annot = MAGMA_OUT / ("annot_%s_35_10" % anc.lower())
    if not Path(str(annot) + ".genes.annot").exists():
        rc, o, e = sh(NICE + [MAGMA, "--annotate", "window=35,10", "--snp-loc", str(bfile) + ".bim", "--gene-loc", GENE_LOC_B37,
                              "--out", annot], timeout=3600)
        say("  [ANNOT] %s rc=%d" % (anc, rc))
    c, pfile = prepare(kf, trait, bfile)
    if pfile is None:
        return None
    pv = [str(pfile), "ncol=N"] if c["n"] else [str(pfile)]
    nn = [] if c["n"] else ["N=100000"]
    if not c["n"]:
        say("  [MAGMA_GENE] WARNING %s has no per-SNP N column; N=100,000 used" % trait)
    t0 = time.time()
    say("  [MAGMA_GENE] %s running (20-40 min, quiet)" % trait)
    rc, o, e = sh(NICE + [MAGMA, "--bfile", bfile, "--pval"] + pv + nn + ["--gene-annot", str(annot) + ".genes.annot",
                  "--out", MAGMA_OUT / (trait + "_35_10")], timeout=8 * 3600)
    say("  [MAGMA_GENE] %s rc=%d %.0f min%s" % (trait, rc, (time.time() - t0) / 60, "" if raw.exists() else " " + (o + e).strip()[-300:]))
    if raw.exists():
        pfile.unlink()
        lg = Path(str(MAGMA_OUT / (trait + "_35_10")) + ".log")
        if lg.exists():
            for line in lg.read_text(errors="ignore").splitlines():
                if re.search(r"SNPs|genes|sample size", line, re.I) and len(line) < 160:
                    say("    [MAGMA_LOG] " + line.strip())
        return raw
    return None


# ----------------------------------------------------------------------------
# N2 GSE266265 families and tests
# ----------------------------------------------------------------------------
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


def magma_covar(raw, cf, tag, cond):
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


def group_of(label):
    for g, rx in GROUPS:
        if re.search(rx, str(label)):
            return g
    return "immune / other"


def build(DE, CPM, sym2ent):
    meta, cols = {}, {}
    uni = sorted({g for g in DE["gene"].astype(str) if g in sym2ent})
    uset = set(uni)
    for g, d in DE.groupby("group", sort=False):
        d = d.assign(gene=d["gene"].astype(str))
        sig = d["padj"] < PADJ
        sel = {"ABS": d.loc[(d["lfc"].abs() > LFC) & sig, "gene"], "UP": d.loc[(d["lfc"] > LFC) & sig, "gene"],
               "DOWN": d.loc[(d["lfc"] < -LFC) & sig, "gene"],
               "ABSN": d.reindex(d["score"].abs().sort_values(ascending=False).index)["gene"].head(N_MATCH),
               "UPN": d.sort_values("score", ascending=False)["gene"].head(N_MATCH),
               "DOWNN": d.sort_values("score", ascending=True)["gene"].head(N_MATCH)}
        for fam, genes in sel.items():
            gl = [x for x in genes if x in uset]
            if len(gl) >= MIN_SET:
                k = "V%04d" % (len(meta) + 1)
                meta[k], cols[k] = (fam, g), gl
    T = pd.DataFrame({k: np.isin(uni, gl).astype(int) for k, gl in cols.items()}, index=[sym2ent[g] for g in uni])
    T = T[~T.index.duplicated()]
    C = CPM.copy()
    C.index = C.index.astype(str)
    C = C.groupby(level=0).mean()
    C = C.loc[[g for g in C.index if g in sym2ent]]
    S = C.div(C.sum(axis=1), axis=0)
    spec, sd = {}, {}
    for i, c in enumerate(S.columns):
        k = "S%03d" % (i + 1)
        spec[k], sd[k] = c, S[c].values
    sd["Average"] = np.log2(C.mean(axis=1).values + 1.0)
    Ts = pd.DataFrame(sd, index=[sym2ent[g] for g in C.index])
    Ts = Ts[~Ts.index.duplicated()]
    return T, Ts, meta, spec


def ranked(res, meta, spec, fam):
    it = ([(spec[k], v) for k, v in res.items() if k in spec] if fam == "SPEC"
          else [(meta[k][1], v) for k, v in res.items() if k in meta and meta[k][0] == fam])
    it.sort(key=lambda x: x[1]["p"])
    q = bh([v["p"] for _, v in it]) if it else []
    return [(l, v, qq) for (l, v), qq in zip(it, q)]


def part_n2(traits, sym2ent):
    say("")
    say("-" * 72)
    say("[SECTION] N2 GSE266265 (60 donors): seven families per ancestry (top 3: p / BH q)")
    f_de, f_cpm = DE_DIR / "de_all.csv.gz", DE_DIR / "mean_cpm_by_type.csv.gz"
    if not f_de.exists():
        say("[N2] v050 DE cache missing")
        return
    DE, CPM = pd.read_csv(f_de), pd.read_csv(f_cpm, index_col=0)
    T, Ts, meta, spec = build(DE, CPM, sym2ent)
    cf, cfs = WORK / "covar_sets.txt", WORK / "covar_spec.txt"
    T.to_csv(cf, sep="\t", index_label="GENE")
    Ts.to_csv(cfs, sep="\t", float_format="%.6g", index_label="GENE")
    res = {}
    for t in traits:
        raw = MAGMA_OUT / (t + "_35_10.genes.raw")
        r0 = magma_covar(raw, cf, "sets_%s" % t, [])
        r1 = magma_covar(raw, cfs, "spec_%s" % t, ["Average"])
        if r0 and r1:
            res[t] = {fam: ranked(r1 if fam == "SPEC" else r0, meta, spec, fam) for fam in FAMILIES}
    rows = []
    for t, F in res.items():
        say("  %s" % t)
        for fam in FAMILIES:
            say("    %-6s %s" % (fam, "; ".join("%s %.0e/%.2f" % (l[:20], v["p"], q) for l, v, q in F[fam][:3])))
            for i, (l, v, q) in enumerate(F[fam]):
                rows.append({"trait": t, "family": fam, "rank": i + 1, "label": l, "group": group_of(l), "p": v["p"], "q": q,
                             "beta_std": v["beta_std"], "beta": v["beta"], "se": v["se"]})
    A = pd.DataFrame(rows)
    A.to_csv(OUT / "n2_rankings.csv", index=False)
    say("")
    say("-" * 72)
    say("[SECTION] N2 verdicts")
    zs = []
    for t in res:
        d = A[(A["trait"] == t) & (A["family"] == "UPN")]
        e = d[d["label"] == KEY_TYPE]
        br = {}
        for _, r in d.iterrows():
            br.setdefault(r["group"], (int(r["rank"]), r["label"], r["p"]))
        rs = br.get("endometrial stroma", (99,))[0]
        oth = [br[g][0] for g in ("perivascular / MSC-like", "mesothelial-like", "epithelium") if g in br]
        origin = bool(oth) and all(rs < x for x in oth)
        pe = float(e["p"].iloc[0]) if len(e) else np.nan
        say("  [N2] %-11s UPN %s p %.1e (rank %s of %d) | best stromal type %s rank %d | perivascular / mesothelial-like / epithelial best ranks %s -> origin %s" % (
            t, KEY_TYPE, pe, int(e["rank"].iloc[0]) if len(e) else "-", len(d), br.get("endometrial stroma", ("-", "-"))[1], rs,
            ", ".join(str(x) for x in oth), "yes" if origin else "no"))
        SUM["results"][t] = {"eF3_UPN_p": pe, "origin": origin, "stroma_rank": rs}
        if t != REF_TRAIT and np.isfinite(pe):
            from scipy.stats import norm
            zs.append(norm.isf(pe))
    eas = SUM["results"].get("KOLLER_EAS", {})
    v1 = "NOT EVALUABLE" if not eas else ("SUPPORTED" if eas.get("eF3_UPN_p", 1) < 0.05 else "NOT SUPPORTED")
    v2 = "NOT EVALUABLE" if not eas else ("SUPPORTED" if eas.get("origin") else "NOT SUPPORTED")
    say("[H_NE_EAS_VERDICT] KOLLER_EAS: %s size-matched up-only enrichment p < 0.05 -> %s" % (KEY_TYPE, v1))
    say("[H_NE_ORIGIN_VERDICT] KOLLER_EAS: best endometrial-stroma type above perivascular, mesothelial-like and epithelial types -> %s" % v2)
    if zs:
        from scipy.stats import norm
        zc = sum(zs) / math.sqrt(len(zs))
        say("[N2] Stouffer combination of the non-European %s UPN p-values (%d ancestries): z %.2f, p %.1e" % (KEY_TYPE, len(zs), zc, norm.sf(zc)))
        SUM["results"]["stouffer"] = {"z": zc, "p": float(norm.sf(zc)), "k": len(zs)}


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.67 | non-European replication of the stromal signal (Koller EAS / AFR / AMR) | " + TS)
    say("=" * 72)
    files = part_n0()
    say("")
    say("-" * 72)
    say("[SECTION] N1 reference panels and MAGMA gene analysis per ancestry")
    links = magma_links()
    say("  [REF] MAGMA page links: %s" % (", ".join("%s %s" % (k, v[-60:]) for k, v in links.items()) or "none found (fallback URLs used)"))
    traits = [REF_TRAIT] if (MAGMA_OUT / (REF_TRAIT + "_35_10.genes.raw")).exists() else []
    for a in ANC:
        if gene_analysis(a, files, links) is not None:
            traits.append("KOLLER_" + a)
    if len(traits) <= 1:
        say("[FATAL] no non-European gene analysis available")
        _fh.close()
        return
    part_n2(traits, gene_loc())
    (OUT / "summary.json").write_text(json.dumps(SUM, indent=1, default=str), encoding="utf-8")
    say("")
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V067_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
