#!/usr/bin/env python3
# Paper 9 | v0.0.19 | (1) G6: GTEx v8 tissue-level MAGMA gene-property ("outside the uterus" arm)
#                    (2) G1-S post hoc diagnostics of the G1 eutopic tilt: T2D control (S1), conditioning on
#                        gnomAD LOEUF + GTEx expression breadth (S2), on control-trait gene z (S3), lymphocyte
#                        label composition (S4).
# Rules: CHANGELOG v0.2.10 (fixed before this run). Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single
# thread, nice 10). Expected time: T2D gene analysis 20-40 min (once), everything else ~10-20 min.

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
FINNGEN_DIR = RAW / "finngen"
FINNGEN_MANIFEST = FINNGEN_DIR / "finngen_R13_manifest.tsv"
COVAR_G1 = PROC / "magma_covar" / "g1"
PB_GROUPS = PROC / "gse213216_pseudobulk_clusters" / "pseudobulk_groups.csv"

GTEX_URLS = [
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
    "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
]
GTEX_FILE = EXT / "gtex" / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
GNOMAD_URLS = [
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/v4.1/constraint/gnomad.v4.1.constraint_metrics.tsv",
]
GNOMAD_DIR = EXT / "gnomad"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "g6_g1s" / ("p9_v019_" + TS)
WORK = INTERIM / "magma_v019"
LOG = REPO / "logs" / ("p9_v019_tissue_g1s_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Pre-registered (CHANGELOG v0.2.10)
G6_PRIMARY = "GCST90205183"
G6_TRAITS = ["GCST90205183", "ENDO_ALL", "ENDO_OVARY", "ENDO_PERITONEUM", "ENDO_DEEP", "LEIOMYOMA", "HYPERTENSION", "T2D"]
REPRO = ["Uterus", "Ovary", "Vagina", "Fallopian Tube", "Cervix - Ectocervix", "Cervix - Endocervix"]
ARTERY = ["Artery - Aorta", "Artery - Coronary", "Artery - Tibial"]
LINEAGES = ["stromal", "epithelial", "endothelial", "myeloid", "lymphocyte"]
G1S_TRAITS = ["GCST90205183", "HYPERTENSION", "T2D", "LEIOMYOMA", "ENDO_ALL"]

for d in (OUT, WORK, LOG.parent, GTEX_FILE.parent, GNOMAD_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V019_TISSUE_G1S", "timestamp_utc": TS, "gates": {}, "results": {}}


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


def download(urls, dst, label):
    dst = Path(dst)
    if dst.exists() and dst.stat().st_size > 1000:
        say("[DL] %-34s SKIP_PRESENT %.1f MB" % (label, dst.stat().st_size / 1e6))
        return True
    for url in urls:
        part = dst.with_name(dst.name + ".part")
        rc, o, e = sh(["curl", "-L", "-sS", "--fail", "--retry", "3", "--connect-timeout", "30", "-A", "paper9",
                       "-o", part, url], timeout=1800)
        if rc == 0 and part.exists() and part.stat().st_size > 1000:
            part.rename(dst)
            say("[DL] %-34s OK %.1f MB  (%s)" % (label, dst.stat().st_size / 1e6, url.split("/")[2]))
            return True
        say("[DL] %-34s FAILED rc=%s %s" % (label, rc, e.strip()[:120]))
        if part.exists():
            part.unlink()
    return False


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
    # Gene-property analysis, direction=greater, each non-conditioned covariate tested marginally
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    if not raw.exists():
        return None, "no genes.raw"
    outp = WORK / ("%s__%s" % (trait, tag))
    gsa = Path(str(outp) + ".gsa.out")
    err = ""
    for key in ("condition-hide=", "condition="):
        if gsa.exists():
            gsa.unlink()
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--gene-covar", covf, "--model", key + cond,
                              "direction=greater", "--out", outp], timeout=3 * 3600)
        if rc == 0 and gsa.exists():
            rows = [l for l in open(gsa) if not l.startswith("#") and l.strip()]
            hdr = rows[0].split()
            res = {}
            for l in rows[1:]:
                d = dict(zip(hdr, l.split()))
                name = d.get("FULL_NAME", d.get("VARIABLE"))
                try:
                    res[name] = {"beta": float(d["BETA"]), "se": float(d["SE"]), "p": float(d["P"]),
                                 "ngenes": int(d["NGENES"])}
                except (KeyError, ValueError):
                    continue
            return res, ""
        err = (o + e).strip()[-300:]
    say("    [MAGMA_FAIL] %s %s %s" % (trait, tag, err))
    return None, err


# ----------------------------------------------------------------------------
# 0. T2D gene analysis (second unrelated control; same settings as v013)
# ----------------------------------------------------------------------------
def manifest_counts():
    out = {}
    lines = FINNGEN_MANIFEST.read_text().splitlines()
    ix = {h: i for i, h in enumerate(lines[0].split("\t"))}
    for line in lines[1:]:
        t = line.split("\t")
        try:
            out[t[0]] = (int(t[ix["num_cases"]]), int(t[ix["num_controls"]]))
        except (KeyError, ValueError, IndexError):
            continue
    return out


def run_gene_analysis(trait, code):
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    if raw.exists():
        say("[MAGMA_GENE] %s present" % trait)
        return True
    src = next((p for p in FINNGEN_DIR.glob("*" + code + ".gz") if p.name.endswith("_" + code + ".gz")), None)
    if src is None:
        say("[MAGMA_GENE] %s source missing in %s" % (trait, FINNGEN_DIR))
        return False
    pfile = MAGMA_OUT / (trait + ".pval.txt")
    if not pfile.exists():
        n = 0
        tmp = pfile.with_name(pfile.name + ".part")
        with gzip.open(src, "rt") as f, open(tmp, "w") as o:
            hdr = f.readline().rstrip("\n").split("\t")
            ix = {h: i for i, h in enumerate(hdr)}
            o.write("SNP\tP\n")
            for line in f:
                t = line.rstrip("\n").split("\t")
                rs = t[ix["rsids"]].split(",")[0]
                if rs.startswith("rs"):
                    o.write(rs + "\t" + t[ix["pval"]] + "\n")
                    n += 1
        tmp.rename(pfile)
        say("[MAGMA_PFILE] %s %d SNPs" % (trait, n))
    counts = manifest_counts()
    ncas, ncon = counts[code]
    t0 = time.time()
    say("[MAGMA_GENE] %s running (N=%d; 20-40 min, quiet)" % (trait, ncas + ncon))
    rc, o, e = sh(NICE + [MAGMA, "--bfile", G1000, "--pval", pfile, "N=%d" % (ncas + ncon),
                          "--gene-annot", ANNOT, "--out", MAGMA_OUT / (trait + "_35_10")], timeout=6 * 3600)
    say("[MAGMA_GENE] %s rc=%d %.0f min" % (trait, rc, (time.time() - t0) / 60))
    return raw.exists()


# ----------------------------------------------------------------------------
# 1. Gene annotations: GTEx tissues, tau, LOEUF
# ----------------------------------------------------------------------------
def load_gtex(sym2):
    df = pd.read_csv(GTEX_FILE, sep="\t", skiprows=2, compression="gzip")
    tissues = [c for c in df.columns[2:] if not c.startswith("Cells")]
    df = df[df["Description"].isin(sym2.keys())].copy()
    tpm = df[tissues].astype(float)
    df["_m"] = tpm.mean(axis=1)
    df = df.sort_values("_m", ascending=False).drop_duplicates("Description")
    tpm = df[tissues].astype(float)
    keep = (tpm.values >= 1.0).any(axis=1)
    L = np.log2(tpm[keep] + 1.0)
    L.index = df.loc[keep, "Description"].values
    mx = L.max(axis=1).values
    tau = ((1.0 - L.values / mx[:, None]).sum(axis=1)) / (L.shape[1] - 1)
    say("[GTEX] tissues=%d (cell lines excluded) genes kept=%d (TPM>=1 in >=1 tissue, unique NCBI37.3 symbol)" % (
        len(tissues), len(L)))
    miss = [t for t in REPRO + ARTERY if t not in tissues]
    if miss:
        say("[GTEX_WARN] tissue names not found: " + ", ".join(miss))
    return L, pd.Series(tau, index=L.index), tissues


def load_loeuf(sym2):
    for url in GNOMAD_URLS:
        dst = GNOMAD_DIR / url.split("/")[-1]
        if not download([url], dst, "gnomAD constraint " + ("v2.1.1" if "2.1.1" in url else "v4.1")):
            continue
        df = pd.read_csv(dst, sep="\t", compression="gzip" if dst.suffix in (".bgz", ".gz") else None, low_memory=False)
        col = next((c for c in ["oe_lof_upper", "lof.oe_ci.upper"] if c in df.columns), None)
        if col is None or "gene" not in df.columns:
            say("[GNOMAD] unexpected columns in %s" % dst.name)
            continue
        for flag in ("canonical", "mane_select"):
            if flag in df.columns:
                df = df[df[flag].astype(str).str.lower() == "true"]
                break
        s = pd.to_numeric(df[col], errors="coerce")
        s.index = df["gene"].values
        s = s.dropna()
        s = s[~s.index.duplicated(keep="first")]
        s = s[s.index.isin(sym2.keys())]
        say("[GNOMAD] %s LOEUF genes mapped=%d (median %.2f)" % (dst.name, len(s), s.median()))
        return s
    return None


# ----------------------------------------------------------------------------
# 2. G6: tissue-level gene-property
# ----------------------------------------------------------------------------
def g6(L, tissues, sym2):
    say("")
    say("-" * 72)
    say("[SECTION] G6 (registered v0.2.10): GTEx v8 tissue specificity, MAGMA gene-property | primary " + G6_PRIMARY)
    keys = ["T%02d" % (i + 1) for i in range(len(tissues))]
    lab = dict(zip(keys, tissues))
    cov = pd.DataFrame(L[tissues].values, columns=keys)
    cov["Average"] = L[tissues].mean(axis=1).values
    cov.insert(0, "GENE", [sym2[s] for s in L.index])
    covf = WORK / "G6_gtex_v8_52.covar.txt"
    cov.to_csv(covf, sep="\t", index=False, float_format="%.5g")
    pd.DataFrame({"key": keys, "tissue": tissues}).to_csv(OUT / "G6_tissue_keys.csv", index=False)
    allres = []
    verdicts = {}
    for trait in G6_TRAITS:
        res, err = magma_gp(trait, covf, "Average", "G6")
        if not res:
            say("[G6] %s no result (%s)" % (trait, err[:120]))
            continue
        df = pd.DataFrame([{"key": k, "tissue": lab[k], **res[k]} for k in keys if k in res])
        df["q_bh"] = bh(df["p"].values)
        df = df.sort_values("p").reset_index(drop=True)
        df["rank"] = np.arange(1, len(df) + 1)
        rmask = df["tissue"].isin(REPRO)
        df["holm_repro"] = np.nan
        df.loc[rmask, "holm_repro"] = holm(df.loc[rmask, "p"].values)
        df["trait"] = trait
        allres.append(df)
        say("[G6_TOP] %s (genes=%d)" % (trait, int(df["ngenes"].max())))
        for _, r in df.head(5).iterrows():
            say("    %2d %-40s beta=%7.4f p=%.2e q=%.3f" % (r["rank"], r["tissue"][:40], r["beta"], r["p"], r["q_bh"]))
        rep = df[rmask].sort_values("p")
        say("    [REPRO] " + " | ".join("%s r%d p=%.3f H=%.2f" % (t, k, p, h) for t, k, p, h in
                                       zip(rep["tissue"], rep["rank"], rep["p"], rep["holm_repro"])))
        top3 = list(df.head(3)["tissue"])
        if trait == G6_PRIMARY:
            if (rep["holm_repro"] < 0.05).any():
                v = "REPRODUCTIVE: " + ", ".join(rep.loc[rep["holm_repro"] < 0.05, "tissue"])
            elif ((~rmask) & (df["q_bh"] < 0.05)).any():
                v = "NON-REPRODUCTIVE: " + ", ".join(df.loc[(~rmask) & (df["q_bh"] < 0.05), "tissue"])
            else:
                v = "NO TISSUE CONCENTRATION (min Holm_R=%.2f, min q=%.2f)" % (rep["holm_repro"].min(), df["q_bh"].min())
            verdicts["G6"] = v
            say("[G6_VERDICT] " + v)
        if trait == "HYPERTENSION":
            ok = any(t in ARTERY for t in top3)
            verdicts["PC_HYPERTENSION_artery_top3"] = ok
            say("[G6_POSCTRL] HYPERTENSION artery in top 3: %s (%s)" % ("PASS" if ok else "FAIL", "; ".join(top3)))
        if trait == "LEIOMYOMA":
            ok = "Uterus" in top3
            verdicts["PC_LEIOMYOMA_uterus_top3"] = ok
            say("[G6_POSCTRL] LEIOMYOMA uterus in top 3: %s (%s)" % ("PASS" if ok else "FAIL", "; ".join(top3)))
    if allres:
        pd.concat(allres, ignore_index=True).to_csv(OUT / "G6_tissue_results.csv", index=False)
    SUM["results"]["G6"] = verdicts


# ----------------------------------------------------------------------------
# 3. G1-S post hoc diagnostics
# ----------------------------------------------------------------------------
def contrast_row(trait, res, tag, lin, model):
    if not res or "D" not in res or "DN" not in res:
        return None
    return {"trait": trait, "contrast": tag, "lineage": lin, "model": model, "beta": res["D"]["beta"],
            "se": res["D"]["se"], "p_greater": res["D"]["p"], "p_less": res["DN"]["p"], "ngenes": res["D"]["ngenes"]}


def print_block(df, label):
    for _, r in df.iterrows():
        say("    [%s] %-13s %-11s %-9s beta=%8.4f se=%7.4f p_less=%.2e Holm=%.3f | p_greater=%.2e genes=%d" % (
            label, r["trait"][:13], r["lineage"], r["model"], r["beta"], r["se"], r["p_less"], r["holm_less"],
            r["p_greater"], r["ngenes"]))


def print_pairs(df, label, key, m_a, m_b):
    # One line per (key, lineage): model A -> model B
    for (k, lin), g in df.groupby([key, "lineage"]):
        a = g[g["model"] == m_a]
        b = g[g["model"] == m_b]
        if not len(a) or not len(b):
            continue
        a, b = a.iloc[0], b.iloc[0]
        say("    [%s] %-13s %-11s beta %8.4f -> %8.4f | p_less %.2e -> %.2e (Holm %.3f -> %.3f) genes=%d" % (
            label, str(k)[:13], lin, a["beta"], b["beta"], a["p_less"], b["p_less"], a["holm_less"], b["holm_less"], b["ngenes"]))


def add_holm(df):
    df = df.copy()
    df["holm_less"] = np.nan
    df["holm_greater"] = np.nan
    for _, idx in df.groupby(["trait", "contrast", "model"]).groups.items():
        df.loc[idx, "holm_less"] = holm(df.loc[idx, "p_less"].values)
        df.loc[idx, "holm_greater"] = holm(df.loc[idx, "p_greater"].values)
    return df


def gene_z(trait):
    f = MAGMA_OUT / (trait + "_35_10.genes.out")
    if not f.exists():
        return None
    g = pd.read_csv(f, sep=r"\s+")
    return pd.Series(g["ZSTAT"].values, index=g["GENE"].astype(str))


def g1s(sym2, tau, loeuf):
    say("")
    say("-" * 72)
    say("[SECTION] G1-S (post hoc, descriptive; CHANGELOG v0.2.10): is the G1 eutopic tilt endometriosis-specific?")
    rows = []
    covs = {}
    for tag in ["LES_vs_EU", "OMA_vs_EU"]:
        for lin in LINEAGES:
            f = COVAR_G1 / ("%s__%s.covar.txt" % (tag, lin))
            if f.exists():
                c = pd.read_csv(f, sep="\t")
                c["GENE"] = c["GENE"].astype(str)
                covs[(tag, lin)] = (f, c)
    say("[G1S_COVARIATES] %d contrast x lineage files from v018" % len(covs))
    if not covs:
        return

    # S1: T2D as second unrelated control (+ the v018 traits refitted for one table)
    say("[S1] AVG-only model (as v018), traits: " + ", ".join(G1S_TRAITS))
    for (tag, lin), (f, c) in covs.items():
        for trait in G1S_TRAITS:
            res, _ = magma_gp(trait, f, "AVG", "S1_%s_%s" % (tag, lin))
            r = contrast_row(trait, res, tag, lin, "AVG")
            if r:
                rows.append(r)
    s1 = add_holm(pd.DataFrame(rows))
    for tag in ["LES_vs_EU", "OMA_vs_EU"]:
        say("  [S1] contrast %s (T2D new; other traits refitted, compact)" % tag)
        print_block(s1[(s1["contrast"] == tag) & (s1["trait"] == "T2D")].sort_values("lineage"), "S1")
        for trait in [x for x in G1S_TRAITS if x != "T2D"]:
            g = s1[(s1["contrast"] == tag) & (s1["trait"] == trait)].set_index("lineage")
            say("    [S1_REFIT] %-13s " % trait[:13] + " ".join("%s:%+.3f(H%.2f)" % (l[:5], g.loc[l, "beta"], g.loc[l, "holm_less"])
                                                          for l in LINEAGES if l in g.index))

    # S2: conditioning on LOEUF + tau (same genes for both models)
    rows2 = []
    ent_tau = pd.Series(tau.values, index=[sym2[s] for s in tau.index])
    ent_loeuf = pd.Series(loeuf.values, index=[sym2[s] for s in loeuf.index]) if loeuf is not None else None
    if ent_loeuf is None:
        say("[S2] SKIP (no LOEUF)")
    else:
        say("[S2] conditioning AVG vs AVG+LOEUF+TAU on the same genes (lesion vs eutopic)")
        for lin in LINEAGES:
            if ("LES_vs_EU", lin) not in covs:
                continue
            f, c = covs[("LES_vs_EU", lin)]
            c2 = c.copy()
            c2["LOEUF"] = c2["GENE"].map(ent_loeuf)
            c2["TAU"] = c2["GENE"].map(ent_tau)
            c2 = c2.dropna()
            f2 = WORK / ("S2_LES_vs_EU__%s.covar.txt" % lin)
            c2.to_csv(f2, sep="\t", index=False, float_format="%.6g")
            rho_l = c2["D"].corr(c2["LOEUF"], method="spearman")
            rho_t = c2["D"].corr(c2["TAU"], method="spearman")
            say("    [S2_GENES] %-11s genes=%d (of %d)  spearman D~LOEUF=%.3f  D~TAU=%.3f" % (lin, len(c2), len(c), rho_l, rho_t))
            for trait in G1S_TRAITS:
                for model, cond in (("AVG_same", "AVG"), ("AVG+L+T", "AVG,LOEUF,TAU")):
                    res, _ = magma_gp(trait, f2, cond, "S2_%s_%s" % (lin, model.replace("+", "")))
                    r = contrast_row(trait, res, "LES_vs_EU", lin, model)
                    if r:
                        rows2.append(r)
        if rows2:
            s2 = add_holm(pd.DataFrame(rows2))
            say("  [S2] AVG_same -> AVG+LOEUF+TAU")
            print_pairs(s2, "S2", "trait", "AVG_same", "AVG+L+T")
            # Fixed interpretation (v0.2.10): lymphocyte, GCST and HYPERTENSION both shrink >= 50% and p_less >= 0.05
            ok = []
            for trait in ["GCST90205183", "HYPERTENSION"]:
                a = s2[(s2["trait"] == trait) & (s2["lineage"] == "lymphocyte") & (s2["model"] == "AVG_same")]
                b = s2[(s2["trait"] == trait) & (s2["lineage"] == "lymphocyte") & (s2["model"] == "AVG+L+T")]
                if len(a) and len(b):
                    shrink = 1.0 - b["beta"].iloc[0] / a["beta"].iloc[0] if a["beta"].iloc[0] != 0 else np.nan
                    say("    [S2_SHRINK] %-13s lymphocyte beta %.4f -> %.4f (shrink %.0f%%), p_less %.3g -> %.3g" % (
                        trait, a["beta"].iloc[0], b["beta"].iloc[0], 100 * shrink, a["p_less"].iloc[0], b["p_less"].iloc[0]))
                    ok.append(bool(shrink >= 0.5 and b["p_less"].iloc[0] >= 0.05))
                else:
                    ok.append(False)
            v = "SUPPORTED" if len(ok) == 2 and all(ok) else "NOT SUPPORTED"
            say("[S2_VERDICT] generic gene-property explanation (LOEUF + breadth): " + v)
            SUM["results"]["G1S_S2"] = v
            rows += rows2

    # S3: GCST conditional on a control trait's gene z
    rows3 = []
    say("[S3] GCST90205183 conditional on control-trait gene z (lesion vs eutopic)")
    for ctrl in ["HYPERTENSION", "T2D"]:
        z = gene_z(ctrl)
        if z is None:
            say("    [S3] %s gene z missing" % ctrl)
            continue
        for lin in LINEAGES:
            if ("LES_vs_EU", lin) not in covs:
                continue
            f, c = covs[("LES_vs_EU", lin)]
            c3 = c.copy()
            c3["ZC"] = c3["GENE"].map(z)
            c3 = c3.dropna()
            f3 = WORK / ("S3_LES_vs_EU__%s__%s.covar.txt" % (lin, ctrl))
            c3.to_csv(f3, sep="\t", index=False, float_format="%.6g")
            for model, cond in (("AVG|" + ctrl[:3], "AVG"), ("AVG+z" + ctrl[:3], "AVG,ZC")):
                res, _ = magma_gp("GCST90205183", f3, cond, "S3_%s_%s_%s" % (lin, ctrl, "z" if "+" in model else "avg"))
                r = contrast_row("GCST90205183", res, "LES_vs_EU", lin, model)
                if r:
                    r["control"] = ctrl
                    rows3.append(r)
    if rows3:
        s3 = add_holm(pd.DataFrame(rows3))
        for ctrl in ["HYPERTENSION", "T2D"]:
            sub = s3[s3["control"] == ctrl]
            if len(sub):
                say("  [S3] GCST90205183: AVG -> AVG + %s gene z" % ctrl)
                print_pairs(sub, "S3", "control", "AVG|" + ctrl[:3], "AVG+z" + ctrl[:3])
        rows += rows3

    # S4: lymphocyte label composition per context
    if PB_GROUPS.exists():
        G = pd.read_csv(PB_GROUPS)
        G = G[(G["agree"]) & (G["lineage_ref"] == "lymphocyte")]
        say("[S4] lymphocyte lineage composition (agreeing clusters), top labels per context")
        for ctx in ["eutopic", "peritoneal_lesion", "endometrioma"]:
            t = G[G["context"] == ctx].groupby("label")["n_cells"].sum().sort_values(ascending=False)
            tot = t.sum()
            say("    [S4] %-18s cells=%d  %s" % (ctx, tot, " | ".join("%s %.0f%%" % (k, 100 * v / tot) for k, v in t.head(6).items())))
    pd.DataFrame(rows).to_csv(OUT / "G1S_results.csv", index=False)


# ----------------------------------------------------------------------------
def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.19 | G6 tissue level (GTEx v8) + G1-S specificity diagnostics | " + TS)
    say("=" * 72)
    t0 = time.time()
    sym2 = entrez_unique()
    ok_t2d = run_gene_analysis("T2D", "T2D")
    SUM["gates"]["T2D_GENES_RAW"] = ok_t2d
    ok_gtex = download(GTEX_URLS, GTEX_FILE, "GTEx v8 gene median TPM")
    SUM["gates"]["GTEX"] = ok_gtex
    if not ok_gtex:
        say("[FATAL] GTEx file not available")
    else:
        L, tau, tissues = load_gtex(sym2)
        loeuf = load_loeuf(sym2)
        SUM["gates"]["LOEUF"] = loeuf is not None
        try:
            g6(L, tissues, sym2)
        except Exception as e:
            import traceback
            say("[STEP_ERROR] g6 %s: %s" % (type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1200:])
        try:
            g1s(sym2, tau, loeuf)
        except Exception as e:
            import traceback
            say("[STEP_ERROR] g1s %s: %s" % (type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1200:])
    with open(OUT / "p9_v019_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    for k, v in SUM["gates"].items():
        say("[GATE] %-24s %s" % (k, v))
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] G6 rules and G1-S interpretation fixed in CHANGELOG v0.2.10 before this run")
    say("[DONE] PAPER9_V019_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
