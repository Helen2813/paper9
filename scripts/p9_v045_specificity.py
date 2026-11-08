#!/usr/bin/env python3
# Paper 9 | v0.0.45 | (rules CHANGELOG v0.2.35)
#   C1  post hoc label check of ambiguous GTEx v9 cell types (author annotation, tissue, canonical markers)
#   C2  trait-specificity of each marker family: cross-trait rank correlation (HECA v043 and GTEx v044 rankings)
#   C3  |log2FC| sets in GTEx: depleted-gene share; do |log2FC| rankings track down-set or up-set rankings?
#   C4  conditioning on generic gene properties (mean expression, expression breadth, LoF constraint):
#       MAGMA gene-property analysis of binary set indicators, both atlases, six GWAS
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). ~20-40 min; gnomAD constraint (~5 MB) once.

import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
try:
    os.nice(10)
except Exception:
    pass

import datetime as dt
import itertools
import json
import re
import subprocess
import time
import warnings
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RAW = REPO / "data" / "raw"
INTERIM = REPO / "data" / "interim"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
GTEX_DIR = RAW / "cellxgene" / "gtex_eraslan2022"
GTEX_DE = PROC / "genesets" / "gtex_v044_de"
HECA_DE = PROC / "genesets" / "heca_koller_v041_de"
GNOMAD_URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
GNOMAD_FILE = RAW / "gnomad" / "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "specificity" / ("p9_v045_" + TS)
WORK = INTERIM / "magma_v045"
LOG = REPO / "logs" / ("p9_v045_specificity_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.35)
SEED = 20261108
N_BOOT = 1000
LFC, PADJ, MIN_SET, N_MATCH = 1.2, 0.01, 20, 400
FAMILIES = ["ABS", "UP", "DOWN", "ABSN", "UPN", "DOWNN", "SPEC"]
TRAITS = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL", "LEIOMYOMA", "HYPERTENSION", "T2D"]
PAIR_TRAITS = ["KOLLER_EUR", "LEIOMYOMA", "HYPERTENSION", "T2D"]
DOWN_T = ["DOWN", "DOWNN"]
UP_T = ["UP", "UPN", "SPEC"]
C1_TYPES = ["contractile cell", "smooth muscle cell", "enteric smooth muscle cell", "fibroblast",
            "basal-myoepithelial cell of mammary gland", "professional antigen presenting cell"]
PANELS = {"pericyte": ["RGS5", "KCNJ8", "ABCC9", "PDGFRB", "NOTCH3"], "vascular SMC": ["MYH11", "ACTA2", "CNN1", "TAGLN", "DES"],
          "myoepithelial": ["KRT5", "KRT14", "KRT17", "OXTR"], "endothelial": ["PECAM1", "VWF", "CDH5"],
          "fibroblast": ["DCN", "LUM", "PDGFRA"], "striated muscle": ["TNNT2", "ACTA1", "MYH7", "TTN"]}
MURAL_RE = r"pericyte|smc|smooth|vascular|mural"
EXPECT = {"GTEx": {"LEIOMYOMA": r"smooth muscle", "HYPERTENSION": r"smooth muscle|endothel|pericyte|vascular"},
          "HECA": {"LEIOMYOMA": r"^(uSMCs|mPV|ePV)", "HYPERTENSION": r"^(Arterial|Venous|ePV|mPV|uSMCs)"}}
KEY_HECA = [("ABS", "Venous/CTL"), ("ABS", "Arterial/CTL"), ("UP", "dStromal_mid/ENDO"), ("UPN", "dStromal_mid/ENDO"),
            ("DOWN", "Luminal/ENDO"), ("DOWN", "Ciliated/ENDO"), ("SPEC", "dStromal_mid/ENDO")]
STD_OBS = {"cell_type", "tissue", "assay", "disease", "sex", "organism", "development_stage", "self_reported_ethnicity",
           "donor_id", "suspension_type", "is_primary_data", "tissue_type", "observation_joinid", "ethnicity"}

for d in (OUT, WORK, LOG.parent):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V045_SPECIFICITY", "timestamp_utc": TS, "results": {}}
ERR = []


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


def magma_covar(raw, cf, tag, cond):
    outp = WORK / tag
    gsa = Path(str(outp) + ".gsa.out")
    model = (["condition-hide=" + ",".join(cond)] if cond else []) + ["direction=greater"]
    if gsa.exists():
        gsa.unlink()
    rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--gene-covar", cf, "--model"] + model + ["--out", outp], timeout=7200)
    if gsa.exists():
        return parse_gsa(gsa)
    say("    [MAGMA_FAIL] %s %s" % (tag, (o + e).strip()[-300:]))
    return None


def dec(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def obs_col(f, name):
    g = f["obs"][name]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = np.array([dec(c) for c in g["categories"][:]], dtype=object)
        codes = g["codes"][:]
        out = np.empty(len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        return out
    v = g[:]
    return np.array([dec(x) for x in v], dtype=object) if v.dtype.kind in "SO" else v


def top3(v):
    s = pd.Series(v).value_counts(normalize=True)
    return ", ".join("%s %.0f%%" % (k, 100 * x) for k, x in s.head(3).items())


def latest(base, pattern):
    fs = sorted(base.glob(pattern))
    return fs[-1] if fs else None


def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return np.nan
    return float(pd.Series(x[m]).rank().corr(pd.Series(y[m]).rank()))


def short(label, n=26):
    s = str(label).replace("SOX9_functionalis_", "SOX9f_").replace("dStromal_", "dS_")
    return s[:n]


# ----------------------------------------------------------------------------
# C1. Ambiguous GTEx labels (post hoc)
# ----------------------------------------------------------------------------
def c1():
    say("")
    say("-" * 72)
    say("[SECTION] C1  post hoc label check of ambiguous GTEx v9 cell types")
    fs = [p for p in sorted(GTEX_DIR.glob("gtex_*.h5ad")) if p.stat().st_size > 0]
    mural = None
    if fs:
        f = h5py.File(fs[-1], "r")
        cols = [k for k in f["obs"].keys() if k not in ("_index", "__categories")]
        say("[C1] obs columns (%d): %s" % (len(cols), ", ".join(cols)[:800]))
        ct = obs_col(f, "cell_type")
        tissue = obs_col(f, "tissue") if "tissue" in f["obs"] else np.full(len(ct), "-", dtype=object)
        cand = []
        for c in cols:
            if c in STD_OBS or c.endswith("_ontology_term_id"):
                continue
            try:
                v = obs_col(f, c)
            except Exception:
                continue
            if getattr(v, "dtype", None) is not None and v.dtype == object:
                nun = len(set(v.tolist()))
                if 2 <= nun <= 400:
                    cand.append((c, v))
        f.close()
        for t in C1_TYPES:
            m = ct == t
            if not m.any():
                continue
            say("  [C1] %s: %d nuclei; tissue %s" % (t, int(m.sum()), top3(tissue[m])))
            for c, v in cand:
                say("       %s: %s" % (c, top3(v[m])))
            if t == "contractile cell":
                shares = []
                for c, v in cand:
                    if re.search(r"cell.?type|annotation|granular|broad", c, re.I):
                        s = pd.Series(v[m]).value_counts(normalize=True)
                        share = float(s[[k for k in s.index if re.search(MURAL_RE, str(k), re.I)]].sum())
                        shares.append(share)
                        say("       -> %s: mural-matching annotation share %.2f" % (c, share))
                if shares:
                    mural = all(x >= 0.5 for x in shares)
    else:
        say("[C1] GTEx h5ad not found")
    cf = GTEX_DE / "mean_cpm_by_type.csv.gz"
    df = GTEX_DE / "de_all.csv.gz"
    panel_vals = {}
    if cf.exists() and df.exists():
        CPM = pd.read_csv(cf, index_col=0)
        DE = pd.read_csv(df)
        for t in C1_TYPES:
            if t not in CPM.columns:
                continue
            pv = {k: float(CPM.loc[[g for g in gl if g in CPM.index], t].mean()) if any(g in CPM.index for g in gl) else np.nan for k, gl in PANELS.items()}
            panel_vals[t] = pv
            topg = DE[DE["group"] == t].sort_values("score", ascending=False)["gene"].head(12).astype(str)
            say("  [C1] %s markers (mean CPM): %s | top up genes: %s" % (t, ", ".join("%s %.0f" % (k, v) for k, v in pv.items()), ", ".join(topg)))
    if mural is None and "contractile cell" in panel_vals:
        pv = panel_vals["contractile cell"]
        mur = np.nanmax([pv["pericyte"], pv["vascular SMC"]])
        mural = bool(mur > 2 * np.nan_to_num(pv["myoepithelial"]) and mur > 2 * np.nan_to_num(pv["striated muscle"]))
        say("  [C1] no author annotation column; marker rule -> %s" % ("mural" if mural else "not mural"))
    say("[C1_LABEL] contractile cell -> %s (post hoc)" % ({True: "vascular mural (pericyte / smooth muscle)", False: "not vascular mural", None: "undetermined"}[mural]))
    SUM["results"]["c1_contractile_mural"] = mural
    return bool(mural)


# ----------------------------------------------------------------------------
# C2 / C3. Cross-trait sharing and |log2FC| composition (from v043 / v044 rankings)
# ----------------------------------------------------------------------------
def rankings():
    R = {}
    for atlas, base, pat in (("HECA", REPO / "results" / "benchmark", "p9_v043_*/e1_real_family_rankings.csv"),
                             ("GTEx", REPO / "results" / "external_atlas", "p9_v044_*/b_gtex_family_rankings.csv")):
        f = latest(base, pat)
        if f is None:
            say("[C2] %s rankings missing" % atlas)
            continue
        R[atlas] = pd.read_csv(f)
        say("[C2] %s rankings from %s" % (atlas, f.parent.name))
    return R


def lp_matrix(D, fam, traits, col="p"):
    d = D[(D["family"] == fam) & (D["trait"].isin(traits))]
    M = d.pivot_table(index="label", columns="trait", values=col, aggfunc="min")
    return -np.log10(M.clip(lower=1e-300))


def mean_cross(M, traits):
    vals = []
    for a, b in itertools.combinations(traits, 2):
        if a in M.columns and b in M.columns:
            vals.append(spearman(M[a].values, M[b].values))
    vals = [v for v in vals if np.isfinite(v)]
    return (float(np.mean(vals)) if vals else np.nan), vals


def cross_trait(D, atlas, tag, col="p"):
    mats = {fam: lp_matrix(D, fam, PAIR_TRAITS, col) for fam in FAMILIES}
    common = sorted(set.intersection(*[set(m.dropna().index) for m in mats.values() if len(m)]))
    rho = {}
    for fam in FAMILIES:
        rho[fam] = mean_cross(mats[fam].loc[common], PAIR_TRAITS)[0]
    rng = np.random.default_rng(SEED)
    idx = np.array(common, dtype=object)
    diffs = []
    for _ in range(N_BOOT):
        lab = idx[rng.integers(0, len(idx), len(idx))]
        vv = {fam: mean_cross(mats[fam].loc[lab].reset_index(drop=True), PAIR_TRAITS)[0] for fam in DOWN_T + UP_T}
        diffs.append(np.nanmean([vv[f] for f in DOWN_T]) - np.nanmean([vv[f] for f in UP_T]))
    d0 = np.nanmean([rho[f] for f in DOWN_T]) - np.nanmean([rho[f] for f in UP_T])
    lo, hi = np.nanpercentile(diffs, [2.5, 97.5])
    say("  [%s] %s (%d labels, %d trait pairs): mean cross-trait Spearman rho: %s" % (tag, atlas, len(common), 6, "; ".join("%s %.2f" % (f, rho[f]) for f in FAMILIES)))
    say("  [%s] %s down-type minus up-type: %.2f (bootstrap 95%% CI %.2f to %.2f)" % (tag, atlas, d0, lo, hi))
    return {"rho": rho, "diff": float(d0), "ci": [float(lo), float(hi)], "n_labels": len(common)}


def c2(R):
    say("")
    say("-" * 72)
    say("[SECTION] C2  trait-specificity of marker families (traits %s)" % ", ".join(PAIR_TRAITS))
    res = {}
    for atlas, D in R.items():
        res[atlas] = cross_trait(D, atlas, "C2")
        for fam in ("ABS", "UP", "UPN", "SPEC", "DOWN"):
            d = D[D["family"] == fam]
            k = d[d["trait"] == "KOLLER_EUR"].sort_values("p").head(3)
            parts = []
            for _, r in k.iterrows():
                rk = []
                for t in PAIR_TRAITS[1:]:
                    e = d[d["trait"] == t].sort_values("p").reset_index(drop=True)
                    hit = e.index[e["label"] == r["label"]]
                    rk.append("%s %s" % (t[:4], int(hit[0]) + 1 if len(hit) else "-"))
                parts.append("%s (%s)" % (short(r["label"]), ", ".join(rk)))
            say("    [C2_TOP] %s %-4s KOLLER top 3 and their ranks in comparators: %s" % (atlas, fam, "; ".join(parts)))
    ok = len(res) == 2 and all(v["diff"] > 0 and v["ci"][0] > 0 for v in res.values())
    say("[H_NONSPEC_VERDICT] down-type families more shared across unrelated traits than up-type families in both atlases -> %s" % (
        "SUPPORTED" if ok else "NOT SUPPORTED"))
    SUM["results"]["c2"] = res
    SUM["results"]["h_nonspec"] = ok


def c3(R):
    say("")
    say("-" * 72)
    say("[SECTION] C3  |log2FC| sets: composition and which direction their rankings follow")
    med = {}
    for atlas, fn in (("HECA", HECA_DE / "de_all.csv.gz"), ("GTEx", GTEX_DE / "de_all.csv.gz")):
        if not fn.exists():
            continue
        DE = pd.read_csv(fn)
        fr = []
        for g, d in DE.groupby("group"):
            s = d[(d["lfc"].abs() > LFC) & (d["padj"] < PADJ)]
            if len(s) >= MIN_SET:
                fr.append((g, float((s["lfc"] < 0).mean()), len(s)))
        F = pd.DataFrame(fr, columns=["group", "frac_down", "n"])
        F.to_csv(OUT / ("c3_abs_composition_%s.csv" % atlas), index=False)
        med[atlas] = float(F["frac_down"].median())
        say("  [C3] %s |log2FC| sets: depleted-gene share median %.2f (IQR %.2f-%.2f); %d of %d sets > 0.5; median size %d" % (
            atlas, med[atlas], F["frac_down"].quantile(0.25), F["frac_down"].quantile(0.75), int((F["frac_down"] > 0.5).sum()), len(F), int(F["n"].median())))
    wins = {}
    for atlas, D in R.items():
        w, parts = 0, []
        for t in TRAITS:
            x = D[D["trait"] == t].pivot_table(index="label", columns="family", values="p", aggfunc="min")
            if not {"ABS", "DOWN", "UP"} <= set(x.columns):
                continue
            x = -np.log10(x.clip(lower=1e-300))
            rd, ru = spearman(x["ABS"].values, x["DOWN"].values), spearman(x["ABS"].values, x["UP"].values)
            w += int(rd > ru)
            parts.append("%s %.2f/%.2f" % (t, rd, ru))
        wins[atlas] = w
        say("  [C3] %s rho(ABS,DOWN) / rho(ABS,UP) across labels: %s -> ABS follows DOWN in %d/%d traits" % (atlas, "; ".join(parts), w, len(parts)))
    ok = med.get("GTEx", 0) > 0.5 and len(wins) == 2 and all(v >= 4 for v in wins.values())
    say("[H_ABSDOWN_VERDICT] GTEx |log2FC| sets mostly depleted genes and |log2FC| rankings follow down-set rankings in both atlases -> %s" % (
        "SUPPORTED" if ok else "NOT SUPPORTED"))
    SUM["results"]["c3"] = {"median_frac_down": med, "abs_follows_down": wins, "supported": ok}


# ----------------------------------------------------------------------------
# C4. Conditioning on generic gene properties
# ----------------------------------------------------------------------------
def gnomad():
    if not GNOMAD_FILE.exists():
        GNOMAD_FILE.parent.mkdir(parents=True, exist_ok=True)
        part = GNOMAD_FILE.with_name(GNOMAD_FILE.name + ".part")
        rc, o, e = sh(["curl", "-L", "-sS", "--fail", "--retry", "3", "-A", "paper9", "-o", part, GNOMAD_URL], timeout=1800)
        if rc == 0 and part.exists() and part.stat().st_size > 1_000_000:
            part.rename(GNOMAD_FILE)
        else:
            say("[C4] gnomAD constraint download failed (%s); conditioning on expression level and breadth only" % (e or o).strip()[:150])
            return None
    try:
        G = pd.read_csv(GNOMAD_FILE, sep="\t", compression="gzip")
        sc = next(c for c in ("gene", "gene_symbol") if c in G.columns)
        lc = next(c for c in ("oe_lof_upper", "LOEUF") if c in G.columns)
        G = pd.to_numeric(G.set_index(sc)[lc], errors="coerce").dropna()
        G = G.groupby(level=0).min()
        say("[C4] gnomAD v2.1.1 LOEUF for %d genes" % len(G))
        return G
    except Exception as ex:
        say("[C4] gnomAD constraint unreadable (%s); expression covariates only" % str(ex)[:150])
        return None


def build(DE, CPM, sym2ent, loeuf):
    meta, sets = {}, {}
    for g, d in DE.groupby("group", sort=False):
        sig = d["padj"] < PADJ
        sel = {"ABS": d.loc[(d["lfc"].abs() > LFC) & sig, "gene"], "UP": d.loc[(d["lfc"] > LFC) & sig, "gene"],
               "DOWN": d.loc[(d["lfc"] < -LFC) & sig, "gene"],
               "ABSN": d.reindex(d["score"].abs().sort_values(ascending=False).index)["gene"].head(N_MATCH),
               "UPN": d.sort_values("score", ascending=False)["gene"].head(N_MATCH),
               "DOWNN": d.sort_values("score", ascending=True)["gene"].head(N_MATCH)}
        for fam, genes in sel.items():
            gl = {x for x in genes.astype(str) if x in sym2ent}
            if len(gl) >= MIN_SET:
                k = "V%04d" % (len(meta) + 1)
                meta[k], sets[k] = (fam, g), gl
    C = CPM.copy()
    C.index = C.index.astype(str)
    C = C.groupby(level=0).mean()
    C = C[(C >= 1.0).any(axis=1)]
    C = C.loc[[g for g in C.index if g in sym2ent]]
    covs = ["Average", "Breadth"]
    cols = {}
    S = C.div(C.sum(axis=1), axis=0)
    spec = {}
    for i, c in enumerate(S.columns):
        k = "S%03d" % (i + 1)
        spec[k] = c
        cols[k] = S[c].values
    cols["Average"] = np.log2(C.mean(axis=1).values + 1.0)
    cols["Breadth"] = (C >= 1.0).mean(axis=1).values
    if loeuf is not None:
        cols["LOEUF"] = loeuf.reindex(C.index).values
        covs.append("LOEUF")
    for k, gl in sets.items():
        cols[k] = np.isin(C.index.values, list(gl)).astype(int)
    T = pd.DataFrame(cols, index=C.index)
    T = T.dropna()
    T.index = [sym2ent[g] for g in T.index]
    T = T[~T.index.duplicated()]
    keep = [k for k in sets if T[k].sum() >= MIN_SET]
    meta = {k: meta[k] for k in keep}
    T = T[list(spec) + covs + keep]
    return T, meta, spec, covs


def ranked(res, meta, spec, fam):
    if not res:
        return []
    if fam == "SPEC":
        items = [(spec[k], p) for k, p in res.items() if k in spec]
    else:
        items = [(meta[k][1], p) for k, p in res.items() if k in meta and meta[k][0] == fam]
    items.sort(key=lambda x: x[1])
    if not items:
        return []
    q = bh(np.array([p for _, p in items]))
    return [(l, p, qq) for (l, p), qq in zip(items, q)]


def c4(R, sym2ent, mural):
    say("")
    say("-" * 72)
    say("[SECTION] C4  conditioning on generic gene properties (MAGMA gene-property analysis, binary set indicators)")
    loeuf = gnomad()
    atl = {}
    for atlas, ddir, cpmf in (("HECA", HECA_DE, "mean_cpm_by_group.csv.gz"), ("GTEx", GTEX_DE, "mean_cpm_by_type.csv.gz")):
        if (ddir / "de_all.csv.gz").exists() and (ddir / cpmf).exists():
            atl[atlas] = (pd.read_csv(ddir / "de_all.csv.gz"), pd.read_csv(ddir / cpmf, index_col=0))
    allrows = []
    verdict = {}
    for atlas, (DE, CPM) in atl.items():
        T, meta, spec, covs = build(DE, CPM, sym2ent, loeuf)
        cf = WORK / ("covar_%s.txt" % atlas)
        T.to_csv(cf, sep="\t", float_format="%.6g", index_label="GENE")
        say("[C4] %s: %d genes; %d set indicators; %d specificity columns; conditioning on %s" % (atlas, len(T), len(meta), len(spec), ", ".join(covs)))
        res = {}
        t0 = time.time()
        for t in TRAITS:
            raw = MAGMA_OUT / (t + "_35_10.genes.raw")
            if not raw.exists():
                say("[C4] %s genes.raw missing" % t)
                continue
            r0 = magma_covar(raw, cf, "%s_%s_r0" % (atlas, t), [])
            r1 = magma_covar(raw, cf, "%s_%s_r1" % (atlas, t), ["Average"])
            r2 = magma_covar(raw, cf, "%s_%s_r2" % (atlas, t), covs)
            if not (r0 and r1 and r2):
                continue
            base = {k: (r1 if k in spec else r0).get(k, np.nan) for k in list(spec) + list(meta)}
            cond = {k: r2.get(k, np.nan) for k in list(spec) + list(meta)}
            res[t] = {}
            for fam in FAMILIES:
                rb, rc = ranked(base, meta, spec, fam), ranked(cond, meta, spec, fam)
                res[t][fam] = (rb, rc)
                pc = {l: (p, q) for l, p, q in rc}
                pa = {l: p for l, p in ((meta[k][1], r1.get(k, np.nan)) for k in meta if meta[k][0] == fam)} if fam != "SPEC" else {}
                for i, (l, p, q) in enumerate(rb):
                    allrows.append({"atlas": atlas, "trait": t, "family": fam, "label": l, "p_base": p, "q_base": q, "rank_base": i + 1,
                                    "p_avg_only": pa.get(l, np.nan), "p_cond": pc.get(l, (np.nan, np.nan))[0], "q_cond": pc.get(l, (np.nan, np.nan))[1],
                                    "rank_cond": [x[0] for x in rc].index(l) + 1 if l in pc else np.nan})
        say("[C4] %s MAGMA runs done (%.0f min)" % (atlas, (time.time() - t0) / 60))
        if not res:
            continue
        A = pd.DataFrame([r for r in allrows if r["atlas"] == atlas])
        if atlas in R:
            P = R[atlas].merge(A, on=["trait", "family", "label"], how="inner")
            P = P[P["family"] != "SPEC"]
            say("  [C4_CHECK] %s unconditioned indicator p vs v043/v044 set-analysis p: Spearman %.3f over %d tests" % (
                atlas, spearman(-np.log10(P["p"].values), -np.log10(P["p_base"].values)), len(P)))
        nb, nc = {}, {}
        for fam in FAMILIES:
            d = A[A["family"] == fam]
            nb[fam] = int((d["q_base"] < 0.05).sum())
            nc[fam] = int((d["q_cond"] < 0.05).sum())
        say("  [C4] %s significant labels (BH q < 0.05, summed over %d traits) base -> conditioned: %s" % (
            atlas, len(res), "; ".join("%s %d->%d" % (f, nb[f], nc[f]) for f in FAMILIES)))
        cb = cross_trait(A.rename(columns={"p_base": "p"})[["trait", "family", "label", "p"]], atlas, "C4 base", "p")
        cc = cross_trait(A.rename(columns={"p_cond": "p"})[["trait", "family", "label", "p"]], atlas, "C4 cond", "p")
        ret_down = sum(nc[f] for f in DOWN_T) / max(1, sum(nb[f] for f in DOWN_T))
        ret_up = sum(nc[f] for f in UP_T) / max(1, sum(nb[f] for f in UP_T))
        say("  [C4] %s retained significant labels: down-type %.2f, up-type %.2f" % (atlas, ret_down, ret_up))
        verdict[atlas] = {"retained_down": ret_down, "retained_up": ret_up, "n_base": nb, "n_cond": nc, "rho_base": cb, "rho_cond": cc}
        for t, rx in EXPECT[atlas].items():
            if t not in res:
                continue
            parts = []
            for fam in FAMILIES:
                rb, rc = res[t][fam]
                for tag, r in (("base", rb), ("cond", rc)):
                    hits = [i + 1 for i, (l, p, q) in enumerate(r) if re.search(rx, str(l), re.I) or (atlas == "GTEx" and mural and l == "contractile cell")]
                    if tag == "cond":
                        parts.append("%s %s->%s" % (fam, hb, hits[0] if hits else "-"))
                    hb = hits[0] if hits else "-"
            say("  [C4_CTRL] %s %-12s rank of expected group base->conditioned: %s" % (atlas, t, "; ".join(parts)))
        for t in ("KOLLER_EUR",):
            if t not in res:
                continue
            for fam in FAMILIES:
                rc = res[t][fam][1]
                say("    [C4_TOP] %s %s %-5s conditioned top 3: %s" % (atlas, t, fam, "; ".join("%s %.1e/%.2f" % (short(l), p, q) for l, p, q in rc[:3])))
        if atlas == "HECA" and "KOLLER_EUR" in res:
            kv = {}
            for fam, lab in KEY_HECA:
                d = A[(A["trait"] == "KOLLER_EUR") & (A["family"] == fam) & (A["label"] == lab)]
                if len(d):
                    r = d.iloc[0]
                    kv[fam + " " + lab] = float(r["p_cond"])
                    say("    [C4_KEY] KOLLER_EUR %-4s %-18s p %.1e (rank %d) -> expression level only %s -> all covariates %.1e (rank %s)" % (
                        fam, lab, r["p_base"], r["rank_base"], ("%.1e" % r["p_avg_only"]) if np.isfinite(r["p_avg_only"]) else "-",
                        r["p_cond"], int(r["rank_cond"]) if np.isfinite(r["rank_cond"]) else "-"))
            ok_d = kv.get("UP dStromal_mid/ENDO", 1) < 0.05 and kv.get("UPN dStromal_mid/ENDO", 1) < 0.05
            say("[H_DSTROMA_VERDICT] up-only and size-matched up-only decidual stroma keep p < 0.05 after conditioning (KOLLER_EUR) -> %s" % (
                "SUPPORTED" if ok_d else "NOT SUPPORTED"))
            SUM["results"]["h_dstroma"] = {"supported": ok_d, **kv}
    pd.DataFrame(allrows).to_csv(OUT / "c4_conditioned.csv", index=False)
    ok = len(verdict) == 2 and all(v["retained_down"] < v["retained_up"] for v in verdict.values())
    say("[H_GENERIC_VERDICT] conditioning removes a larger share of down-type than up-type significance in both atlases -> %s" % (
        "SUPPORTED" if ok else "NOT SUPPORTED"))
    SUM["results"]["c4"] = verdict
    SUM["results"]["h_generic"] = ok


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.45 | label check, trait-specificity of marker families, generic gene-property conditioning | " + TS)
    say("=" * 72)
    t0 = time.time()
    sym2ent = gene_loc()
    state = {"mural": False, "R": {}}

    def s1():
        state["mural"] = c1()

    def s2():
        state["R"] = rankings()
        c2(state["R"])

    def s3():
        c3(state["R"])

    def s4():
        c4(state["R"], sym2ent, state["mural"])

    for name, fn in (("C1", s1), ("C2", s2), ("C3", s3), ("C4", s4)):
        try:
            fn()
        except Exception as e:
            import traceback
            ERR.append(name)
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v045_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.35 before this run (C1 post hoc; C2-C4 formulated after v044)")
    if ERR:
        say("[INCOMPLETE] steps with errors: " + ", ".join(ERR))
    else:
        say("[DONE] PAPER9_V045_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
