#!/usr/bin/env python3
# Paper 9 | v0.0.49 | review-proofing on cached data (rules CHANGELOG v0.2.39)
#   R1 threshold / size sensitivity of the rule flip and of the leiomyoma control (HECA, MAGMA)
#   R3 trait contrast: is the decidual-stroma enrichment stronger for endometriosis than for T2D / hypertension?
#   O2 where the genes of the 57 endometriosis loci are expressed across HECA cell states (heat map)
#   I0 read-only inventory of the local GSE266265 pseudobulk and GSE213216 cell tables (labels, contexts, donors)
#      so that the v050 rules (independent atlas; candidate cells of origin) can be fixed before any test
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). ~15-30 min; no downloads.

import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
try:
    os.nice(10)
except Exception:
    pass

import datetime as dt
import json
import math
import re
import subprocess
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
INTERIM = REPO / "data" / "interim"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
HECA_DE = PROC / "genesets" / "heca_koller_v041_de"
RES = REPO / "results"
PB266 = PROC / "gse266265_pseudobulk"
CELLS213 = PROC / "gse213216_cells_v028"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = RES / "review" / ("p9_v049_" + TS)
WORK = INTERIM / "magma_v049"
LOG = REPO / "logs" / ("p9_v049_review_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.39)
MIN_SET = 20
LFCS = [0.5, 1.0, 1.2, 2.0]
PADJS = [0.05, 0.01]
SIZES = [200, 400, 800]
BASE = (1.2, 0.01, 400)
TRAITS = ["KOLLER_EUR", "GCST90205183", "LEIOMYOMA", "HYPERTENSION", "T2D"]
DECIDUAL_RE = r"^dStromal"
LEIO_RE = r"^(uSMCs|mPV|ePV)"
CONTRAST = [("UPN", "dStromal_mid/ENDO"), ("UP", "dStromal_mid/ENDO"), ("UPN", "dStromal_early/ENDO"), ("UP", "dStromal_early/ENDO"),
            ("SPEC", "dStromal_mid/ENDO"), ("ABS", "Venous/CTL")]
PRIMARY = ("UPN", "dStromal_mid/ENDO")
COMPARE = [("KOLLER_EUR", "T2D"), ("KOLLER_EUR", "HYPERTENSION"), ("KOLLER_EUR", "LEIOMYOMA"), ("GCST90205183", "T2D"),
           ("GCST90205183", "HYPERTENSION")]
COMP_ORDER = ["stroma", "perivascular/SM", "endothelium", "epithelium", "immune", "broad"]

for d in (OUT, WORK, LOG.parent):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V049_REVIEW", "timestamp_utc": TS, "results": {}}
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


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def short(label, n=22):
    return str(label).replace("SOX9_functionalis_", "SOX9f_").replace("dStromal_", "dS_")[:n]


def latest(base, pattern):
    fs = sorted(base.glob(pattern))
    return fs[-1] if fs else None


# ----------------------------------------------------------------------------
# Sets for the grid
# ----------------------------------------------------------------------------
def build(DE, CPM, sym2ent):
    uni = sorted({g for g in DE["gene"].astype(str) if g in sym2ent})
    uset = set(uni)
    meta, cols = {}, {}

    def add(fam, setting, g, genes):
        gl = [x for x in genes if x in uset]
        if len(gl) >= MIN_SET:
            k = "V%05d" % (len(meta) + 1)
            meta[k] = (fam, setting, g)
            cols[k] = gl

    for g, d in DE.groupby("group", sort=False):
        d = d.assign(gene=d["gene"].astype(str))
        for lfc in LFCS:
            for pa in PADJS:
                sig = d["padj"] < pa
                st = "lfc%.1f_p%.2f" % (lfc, pa)
                add("ABS", st, g, d.loc[(d["lfc"].abs() > lfc) & sig, "gene"])
                add("UP", st, g, d.loc[(d["lfc"] > lfc) & sig, "gene"])
                add("DOWN", st, g, d.loc[(d["lfc"] < -lfc) & sig, "gene"])
        for n in SIZES:
            st = "N%d" % n
            add("ABSN", st, g, d.reindex(d["score"].abs().sort_values(ascending=False).index)["gene"].head(n))
            add("UPN", st, g, d.sort_values("score", ascending=False)["gene"].head(n))
            add("DOWNN", st, g, d.sort_values("score", ascending=True)["gene"].head(n))
    data = {k: np.isin(uni, gl).astype(int) for k, gl in cols.items()}
    T = pd.DataFrame(data, index=[sym2ent[g] for g in uni])
    T = T[~T.index.duplicated()]
    # Specificity covariates on expressed genes only (as v043): shares + Average
    C = CPM.copy()
    C.index = C.index.astype(str)
    C = C.groupby(level=0).mean()
    C = C[(C >= 1.0).any(axis=1)]
    C = C.loc[[g for g in C.index if g in sym2ent]]
    S = C.div(C.sum(axis=1), axis=0)
    spec, sd = {}, {}
    for i, c in enumerate(S.columns):
        k = "S%03d" % (i + 1)
        spec[k] = c
        sd[k] = S[c].values
    sd["Average"] = np.log2(C.mean(axis=1).values + 1.0)
    Ts = pd.DataFrame(sd, index=[sym2ent[g] for g in C.index])
    Ts = Ts[~Ts.index.duplicated()]
    return T, Ts, meta, spec


def first_labels(res, meta, fam, setting):
    items = sorted([(meta[k][2], v["p"]) for k, v in res.items() if k in meta and meta[k][0] == fam and meta[k][1] == setting], key=lambda x: x[1])
    return items


# ----------------------------------------------------------------------------
# R1 + R3
# ----------------------------------------------------------------------------
def r1_r3(sym2ent):
    DE = pd.read_csv(HECA_DE / "de_all.csv.gz")
    CPM = pd.read_csv(HECA_DE / "mean_cpm_by_group.csv.gz", index_col=0)
    T, Ts, meta, spec = build(DE, CPM, sym2ent)
    cf = WORK / "covar_grid.txt"
    T.to_csv(cf, sep="\t", index_label="GENE")
    cfs = WORK / "covar_spec.txt"
    Ts.to_csv(cfs, sep="\t", float_format="%.6g", index_label="GENE")
    say("[R1] HECA (v041 DE, variant P cells): %d genes; %d set indicators over %d threshold settings and %d sizes; %d specificity columns" % (
        len(T), len(meta), len(LFCS) * len(PADJS), len(SIZES), len(spec)))
    R0, R1 = {}, {}
    t0 = time.time()
    for t in TRAITS:
        raw = MAGMA_OUT / (t + "_35_10.genes.raw")
        if not raw.exists():
            say("[R1] %s genes.raw missing" % t)
            continue
        r0 = magma_covar(raw, cf, "grid_%s_r0" % t, [])
        r1 = magma_covar(raw, cfs, "spec_%s_r1" % t, ["Average"])
        if r0 and r1:
            R0[t], R1[t] = r0, r1
    say("[R1] MAGMA runs done (%.0f min)" % ((time.time() - t0) / 60))
    rows = []
    for t in R0:
        for k, (fam, st, g) in meta.items():
            if k in R0[t]:
                rows.append({"trait": t, "family": fam, "setting": st, "label": g, **R0[t][k]})
        for k, g in spec.items():
            if k in R1[t]:
                rows.append({"trait": t, "family": "SPEC", "setting": "continuous", "label": g, **R1[t][k]})
    A = pd.DataFrame(rows)
    A.to_csv(OUT / "r1_grid_results.csv", index=False)
    # R1 report
    say("")
    say("-" * 72)
    say("[SECTION] R1 first-ranked HECA label across thresholds (|log2FC| > x, adj. p < y) and set sizes")
    settings = ["lfc%.1f_p%.2f" % (l, p) for l in LFCS for p in PADJS]
    verdict = {}
    for t in ("KOLLER_EUR", "GCST90205183", "LEIOMYOMA", "HYPERTENSION"):
        if t not in R0:
            continue
        say("  %s" % t)
        cnt = {"ABS": 0, "UP": 0, "DOWN": 0}
        for st in settings:
            parts = []
            for fam in ("ABS", "UP", "DOWN"):
                it = first_labels(R0[t], meta, fam, st)
                if not it:
                    parts.append("%s -" % fam)
                    continue
                lab, p = it[0]
                n_sets = len(it)
                rx = DECIDUAL_RE if t in ("KOLLER_EUR", "GCST90205183") else (LEIO_RE if t == "LEIOMYOMA" else r"^(Arterial|Venous|ePV|mPV|uSMCs)")
                cnt[fam] += int(bool(re.match(rx, lab)))
                parts.append("%s %s %.0e (%d sets)" % (fam, short(lab, 18), p, n_sets))
            say("    %-14s %s" % (st, " | ".join(parts)))
        for n in SIZES:
            st = "N%d" % n
            parts = []
            for fam in ("ABSN", "UPN", "DOWNN"):
                it = first_labels(R0[t], meta, fam, st)
                parts.append("%s %s %.0e" % (fam, short(it[0][0], 18), it[0][1]) if it else "%s -" % fam)
            say("    %-14s %s" % (st, " | ".join(parts)))
        verdict[t] = cnt
        say("    [R1] %s settings (of %d) with the expected / decidual group first: ABS %d, UP %d, DOWN %d" % (t, len(settings), cnt["ABS"], cnt["UP"], cnt["DOWN"]))
    k, l = verdict.get("KOLLER_EUR", {}), verdict.get("LEIOMYOMA", {})
    flip = k.get("UP", 0) >= 6 and (len(settings) - k.get("ABS", 0)) >= 6
    ctrl = l.get("UP", 0) >= 6 and l.get("ABS", 0) < l.get("UP", 0)
    say("[H_ROBUST_FLIP_VERDICT] KOLLER_EUR: up-only first = decidual stroma in %d/8 settings; |log2FC| first not decidual in %d/8 -> %s" % (
        k.get("UP", 0), len(settings) - k.get("ABS", 0), "SUPPORTED" if flip else "NOT SUPPORTED"))
    say("[H_ROBUST_CTRL_VERDICT] LEIOMYOMA: smooth-muscle group first with up-only in %d/8, with |log2FC| in %d/8 -> %s" % (
        l.get("UP", 0), l.get("ABS", 0), "SUPPORTED" if ctrl else "NOT SUPPORTED"))
    SUM["results"]["r1"] = {"counts": verdict, "flip": flip, "ctrl": ctrl}
    # R3 report
    say("")
    say("-" * 72)
    say("[SECTION] R3 trait contrast of standardized enrichment (MAGMA BETA_STD; base thresholds |log2FC| > 1.2, adj. p < 0.01, N = 400)")
    base_st = {"ABS": "lfc1.2_p0.01", "UP": "lfc1.2_p0.01", "DOWN": "lfc1.2_p0.01", "UPN": "N400", "ABSN": "N400", "DOWNN": "N400", "SPEC": "continuous"}
    est = {}
    for fam, lab in CONTRAST:
        d = A[(A["family"] == fam) & (A["setting"] == base_st[fam]) & (A["label"] == lab)]
        for _, r in d.iterrows():
            b, bs, se = r["beta"], r["beta_std"], r["se"]
            ses = abs(se * bs / b) if (np.isfinite(b) and b != 0 and np.isfinite(bs)) else np.nan
            est[(fam, lab, r["trait"])] = (bs, ses, r["p"])
        say("  %-5s %-20s %s" % (fam, short(lab, 20), " | ".join("%s %.3f (se %.3f, p %.0e)" % (t[:6], *est[(fam, lab, t)]) for t in TRAITS if (fam, lab, t) in est)))
    out = []
    for fam, lab in CONTRAST:
        for a, b in COMPARE:
            if (fam, lab, a) in est and (fam, lab, b) in est:
                (ba, sa, pa), (bb, sb, pb) = est[(fam, lab, a)], est[(fam, lab, b)]
                z = (ba - bb) / math.sqrt(sa ** 2 + sb ** 2) if np.isfinite(sa) and np.isfinite(sb) and (sa + sb) > 0 else np.nan
                p1 = 0.5 * math.erfc(z / math.sqrt(2)) if np.isfinite(z) else np.nan
                out.append({"family": fam, "label": lab, "trait_a": a, "trait_b": b, "beta_std_a": ba, "beta_std_b": bb, "z_diff": z, "p_one_sided": p1})
    O = pd.DataFrame(out)
    O.to_csv(OUT / "r3_trait_contrast.csv", index=False)
    for (fam, lab), g in O.groupby(["family", "label"], sort=False):
        say("  [R3] %-5s %-20s %s" % (fam, short(lab, 20), " | ".join("%s>%s z %.2f p %.3f" % (r["trait_a"][:6], r["trait_b"][:6], r["z_diff"], r["p_one_sided"]) for _, r in g.iterrows())))
    pr = O[(O["family"] == PRIMARY[0]) & (O["label"] == PRIMARY[1]) & (O["trait_a"] == "KOLLER_EUR") & O["trait_b"].isin(["T2D", "HYPERTENSION"])]
    ok = len(pr) == 2 and bool((pr["p_one_sided"] < 0.05).all())
    say("[H_DS_SPECIFIC_VERDICT] %s %s: KOLLER_EUR enrichment exceeds T2D and hypertension (one-sided p < 0.05 both) -> %s" % (
        PRIMARY[0], PRIMARY[1], "SUPPORTED" if ok else "NOT SUPPORTED"))
    SUM["results"]["r3"] = {"supported": ok, "primary": pr.to_dict("records")}


# ----------------------------------------------------------------------------
# O2 locus-gene expression map
# ----------------------------------------------------------------------------
def o2():
    say("")
    say("-" * 72)
    say("[SECTION] O2 expression of the KOLLER_EUR locus genes across HECA cell states")
    lf = latest(RES / "locus_atlas", "p9_v042_*/l1_loci_compartments.csv")
    if lf is None:
        say("[O2] v042 locus table missing")
        return
    L = pd.read_csv(lf)
    L = L[L["trait"] == "KOLLER_EUR"]
    CPM = pd.read_csv(HECA_DE / "mean_cpm_by_group.csv.gz", index_col=0)
    CPM.index = CPM.index.astype(str)
    CPM = CPM.groupby(level=0).mean()
    lab = pd.Series([c.split("/")[0] for c in CPM.columns], index=CPM.columns)
    M = CPM.T.groupby(lab.values).mean().T
    genes = L[["expressed_gene_used", "compartment", "prioritized_gene", "lead_snp"]].dropna(subset=["expressed_gene_used"])
    genes = genes[genes["expressed_gene_used"].astype(str).isin(M.index)].drop_duplicates("expressed_gene_used")
    genes["order"] = genes["compartment"].map({c: i for i, c in enumerate(COMP_ORDER)}).fillna(len(COMP_ORDER))
    genes = genes.sort_values(["order", "expressed_gene_used"])
    X = np.log2(M.loc[genes["expressed_gene_used"].astype(str)] + 1.0)
    Z = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1).replace(0, np.nan), axis=0).fillna(0.0)

    def grp(c):
        if re.match(r"^(dStromal|eStromal|sHormones|dHormones|Fibroblast)", c):
            return 0
        if re.match(r"^(ePV|mPV|uSMC|Arterial|Venous|Lymphatic)", c):
            return 1
        if re.match(r"^(SOX9|Luminal|preLuminal|Ciliated|preCiliated|Glandular|eHormones|MUC5B|KRT5)", c):
            return 2
        return 3
    cols = sorted(Z.columns, key=lambda c: (grp(c), c))
    Z = Z[cols]
    Z.insert(0, "compartment", genes["compartment"].values)
    Z.to_csv(OUT / "o2_locus_gene_expression_z.csv")
    for _, r in genes[genes["compartment"] == "stroma"].iterrows():
        g = str(r["expressed_gene_used"])
        top = M.loc[g].sort_values(ascending=False).head(3)
        say("  [O2] stromal locus %-10s top states: %s" % (g, ", ".join("%s %.0f" % (k, v) for k, v in top.items())))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 6.5,
                             "pdf.fonttype": 42, "axes.titlecolor": "#0b0b0b"})
        cmap = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#f0efec", "#ec835a", "#b8322f"])
        V = Z.drop(columns=["compartment"]).values
        fig, ax = plt.subplots(figsize=(7.1, 0.13 * V.shape[0] + 1.6))
        ax.imshow(np.clip(V, -2.5, 2.5), cmap=cmap, vmin=-2.5, vmax=2.5, aspect="auto")
        ax.set_yticks(range(V.shape[0]))
        ax.set_yticklabels(["%s  (%s)" % (g, c) for g, c in zip(genes["expressed_gene_used"], genes["compartment"])], color="#52514e", fontsize=5.6)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([short(c, 20) for c in cols], rotation=60, ha="right", color="#52514e", fontsize=5.6)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        comp = list(genes["compartment"])
        for i in range(1, len(comp)):
            if comp[i] != comp[i - 1]:
                ax.axhline(i - 0.5, color="white", linewidth=2)
        gc = [grp(c) for c in cols]
        for j in range(1, len(gc)):
            if gc[j] != gc[j - 1]:
                ax.axvline(j - 0.5, color="white", linewidth=2)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-2.5, 2.5))
        cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01)
        cb.set_label("row z-score of log2(CPM + 1)", color="#52514e")
        cb.outline.set_visible(False)
        ax.set_title("Genes at the KOLLER_EUR genome-wide loci across HECA cell states", loc="left", fontweight="bold", fontsize=7.5)
        for ext in ("png", "pdf"):
            fig.savefig(OUT / ("F7_locus_genes_heca.%s" % ext), dpi=300, bbox_inches="tight")
        plt.close(fig)
        say("[O2] figure F7_locus_genes_heca written (%d genes x %d states)" % V.shape)
    except Exception as e:
        say("[O2] figure skipped: %s" % str(e)[:200])



# ----------------------------------------------------------------------------
# I0 inventory (read-only)
# ----------------------------------------------------------------------------
def inventory():
    say("")
    say("-" * 72)
    say("[SECTION] I0 inventory of local single-cell tables (read-only, for the v050 rules)")
    f = PB266 / "pseudobulk_groups.csv"
    if f.exists():
        G = pd.read_csv(f)
        say("[I0] GSE266265 pseudobulk groups: %d rows; columns: %s" % (len(G), ", ".join(G.columns)))
        for c in ("status", "main", "phase5"):
            if c in G.columns:
                say("[I0]   %s: %s" % (c, ", ".join("%s %d" % (k, v) for k, v in G[c].value_counts().items())))
        if "refined" in G.columns:
            agg = G.groupby("refined").agg(rows=("refined", "size"), donors=("sample", "nunique") if "sample" in G.columns else ("refined", "size"),
                                           cells=("n_cells", "sum") if "n_cells" in G.columns else ("refined", "size"))
            if "status" in G.columns and "sample" in G.columns:
                for st in G["status"].dropna().unique():
                    agg["donors_" + str(st)] = G[G["status"] == st].groupby("refined")["sample"].nunique()
            say("[I0]   refined types (%d):" % len(agg))
            for k, r in agg.sort_values("cells", ascending=False).iterrows():
                say("[I0]     %-28s %s" % (k, "  ".join("%s %s" % (c, int(v) if pd.notna(v) else 0) for c, v in r.items())))
        gf = PB266 / "genes.csv"
        if gf.exists():
            say("[I0]   genes: %d" % len(pd.read_csv(gf)))
    else:
        say("[I0] GSE266265 pseudobulk not found: " + str(f))
    f = CELLS213 / "obs.csv"
    if f.exists():
        O = pd.read_csv(f)
        say("[I0] GSE213216 cells (v028 cache): %d cells; columns: %s" % (len(O), ", ".join(O.columns)))
        for c in ("context", "lineage"):
            if c in O.columns:
                say("[I0]   %s: %s" % (c, ", ".join("%s %d" % (k, v) for k, v in O[c].value_counts().items())))
        if {"label", "context"} <= set(O.columns):
            T = pd.crosstab(O["label"], O["context"])
            if "patient" in O.columns:
                T["patients"] = O.groupby("label")["patient"].nunique()
            if "lineage" in O.columns:
                T["lineage"] = O.groupby("label")["lineage"].agg(lambda x: x.value_counts().index[0])
            say("[I0]   label x context (%d labels; contexts %s):" % (len(T), ", ".join(str(c) for c in T.columns)))
            for k, r in T.iterrows():
                say("[I0]     %-24s %s" % (k, "  ".join(str(v) for v in r.values)))
    else:
        say("[I0] GSE213216 cell cache not found: " + str(f))


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.49 | review-proofing: threshold sensitivity, trait contrast, locus-gene map | " + TS)
    say("=" * 72)
    t0 = time.time()
    sym2ent = gene_loc()
    for name, fn in (("R1_R3", lambda: r1_r3(sym2ent)), ("O2", o2), ("I0", inventory)):
        try:
            fn()
        except Exception as e:
            import traceback
            ERR.append(name)
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v049_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.39 before this run")
    if ERR:
        say("[INCOMPLETE] steps with errors: " + ", ".join(ERR))
    else:
        say("[DONE] PAPER9_V049_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
