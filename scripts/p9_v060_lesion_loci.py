#!/usr/bin/env python3
# Paper 9 | v0.0.60 | endometriosis locus genes inside lesions (rules CHANGELOG v0.2.50, fixed before the run)
#   LES-L1 H-LESION-LOC: do the genes of the KOLLER_EUR loci keep their HECA compartment inside peritoneal lesions and
#      endometriomas? Per sample: log2 ratio of a gene's CPM in its lineage group vs the rest of the sample (GSE213216,
#      v052 marker-based pseudobulk); set mean vs expression-matched random gene sets. Eutopic endometrium and
#      unaffected peritoneum reported as references; endothelial / epithelial / immune / perivascular loci as design controls.
#   LES-L2 descriptive: stromal-locus genes in endometrial-type stroma (eutopic, lesions, endometriomas) vs stroma of
#      unaffected peritoneum (host fibroblasts).
#   S5 figure: locus-gene expression by lineage and tissue context; set statistics vs null.
# Localization of genetic signal only: no PGR-response / decidualization tests (Paper 7), no lesion-vs-eutopic DE claims
# (Paper 8). Reads the v052 cache (data/processed/gse213216_pseudobulk_v052); no MAGMA, no downloads; ~2-5 min.

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
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
PROC = REPO / "data" / "processed"
RES = REPO / "results"
PB_DIR = PROC / "gse213216_pseudobulk_v052"
LOCUS_DIR = RES / "locus_atlas"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = RES / "review" / ("p9_v060_" + TS)
LOG = REPO / "logs" / ("p9_v060_lesion_loci_" + TS + ".log")

# Fixed (CHANGELOG v0.2.50)
TRAIT = "KOLLER_EUR"
CONTEXTS = ["eutopic", "peritoneal_lesion", "endometrioma", "peritoneum_unaffected"]
LESIONS = ["peritoneal_lesion", "endometrioma"]
CTX_NAME = {"eutopic": "eutopic endometrium", "peritoneal_lesion": "peritoneal lesions", "endometrioma": "endometriomas",
            "peritoneum_unaffected": "unaffected peritoneum"}
LINEAGES = ["stroma", "perivascular / smooth muscle", "endothelium", "epithelium", "immune", "mesothelium"]
COMP_TO_LINEAGE = {"stroma": "stroma", "perivascular/SM": "perivascular / smooth muscle", "endothelium": "endothelium",
                   "epithelium": "epithelium", "immune": "immune"}
COMP_ORDER = ["stroma", "perivascular/SM", "endothelium", "epithelium", "immune"]
MIN_CELLS_UNIT = 10       # as v052
MIN_KEYS = 4              # donors (or samples without donor id) per context and lineage
MIN_GENES_SET = 3
MIN_CPM = 1.0
N_NULL = 10000
N_BINS = 10
SEED = 20261111

for d in (OUT, LOG.parent):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V060_LESION_LOCI", "timestamp_utc": TS, "results": {}}


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def latest(base, pattern):
    fs = sorted(base.glob(pattern))
    return fs[-1] if fs else None


def lineage_of(group):
    g = str(group)
    return "stroma" if g.startswith("stroma") else g


def tsf_one(t, df):
    # one-sided upper-tail p of Student t
    from scipy import stats
    return float(stats.t.sf(t, df))


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "-"
    return "%.3f" % p if p >= 0.001 else "%.0e" % p


# ----------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------
def load():
    f_pb, f_u = PB_DIR / "pseudobulk.csv.gz", PB_DIR / "units.csv"
    if not f_pb.exists() or not f_u.exists():
        raise FileNotFoundError("v052 pseudobulk cache missing: " + str(PB_DIR))
    PB = pd.read_csv(f_pb, index_col=0)
    PB.index = PB.index.astype(str)
    PB = PB.groupby(level=0).sum()
    U = pd.read_csv(f_u)
    U = U[U["context"].isin(CONTEXTS)].copy()
    U["lineage"] = U["group"].map(lineage_of)
    U["key"] = [("P:" + str(p)) if pd.notna(p) and str(p) not in ("", "nan", "None") else ("S:" + str(g))
                for p, g in zip(U["patient"], U["gsm"])]
    lf = latest(LOCUS_DIR, "p9_v042_*/l1_loci_compartments.csv")
    if lf is None:
        raise FileNotFoundError("v042 locus table missing (results/locus_atlas/p9_v042_*/l1_loci_compartments.csv)")
    L = pd.read_csv(lf)
    L = L[L["trait"] == TRAIT].copy()
    L["gene"] = L["expressed_gene_used"].astype(str).str.strip()
    L = L[L["compartment"].isin(COMP_ORDER) & (L["gene"] != "") & (L["gene"] != "nan")]
    L = L.drop_duplicates("gene")
    say("[INPUT] pseudobulk %d genes x %d units; locus table %s" % (PB.shape[0], PB.shape[1], lf.relative_to(REPO)))
    say("[INPUT] %s locus genes by HECA compartment: %s" % (TRAIT, "; ".join(
        "%s %d (%s)" % (c, int((L["compartment"] == c).sum()), ",".join(L.loc[L["compartment"] == c, "gene"])) for c in COMP_ORDER)))
    return PB, U, L


def context_tables(PB, U, ctx):
    # per sample: lineage CPM, rest-of-sample CPM; aggregated to donor keys
    T = U[(U["context"] == ctx) & (U["group"] == "__total__")]
    S = U[(U["context"] == ctx) & (U["group"] != "__total__") & (U["n_cells"] >= MIN_CELLS_UNIT)]
    genes = PB.index.values
    tot_cpm = {}
    for _, r in T.iterrows():
        x = PB[r["unit"]].values.astype(float)
        tot_cpm[r["gsm"]] = x / max(x.sum(), 1.0) * 1e6
    lfc, cpm = {}, {}
    for lin in LINEAGES:
        s = S[S["lineage"] == lin]
        per_key_l, per_key_c = {}, {}
        for _, r in s.iterrows():
            tu = T[T["gsm"] == r["gsm"]]
            if tu.empty:
                continue
            x = PB[r["unit"]].values.astype(float)
            rest = PB[tu["unit"].iloc[0]].values.astype(float) - x
            rest = np.maximum(rest, 0.0)
            if rest.sum() <= 0:
                continue
            c_l = x / max(x.sum(), 1.0) * 1e6
            c_r = rest / rest.sum() * 1e6
            per_key_l.setdefault(r["key"], []).append(np.log2(c_l + 1) - np.log2(c_r + 1))
            per_key_c.setdefault(r["key"], []).append(c_l)
        if len(per_key_l) >= 1:
            lfc[lin] = pd.DataFrame({k: np.mean(v, axis=0) for k, v in per_key_l.items()}, index=genes)
            cpm[lin] = pd.DataFrame({k: np.mean(v, axis=0) for k, v in per_key_c.items()}, index=genes)
    mean_tot = pd.Series(np.mean(list(tot_cpm.values()), axis=0), index=genes) if tot_cpm else pd.Series(0.0, index=genes)
    n_samples = len(tot_cpm)
    return lfc, cpm, mean_tot, n_samples


def matched_null(mean_tot, exclude, genes_obs, stat, rng):
    # random gene sets matched on expression decile (mean total CPM), locus genes excluded from the pool
    expr = mean_tot[mean_tot >= MIN_CPM]
    logm = np.log10(expr + 1e-9)
    bins = pd.qcut(logm.rank(method="first"), N_BINS, labels=False)
    pool = {b: np.array([g for g in bins.index[bins == b] if g not in exclude]) for b in range(N_BINS)}
    want = [int(bins.get(g)) for g in genes_obs if g in bins.index]
    vals = stat.reindex(expr.index).values
    idx = {g: i for i, g in enumerate(expr.index)}
    pool_i = {b: np.array([idx[g] for g in pool[b]]) for b in pool}
    if not want:
        return None
    parts = []
    for b in sorted(set(want)):
        k = want.count(b)
        if len(pool_i[b]) < 10 * k:
            return None
        parts.append(pool_i[b][rng.integers(0, len(pool_i[b]), size=(N_NULL, k))])
    sel = np.concatenate(parts, axis=1)
    return np.nanmean(vals[sel], axis=1)


# ----------------------------------------------------------------------------
# L1 locus genes keep their compartment inside lesions
# ----------------------------------------------------------------------------
def l1(PB, U, L):
    say("")
    say("-" * 72)
    say("[SECTION] L1 locus genes by lineage inside lesions (log2 CPM lineage vs rest of sample; donor means)")
    rng = np.random.default_rng(SEED)
    exclude = set(L["gene"])
    gene_rows, set_rows, heat = [], [], {}
    for ctx in CONTEXTS:
        lfc, cpm, mean_tot, n_s = context_tables(PB, U, ctx)
        say("[L1] %-22s samples %d; donor keys per lineage: %s" % (CTX_NAME[ctx], n_s, ", ".join(
            "%s %d" % (lin.split(" ")[0], lfc[lin].shape[1]) for lin in LINEAGES if lin in lfc)))
        for lin in LINEAGES:
            if lin in cpm:
                heat[(ctx, lin)] = cpm[lin].mean(axis=1)
        top_lin = pd.DataFrame({lin: cpm[lin].mean(axis=1) for lin in cpm}) if cpm else pd.DataFrame()
        for comp in COMP_ORDER:
            lin = COMP_TO_LINEAGE[comp]
            G = [g for g in L.loc[L["compartment"] == comp, "gene"] if g in PB.index]
            if lin not in lfc or lfc[lin].shape[1] < MIN_KEYS:
                if G:
                    say("  [L1] %-22s %-16s not evaluable (%s donor keys with >= %d cells)" % (
                        CTX_NAME[ctx], comp, lfc[lin].shape[1] if lin in lfc else 0, MIN_CELLS_UNIT))
                continue
            M = lfc[lin]
            n = M.shape[1]
            Gx = [g for g in G if mean_tot.get(g, 0) >= MIN_CPM]
            for g in G:
                v = M.loc[g].values
                m, sd = float(np.mean(v)), float(np.std(v, ddof=1)) if n > 1 else float("nan")
                t = m / max(sd / math.sqrt(n), 1e-9) if np.isfinite(sd) else float("nan")
                tl = top_lin.loc[g].idxmax() if (not top_lin.empty and top_lin.loc[g].max() > 0) else "-"
                gene_rows.append({"context": ctx, "compartment": comp, "lineage": lin, "gene": g, "expressed": g in Gx,
                                  "n_keys": n, "mean_log2_ratio": m, "t": t, "p_one_sided": tsf_one(t, n - 1) if np.isfinite(t) else np.nan,
                                  "lineage_cpm": float(cpm[lin].loc[g].mean()), "top_lineage": tl})
            if len(Gx) < MIN_GENES_SET:
                say("  [L1] %-22s %-16s fewer than %d expressed genes (%d)" % (CTX_NAME[ctx], comp, MIN_GENES_SET, len(Gx)))
                continue
            stat = M.mean(axis=1)
            obs = float(stat.loc[Gx].mean())
            null = matched_null(mean_tot, exclude, Gx, stat, rng)
            if null is None:
                say("  [L1] %-22s %-16s matched pool too small" % (CTX_NAME[ctx], comp))
                continue
            p = (1 + int((null >= obs).sum())) / (N_NULL + 1)
            n_top = sum(1 for g in Gx if (not top_lin.empty and top_lin.loc[g].idxmax() == lin))
            gl = [r for r in gene_rows if r["context"] == ctx and r["compartment"] == comp and r["expressed"]]
            n_sig = sum(1 for r in gl if np.isfinite(r["p_one_sided"]) and r["p_one_sided"] < 0.05)
            set_rows.append({"context": ctx, "compartment": comp, "lineage": lin, "n_genes": len(Gx), "n_keys": n, "observed": obs,
                             "null_median": float(np.median(null)), "null_lo": float(np.percentile(null, 2.5)),
                             "null_hi": float(np.percentile(null, 97.5)), "null_p95": float(np.percentile(null, 95)),
                             "p_perm": p, "n_top_lineage": n_top, "n_gene_p05": n_sig})
            say("  [L1] %-22s %-16s genes %d, donors %2d: mean log2 ratio %+.2f (null median %+.2f, 95th %+.2f) perm p %s | top lineage %d/%d | gene p<0.05 %d/%d" % (
                CTX_NAME[ctx], comp, len(Gx), n, obs, np.median(null), np.percentile(null, 95), fmt_p(p), n_top, len(Gx), n_sig, len(Gx)))
    G = pd.DataFrame(gene_rows)
    S = pd.DataFrame(set_rows)
    G.to_csv(OUT / "l1_locus_genes_by_context.csv", index=False)
    S.to_csv(OUT / "l1_set_tests.csv", index=False)
    st = G[(G["compartment"] == "stroma")]
    for ctx in CONTEXTS:
        s = st[st["context"] == ctx]
        if len(s):
            say("  [L1] stromal-locus genes, %-22s: %s" % (CTX_NAME[ctx], ", ".join(
                "%s %+.1f%s" % (r["gene"], r["mean_log2_ratio"], "" if r["top_lineage"] == "stroma" else " (top %s)" % str(r["top_lineage"]).split(" ")[0])
                for _, r in s.iterrows())))
    ok = []
    for ctx in LESIONS:
        r = S[(S["context"] == ctx) & (S["compartment"] == "stroma")]
        ok.append(bool(len(r) and r["p_perm"].iloc[0] < 0.05))
    verdict = "SUPPORTED" if all(ok) and len(ok) == 2 else "NOT SUPPORTED"
    ctrl = S[S["context"].isin(LESIONS) & (S["compartment"] != "stroma")]
    say("[H_LESION_LOC_VERDICT] stromal-locus genes stroma-enriched beyond matched random genes in peritoneal lesions (%s) and endometriomas (%s) -> %s" % (
        "yes" if ok[0] else "no", "yes" if ok[1] else "no", verdict))
    if len(ctrl):
        say("[L1_DESIGN_CONTROL] other compartments in lesions (perm p < 0.05): %s" % "; ".join(
            "%s %s %s" % (CTX_NAME[r["context"]], r["compartment"], fmt_p(r["p_perm"])) for _, r in ctrl.iterrows()))
    SUM["results"]["L1"] = {"verdict": verdict, "sets": S.to_dict(orient="records")}
    return G, S, heat


# ----------------------------------------------------------------------------
# L2 endometrial-type stroma vs host peritoneal stroma (descriptive)
# ----------------------------------------------------------------------------
def l2(PB, U, L):
    say("")
    say("-" * 72)
    say("[SECTION] L2 stromal-locus genes: endometrial-type stroma vs stroma of unaffected peritoneum (descriptive)")
    rng = np.random.default_rng(SEED + 1)
    tabs = {}
    for ctx in CONTEXTS:
        _, cpm, mean_tot, _ = context_tables(PB, U, ctx)
        if "stroma" in cpm:
            tabs[ctx] = np.log2(cpm["stroma"] + 1)
    ref = "peritoneum_unaffected"
    if ref not in tabs or tabs[ref].shape[1] < MIN_KEYS:
        say("[L2] unaffected-peritoneum stroma not evaluable")
        return pd.DataFrame()
    G = [g for g in L.loc[L["compartment"] == "stroma", "gene"] if g in PB.index]
    exclude = set(L["gene"])
    rows = []
    B = tabs[ref]
    for ctx in ["eutopic"] + LESIONS:
        if ctx not in tabs or tabs[ctx].shape[1] < MIN_KEYS:
            continue
        A = tabs[ctx]
        shared = sorted(set(A.columns) & set(B.columns))
        for g in G:
            a, b = A.loc[g].values, B.loc[g].values
            d = float(a.mean() - b.mean())
            se = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
            rows.append({"context": ctx, "gene": g, "n_ctx": len(a), "n_ref": len(b), "diff_log2": d,
                         "ci_lo": d - 1.96 * se, "ci_hi": d + 1.96 * se, "shared_donors": len(shared)})
        diff = A.mean(axis=1) - B.mean(axis=1)
        mean_all = (2 ** A - 1).mean(axis=1) * 0.5 + (2 ** B - 1).mean(axis=1) * 0.5
        Gx = [g for g in G if mean_all.get(g, 0) >= MIN_CPM]
        obs = float(diff.loc[Gx].mean())
        null = matched_null(mean_all, exclude, Gx, diff, rng)
        p = (1 + int((null >= obs).sum())) / (N_NULL + 1) if null is not None else np.nan
        say("  [L2] %-20s (%d donors) vs unaffected peritoneum (%d): stromal-locus mean diff %+.2f log2 (matched null median %+.2f) p %s; shared donors %d" % (
            CTX_NAME[ctx], A.shape[1], B.shape[1], obs, float(np.median(null)) if null is not None else float("nan"), fmt_p(p), len(shared)))
        say("       per gene: " + ", ".join("%s %+.1f" % (r["gene"], r["diff_log2"]) for r in rows if r["context"] == ctx))
        SUM["results"].setdefault("L2", {})[ctx] = {"mean_diff": obs, "p_matched": p, "n_genes": len(Gx)}
    D = pd.DataFrame(rows)
    D.to_csv(OUT / "l2_host_stroma_comparison.csv", index=False)
    return D


# ----------------------------------------------------------------------------
# S5 figure
# ----------------------------------------------------------------------------
def s5(G, S, heat, L):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#898781"
    INK, INK2, MUTED, GRID, AXIS, LIGHT = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#d9d8d2"
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7,
        "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6.5,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
        "text.color": INK, "axes.titlecolor": INK, "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
    cmap = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#f0efec", "#ec835a", "#b8322f"])
    LIN_SHORT = {"stroma": "Str", "perivascular / smooth muscle": "PV/SM", "endothelium": "Endo", "epithelium": "Epi",
                 "immune": "Imm", "mesothelium": "Meso"}
    CTX_SHORT = {"eutopic": "Eutopic", "peritoneal_lesion": "Peritoneal lesion", "endometrioma": "Endometrioma",
                 "peritoneum_unaffected": "Unaffected peritoneum"}
    CTX_HEAD = {"eutopic": "Eutopic", "peritoneal_lesion": "Peritoneal\nlesion", "endometrioma": "Endo-\nmetrioma",
                "peritoneum_unaffected": "Unaffected\nperitoneum"}
    COMP_SHORT = {"stroma": "stroma", "perivascular/SM": "PV/SM", "endothelium": "endoth.", "epithelium": "epith.", "immune": "imm."}
    cols = [(c, l) for c in CONTEXTS for l in LINEAGES if (c, l) in heat]
    genes = []
    for comp in COMP_ORDER:
        genes += [g for g in L.loc[L["compartment"] == comp, "gene"] if all(g in heat[k].index for k in cols)]
    if not cols or not genes:
        say("[S5] nothing to plot")
        return
    H = pd.DataFrame({k: heat[k].reindex(genes) for k in cols})
    X = np.log2(H + 1.0)
    Z = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    comp_of = L.set_index("gene")["compartment"].to_dict()
    nrow = len(genes)
    fig = plt.figure(figsize=(7.2, max(3.8, 0.115 * nrow + 1.6)))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[1.0, 0.025], wspace=0.95, hspace=0.42)
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(Z.values, aspect="auto", cmap=cmap, vmin=-2.5, vmax=2.5, interpolation="nearest")
    ax.set_yticks(range(nrow))
    ax.set_yticklabels(genes, fontsize=5.4)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([LIN_SHORT[l] for _, l in cols], rotation=90, fontsize=5.2)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    ax.tick_params(length=0)
    xs = 0
    for c in CONTEXTS:
        k = sum(1 for cc, _ in cols if cc == c)
        if k:
            ax.text(xs + (k - 1) / 2, -0.8, CTX_HEAD[c], ha="center", va="bottom", fontsize=5.6, color=INK, linespacing=0.95)
            if xs > 0:
                ax.axvline(xs - 0.5, color="white", linewidth=1.8)
            xs += k
    prev = None
    for i, g in enumerate(genes):
        c = comp_of.get(g)
        if prev is not None and c != prev:
            ax.axhline(i - 0.5, color="white", linewidth=1.8)
        prev = c
    for comp in COMP_ORDER:
        ii = [i for i, g in enumerate(genes) if comp_of.get(g) == comp]
        if ii:
            ax.text(len(cols) - 0.1, np.mean(ii), COMP_SHORT[comp], va="center", ha="left", fontsize=5.6, color=INK2, rotation=90)
    ax.set_title("a  Locus genes by lineage and tissue", loc="left", fontweight="bold", pad=20)
    cax = fig.add_subplot(gs[1, 0])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=5.5, length=0)
    cb.set_label("row z-score of log2 CPM", fontsize=5.8)
    ax2 = fig.add_subplot(gs[0, 1])
    R = S.copy()
    R["order_c"] = R["compartment"].map({c: i for i, c in enumerate(COMP_ORDER)})
    R["order_x"] = R["context"].map({c: i for i, c in enumerate(CONTEXTS)})
    R = R.sort_values(["order_c", "order_x"]).reset_index(drop=True)
    y = np.arange(len(R))[::-1]
    for yi, (_, r) in zip(y, R.iterrows()):
        ax2.plot([r["null_lo"], r["null_hi"]], [yi, yi], color=LIGHT, linewidth=4.0, solid_capstyle="butt")
        ax2.plot([r["null_median"]], [yi], marker="|", color=GRAY, markersize=6)
        col = ORANGE if r["p_perm"] < 0.05 else INK2
        ax2.plot([r["observed"]], [yi], marker="o", color=col, markersize=3.6)
        ax2.text(max(r["observed"], r["null_hi"]) + 0.15, yi, fmt_p(r["p_perm"]), va="center", fontsize=5.4, color=INK2)
    ax2.set_yticks(y)
    ax2.set_yticklabels(["%s loci · %s" % (r["compartment"].replace("perivascular/SM", "PV/SM"), CTX_SHORT[r["context"]]) for _, r in R.iterrows()], fontsize=5.6)
    ax2.axvline(0, color=AXIS, linewidth=0.6)
    ax2.spines["left"].set_visible(False)
    ax2.tick_params(axis="y", length=0)
    ax2.grid(axis="x", color=GRID, linewidth=0.5)
    ax2.set_axisbelow(True)
    ax2.set_xlabel("mean log2 ratio, own lineage vs rest of sample")
    ax2.set_title("b  Own-lineage enrichment\n    vs matched random genes", loc="left", fontweight="bold")
    data = Z.copy()
    data.insert(0, "compartment", [comp_of.get(g) for g in genes])
    data.columns = ["compartment"] + ["%s|%s" % k for k in cols]
    for ext in ("png", "pdf"):
        fig.savefig(OUT / ("S5_lesion_locus_genes.%s" % ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    data.to_csv(OUT / "S5_lesion_locus_genes_data.csv")
    R.to_csv(OUT / "S5b_set_tests_data.csv", index=False)
    leg = ("Genes of the %s risk loci inside endometriosis lesions (GSE213216, Fonseca et al. 2023; v052 marker-based lineage "
           "pseudobulk, >= %d cells per sample and lineage, donor means). (a) Expression of each locus gene (HECA compartment from the "
           "locus atlas) across lineages in eutopic endometrium, peritoneal lesions, endometriomas and unaffected peritoneum; row "
           "z-score of log2 CPM. (b) Mean log2 ratio of the locus genes in their own lineage versus the rest of the same sample "
           "(dot; orange, permutation p < 0.05), with the median (tick) and 95%% range (bar) of %d random gene sets matched on "
           "expression decile; one-sided permutation p at right." % (TRAIT, MIN_CELLS_UNIT, N_NULL))
    (OUT / "S5_legend.md").write_text("## S5_lesion_locus_genes\n\n" + leg + "\n", encoding="utf-8")
    say("[FIG] S5_lesion_locus_genes written (%d genes x %d columns)" % (nrow, len(cols)))


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.60 | endometriosis locus genes inside lesions (GSE213216) | " + TS)
    say("=" * 72)
    try:
        PB, U, L = load()
    except FileNotFoundError as e:
        say("[FATAL] " + str(e))
        _fh.close()
        return
    G, S, heat = l1(PB, U, L)
    l2(PB, U, L)
    try:
        s5(G, S, heat, L)
    except Exception as e:
        say("[S5] figure failed: %s: %s" % (type(e).__name__, e))
    (OUT / "summary.json").write_text(json.dumps(SUM, indent=1, default=str), encoding="utf-8")
    say("")
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V060_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
