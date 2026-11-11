#!/usr/bin/env python3
# Paper 9 | v0.0.56 | figures for the review audits (v054, v055) and the final figure set (no new analyses)
#   F9  donor structure and donor stability: largest-donor share of every HECA status x label group (v055 E1);
#       first-ranked group across leave-one-donor-out runs in GSE266265 (v054 B1) and three HECA designs (v054 B2, v055 E2, E3)
#   S2  trait contrasts with GWAS covariance (v054 A2): BETA_STD +- 95% CI per trait, covariance-aware p
#   S3  cross-atlas locus test (v054 C): random-subset null vs the observed loss of significance, both directions
#   S4  benchmark intervals (v054 D): top-1 exact recovery with Wilson 95% CI under specificity truths
#   Collect: F1-F8, S1 (final_v053) + F9, S2-S4 into results/paper_figures/final_v056_<TS> with combined legends.
# Writes PNG (300 dpi), PDF and the plotted values (CSV) per figure. ~1 min, no network.

import datetime as dt
import math
import re
import shutil
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RES = REPO / "results"
INTERIM = REPO / "data" / "interim"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = RES / "paper_figures" / ("p9_v056_" + TS)
FINAL = RES / "paper_figures" / ("final_v056_" + TS)
LOG = REPO / "logs" / ("p9_v056_figures_" + TS + ".log")
for d in (OUT, FINAL, LOG.parent):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")

# Palette as v048 (validated)
BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#898781"
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
LIGHT = "#d9d8d2"
FAMS = ["ABS", "UP", "DOWN", "ABSN", "UPN", "DOWNN", "SPEC"]
FAM_NAME = {"ABS": "|log2FC|", "UP": "up-only", "DOWN": "down-only", "ABSN": "|log2FC|, N=400", "UPN": "up-only, N=400",
            "DOWNN": "down-only, N=400", "SPEC": "specificity"}
STROMA_H = r"^(dStromal|eStromal|sHormones|dHormones|Fibroblast)"
VASC_H = r"^(Arterial|Venous|Lymphatic|ePV|mPV|uSMC)"
TRAIT_NAME = {"KOLLER_EUR": "KOLLER_EUR", "GCST90205183": "GCST90205183", "ENDO_ALL": "FinnGen ENDO_ALL", "LEIOMYOMA": "Leiomyoma",
              "HYPERTENSION": "Hypertension", "T2D": "Type 2 diabetes"}
ENDO_TRAITS = ("KOLLER_EUR", "GCST90205183", "ENDO_ALL")
LEG = []

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7,
    "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "axes.titlecolor": INK, "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 0})


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def latest(base, pattern):
    fs = sorted(base.glob(pattern))
    return fs[-1] if fs else None


def need(base, pattern, what):
    f = latest(base, pattern)
    if f is None:
        raise FileNotFoundError("%s (%s/%s)" % (what, base.name, pattern))
    say("  input %-30s %s" % (what, f.relative_to(REPO)))
    return f


def nice(label, n=30):
    s = str(label).replace("SOX9_functionalis_", "SOX9 functionalis ").replace("dStromal_", "dStromal ").replace("eStromal_", "eStromal ")
    s = s.replace("Glandular_secretory", "Glandular secretory").replace("_", " ")
    return s if len(s) <= n else s[:n - 1] + "…"


def style(ax, xgrid=True):
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    if xgrid:
        ax.grid(axis="x", color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)


def save(fig, name, data, legend):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / ("%s.%s" % (name, ext)), dpi=300, bbox_inches="tight")
    plt.close(fig)
    if data is not None:
        data.to_csv(OUT / ("%s_data.csv" % name), index=False)
    LEG.append("## %s\n\n%s\n" % (name, legend))
    say("[FIG] %s written" % name)


def fmt_p(p):
    if not np.isfinite(p):
        return "-"
    return "%.3f" % p if p >= 0.001 else "%.0e" % p


# ----------------------------------------------------------------------------
# F9 donor structure and donor stability
# ----------------------------------------------------------------------------
def cls_heca(label):
    if re.search(STROMA_H, str(label)):
        return "stroma"
    if re.search(VASC_H, str(label)):
        return "perivascular / vascular"
    return "other"


def cls_gse(group):
    if group == "endometrial stroma":
        return "stroma"
    if group in ("perivascular / MSC-like", "endothelium"):
        return "perivascular / vascular"
    return "other"


def f9():
    e1 = need(RES / "review", "p9_v055_*/e1_heca_donor_composition.csv", "v055 E1 composition")
    b1 = need(RES / "review", "p9_v054_*/b1_gse266265_loo.csv", "v054 B1 GSE266265 LOO")
    b2 = need(RES / "review", "p9_v054_*/b2_b3_heca_loo_pseudobulk.csv", "v054 B2 HECA LOO")
    e2 = need(RES / "review", "p9_v055_*/e2_balanced_pooled_loo.csv", "v055 E2 balanced")
    e3 = need(RES / "review", "p9_v055_*/e3_pooled_pseudobulk_loo.csv", "v055 E3 pooled pseudobulk")
    E = pd.read_csv(e1)
    E = E[E["variant"] == "P"].sort_values("top_share", ascending=True).reset_index(drop=True)
    rows = []
    B1 = pd.read_csv(b1)
    for t in ENDO_TRAITS:
        g = B1[(B1["trait"] == t) & (B1["left_out"] != "FULL")]
        if len(g):
            c = g["first_group"].map(cls_gse).value_counts()
            rows.append({"design": "GSE266265, donor-paired DE", "trait": t, "n": len(g), "need": math.ceil(0.9 * len(g)) if t != "ENDO_ALL" else np.nan,
                         **{k: int(c.get(k, 0)) for k in ("stroma", "perivascular / vascular", "other")}})
    for f, design, fam in ((b2, "HECA status-split, cell-level", "UPN"), (e2, "HECA pooled, donor-balanced", "UPN"), (e3, "HECA pooled, pseudobulk", "UPN")):
        D = pd.read_csv(f)
        g = D[(D["trait"] == "KOLLER_EUR") & (D["family"] == fam) & (D["left_out"] != "FULL") & (D["left_out"] != "-")]
        if len(g):
            c = g["first"].map(cls_heca).value_counts()
            rows.append({"design": design, "trait": "KOLLER_EUR", "n": len(g), "need": len(g) - 1,
                         **{k: int(c.get(k, 0)) for k in ("stroma", "perivascular / vascular", "other")}})
    L = pd.DataFrame(rows)
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.75)
    ax = fig.add_subplot(gs[0, 0])
    yy = np.arange(len(E))
    col = [ORANGE if re.search(STROMA_H, str(l)) else GRAY for l in E["label"]]
    ax.hlines(yy, 0, E["top_share"], color=GRID, linewidth=0.8, zorder=1)
    ax.scatter(E["top_share"], yy, s=16, c=col, edgecolor="white", linewidth=0.5, zorder=3)
    for x in (0.5, 0.8):
        ax.axvline(x, color=AXIS, linestyle=(0, (2, 2)), linewidth=0.6, zorder=0)
    ax.set_yticks(yy)
    ax.set_yticklabels([nice(g, 28) for g in E["group"]], fontsize=5.2, color=INK2)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(-0.8, len(E) - 0.2)
    ax.set_xlabel("share of cells from the largest donor")
    ax.set_title("a  HECA groups (cell type x status)", loc="left", fontweight="bold")
    style(ax)
    ax.legend(handles=[Line2D([], [], marker="o", linestyle="", markersize=4, color=ORANGE, label="stromal state"),
                       Line2D([], [], marker="o", linestyle="", markersize=4, color=GRAY, label="other state")],
              loc="lower right", frameon=False, handletextpad=0.3)
    ax2 = fig.add_subplot(gs[0, 1])
    cols = [("stroma", ORANGE), ("perivascular / vascular", BLUE), ("other", LIGHT)]
    yy2 = np.arange(len(L))[::-1]
    for yi, (_, r) in zip(yy2, L.iterrows()):
        left = 0.0
        for k, c in cols:
            w = r[k] / r["n"]
            if w > 0:
                ax2.barh(yi, w - (0.006 if left + w < 0.999 else 0), left=left, height=0.56, color=c, edgecolor="none")
            left += w
        if np.isfinite(r["need"]):
            xn = r["need"] / r["n"]
            ax2.plot([xn, xn], [yi - 0.38, yi + 0.38], color=INK, linewidth=1.0)
        ax2.text(1.02, yi, "%d/%d" % (r["stroma"], r["n"]), va="center", fontsize=6.2, color=INK2)
    ax2.set_yticks(yy2)
    ax2.set_yticklabels(["%s\n%s" % (r["design"], TRAIT_NAME.get(r["trait"], r["trait"])) for _, r in L.iterrows()], fontsize=6, color=INK2)
    ax2.set_xlim(0, 1.12)
    ax2.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax2.set_xticklabels(["0", "25", "50", "75", "100%"])
    ax2.set_xlabel("leave-one-donor-out runs by first-ranked group (up-only, N = 400)")
    ax2.set_title("b  Donor stability of the first-ranked group", loc="left", fontweight="bold")
    style(ax2)
    ax2.legend(handles=[Patch(color=c, label=k) for k, c in cols] + [Line2D([], [], color=INK, linewidth=1.0, label="pre-registered threshold")],
               loc="lower center", bbox_to_anchor=(0.45, 1.06), ncol=2, frameon=False, handletextpad=0.4, columnspacing=1.0)
    D = pd.concat([E.assign(panel="a"), L.assign(panel="b")], ignore_index=True, sort=False)
    nhi = int((E["top_share"] > 0.8).sum())
    leg = ("(a) Share of cells contributed by the largest donor in each HECA cell type x disease-status group with >= 30 cells (variant P "
           "cells, the design of Koller et al.; %d of %d groups > 80%%, median effective donors %.1f). Stromal states in orange. "
           "(b) Leave-one-donor-out stability of the first-ranked group (MAGMA, size-matched up-only sets): GSE266265 donor-paired DE "
           "(60 donors; three endometriosis GWAS) and three HECA designs for KOLLER_EUR: status-split labels with cell-level DE (7 "
           "endometriosis donors left out; the pre-registered rule there concerned decidual states, 6/7), labels pooled over status with "
           "<= 50 cells per donor and type (14 donors), and donor-paired pseudobulk DE with pooled labels (14 donors). Right: runs with "
           "a stromal state first. Black ticks: pre-registered thresholds (90%% of runs, or all runs but one)." % (nhi, len(E), E["effective_donors"].median()))
    save(fig, "F9_donor_structure", D, leg)


# ----------------------------------------------------------------------------
# S2 trait contrasts with GWAS covariance
# ----------------------------------------------------------------------------
def s2():
    f = need(RES / "review", "p9_v054_*/a2_trait_contrasts_rotation.csv", "v054 A2 contrasts")
    C = pd.read_csv(f)
    panels = [("HECA", "UPN", "dStromal_mid/ENDO", "a  HECA dStromal mid / ENDO"), ("GSE", "UPN", "eF3", "b  GSE266265 eF3")]
    traits = ["KOLLER_EUR", "GCST90205183", "LEIOMYOMA", "HYPERTENSION", "T2D"]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.3), sharex=True)
    rows = []
    for ax, (atlas, fam, lab, title) in zip(axes, panels):
        g = C[(C["atlas"] == atlas) & (C["family"] == fam) & (C["label"] == lab)]
        est = {}
        for _, r in g.iterrows():
            est[r["trait_a"]] = (r["beta_std_a"], r["se_a"])
            est[r["trait_b"]] = (r["beta_std_b"], r["se_b"])
        yy = np.arange(len(traits))[::-1]
        for yi, t in zip(yy, traits):
            if t not in est:
                continue
            b, se = est[t]
            c = BLUE if t in ENDO_TRAITS else GRAY
            ax.plot([b - 1.96 * se, b + 1.96 * se], [yi, yi], color=c, linewidth=1.4, solid_capstyle="round", zorder=2)
            ax.scatter([b], [yi], s=26, color=c, edgecolor="white", linewidth=0.8, zorder=3)
            k = g[(g["trait_a"] == "KOLLER_EUR") & (g["trait_b"] == t)]
            if len(k):
                ax.text(0.078, yi, "p %s" % fmt_p(float(k["p_covariance"].iloc[0])), va="center", fontsize=6.2, color=INK2)
            rows.append({"panel": title, "trait": t, "beta_std": b, "se": se,
                         "p_cov_koller_vs": float(k["p_covariance"].iloc[0]) if len(k) else np.nan,
                         "p_rot_koller_vs": float(k["p_rot_empirical"].iloc[0]) if len(k) else np.nan})
        ax.axvline(0, color=AXIS, linewidth=0.6)
        ax.text(0.078, yy[2] + 0.55, "KOLLER_EUR > trait", fontsize=6, color=MUTED, va="bottom")
        ax.set_yticks(yy)
        ax.set_yticklabels([TRAIT_NAME[t] for t in traits], color=INK2)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xlabel("standardized enrichment (BETA_STD, 95% CI)")
        ax.set_xlim(-0.02, 0.10)
        ax.set_xticks([0, 0.02, 0.04, 0.06])
        style(ax)
    fig.subplots_adjust(wspace=0.5)
    D = pd.DataFrame(rows)
    leg = ("MAGMA gene-property enrichment of the size-matched up-only sets (N = 400) per GWAS; whiskers 95% CI from the MAGMA standard "
           "error. Endometriosis GWAS in blue, control traits in gray. Right column: one-sided p that KOLLER_EUR exceeds the control, with the "
           "covariance of the two estimates taken from 1,000 paired genome rotations (correlation -0.01 to 0.14); the empirical rotation "
           "p was smaller in every case.")
    save(fig, "S2_trait_contrasts", D, leg)


# ----------------------------------------------------------------------------
# S3 cross-atlas locus test
# ----------------------------------------------------------------------------
def read_gsa_p(f):
    rows = [l for l in open(f) if not l.startswith("#") and l.strip()]
    hdr = rows[0].split()
    out = {}
    for l in rows[1:]:
        d = dict(zip(hdr, l.split()))
        try:
            out[d.get("FULL_NAME", d.get("VARIABLE"))] = float(d["P"])
        except (KeyError, ValueError):
            continue
    return out


def s3():
    fc = need(RES / "review", "p9_v054_*/c_cross_atlas_removal.csv", "v054 C summary")
    S = pd.read_csv(fc).set_index("test")
    work = INTERIM / "magma_v054"
    panels = [("a", "c_heca_dsmid.gsa.out", "a  HECA dStromal mid / ENDO set (up-only)\n    minus loci stromal in GSE266265"),
              ("b", "c_gse_ef3.gsa.out", "b  GSE266265 eF3 set (up-only, N = 400)\n    minus loci stromal in HECA")]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.4), sharey=True)
    rows = []
    for ax, (key, fn, title) in zip(axes, panels):
        f = work / fn
        if not f.exists() or key not in S.index:
            ax.set_axis_off()
            say("  [S3] %s missing" % fn)
            continue
        p = read_gsa_p(f)
        base = p["BASE"]
        null = np.array([math.log10(v) - math.log10(base) for k, v in p.items() if re.match(r"^R\d+$", k)])
        obs = math.log10(p["DROPA"]) - math.log10(base)
        bins = np.linspace(min(null.min(), 0) - 0.05, max(null.max(), obs) + 0.1, 28)
        ax.hist(null, bins=bins, color=LIGHT, edgecolor="white", linewidth=0.6)
        ax.axvline(obs, color=ORANGE, linewidth=1.6)
        ax.text(obs - 0.03, ax.get_ylim()[1] * 0.92, "observed +%.2f\nperm p %s" % (obs, fmt_p(float(S.loc[key, "p_perm"]))), color=INK, fontsize=6.2,
                va="top", ha="right")
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xlabel("loss of significance after removal (Δ log10 p)")
        ax.grid(axis="y", color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)
        rows += [{"panel": key, "kind": "random", "delta": x} for x in null] + [{"panel": key, "kind": "observed", "delta": obs}]
        say("  [S3] %s: %d random subsets, observed %.2f, %d >= observed" % (key, len(null), obs, int((null >= obs).sum())))
    axes[0].set_ylabel("random locus subsets")
    fig.subplots_adjust(wspace=0.2)
    leg = ("Cross-atlas locus test (KOLLER_EUR, MAGMA, genes within 500 kb of each removed lead SNP dropped from the set). (a) The 8 loci "
           "classified as stromal with GSE266265 expression removed from the HECA up-only dStromal_mid/ENDO set; (b) the 9 loci classified "
           "as stromal with HECA expression removed from the GSE266265 size-matched up-only eF3 set. Gray: 200 random equal-size subsets "
           "of the other loci; orange: observed. Classification and enrichment come from different atlases (agreement 82%, kappa 0.77).")
    save(fig, "S3_cross_atlas_loci", pd.DataFrame(rows), leg)


# ----------------------------------------------------------------------------
# S4 benchmark intervals
# ----------------------------------------------------------------------------
def s4():
    f = need(RES / "review", "p9_v054_*/d_benchmark_rates_ci.csv", "v054 D intervals")
    W = pd.read_csv(f)
    W = W[W["metric"] == "top1_exact"]
    scen = [("POS_dStromal_mid", "a  decidual stroma (mid)"), ("POS_Venous", "b  venous endothelium"), ("POS_ePV_2", "c  perivascular ePV_2")]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.2), sharey=True)
    fams = [x for x in ("ABS", "ABSN", "UP", "UPN", "DOWN", "DOWNN", "SPEC")]
    yy = np.arange(len(fams))[::-1]
    for ax, (sc, title) in zip(axes, scen):
        g = W[W["scenario"] == sc].set_index("family")
        for yi, fa in zip(yy, fams):
            if fa not in g.index:
                continue
            r = g.loc[fa]
            c = ORANGE if fa == "ABS" else BLUE
            ax.plot([r["wilson_lo"], r["wilson_hi"]], [yi, yi], color=c, linewidth=1.4, solid_capstyle="round", zorder=2)
            ax.scatter([r["rate"]], [yi], s=24, color=c, edgecolor="white", linewidth=0.8, zorder=3)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xlim(-0.02, 1.02)
        ax.set_xlabel("top-1 exact recovery (95% CI)")
        style(ax)
    axes[0].set_yticks(yy)
    axes[0].set_yticklabels([FAM_NAME[x] for x in fams], color=INK2)
    fig.subplots_adjust(wspace=0.15)
    D = W[W["scenario"].isin([s for s, _ in scen])]
    leg = ("Known-source benchmark (v041 replicates, 40 per truth): share of replicates in which the planted state is ranked first, with "
           "Wilson 95% intervals, under specificity truths (genes specific to the named state). Koller et al.'s thresholded |log2FC| rule in orange. Up-only beat |log2FC| for "
           "all three truths (paired McNemar p < 0.001, < 0.001 and 0.039 for venous, ePV_2 and decidual stroma).")
    save(fig, "S4_benchmark_intervals", D, leg)


# ----------------------------------------------------------------------------
def collect():
    src = latest(RES / "paper_figures", "final_v053_*/figure_legends.md")
    legends = []
    if src is not None:
        for f in sorted(src.parent.glob("*")):
            if f.suffix in (".png", ".pdf", ".csv"):
                shutil.copy2(f, FINAL / f.name)
        legends.append(src.read_text().split("\n", 2)[-1])
        say("[COPY] F1-F8, S1 from %s" % src.parent.name)
    else:
        say("[MISS] final_v053 figure set not found")
    for f in sorted(OUT.glob("*")):
        if f.suffix in (".png", ".pdf", ".csv"):
            shutil.copy2(f, FINAL / f.name)
    legends.append("\n".join(LEG))
    (FINAL / "figure_legends.md").write_text("# Paper 9 figure legends (v056 collection)\n\n" + "\n".join(legends))
    (FINAL / "MANIFEST.txt").write_text("F1-F8, S1: %s\nF9, S2-S4: %s\n" % (src.parent.name if src else "-", OUT.name))


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.56 | audit figures (F9, S2-S4) and final figure set | " + TS)
    say("=" * 72)
    bad = []
    for name, fn in (("F9", f9), ("S2", s2), ("S3", s3), ("S4", s4)):
        try:
            fn()
        except Exception as e:
            import traceback
            bad.append(name)
            say("[FIG_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1200:])
            plt.close("all")
    collect()
    say("")
    say("[OUTPUT] " + str(FINAL))
    say("[LOG] " + str(LOG))
    say("[INCOMPLETE] figures with errors: " + ", ".join(bad) if bad else "[DONE] PAPER9_V056_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
