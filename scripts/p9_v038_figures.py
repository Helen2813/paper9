#!/usr/bin/env python3
# Paper 9 | v0.0.38 | draft paper figures for the Koller-pipeline / marker-rule results (no new analysis):
#   F1 rule flip (S-LDSC |log2FC| vs up-only vs down-only), F2 sign decomposition, F3 joint resolution,
#   F4 known-source simulations (recovery and detection), F5 composition of |log2FC| sets.
# Reads the newest results of v031, v034, v036, v037; writes PNG (300 dpi) + PDF + source-data CSV per figure.
# Runs in ~/Desktop/paper9/.venv; a few seconds, no heavy computation.

import datetime as dt
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RES = REPO / "results" / "koller_repro"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "paper_figures" / ("p9_v038_" + TS)
LOG = REPO / "logs" / ("p9_v038_figures_" + TS + ".log")
OUT.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")

TRAIT = "KOLLER_EUR"
# Colors (validated: categorical slots 1-3 pass all-pairs CVD / normal-vision checks on the light surface;
# slot 3 is below 3:1 contrast, so every figure carries a legend and direct labels)
COL = {"dStromal": "#2a78d6", "endothelial": "#eb6834", "perivascular/SM": "#1baf7a", "other": "#a3a29c"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SEQ = ["#f4f8fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
VASC = ["Arterial", "Venous", "ePV_1a", "ePV_1b", "ePV_2", "mPV", "uSMCs"]

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6, "legend.frameon": False,
                     "font.family": "DejaVu Sans", "pdf.fonttype": 42, "figure.facecolor": "white", "savefig.facecolor": "white"})


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def newest(pattern, fname):
    dirs = sorted(RES.glob(pattern))
    for d in reversed(dirs):
        if (d / fname).exists():
            return d / fname
    return None


def group_of(label):
    b = str(label).split("/")[0]
    if b.startswith("dStromal"):
        return "dStromal"
    if b in ("Arterial", "Venous"):
        return "endothelial"
    if b in VASC:
        return "perivascular/SM"
    return "other"


def short(label, n=24):
    s = str(label).replace("SOX9_functionalis_", "SOX9f_").replace("dStromal_", "dS_").replace("Glandular_secretory", "Gland_sec")
    return s[:n]


def nl10(p):
    return -np.log10(np.clip(np.asarray(p, dtype=float), 1e-300, 1))


def save(fig, name, data=None):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / ("%s.%s" % (name, ext)), dpi=300, bbox_inches="tight")
    plt.close(fig)
    if data is not None:
        data.to_csv(OUT / ("%s_source_data.csv" % name), index=False)
    say("[FIG] %s" % name)


def legend_groups(ax, loc="lower right"):
    h = [Line2D([0], [0], marker="o", color="none", markerfacecolor=COL[g], markeredgecolor=COL[g], markersize=5, label=g) for g in COL]
    ax.legend(handles=h, loc=loc, fontsize=7, handletextpad=0.3)


# ----------------------------------------------------------------------------
# F1 rule flip
# ----------------------------------------------------------------------------
def f1():
    fams = {}
    for mode in ("ABS", "UP", "DOWN"):
        f = newest("p9_v031_*", "sldsc_P_%s.csv" % mode)
        if f is None:
            say("[F1] sldsc_P_%s.csv not found; skipped" % mode)
            return
        d = pd.read_csv(f)
        fams[mode] = d[d["trait"] == TRAIT]
    rows = []
    for mode, d in fams.items():
        d = d.copy()
        d["base"] = d["label"].astype(str).str.split("/").str[0]
        b = d.sort_values("Coefficient_P_value").groupby("base", as_index=False).first()
        for _, r in b.iterrows():
            rows.append({"mode": mode, "cell_type": r["base"], "best_status": str(r["label"]).split("/")[-1], "p": r["Coefficient_P_value"], "q_bh": r["q_bh"]})
    D = pd.DataFrame(rows)
    order = D[D["mode"] == "ABS"].sort_values("p")["cell_type"].tolist()
    order += [c for c in D["cell_type"].unique() if c not in order]
    n = len(order)
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 0.16 * n + 1.2), sharey=True)
    titles = {"ABS": "|log2FC| > 1.2 (Koller rule)", "UP": "up-regulated only", "DOWN": "down-regulated only"}
    xmax = max(3.5, float(nl10(D["p"]).max()) * 1.08)
    for ax, mode in zip(axes, ("ABS", "UP", "DOWN")):
        d = D[D["mode"] == mode].set_index("cell_type").reindex(order)
        y = np.arange(n)[::-1]
        x = nl10(d["p"].values)
        g = [group_of(c) for c in order]
        ax.axvline(-math.log10(0.05), color=INK2, lw=0.6, ls=(0, (3, 3)), zorder=1)
        ax.hlines(y, 0, x, color=GRID, lw=1.0, zorder=1)
        for yi, xi, gi, qi in zip(y, x, g, d["q_bh"].values):
            filled = bool(qi < 0.05)
            ax.scatter([xi], [yi], s=22, color=COL[gi] if filled else "white", edgecolor=COL[gi], linewidth=1.2, zorder=3)
        ax.set_title(titles[mode], loc="left", color=INK)
        ax.set_xlim(0, xmax)
        ax.set_xlabel("-log10 p (S-LDSC, %s)" % TRAIT)
        ax.grid(axis="x", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
    axes[0].set_yticks(np.arange(n)[::-1])
    axes[0].set_yticklabels([short(c) for c in order])
    h = [Line2D([0], [0], marker="o", color="none", markerfacecolor=COL[g], markeredgecolor=COL[g], markersize=5, label=g) for g in COL]
    h += [Line2D([0], [0], marker="o", color="none", markerfacecolor=INK2, markeredgecolor=INK2, markersize=5, label="BH q < 0.05"),
          Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=INK2, markersize=5, label="q >= 0.05"),
          Line2D([0], [0], color=INK2, lw=0.6, ls=(0, (3, 3)), label="p = 0.05")]
    fig.suptitle("The prioritized cell type depends on the marker-gene selection rule (HECA, Koller pipeline; best status per type)", x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.legend(handles=h, loc="lower center", ncol=7, fontsize=7, handletextpad=0.3, columnspacing=1.2, bbox_to_anchor=(0.5, 0.0))
    save(fig, "F1_rule_flip", D)


# ----------------------------------------------------------------------------
# F2 sign decomposition
# ----------------------------------------------------------------------------
def f2():
    f = newest("p9_v034_*", "k7_sign_decomposition.csv")
    if f is None:
        say("[F2] k7 file not found; skipped")
        return
    d = pd.read_csv(f)
    d = d[(d["trait"] == TRAIT) & (d["class"] != "ABS not nominal")].copy()
    if not len(d):
        say("[F2] no ABS-nominal types; skipped")
        return
    d["x"] = nl10(d["p_up_given_down"])
    d["y"] = nl10(d["p_down_given_up"])
    d["group"] = d["type"].map(group_of)
    t = -math.log10(0.05)
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    lim = max(3.2, float(max(d["x"].max(), d["y"].max())) * 1.1)
    ax.axvline(t, color=INK2, lw=0.6, ls=(0, (3, 3)))
    ax.axhline(t, color=INK2, lw=0.6, ls=(0, (3, 3)))
    for gname, g in d.groupby("group"):
        ax.scatter(g["x"], g["y"], s=34, color=COL[gname], edgecolor="white", linewidth=0.8, zorder=3, label=gname)
    lab = d[(d["x"] > t) | (d["y"] > t)].sort_values("p_abs").head(12)
    for _, r in lab.iterrows():
        ax.annotate(short(r["type"], 18), (r["x"], r["y"]), xytext=(4, 2), textcoords="offset points", fontsize=6.3, color=INK2)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("-log10 p (UP set | DOWN set)")
    ax.set_ylabel("-log10 p (DOWN set | UP set)")
    cnt = d["class"].value_counts().to_dict()
    ax.set_title("%s: %d nominal |log2FC| annotations — %d down-carried, %d up-carried\nupper left: signal from genes depleted in the named type; lower right: from its up-regulated genes" % (
        TRAIT, len(d), cnt.get("DOWN-carried", 0), cnt.get("UP-carried", 0)), loc="left", fontsize=7.5)
    ax.legend(loc="upper right", fontsize=7, handletextpad=0.2)
    ax.grid(color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    fig.tight_layout()
    save(fig, "F2_sign_decomposition", d)


# ----------------------------------------------------------------------------
# F3 joint resolution heatmaps
# ----------------------------------------------------------------------------
def heat(ax, M, title, rows, cols):
    from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
    cmap = LinearSegmentedColormap.from_list("seq", SEQ)
    bounds = [0, -math.log10(0.5), -math.log10(0.2), -math.log10(0.1), -math.log10(0.05), -math.log10(0.01), -math.log10(0.001), 6]
    norm = BoundaryNorm(bounds, cmap.N)
    im = ax.imshow(nl10(M), cmap=cmap, norm=norm, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if np.isfinite(v):
                txt = "<.001" if v < 0.001 else (("%.3f" % v)[1:] if v < 0.01 else ("%.2f" % v)[1:] if v < 1 else "1")
                ax.text(j, i, txt, ha="center", va="center", fontsize=5.4, color="white" if v < 0.01 else INK)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([short(c, 14) for c in cols], rotation=45, ha="right")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([short(r, 16) for r in rows])
    ax.set_title(title, loc="left", fontsize=8)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    return im


def f3():
    fa = newest("p9_v031_*", "k4_joint_P_ABS.csv")
    fu = newest("p9_v034_*", "k8_dstromal_vs_vascular_up.csv")
    if fa is None or fu is None:
        say("[F3] joint files not found; skipped")
        return
    A = pd.read_csv(fa)
    A = A[A["trait"] == TRAIT]
    U = pd.read_csv(fu)
    U = U[U["trait"] == TRAIT]
    ra, ca = list(dict.fromkeys(A["focus"])), list(dict.fromkeys(A["other"]))
    ru, cu = list(dict.fromkeys(U["dstromal"])), list(dict.fromkeys(U["other"]))
    pa1 = A.pivot(index="focus", columns="other", values="p_focus_given_other").reindex(index=ra, columns=ca).values
    pa2 = A.pivot(index="focus", columns="other", values="p_other_given_focus").reindex(index=ra, columns=ca).values
    pu1 = U.pivot(index="dstromal", columns="other", values="p_ds_given_other").reindex(index=ru, columns=cu).values
    pu2 = U.pivot(index="dstromal", columns="other", values="p_other_given_ds").reindex(index=ru, columns=cu).values
    fig, axes = plt.subplots(2, 2, figsize=(7.6, 6.4), layout="constrained",
                             gridspec_kw={"height_ratios": [len(ra), max(len(ru), 4)], "hspace": 0.12})
    im = heat(axes[0, 0], pa1, "a  |log2FC| sets: p(focus | other)", ra, ca)
    heat(axes[0, 1], pa2, "b  |log2FC| sets: p(other | focus)", ra, ca)
    heat(axes[1, 0], pu1, "c  up-only sets: p(decidual | vascular)", ru, cu)
    heat(axes[1, 1], pu2, "d  up-only sets: p(vascular | decidual)", ru, cu)
    for ax in (axes[0, 1], axes[1, 1]):
        ax.set_yticklabels([])
    pt = [0.5, 0.2, 0.1, 0.05, 0.01, 0.001]
    cb = fig.colorbar(im, ax=axes, shrink=0.5, pad=0.02, ticks=[-math.log10(x) for x in pt])
    cb.ax.set_yticklabels([str(x) for x in pt])
    cb.set_label("p (joint S-LDSC; darker = survives)")
    cb.outline.set_visible(False)
    fig.suptitle("%s joint S-LDSC: which annotation survives conditioning on the other" % TRAIT, x=0.01, ha="left", fontsize=9)
    src = pd.concat([A.assign(panel="ab"), U.rename(columns={"dstromal": "focus", "p_ds": "p_focus", "p_ds_given_other": "p_focus_given_other",
                                                             "p_other_given_ds": "p_other_given_focus"}).assign(panel="cd")], ignore_index=True)
    save(fig, "F3_joint_resolution", src)


# ----------------------------------------------------------------------------
# F4 known-source simulations
# ----------------------------------------------------------------------------
def f4():
    f36 = newest("p9_v036_*", "s_abs_rotation.csv")
    if f36 is None:
        say("[F4] v036 simulation file not found; skipped")
        return
    R = pd.read_csv(f36)
    truths = [("dStromal_mid", 2.0, "dStromal", "decidual stroma x2"), ("dStromal_mid", 4.0, "dStromal", "decidual stroma x4"),
              ("Venous", 2.0, "endothelial", "venous endothelium x2"), ("ePV_2", 2.0, "perivascular/SM", "perivascular ePV_2 x2")]
    rows = []
    for name, st, grp, lab in truths:
        for mode in ("ABS", "UP"):
            g = R[(R["truth"] == name) & (R["strength"] == st) & (R["mode"] == mode)]
            n0 = R[(R["truth"] == "null") & (R["mode"] == mode)]
            if not len(g):
                continue
            rows.append({"truth": lab, "group": grp, "rule": mode, "n": len(g), "top_correct": float((g["top_group"] == grp).mean()),
                         "top_correct_null": float((n0["top_group"] == grp).mean()) if len(n0) else np.nan,
                         "detected_any_q05": float(g["any_q05"].mean()), "detected_null": float(n0["any_q05"].mean()) if len(n0) else np.nan})
    D = pd.DataFrame(rows)
    if not len(D):
        say("[F4] no scenarios; skipped")
        return
    labs = list(dict.fromkeys(D["truth"]))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
    w = 0.36
    for ax, col, ncol, title in ((axes[0], "top_correct", "top_correct_null", "a  true cell group ranked first"),
                                 (axes[1], "detected_any_q05", "detected_null", "b  any annotation BH q < 0.05")):
        for i, lab in enumerate(labs):
            for k, mode in enumerate(("ABS", "UP")):
                r = D[(D["truth"] == lab) & (D["rule"] == mode)]
                if not len(r):
                    continue
                r = r.iloc[0]
                x = i + (k - 0.5) * (w + 0.04)
                c = COL[r["group"]]
                ax.bar(x, r[col], width=w, color=c if mode == "UP" else "white", edgecolor=c, linewidth=1.2,
                       hatch=None if mode == "UP" else "////", zorder=2)
                ax.hlines(r[ncol], x - w / 2, x + w / 2, color=INK, lw=1.0, zorder=3)
                ax.text(x, r[col] + 0.02, "%.0f" % (100 * r[col]), ha="center", va="bottom", fontsize=6.3, color=INK2)
        ax.set_xticks(range(len(labs)))
        ax.set_xticklabels(labs, rotation=20, ha="right")
        ax.set_ylim(0, 1.08)
        ax.set_title(title, loc="left", fontsize=8)
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("share of 60 replicates")
    h = [Patch(facecolor="white", edgecolor=INK2, hatch="////", label="|log2FC| rule"), Patch(facecolor=INK2, edgecolor=INK2, label="up-only rule"),
         Line2D([0], [0], color=INK, lw=1.0, label="same quantity under the null")]
    axes[1].legend(handles=h, loc="upper right", fontsize=6.5)
    fig.suptitle("Known-source simulations (rotation design; real %s z-scores)" % TRAIT, x=0.01, ha="left", fontsize=9)
    fig.tight_layout()
    save(fig, "F4_simulation_recovery", D)


# ----------------------------------------------------------------------------
# F5 composition of |log2FC| sets
# ----------------------------------------------------------------------------
def f5():
    f = newest("p9_v031_*", "k2_abs_composition.csv")
    if f is None:
        say("[F5] composition file not found; skipped")
        return
    d = pd.read_csv(f)
    d["group"] = d["label"].map(group_of)
    d["pct"] = 100 * d["down_share"]
    order = d.groupby("label")["pct"].mean().sort_values().index.tolist()
    n = len(order)
    fig, ax = plt.subplots(figsize=(4.6, 0.17 * n + 1.3))
    ypos = {l: i for i, l in enumerate(order)}
    for status, mk in (("CTL", "o"), ("ENDO", "^")):
        s_ = d[d["status"] == status]
        for _, r in s_.iterrows():
            ax.scatter(r["pct"], ypos[r["label"]], marker=mk, s=26, color=COL[r["group"]], edgecolor="white", linewidth=0.6, zorder=3)
    for l in order:
        v = d.loc[d["label"] == l, "pct"]
        ax.hlines(ypos[l], v.min(), v.max(), color=GRID, lw=1.2, zorder=1)
    med = float(d["pct"].median())
    ax.axvline(50, color=INK2, lw=0.6, ls=(0, (3, 3)))
    ax.axvline(med, color=INK, lw=0.8)
    ax.text(med + 1, n - 0.4, "median %.0f%%" % med, fontsize=6.5, color=INK, va="top")
    ax.set_yticks(range(n))
    ax.set_yticklabels([short(l) for l in order])
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlabel("% of genes in the |log2FC| set depleted in the named cell type")
    ax.set_title("|log2FC| marker sets are mostly depleted genes", loc="left", fontsize=8)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    h = [Line2D([0], [0], marker="o", color="none", markerfacecolor=COL[g], markeredgecolor=COL[g], markersize=5, label=g) for g in COL]
    h += [Line2D([0], [0], marker="o", color="none", markerfacecolor=INK2, markersize=5, label="control donors"),
          Line2D([0], [0], marker="^", color="none", markerfacecolor=INK2, markersize=5, label="endometriosis donors")]
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.legend(handles=h, loc="lower center", ncol=3, fontsize=6.5, handletextpad=0.3, columnspacing=1.0, bbox_to_anchor=(0.5, 0.0))
    save(fig, "F5_abs_composition", d)


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.38 | draft figures (Koller pipeline, marker rules, simulations) | " + TS)
    say("=" * 72)
    for name, fn in (("F1", f1), ("F2", f2), ("F3", f3), ("F4", f4), ("F5", f5)):
        try:
            fn()
        except Exception as e:
            import traceback
            say("[FIG_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1200:])
    say("")
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V038_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
