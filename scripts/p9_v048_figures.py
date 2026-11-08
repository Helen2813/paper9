#!/usr/bin/env python3
# Paper 9 | v0.0.48 | manuscript figures from finished runs (no new analyses)
#   F1 rule flip in Koller's exact cohort (S-LDSC, v047)
#   F2 mechanism: sign decomposition (v047 K7), composition (v047 K2), masking of the true source (v047, v046)
#   F3 joint resolution, decidual / SOX9-functionalis sets given the top non-focus sets (v047 K4)
#   F4 known-source benchmark, top-1 exact recovery (v041)
#   F5 positive controls across frameworks and atlases: rank of the expected group (v047, v043, v044)
#   F6 stromal loci vs the first-ranked sets (v044 part A)
#   S1 conditioning on generic gene properties, first-ranked label per family (v045 C4)
# Writes PNG (300 dpi), PDF and the plotted values (CSV) per figure plus figure_legends.md. ~1 min, no network.

import datetime as dt
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RES = REPO / "results"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = RES / "paper_figures" / ("p9_v048_" + TS)
LOG = REPO / "logs" / ("p9_v048_figures_" + TS + ".log")
OUT.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")

# Palette (validated: slots 1-3 all-pairs, light surface; aqua < 3:1 -> every mark is also named on an axis)
BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#898781"
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95",
       "#104281", "#0d366b"]
GROUP_COLOR = {"endothelial / perivascular / smooth muscle": BLUE, "stromal": ORANGE, "epithelial": AQUA, "other": GRAY}
FAMS = ["ABS", "UP", "DOWN", "ABSN", "UPN", "DOWNN", "SPEC"]
FAM_NAME = {"ABS": "|log2FC|", "UP": "up-only", "DOWN": "down-only", "ABSN": "|log2FC|, N=400", "UPN": "up-only, N=400",
            "DOWNN": "down-only, N=400", "SPEC": "specificity"}
LEIO_RE = r"^(uSMCs|mPV)/"
VASC_RE = r"^(Arterial|Venous|mPV|ePV_1a|ePV_1b|ePV_2|uSMCs)/"
EXPECT = {"HECA": {"LEIOMYOMA": r"^(uSMCs|mPV|ePV)", "HYPERTENSION": r"^(Arterial|Venous|ePV|mPV|uSMCs)"},
          "GTEx": {"LEIOMYOMA": r"smooth muscle", "HYPERTENSION": r"smooth muscle|endothel|pericyte|vascular"}}
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
    say("  input %-28s %s" % (what, f.relative_to(RES)))
    return f


def nlp(p):
    p = np.asarray(p, dtype=float)
    return -np.log10(np.clip(p, 1e-300, 1.0))


def group_of(label):
    b = str(label).split("/")[0]
    if re.search(r"^(Arterial|Venous|Lymphatic|ePV|mPV|uSMC)|endothel|pericyte|smooth muscle|contractile", b, re.I):
        return "endothelial / perivascular / smooth muscle"
    if re.search(r"^(dStromal|eStromal|sHormones|dHormones|Fibroblast)|fibroblast|stromal", b, re.I):
        return "stromal"
    if re.search(r"^(SOX9|Luminal|preLuminal|Ciliated|preCiliated|Glandular|Gland|eHormones|MUC5B|KRT5)|epitheli", b, re.I):
        return "epithelial"
    return "other"


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


def cell_gaps(ax, nrow, ncol):
    for i in range(1, nrow):
        ax.axhline(i - 0.5, color="white", linewidth=1.5)
    for j in range(1, ncol):
        ax.axvline(j - 0.5, color="white", linewidth=1.5)


def group_legend(fig, y=1.0, extra=None):
    h = [Line2D([], [], marker="o", linestyle="", markersize=5, markerfacecolor=c, markeredgecolor=c, label=g) for g, c in GROUP_COLOR.items()]
    if extra:
        h += extra
    fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, y), ncol=len(h), frameon=False, handletextpad=0.3, columnspacing=1.2)


def v047_dir():
    f = need(RES / "koller_repro", "p9_v047_*/sldsc_H_ABS.csv", "v047 S-LDSC (exact cohort)")
    return f.parent


# ----------------------------------------------------------------------------
def f1():
    d = v047_dir()
    R = {m: pd.read_csv(d / ("sldsc_H_%s.csv" % m)) for m in ("ABS", "UP", "DOWN") if (d / ("sldsc_H_%s.csv" % m)).exists()}
    fig, axs = plt.subplots(1, 3, figsize=(7.1, 3.0), sharex=True)
    rows, xmax = [], 0
    for ax, m in zip(axs, ("ABS", "UP", "DOWN")):
        if m not in R:
            ax.set_visible(False)
            continue
        k = R[m][R[m]["trait"] == "KOLLER_EUR"].sort_values("Coefficient_P_value").head(12).reset_index(drop=True)
        y = np.arange(len(k))[::-1]
        x = nlp(k["Coefficient_P_value"])
        xmax = max(xmax, float(x.max()))
        for yi, xi, lab, q in zip(y, x, k["label"], k["q_bh"]):
            c = GROUP_COLOR[group_of(lab)]
            ax.plot([0, xi], [yi, yi], color=GRID, linewidth=0.8, zorder=1)
            ax.scatter([xi], [yi], s=26, color=c if q < 0.05 else "white", edgecolor=c, linewidth=1.2, zorder=3)
            rows.append({"rule": m, "label": lab, "p": 10 ** -xi, "q_bh": q, "group": group_of(lab)})
        ax.set_yticks(y)
        ax.set_yticklabels([nice(l, 24) for l in k["label"]], color=INK2)
        ax.axvline(-math.log10(0.05), color=AXIS, linestyle=(0, (2, 2)), linewidth=0.6)
        ax.set_title({"ABS": "|log2FC| (Koller rule)", "UP": "Up-only", "DOWN": "Down-only"}[m], loc="left", fontweight="bold")
        ax.set_xlabel("−log10 p (S-LDSC)")
        style(ax)
    axs[0].set_xlim(0, xmax * 1.08 + 0.1)
    group_legend(fig, 1.0, [Line2D([], [], marker="o", linestyle="", markersize=5, markerfacecolor="white", markeredgecolor=INK2, label="BH q ≥ 0.05 (open)")])
    fig.tight_layout(w_pad=2.2)
    D = pd.DataFrame(rows)
    top = {m: D[D["rule"] == m].iloc[0] for m in ("ABS", "UP", "DOWN") if (D["rule"] == m).any()}
    leg = ("Koller et al.'s exact cohort (22 endometriosis / 13 control library-individual pairs, hormone-treated samples included) "
           "analysed with their S-LDSC pipeline (KOLLER_EUR GWAS). The 12 most significant annotations per marker-selection rule; "
           "filled = BH q < 0.05 within the rule; dashed line p = 0.05. First-ranked annotation: " +
           "; ".join("%s %s (p %.1e, q %.2f)" % ({"ABS": "|log2FC|", "UP": "up-only", "DOWN": "down-only"}[m], r["label"], r["p"], r["q_bh"]) for m, r in top.items()) + ".")
    save(fig, "F1_rule_flip_exact_cohort", D, leg)


# ----------------------------------------------------------------------------
def f2():
    d = v047_dir()
    K = pd.read_csv(d / "k7_sign_decomposition_H.csv")
    K = K[K["trait"] == "KOLLER_EUR"].copy()
    C = pd.read_csv(d / "k2_abs_composition.csv") if (d / "k2_abs_composition.csv").exists() else None
    A = pd.read_csv(d / "sldsc_H_ABS.csv")
    fig = plt.figure(figsize=(7.1, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 0.8], width_ratios=[1.0, 1.0], hspace=0.5, wspace=0.95, top=0.9)
    # (a) sign decomposition
    ax = fig.add_subplot(gs[0, 0])
    x, y = nlp(K["p_up_given_down"]), nlp(K["p_down_given_up"])
    nom = (K["p_abs"] < 0.05).values
    for xi, yi, lab, nm in zip(x, y, K["type"], nom):
        c = GROUP_COLOR[group_of(lab)]
        ax.scatter([xi], [yi], s=22, color=c if nm else "white", edgecolor=c, linewidth=1.1, zorder=3)
    t = -math.log10(0.05)
    ax.axvline(t, color=AXIS, linestyle=(0, (2, 2)), linewidth=0.6)
    ax.axhline(t, color=AXIS, linestyle=(0, (2, 2)), linewidth=0.6)
    lim = max(float(np.nanmax(x)), float(np.nanmax(y))) * 1.1 + 0.2
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.text(0.03, 0.97, "DOWN-carried", transform=ax.transAxes, va="top", color=MUTED, fontsize=6)
    ax.text(0.97, 0.97, "UP-carried →", transform=ax.transAxes, ha="right", va="top", color=MUTED, fontsize=6)
    lab_set = list(K.sort_values("p_abs")["type"].head(3)) + [l for l in K["type"] if re.match(r"^dStromal_(mid|early)/ENDO$", str(l))]
    for lab in dict.fromkeys(lab_set):
        r = K[K["type"] == lab].iloc[0]
        px, py = nlp(r["p_up_given_down"]), nlp(r["p_down_given_up"])
        right = px > 0.6 * lim
        ax.annotate(nice(lab, 20), (px, py), xytext=(-4 if right else 4, 3), textcoords="offset points", ha="right" if right else "left",
                    fontsize=5.8, color=INK2)
    ax.set_xlabel("−log10 p, up part | down part")
    ax.set_ylabel("−log10 p, down part | up part")
    ax.set_title("a  Sign decomposition of |log2FC| annotations", loc="left", fontweight="bold")
    ax.grid(color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    cnt = K.loc[K["class"] != "ABS not nominal", "class"].value_counts().to_dict()
    # (b) composition
    ax = fig.add_subplot(gs[0, 1])
    rows_c = []
    if C is not None:
        k = A[A["trait"] == "KOLLER_EUR"].sort_values("Coefficient_P_value")
        keys = list(k["Name"].head(8)) + [n for n, l in zip(k["Name"], k["label"]) if re.match(r"^dStromal_(mid|early)/ENDO$", str(l))]
        keys = list(dict.fromkeys(keys))
        Ci = C.set_index("set_key")
        lab = dict(zip(k["Name"], k["label"]))
        vals = [(lab[n], float(Ci.loc[n, "down_share"]), float(k.set_index("Name").loc[n, "Coefficient_P_value"])) for n in keys if n in Ci.index]
        yy = np.arange(len(vals))[::-1]
        for yi, (l, s, p) in zip(yy, vals):
            ax.barh(yi, s * 100, height=0.62, color=GROUP_COLOR[group_of(l)], edgecolor="white", linewidth=0.8)
            rows_c.append({"label": l, "depleted_share": s, "p_abs": p})
        ax.set_yticks(yy)
        ax.set_yticklabels(["%s (p %.0e)" % (nice(l, 20), p) for l, s, p in vals], color=INK2)
        ax.axvline(50, color=AXIS, linestyle=(0, (2, 2)), linewidth=0.6)
        ax.set_xlim(0, 100)
        ax.set_xlabel("genes depleted in the cell type (% of the |log2FC| set)")
        style(ax)
    ax.set_title("b  Composition of |log2FC| sets", loc="left", fontweight="bold")
    # (c) masking
    ax = fig.add_subplot(gs[1, :])
    rows_m = []
    r = K[K["type"].astype(str).str.match(r"^dStromal_mid/ENDO$")]
    if len(r):
        r = r.iloc[0]
        rows_m.append(("Endometriosis · decidual stroma (mid, ENDO) · S-LDSC, exact cohort", r["p_abs"], r["p_up"]))
    f46 = latest(RES / "specificity", "p9_v046_*/d_sign_decomposition.csv")
    if f46 is not None:
        say("  input %-28s %s" % ("v046 sign decomposition", f46.relative_to(RES)))
        M = pd.read_csv(f46)
        for atl, tr, lab, txt in (("HECA", "KOLLER_EUR", "dStromal_mid/ENDO", "Endometriosis · decidual stroma (mid, ENDO) · MAGMA, HECA"),
                                  ("HECA", "LEIOMYOMA", "uSMCs/CTL", "Leiomyoma · uterine smooth muscle (CTL) · MAGMA, HECA")):
            g = M[(M["atlas"] == atl) & (M["trait"] == tr) & (M["label"] == lab)]
            if len(g):
                rows_m.append((txt, g.iloc[0]["p_abs"], g.iloc[0]["p_up"]))
    yy = np.arange(len(rows_m))[::-1]
    for yi, (txt, pa, pu) in zip(yy, rows_m):
        xa, xu = nlp(pa), nlp(pu)
        ax.plot([xa, xu], [yi, yi], color=AXIS, linewidth=1.2, zorder=1)
        ax.scatter([xa], [yi], s=30, color="white", edgecolor=INK2, linewidth=1.2, zorder=3)
        ax.scatter([xu], [yi], s=30, color=INK, edgecolor=INK, zorder=3)
        ax.text(xu + 0.12, yi, "p %.0e" % pu, va="center", fontsize=6, color=INK2)
        ax.text(xa - 0.12, yi, "p %.2g" % pa, va="center", ha="right", fontsize=6, color=INK2)
        ax.text(0.0, yi + 0.36, txt, transform=ax.get_yaxis_transform(), fontsize=6.3, color=INK2, va="bottom")
    ax.set_yticks([])
    ax.set_ylim(-0.5, len(rows_m) - 0.1)
    ax.axvline(-math.log10(0.05), color=AXIS, linestyle=(0, (2, 2)), linewidth=0.6)
    ax.set_xlim(-0.6, max([nlp(b) for t, a, b in rows_m] + [3]) + 1.0)
    ax.set_xlabel("−log10 p")
    ax.legend(handles=[Line2D([], [], marker="o", linestyle="", markersize=5, markerfacecolor="white", markeredgecolor=INK2, label="|log2FC| set"),
                       Line2D([], [], marker="o", linestyle="", markersize=5, markerfacecolor=INK, markeredgecolor=INK, label="its up part")],
              loc="lower right", bbox_to_anchor=(1.0, 1.0), frameon=False, ncol=2)
    ax.set_title("c  Masking of the true source by its own depleted genes", loc="left", fontweight="bold")
    style(ax)
    group_legend(fig, 0.935, [Line2D([], [], marker="o", linestyle="", markersize=5, markerfacecolor="white", markeredgecolor=INK2, label="|log2FC| p ≥ 0.05 (open)")])
    D = pd.concat([K.assign(panel="a"), pd.DataFrame(rows_c).assign(panel="b"),
                   pd.DataFrame(rows_m, columns=["row", "p_abs", "p_up"]).assign(panel="c")], ignore_index=True)
    leg = ("Exact Koller cohort, KOLLER_EUR. (a) Joint S-LDSC of each cell type's up part and down part of its |log2FC| set "
           "(K7 rule; dashed lines p = 0.05); filled = |log2FC| set nominally significant. Classes among the %d nominal annotations: %s. "
           "(b) Share of depleted genes in the |log2FC| sets of the top-ranked annotations and of decidual stroma (dashed line 50%%). "
           "(c) |log2FC| set vs its up part for the true source of a trait: decidual stroma for endometriosis (S-LDSC exact cohort; MAGMA HECA) "
           "and uterine smooth muscle for leiomyoma (MAGMA HECA). In GTEx the smooth-muscle |log2FC| set was not masked in this way "
           "(up part p 0.021, see Results)." % (sum(cnt.values()), ", ".join("%s %d" % (k2, v) for k2, v in sorted(cnt.items()))))
    save(fig, "F2_mechanism_sign_composition_masking", D, leg)


# ----------------------------------------------------------------------------
def f3():
    d = v047_dir()
    fig, axs = plt.subplots(1, 2, figsize=(7.1, 3.1), gridspec_kw={"wspace": 0.85})
    cmap = ListedColormap(SEQ[1:])
    rows, legs = [], []
    for ax, tag, title in zip(axs, ("H_ABS", "H_UP"), ("a  |log2FC| sets", "b  Up-only sets")):
        f = d / ("k4_joint_%s.csv" % tag)
        if not f.exists():
            ax.set_visible(False)
            continue
        J = pd.read_csv(f)
        J = J[J["trait"] == "KOLLER_EUR"]
        foc = J.groupby("focus")["p_focus"].first().sort_values().index.tolist()
        oth = list(dict.fromkeys(J["other"]))
        M = J.pivot_table(index="focus", columns="other", values="p_focus_given_other").reindex(index=foc, columns=oth)
        V = np.clip(nlp(M.values), 0, 3)
        ax.imshow(V, cmap=cmap, vmin=0, vmax=3, aspect="auto")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                p = M.values[i, j]
                if np.isfinite(p):
                    ax.text(j, i, "%.2f" % p if p >= 0.01 else "%.0e" % p, ha="center", va="center", fontsize=5.6,
                            color="white" if V[i, j] > 1.6 else INK, fontweight="bold" if p < 0.05 else "normal")
        ax.set_xticks(range(len(oth)))
        ax.set_xticklabels([nice(o, 20) for o in oth], rotation=40, ha="right", color=INK2)
        ax.set_yticks(range(len(foc)))
        ax.set_yticklabels([nice(l, 24) for l in foc], color=INK2)
        cell_gaps(ax, M.shape[0], M.shape[1])
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(title + ": p(focus | other)", loc="left", fontweight="bold")
        ok_f = [f_ for f_, g in J.groupby("focus") if (g["p_focus_given_other"] < 0.05).all()]
        ok_o = [o for o, g in J.groupby("other") if (g["p_other_given_focus"] < 0.05).all()]
        legs.append("%s: focus sets surviving every non-focus set: %s; non-focus sets surviving every focus set: %s" % (
            title[3:], ", ".join(ok_f) or "none", ", ".join(ok_o) or "none"))
        rows.append(J.assign(family=tag))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 3))
    cb = fig.colorbar(sm, ax=axs, fraction=0.025, pad=0.02)
    cb.set_label("−log10 p (capped at 3)", color=INK2)
    cb.outline.set_visible(False)
    leg = ("Joint S-LDSC in the exact Koller cohort (KOLLER_EUR): each decidual-stroma / SOX9-functionalis set (rows) fitted together "
           "with each of the six top-ranked non-focus sets of the same rule (columns); cell = p of the focus set given the other "
           "(bold: p < 0.05). " + "; ".join(legs) + ".")
    save(fig, "F3_joint_resolution_exact_cohort", pd.concat(rows, ignore_index=True) if rows else None, leg)


# ----------------------------------------------------------------------------
def f4():
    f = need(RES / "benchmark", "p9_v041_*/recovery_table.csv", "v041 recovery table")
    T = pd.read_csv(f)
    kinds = [k for k in ("POS", "NEG", "EXPR") if k in set(T["kind"])]
    rows_ = [(k, t) for k in kinds for t in T.loc[T["kind"] == k, "target"].drop_duplicates()]
    fams = [x for x in ("ABS", "ABSN", "UP", "UPN", "DOWN", "DOWNN", "SPEC") if x in set(T["family"])]
    M = np.full((len(rows_), len(fams)), np.nan)
    for i, (k, t) in enumerate(rows_):
        for j, fa in enumerate(fams):
            v = T[(T["kind"] == k) & (T["target"] == t) & (T["family"] == fa)]["top1_exact"]
            if len(v):
                M[i, j] = float(v.iloc[0])
    fig, ax = plt.subplots(figsize=(4.6, 0.32 * len(rows_) + 1.1))
    cmap = LinearSegmentedColormap.from_list("seq", SEQ)
    ax.imshow(M, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isfinite(M[i, j]):
                ax.text(j, i, "%.2f" % M[i, j], ha="center", va="center", fontsize=6, color="white" if M[i, j] > 0.55 else INK)
    kn = {"POS": "specific to", "NEG": "depleted in", "EXPR": "expressed in"}
    ax.set_yticks(range(len(rows_)))
    ax.set_yticklabels(["%s %s" % (kn[k], nice(t, 20)) for k, t in rows_], color=INK2)
    ax.set_xticks(range(len(fams)))
    ax.set_xticklabels([FAM_NAME[x] for x in fams], rotation=35, ha="right", color=INK2)
    ax.tick_params(length=0)
    cell_gaps(ax, M.shape[0], M.shape[1])
    for s in ax.spines.values():
        s.set_visible(False)
    for i in range(1, len(kinds)):
        ax.axhline(sum(1 for k, t in rows_ if k in kinds[:i]) - 0.5, color="white", linewidth=2)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("top-1 exact recovery", color=INK2)
    cb.outline.set_visible(False)
    ax.set_title("Known-source benchmark (rotation design)", loc="left", fontweight="bold")
    D = pd.DataFrame(M, columns=fams).assign(kind=[k for k, t in rows_], target=[t for k, t in rows_])
    best = [fams[int(np.nanargmax(M[i]))] for i in range(M.shape[0]) if np.isfinite(M[i]).any()]
    leg = ("Share of replicates in which the planted source cell state is ranked first (MAGMA, HECA, KOLLER_EUR z-scores with rotated "
           "annotation layer) for seven marker-selection families. Truths: genes specific to, depleted in, or highly expressed in the "
           "target state. Best family per truth: %s. The thresholded |log2FC| rule is never best." % ", ".join(
               "%s %s -> %s" % (kn[k], t, b) for (k, t), b in zip(rows_, best)))
    save(fig, "F4_benchmark_recovery", D, leg)


# ----------------------------------------------------------------------------
def rank_expected(D, trait, fam, rx, pcol, qcol, lab="label", extra=None):
    d = D[(D["trait"] == trait) & (D["family"] == fam)].sort_values(pcol).reset_index(drop=True)
    if not len(d):
        return None, None
    m = d[lab].astype(str).str.contains(rx, regex=True, flags=re.I)
    if extra:
        m = m | d[lab].astype(str).str.match(extra)
    hit = d.index[m]
    if not len(hit):
        return None, None
    return int(hit[0]) + 1, float(d.loc[hit[0], qcol])


def f5():
    rows = []
    d = v047_dir()
    S = []
    for m in ("ABS", "UP", "DOWN"):
        f = d / ("sldsc_H_%s.csv" % m)
        if f.exists():
            S.append(pd.read_csv(f).assign(family=m))
    if S:
        S = pd.concat(S)
        for t, rx in (("LEIOMYOMA", LEIO_RE), ("HYPERTENSION", VASC_RE)):
            rows.append(("S-LDSC · HECA exact cohort · " + t.lower(), {fa: rank_expected(S, t, fa, rx, "Coefficient_P_value", "q_bh") for fa in ("ABS", "UP", "DOWN")}))
    f43 = need(RES / "benchmark", "p9_v043_*/e1_real_family_rankings.csv", "v043 HECA rankings")
    H = pd.read_csv(f43)
    for t in ("LEIOMYOMA", "HYPERTENSION"):
        rows.append(("MAGMA · HECA · " + t.lower(), {fa: rank_expected(H, t, fa, EXPECT["HECA"][t], "p", "q") for fa in FAMS}))
    f44 = need(RES / "external_atlas", "p9_v044_*/b_gtex_family_rankings.csv", "v044 GTEx rankings")
    G = pd.read_csv(f44)
    rows.append(("MAGMA · GTEx v9 · leiomyoma", {fa: rank_expected(G, "LEIOMYOMA", fa, EXPECT["GTEx"]["LEIOMYOMA"], "p", "q") for fa in FAMS}))
    rows.append(("MAGMA · GTEx v9 · hypertension (pre-registered labels)", {fa: rank_expected(G, "HYPERTENSION", fa, EXPECT["GTEx"]["HYPERTENSION"], "p", "q") for fa in FAMS}))
    rows.append(("MAGMA · GTEx v9 · hypertension (Pericyte/SMC label, post hoc)",
                 {fa: rank_expected(G, "HYPERTENSION", fa, EXPECT["GTEx"]["HYPERTENSION"], "p", "q", extra=r"^contractile cell$") for fa in FAMS}))
    R = np.full((len(rows), len(FAMS)), np.nan)
    Q = np.full((len(rows), len(FAMS)), np.nan)
    for i, (name, v) in enumerate(rows):
        for j, fa in enumerate(FAMS):
            rk, q = v.get(fa, (None, None)) or (None, None)
            if rk is not None:
                R[i, j], Q[i, j] = rk, q
    cat = np.where(np.isnan(R), np.nan, np.where(R == 1, 3, np.where(R <= 3, 2, np.where(R <= 10, 1, 0))))
    cmap = ListedColormap(["#f0efec", SEQ[3], SEQ[7], SEQ[11]])
    fig, ax = plt.subplots(figsize=(6.2, 0.36 * len(rows) + 1.2))
    ax.imshow(np.ma.masked_invalid(cat), cmap=cmap, norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], 4), aspect="auto")
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            if np.isfinite(R[i, j]):
                ax.text(j, i, "%d%s" % (R[i, j], "*" if Q[i, j] < 0.05 else ""), ha="center", va="center", fontsize=6.3,
                        color="white" if cat[i, j] >= 2 else INK, fontweight="bold" if R[i, j] == 1 else "normal")
            else:
                ax.text(j, i, "–", ha="center", va="center", fontsize=6.3, color=MUTED)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([n for n, v in rows], color=INK2)
    ax.set_xticks(range(len(FAMS)))
    ax.set_xticklabels([FAM_NAME[x] for x in FAMS], rotation=35, ha="right", color=INK2)
    ax.tick_params(length=0)
    cell_gaps(ax, R.shape[0], R.shape[1])
    for s in ax.spines.values():
        s.set_visible(False)
    hs = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (SEQ[11], SEQ[7], SEQ[3], "#f0efec")]
    ax.legend(hs, ["rank 1", "rank 2-3", "rank 4-10", "rank > 10"], loc="upper left", bbox_to_anchor=(0, -0.26), ncol=4, frameon=False)
    ax.set_title("Positive controls: rank of the expected cell group", loc="left", fontweight="bold")
    D = pd.DataFrame(R, columns=FAMS).assign(row=[n for n, v in rows])
    D = D.merge(pd.DataFrame(Q, columns=[f + "_q" for f in FAMS]).assign(row=[n for n, v in rows]), on="row")
    leg = ("Rank of the best-ranked expected group (leiomyoma: uterine / visceral smooth muscle; hypertension: vascular, perivascular "
           "or smooth-muscle cells) per marker-selection family; * BH q < 0.05 within family; – not run (S-LDSC was run for the "
           "three thresholded rules only). The thresholded |log2FC| rule misses the leiomyoma control in all three settings "
           "(ranks %s). Expected-group patterns were fixed before the runs; the last row uses the post hoc Pericyte/SMC label check "
           "(v045). S-LDSC exact-cohort signals for the controls are nominal only (q ≥ 0.05)." % ", ".join(
               "%d" % R[i, 0] for i, (n, v) in enumerate(rows) if "leiomyoma" in n and np.isfinite(R[i, 0])))
    save(fig, "F5_positive_controls", D, leg)


# ----------------------------------------------------------------------------
def f6():
    f = need(RES / "external_atlas", "p9_v044_*/a_stromal_loci_vs_sets.csv", "v044 stromal loci test")
    A = pd.read_csv(f)
    order = [("ABS", "Venous/CTL"), ("ABS", "Arterial/CTL"), ("ABS", "Venous/ENDO"), ("DOWN", "Luminal/ENDO"), ("DOWN", "Ciliated/ENDO"),
             ("UP", "dStromal_mid/ENDO")]
    A["key"] = list(zip(A["family"], A["set"]))
    A = A.set_index("key").reindex([o for o in order if o in set(A["key"])]).reset_index()
    fig, ax = plt.subplots(figsize=(4.8, 2.6))
    yy = np.arange(len(A))[::-1]
    for yi, (_, r) in zip(yy, A.iterrows()):
        ax.plot([r["null_median"], r["null_p95"]], [yi, yi], color=GRID, linewidth=5, solid_capstyle="butt", zorder=1)
        ax.plot([r["null_median"]] * 2, [yi - 0.18, yi + 0.18], color=MUTED, linewidth=0.8, zorder=2)
        ax.scatter([r["delta_endothelium"]], [yi], s=28, color=BLUE, edgecolor="white", linewidth=0.6, zorder=3)
        ax.scatter([r["delta_stroma"]], [yi], s=28, color=ORANGE, edgecolor="white", linewidth=0.6, zorder=4)
        ax.text(max(r["delta_stroma"], r["null_p95"], r["delta_endothelium"]) + 0.08, yi, "perm p %.3f" % r["p_perm"], va="center", fontsize=6, color=INK2)
    ax.set_yticks(yy)
    ax.set_yticklabels(["%s %s" % (FAM_NAME.get(r["family"], r["family"]), nice(r["set"], 20)) for _, r in A.iterrows()], color=INK2)
    ax.axvline(0, color=AXIS, linewidth=0.6)
    ax.set_xlabel("loss of significance after removing loci (Δ log10 p)")
    ax.set_xlim(min(0, float(A[["delta_stroma", "delta_endothelium"]].min().min())) - 0.1, float(A[["delta_stroma", "null_p95", "delta_endothelium"]].max().max()) + 0.7)
    style(ax)
    ax.legend(handles=[Line2D([], [], marker="o", linestyle="", markersize=5, color=ORANGE, label="9 stromal loci removed"),
                       Line2D([], [], marker="o", linestyle="", markersize=5, color=BLUE, label="8 endothelial loci removed"),
                       Line2D([], [], color=GRID, linewidth=5, label="random 9-locus subsets (median-95th)")],
              loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=3, frameon=False, handletextpad=0.8, columnspacing=1.2)
    leg = ("KOLLER_EUR, MAGMA, HECA marker sets. Increase in log10 p of each first-ranked set after removing the genes within 500 kb "
           "of the nine stromal-assigned genome-wide loci (WNT4, GREB1, IGF1 region, TBX3, WT1, TSHZ3, HOXA9-11, OSR2 region, SFRP1) or "
           "of the eight endothelial-assigned loci, against 200 random 9-locus subsets of the non-stromal loci. The same stromal loci "
           "carry the up-only decidual-stroma signal, the down-only luminal / ciliated signals and the primary |log2FC| venous hit "
           "(H-STRVEN, formulated after v042).")
    save(fig, "F6_stromal_loci", A.drop(columns=["key"]), leg)


# ----------------------------------------------------------------------------
def s1():
    f = need(RES / "specificity", "p9_v045_*/c4_conditioned.csv", "v045 conditioning")
    C = pd.read_csv(f)
    C = C[(C["atlas"] == "HECA") & (C["trait"] == "KOLLER_EUR") & (C["rank_base"] == 1)]
    C = C.set_index("family").reindex([x for x in FAMS if x in set(C["family"])]).reset_index()
    fig, ax = plt.subplots(figsize=(4.8, 2.4))
    yy = np.arange(len(C))[::-1]
    for yi, (_, r) in zip(yy, C.iterrows()):
        xb, xc = nlp(r["p_base"]), nlp(r["p_cond"])
        ax.plot([xc, xb], [yi, yi], color=AXIS, linewidth=1.2, zorder=1)
        ax.scatter([xb], [yi], s=28, color="white", edgecolor=INK2, linewidth=1.2, zorder=3)
        ax.scatter([xc], [yi], s=28, color=INK, zorder=3)
        ax.text(max(xb, xc) + 0.15, yi, "rank after: %s" % ("%d" % r["rank_cond"] if np.isfinite(r["rank_cond"]) else "-"), va="center", fontsize=6, color=INK2)
    ax.set_yticks(yy)
    ax.set_yticklabels(["%s · %s" % (FAM_NAME[r["family"]], nice(r["label"], 20)) for _, r in C.iterrows()], color=INK2)
    ax.axvline(-math.log10(0.05), color=AXIS, linestyle=(0, (2, 2)), linewidth=0.6)
    ax.set_xlabel("−log10 p")
    ax.set_xlim(0, float(nlp(C["p_base"]).max()) + 1.6)
    style(ax)
    ax.legend(handles=[Line2D([], [], marker="o", linestyle="", markersize=5, markerfacecolor="white", markeredgecolor=INK2, label="base"),
                       Line2D([], [], marker="o", linestyle="", markersize=5, color=INK, label="conditioned")],
              loc="lower center", bbox_to_anchor=(0.4, 1.0), ncol=2, frameon=False)
    leg = ("First-ranked HECA label per family for KOLLER_EUR before and after conditioning on mean expression, expression breadth and "
           "gnomAD LOEUF (MAGMA gene-property analysis of binary set indicators); rank after conditioning shown on the right.")
    save(fig, "S1_conditioning_first_labels", C, leg)


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.48 | manuscript figures | " + TS)
    say("=" * 72)
    bad = []
    for name, fn in (("F1", f1), ("F2", f2), ("F3", f3), ("F4", f4), ("F5", f5), ("F6", f6), ("S1", s1)):
        try:
            fn()
        except Exception as e:
            import traceback
            bad.append(name)
            say("[FIG_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1200:])
            plt.close("all")
    (OUT / "figure_legends.md").write_text("# Paper 9 figure legends (auto-generated from the plotted data, v048)\n\n" + "\n".join(LEG))
    say("")
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    if bad:
        say("[INCOMPLETE] figures with errors: " + ", ".join(bad))
    else:
        say("[DONE] PAPER9_V048_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
