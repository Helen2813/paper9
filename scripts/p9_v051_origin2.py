#!/usr/bin/env python3
# Paper 9 | v0.0.51 | candidate-origin follow-up on GSE266265 (rules CHANGELOG v0.2.41)
#   O3 identity of the first-ranked stromal types (eF3, dS3): top up genes and a fixed marker panel
#   O4 joint MAGMA: stromal type (UPN) given each perivascular / mesothelial UPN set, and the reverse
#   F8 figure: best rank of each candidate-origin group per trait (v050 rankings)
#   I1 read-only feasibility check for real peritoneal mesothelium in the raw GSE213216 extracts
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). ~10-20 min; no downloads.

import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
try:
    os.nice(10)
except Exception:
    pass

import datetime as dt
import gzip
import json
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
RES = REPO / "results"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
DE_DIR = PROC / "genesets" / "gse266265_v050_de"
X_DIR = INTERIM / "gse213216_extract"
CELLS213 = PROC / "gse213216_cells_v028"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = RES / "review" / ("p9_v051_" + TS)
WORK = INTERIM / "magma_v051"
LOG = REPO / "logs" / ("p9_v051_origin2_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.41)
N_MATCH, MIN_SET = 400, 20
FOCUS = ["eF3", "dS3"]
OTHERS = ["Prv-MYH11", "VSMC", "Prv-CCL19", "Prv-STEAP4", "mesothelial"]
TRAITS_J = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL"]
STROMA_TYPES = ["eF1", "eF2", "eF3", "eF4-CXCL14", "dS1-myofibroblast", "dS2", "dS3", "fib C7", "fib C7 secretory"]
PANEL_TYPES = STROMA_TYPES + ["Prv-MYH11", "Prv-STEAP4", "VSMC", "mesothelial", "glandular"]
PANEL = {"hormone receptors": ["ESR1", "PGR", "AR"], "proliferation": ["MKI67", "TOP2A"], "decidualization": ["IGFBP1", "PRL", "FOXO1"],
         "fibroblast": ["PDGFRA", "DCN", "MMP11", "SFRP4"], "perivascular / MSC": ["SUSD2", "MCAM", "PDGFRB", "RGS5"],
         "stromal-locus genes": ["WNT4", "GREB1", "WT1", "HOXA10", "HOXA11", "SFRP1", "TBX3", "IGF1"], "mesothelial": ["MSLN", "UPK3B", "CALB2"]}
GROUPS = [("endometrial stroma", r"^(eF\d|dS\d|fib C7)"), ("perivascular / MSC-like", r"^(Prv|VSMC)"), ("mesothelium", r"^mesothelial"),
          ("epithelium", r"^(glandular|lumenal|ciliated|mid-secretory|MUC5B|TP63)"), ("endothelium", r"^(EC-|LEC)")]
MESO_GENES = ["MSLN", "UPK3B", "CALB2", "ITLN1", "PRG4"]

for d in (OUT, WORK, LOG.parent):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V051_ORIGIN2", "timestamp_utc": TS, "results": {}}
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


def fp(p):
    return "%.0e" % p if p < 0.001 else "%.3f" % p


# ----------------------------------------------------------------------------
def o3():
    say("")
    say("-" * 72)
    say("[SECTION] O3 identity of the first-ranked stromal types (GSE266265, donor-paired DE of v050)")
    DE = pd.read_csv(DE_DIR / "de_all.csv.gz")
    CPM = pd.read_csv(DE_DIR / "mean_cpm_by_type.csv.gz", index_col=0)
    CPM.index = CPM.index.astype(str)
    for ty in FOCUS:
        d = DE[DE["group"] == ty].sort_values("score", ascending=False)
        say("  [O3] %s top up genes (t): %s" % (ty, ", ".join("%s %.0f" % (g, s) for g, s in zip(d["gene"].head(25), d["score"].head(25)))))
        st = DE[DE["group"].isin([t for t in STROMA_TYPES if t != ty])]
        say("  [O3] %s n donors %d" % (ty, int(d["n_donors"].iloc[0]) if len(d) else 0))
    cols = [t for t in PANEL_TYPES if t in CPM.columns]
    rows = []
    say("  [O3] marker panel, mean CPM (columns: %s)" % ", ".join(cols))
    for grp, genes in PANEL.items():
        for g in genes:
            if g in CPM.index:
                v = CPM.loc[g, cols]
                rows.append({"panel": grp, "gene": g, **v.to_dict()})
                say("    %-20s %-7s %s" % (grp[:20], g, " ".join("%7.0f" % x for x in v.values)))
    pd.DataFrame(rows).to_csv(OUT / "o3_marker_panel.csv", index=False)


# ----------------------------------------------------------------------------
def o4(sym2ent):
    say("")
    say("-" * 72)
    say("[SECTION] O4 joint MAGMA on GSE266265 UPN sets: stromal type given perivascular / mesothelial type, and the reverse")
    DE = pd.read_csv(DE_DIR / "de_all.csv.gz")
    uni = sorted({g for g in DE["gene"].astype(str) if g in sym2ent})
    ent = [sym2ent[g] for g in uni]
    sets = {}
    for ty in FOCUS + OTHERS:
        d = DE[DE["group"] == ty]
        if len(d):
            gl = [g for g in d.sort_values("score", ascending=False)["gene"].astype(str).head(N_MATCH) if g in sym2ent]
            if len(gl) >= MIN_SET:
                sets[ty] = gl
    say("[O4] UPN sets: %s" % ", ".join("%s %d" % (k, len(v)) for k, v in sets.items()))
    rows = []
    for t in TRAITS_J:
        raw = MAGMA_OUT / (t + "_35_10.genes.raw")
        if not raw.exists():
            say("[O4] %s genes.raw missing" % t)
            continue
        say("  %s" % t)
        for a in [f for f in FOCUS if f in sets]:
            cells = []
            for b in [o for o in OTHERS if o in sets]:
                cf = WORK / "pair.txt"
                T = pd.DataFrame({"A": np.isin(uni, sets[a]).astype(int), "B": np.isin(uni, sets[b]).astype(int)}, index=ent)
                T[~T.index.duplicated()].to_csv(cf, sep="\t", index_label="GENE")
                m = magma_covar(raw, cf, "%s_m" % t, []) or {}
                ra = magma_covar(raw, cf, "%s_a" % t, ["B"]) or {}
                rb = magma_covar(raw, cf, "%s_b" % t, ["A"]) or {}
                r = {"trait": t, "focus": a, "other": b, "p_focus": m.get("A", np.nan), "p_other": m.get("B", np.nan),
                     "p_focus_given_other": ra.get("A", np.nan), "p_other_given_focus": rb.get("B", np.nan)}
                rows.append(r)
                cells.append("%s %s/%s" % (b, fp(r["p_focus_given_other"]), fp(r["p_other_given_focus"])))
            say("    %-4s (p %s) | given / reverse: %s" % (a, fp(rows[-1]["p_focus"]) if rows else "-", " | ".join(cells)))
    J = pd.DataFrame(rows)
    J.to_csv(OUT / "o4_joint_stroma_vs_perivascular.csv", index=False)
    k = J[(J["trait"] == "KOLLER_EUR") & (J["focus"] == "eF3")]
    ok = len(k) > 0 and bool((k["p_focus_given_other"] < 0.05).all())
    rev = list(k.loc[k["p_other_given_focus"] < 0.05, "other"]) if len(k) else []
    say("[H_STROMA_SEP_VERDICT] KOLLER_EUR eF3 (UPN) keeps p < 0.05 given every perivascular / mesothelial UPN set: %s -> %s | others keeping p < 0.05 given eF3: %s" % (
        "YES" if ok else "NO", "SUPPORTED" if ok else "NOT SUPPORTED", ", ".join(rev) or "none"))
    SUM["results"]["h_stroma_sep"] = {"supported": ok, "others_surviving": rev}


# ----------------------------------------------------------------------------
def f8():
    say("")
    say("-" * 72)
    say("[SECTION] F8 candidate-origin figure (v050 rankings, UPN)")
    fs = sorted((RES / "review").glob("p9_v050_*/r2_gse266265_family_rankings.csv"))
    if not fs:
        say("[F8] v050 rankings missing")
        return
    A = pd.read_csv(fs[-1])
    traits = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL", "LEIOMYOMA", "HYPERTENSION", "T2D"]
    groups = [g for g, _ in GROUPS]
    R = np.full((len(traits), len(groups)), np.nan)
    L = [["" for _ in groups] for _ in traits]
    for i, t in enumerate(traits):
        d = A[(A["trait"] == t) & (A["family"] == "UPN")].sort_values("rank")
        for j, g in enumerate(groups):
            h = d[d["group"] == g]
            if len(h):
                R[i, j] = h["rank"].iloc[0]
                L[i][j] = str(h["label"].iloc[0])
    D = pd.DataFrame(R, columns=groups).assign(trait=traits)
    D.to_csv(OUT / "F8_candidate_origin_data.csv", index=False)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap, BoundaryNorm
        plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7, "pdf.fonttype": 42})
        cat = np.where(np.isnan(R), np.nan, np.where(R == 1, 3, np.where(R <= 3, 2, np.where(R <= 10, 1, 0))))
        cmap = ListedColormap(["#f0efec", "#86b6ef", "#2a78d6", "#104281"])
        fig, ax = plt.subplots(figsize=(5.6, 2.9))
        ax.imshow(np.ma.masked_invalid(cat), cmap=cmap, norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], 4), aspect="auto")
        for i in range(R.shape[0]):
            for j in range(R.shape[1]):
                if np.isfinite(R[i, j]):
                    ax.text(j, i - 0.13, "%d" % R[i, j], ha="center", va="center", fontsize=7, fontweight="bold" if R[i, j] == 1 else "normal",
                            color="white" if cat[i, j] >= 2 else "#0b0b0b")
                    ax.text(j, i + 0.22, L[i][j][:14], ha="center", va="center", fontsize=5.2, color="white" if cat[i, j] >= 2 else "#52514e")
        for i in range(1, R.shape[0]):
            ax.axhline(i - 0.5, color="white", linewidth=1.5)
        for j in range(1, R.shape[1]):
            ax.axvline(j - 0.5, color="white", linewidth=1.5)
        ax.axhline(2.5, color="white", linewidth=4)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, rotation=20, ha="right", color="#52514e")
        ax.set_yticks(range(len(traits)))
        ax.set_yticklabels(traits, color="#52514e")
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title("Candidate cells of origin: best rank per group (GSE266265, up-only N = 400)", loc="left", fontweight="bold", fontsize=7.5)
        hs = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ("#104281", "#2a78d6", "#86b6ef", "#f0efec")]
        ax.legend(hs, ["rank 1", "rank 2-3", "rank 4-10", "rank > 10"], loc="upper left", bbox_to_anchor=(0, -0.28), ncol=4, frameon=False)
        for ext in ("png", "pdf"):
            fig.savefig(OUT / ("F8_candidate_origin.%s" % ext), dpi=300, bbox_inches="tight")
        plt.close(fig)
        say("[F8] figure written")
    except Exception as e:
        say("[F8] figure skipped: %s" % str(e)[:200])


# ----------------------------------------------------------------------------
def i1():
    say("")
    say("-" * 72)
    say("[SECTION] I1 raw GSE213216 extracts: is real peritoneal mesothelium available? (read-only)")
    if not X_DIR.exists():
        say("[I1] raw extract folder not found: %s" % X_DIR)
        return
    dirs = sorted(p for p in X_DIR.iterdir() if p.is_dir())
    say("[I1] %d sample folders" % len(dirs))
    ctx = {}
    f = CELLS213 / "obs.csv"
    if f.exists():
        O = pd.read_csv(f, usecols=["gsm", "context"])
        ctx = O.drop_duplicates("gsm").set_index("gsm")["context"].to_dict()
    say("[I1] folders by context: %s" % ", ".join("%s %d" % (k, v) for k, v in pd.Series([ctx.get(d.name, "not in cache") for d in dirs]).value_counts().items()))
    pick = [d for d in dirs if ctx.get(d.name) in ("peritoneum_unaffected", "peritoneal_lesion")][:4]
    try:
        import scipy.io
        import scipy.sparse as sp
    except Exception as e:
        say("[I1] scipy missing: %s" % e)
        return
    for d in pick:
        mtx = sorted(d.rglob("*matrix.mtx*"))
        feat = sorted(list(d.rglob("*features.tsv*")) + list(d.rglob("*genes.tsv*")))
        if not mtx or not feat:
            say("[I1] %s: no 10x matrix files (%s)" % (d.name, ", ".join(p.name for p in list(d.rglob("*"))[:5])))
            continue
        op = gzip.open if str(mtx[0]).endswith(".gz") else open
        with op(mtx[0], "rb") as fh:
            M = sp.csr_matrix(scipy.io.mmread(fh))
        opf = gzip.open if str(feat[0]).endswith(".gz") else open
        with opf(feat[0], "rt") as fh:
            names = [l.rstrip("\n").split("\t")[1] if "\t" in l else l.strip() for l in fh]
        if M.shape[0] != len(names):
            M = M.T.tocsr()
        idx = {g: names.index(g) for g in MESO_GENES if g in names}
        n_cells = M.shape[1]
        pos = None
        for g in ("MSLN", "UPK3B"):
            if g in idx:
                v = np.asarray(M[idx[g], :].todense()).ravel() > 0
                pos = v if pos is None else (pos & v)
        say("[I1] %s (%s): %d barcodes; MSLN+ & UPK3B+ barcodes %s; marker detection %s" % (
            d.name, ctx.get(d.name), n_cells, int(pos.sum()) if pos is not None else "-",
            ", ".join("%s %d" % (g, int((np.asarray(M[i, :].todense()).ravel() > 0).sum())) for g, i in idx.items())))


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.51 | candidate-origin follow-up: stromal identity, stroma vs perivascular, mesothelium feasibility | " + TS)
    say("=" * 72)
    t0 = time.time()
    sym2ent = gene_loc()
    for name, fn in (("O3", o3), ("O4", lambda: o4(sym2ent)), ("F8", f8), ("I1", i1)):
        try:
            fn()
        except Exception as e:
            import traceback
            ERR.append(name)
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1200:])
    with open(OUT / "p9_v051_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.41 before this run")
    if ERR:
        say("[INCOMPLETE] steps with errors: " + ", ".join(ERR))
    else:
        say("[DONE] PAPER9_V051_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
