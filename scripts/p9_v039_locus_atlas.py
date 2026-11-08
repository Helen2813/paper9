#!/usr/bin/env python3
# Paper 9 | v0.0.39 | Locus atlas (rules CHANGELOG v0.2.28): where in the endometrium do individual endometriosis
#   risk loci act, and which loci drive the aggregate cell-type rankings?
#   L1 loci (distance clumping, p < 5e-8) and their candidate genes mapped to endometrial compartments (HECA);
#   L2 leave-one-locus-out influence on the MAGMA gene-set rankings (|log2FC|, up-only, down-only sets of v031);
#   L3 loci shared with leiomyoma / T2D / hypertension (descriptive).
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). ~20-40 min.

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
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
from pathlib import Path

import numpy as np
import pandas as pd

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
INTERIM = REPO / "data" / "interim"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
G1000_BIM = TOOLS / "magma" / "ref_eur" / "g1000_eur.bim"
MAGMA_OUT = INTERIM / "magma"
HECA_GS = PROC / "genesets" / "heca_v014"
GS_P = PROC / "genesets" / "heca_koller_v031" / "P"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "locus_atlas" / ("p9_v039_" + TS)
WORK = INTERIM / "magma_v039"
LOG = REPO / "logs" / ("p9_v039_locus_atlas_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.28)
GWS = 5e-8
CLUMP_BP = 500_000
GENE_WIN_BP = 250_000
DROP_WIN_BP = 500_000
PCT_ASSIGN = 0.90   # a gene is assigned to the compartment where its specificity share ranks highest across genes, if >= 90th percentile
SHARED_P = 1e-5
PRIMARY = "KOLLER_EUR"
TRAITS = ["KOLLER_EUR", "GCST90205183", "LEIOMYOMA", "HYPERTENSION", "T2D"]
SHARE_TRAITS = ["LEIOMYOMA", "T2D", "HYPERTENSION"]
COMPARTMENTS = [("stroma", r"Stromal|Fibroblast"),
                ("epithelium", r"SOX9|Lumin|Gland|MUC5B|KRT5|iliated|HOXA13"),
                ("endothelium", r"^(Arterial|Venous|Lymphatic)"),
                ("perivascular/SM", r"^(ePV|mPV|uSMC)"),
                ("immune", r"^(uM|uNK|NK|T_|T$|B_|B$|Mast|DC|cDC|pDC|Macro|Mono|Lymphoid|Myeloid|ILC|Plasma|Neutro|Immune|Treg)")]

for d in (OUT, WORK, LOG.parent):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V039_LOCUS_ATLAS", "timestamp_utc": TS, "results": {}}


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
    ent2sym, sym2ent, pos = {}, {}, {}
    dup = set()
    for line in open(GENE_LOC_B37):
        t = line.split()
        if len(t) >= 6 and t[1] in [str(c) for c in range(1, 23)]:
            ent2sym[t[0]] = t[5]
            pos[t[0]] = (int(t[1]), int(t[2]), int(t[3]))
            if t[5] in sym2ent:
                dup.add(t[5])
            sym2ent[t[5]] = t[0]
    for d in dup:
        sym2ent.pop(d, None)
    return ent2sym, sym2ent, pos


# ----------------------------------------------------------------------------
# L1 loci and compartments
# ----------------------------------------------------------------------------
def read_sig(trait, thr):
    f = MAGMA_OUT / (trait + ".pval.txt")
    if not f.exists():
        return None
    out = {}
    with open(f) as fh:
        hdr = fh.readline().split()
        ip = hdr.index("P") if "P" in hdr else 1
        for line in fh:
            t = line.split()
            try:
                p = float(t[ip])
            except (ValueError, IndexError):
                continue
            if p < thr:
                out[t[0]] = p
    return out


def bim_positions(wanted):
    pos = {}
    with open(G1000_BIM) as fh:
        for line in fh:
            t = line.split()
            if t[1] in wanted:
                pos[t[1]] = (int(t[0]), int(t[3]))
    return pos


def clump(sig, pos):
    rows = sorted([(p, s) for s, p in sig.items() if s in pos])
    leads, taken = [], {}
    for p, s in rows:
        c, b = pos[s]
        if any(abs(b - b2) <= CLUMP_BP for b2 in taken.get(c, [])):
            continue
        taken.setdefault(c, []).append(b)
        leads.append({"lead_snp": s, "chr": c, "bp": b, "p": p})
    return pd.DataFrame(leads)


def compartment_table(S):
    comp = {}
    for c in S.columns:
        g = next((name for name, rx in COMPARTMENTS if re.search(rx, c)), "other")
        comp[c] = g
    names = [n for n, _ in COMPARTMENTS] + ["other"]
    share = pd.DataFrame({n: S[[c for c in S.columns if comp[c] == n]].sum(axis=1) if any(comp[c] == n for c in S.columns) else 0.0
                          for n in names}, index=S.index)
    base = share.median(axis=0)
    ratio = share.rank(pct=True)
    for n in names:
        if not any(comp[c] == n for c in S.columns):
            ratio[n] = np.nan
    return comp, share, ratio, base


def assign(ratio_row):
    r = ratio_row.drop(labels=["other"], errors="ignore").dropna()
    if not len(r):
        return "no expression", np.nan
    best = r.idxmax()
    return (best if r[best] >= PCT_ASSIGN else "broad"), float(r[best])


def loci_for(trait, sig, pos, genes_out, ent2sym, ratio, S):
    L = clump(sig, pos)
    if not len(L):
        return L
    g = genes_out
    rows = []
    for _, r in L.iterrows():
        m = g[(g["CHR"] == r["chr"]) & (g["STOP"] >= r["bp"] - GENE_WIN_BP) & (g["START"] <= r["bp"] + GENE_WIN_BP)].sort_values("P")
        syms = [ent2sym.get(str(e), str(e)) for e in m["GENE"]]
        prio = syms[0] if syms else ""
        expr = [s for s in syms if s in ratio.index]
        comp, best_ratio = assign(ratio.loc[expr[0]]) if expr else ("no expression", np.nan)
        top_type = S.loc[expr[0]].idxmax() if expr else ""
        # sensitivity: top-3 expressed candidate genes, mean ratio
        if expr:
            rr = ratio.loc[expr[:3]].mean(axis=0)
            comp3, _ = assign(rr)
        else:
            comp3 = "no expression"
        rows.append({**r.to_dict(), "n_genes": len(syms), "prioritized_gene": prio, "gene_min_p": float(m["P"].iloc[0]) if len(m) else np.nan,
                     "top3_genes": ",".join(syms[:3]), "expressed_gene_used": expr[0] if expr else "", "compartment": comp,
                     "compartment_pct": best_ratio, "top_heca_type": top_type, "compartment_top3": comp3})
    D = pd.DataFrame(rows)
    D.insert(0, "trait", trait)
    return D


def l1(ent2sym, sym2ent):
    say("")
    say("-" * 72)
    say("[SECTION] L1 risk loci and endometrial compartments (HECA heca_v014 specificity)")
    S = pd.read_csv(HECA_GS / "specificity_by_type.csv.gz", index_col=0)
    comp, share, ratio, base = compartment_table(S)
    for n in [n for n, _ in COMPARTMENTS] + ["other"]:
        say("[COMP] %-15s %s" % (n, ", ".join(c for c in S.columns if comp[c] == n) or "-"))
    say("[COMP] median specificity share per compartment: " + "; ".join("%s %.3f" % (k, v) for k, v in base.items() if np.isfinite(v)))
    say("[COMP] assignment: compartment with the highest percentile of its share across genes, if >= %.0f%%; otherwise 'broad'" % (100 * PCT_ASSIGN))
    sigs = {}
    for t in TRAITS:
        s = read_sig(t, GWS)
        if s is None:
            say("[L1] %s: pval file missing; skipped" % t)
            continue
        sigs[t] = s
    wanted = set().union(*[set(s) for s in sigs.values()]) if sigs else set()
    t0 = time.time()
    pos = bim_positions(wanted)
    say("[L1] genome-wide significant SNPs: %s; positions found %d (%.1f min)" % (
        ", ".join("%s %d" % (t, len(s)) for t, s in sigs.items()), len(pos), (time.time() - t0) / 60))
    all_loci = []
    for t, s in sigs.items():
        go = MAGMA_OUT / (t + "_35_10.genes.out")
        if not go.exists():
            say("[L1] %s genes.out missing" % t)
            continue
        g = pd.read_csv(go, sep=r"\s+")
        D = loci_for(t, s, pos, g, ent2sym, ratio, S)
        if not len(D):
            say("[L1] %s: no loci" % t)
            continue
        all_loci.append(D)
        vc = D["compartment"].value_counts()
        a = D[~D["compartment"].isin(["broad", "no expression"])]
        sh_ = a["compartment"].value_counts(normalize=True)
        say("[L1] %-12s loci %3d | assigned %3d | %s | broad %d, no expr %d" % (
            t[:12], len(D), len(a), ", ".join("%s %d (%.0f%%)" % (k, vc[k], 100 * sh_.get(k, 0)) for k in sh_.index),
            int(vc.get("broad", 0)), int(vc.get("no expression", 0))))
    if not all_loci:
        return None, ratio, S
    A = pd.concat(all_loci, ignore_index=True)
    A.to_csv(OUT / "l1_loci_compartments.csv", index=False)
    k = A[A["trait"] == PRIMARY]
    if len(k):
        say("[L1] %s loci (lead p, prioritized gene -> compartment, top HECA type):" % PRIMARY)
        for _, r in k.sort_values("p").iterrows():
            say("    chr%-2d %11d %-12s p=%.1e  %-12s -> %-15s %-22s (top3: %s)" % (
                r["chr"], r["bp"], r["lead_snp"], r["p"], str(r["prioritized_gene"])[:12], r["compartment"], str(r["top_heca_type"])[:22],
                str(r["top3_genes"])[:40]))
        a = k[~k["compartment"].isin(["broad", "no expression"])]
        sh_ = a["compartment"].value_counts(normalize=True)
        # background: the same assignment applied to every expressed gene
        from scipy.stats import binomtest
        bg = pd.Series([assign(ratio.loc[g_])[0] for g_ in ratio.index])
        bgs = bg.value_counts(normalize=True)
        n_expr = int((k["compartment"] != "no expression").sum())
        parts = []
        for cn in [n for n, _ in COMPARTMENTS] + ["broad"]:
            obs = int((k["compartment"] == cn).sum())
            exp_ = float(bgs.get(cn, 0.0))
            pv = binomtest(obs, n_expr, exp_, alternative="greater").pvalue if (n_expr and exp_ > 0) else np.nan
            parts.append("%s %d/%d (genome %.0f%%, p %.3g)" % (cn, obs, n_expr, 100 * exp_, pv))
            SUM["results"].setdefault("l1_enrichment", {})[cn] = {"obs": obs, "n": n_expr, "genome_share": exp_, "p_greater": pv}
        say("[L1_ENRICH] %s loci vs all expressed genes: %s" % (PRIMARY, "; ".join(parts)))
        n_ge10 = int((sh_ >= 0.10).sum())
        top_share = float(sh_.max()) if len(sh_) else np.nan
        multi = (n_ge10 >= 3) and (top_share <= 0.5)
        lei = A[(A["trait"] == "LEIOMYOMA") & ~A["compartment"].isin(["broad", "no expression"])]
        lei_top = float(lei["compartment"].value_counts(normalize=True).max()) if len(lei) else np.nan
        say("[H_MULTI_VERDICT] %s assigned loci %d: compartments with >= 10%% share %d, top share %.2f -> %s | leiomyoma top share %s (more concentrated: %s)" % (
            PRIMARY, len(a), n_ge10, top_share, "SUPPORTED" if multi else "NOT SUPPORTED",
            "%.2f" % lei_top if np.isfinite(lei_top) else "NA", "NA" if not np.isfinite(lei_top) else ("YES" if lei_top > top_share else "NO")))
        SUM["results"]["h_multi"] = {"assigned": len(a), "n_ge10": n_ge10, "top_share": top_share, "supported": multi, "leiomyoma_top_share": lei_top}
    return A, ratio, S


# ----------------------------------------------------------------------------
# L2 leave-one-locus-out influence on MAGMA gene-set rankings
# ----------------------------------------------------------------------------
def magma_sets_run(raw, sf, tag):
    outp = WORK / tag
    gsa = Path(str(outp) + ".gsa.out")
    for model in (["--model", "condition-hide=ALLGENES"], []):
        if gsa.exists():
            gsa.unlink()
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--set-annot", sf] + model + ["--out", outp], timeout=3600)
        if gsa.exists():
            break
    if not gsa.exists():
        say("    [MAGMA_FAIL] %s %s" % (tag, (o + e).strip()[-200:]))
        return None
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


def l2(A, sym2ent, pos_gene):
    say("")
    say("-" * 72)
    say("[SECTION] L2 leave-one-locus-out: which %s loci drive the gene-set rankings (v031 marker sets)" % PRIMARY)
    mf = GS_P / "manifest.csv"
    raw = MAGMA_OUT / (PRIMARY + "_35_10.genes.raw")
    if not mf.exists() or not raw.exists() or A is None:
        say("[L2] inputs missing; skipped")
        return
    man = pd.read_csv(mf)
    lab = dict(zip(man["set_key"], man["label"] + "/" + man["status"]))
    mode = dict(zip(man["set_key"], man["mode"]))
    sets = {}
    for k in man["set_key"]:
        sets[k] = [sym2ent[s] for s in open(GS_P / (k + ".txt")).read().split() if s in sym2ent]

    def write_sets(drop, path):
        with open(path, "w") as o:
            for k, g in sets.items():
                keep = g if k == "ALLGENES" else [x for x in g if x not in drop]
                o.write(k + " " + " ".join(keep) + "\n")

    sf = WORK / "sets_base.txt"
    write_sets(set(), sf)
    base = magma_sets_run(raw, sf, "base")
    if not base:
        return
    tracked = []
    for m in ("ABS", "UP", "DOWN"):
        ks = sorted([k for k in base if mode.get(k) == m], key=lambda k: base[k])[:3]
        tracked += ks
        say("[L2] baseline %-4s top: %s" % (m, "; ".join("%s p=%.1e" % (lab[k], base[k]) for k in ks)))
    K = A[A["trait"] == PRIMARY].reset_index(drop=True)
    rows = []
    t0 = time.time()
    sfl = WORK / "sets_drop.txt"
    for i, r in K.iterrows():
        drop = {e for e, (c, s, en) in pos_gene.items() if c == r["chr"] and en >= r["bp"] - DROP_WIN_BP and s <= r["bp"] + DROP_WIN_BP}
        write_sets(drop, sfl)
        res = magma_sets_run(raw, sfl, "drop")
        if not res:
            continue
        for k in tracked:
            rows.append({"locus": r["lead_snp"], "chr": r["chr"], "bp": r["bp"], "gene": r["prioritized_gene"], "compartment": r["compartment"],
                         "set": lab[k], "mode": mode[k], "p_base": base[k], "p_drop": res.get(k, np.nan),
                         "delta_log10": -math.log10(base[k]) + math.log10(max(res.get(k, np.nan), 1e-300)) if np.isfinite(res.get(k, np.nan)) else np.nan,
                         "n_dropped_genes": len(drop)})
    say("[L2] %d loci x %d tracked sets (%.0f min)" % (len(K), len(tracked), (time.time() - t0) / 60))
    R = pd.DataFrame(rows)
    if not len(R):
        return
    R.to_csv(OUT / "l2_leave_one_locus_out.csv", index=False)
    drivers = {}
    for s_, g in R.groupby("set", sort=False):
        d = g[g["delta_log10"] >= 1.0].sort_values("delta_log10", ascending=False)
        drivers[s_] = d
        say("    [DRIVERS] %-24s %-4s p=%.1e | %d loci raise p >= 10x: %s" % (
            s_[:24], g["mode"].iloc[0], g["p_base"].iloc[0], len(d),
            "; ".join("%s (%s, %s) %.1e" % (r["gene"], r["compartment"], r["locus"], r["p_drop"]) for _, r in d.head(5).iterrows()) or "-"))
    # all drivers of a set removed together
    for s_, d in drivers.items():
        if not len(d):
            continue
        loci = K[K["lead_snp"].isin(d["locus"])]
        drop = set()
        for _, r in loci.iterrows():
            drop |= {e for e, (c, s, en) in pos_gene.items() if c == r["chr"] and en >= r["bp"] - DROP_WIN_BP and s <= r["bp"] + DROP_WIN_BP}
        write_sets(drop, sfl)
        res = magma_sets_run(raw, sfl, "drop_all")
        k = next(k for k in tracked if lab[k] == s_)
        if res:
            say("    [DRIVERS_ALL] %-24s without its %d driver loci: p=%.1e (base %.1e)" % (s_[:24], len(d), res.get(k, np.nan), base[k]))
    SUM["results"]["l2_drivers"] = {s_: list(d["gene"]) for s_, d in drivers.items()}


# ----------------------------------------------------------------------------
# L3 loci shared with other traits (descriptive)
# ----------------------------------------------------------------------------
def l3(A):
    say("")
    say("-" * 72)
    say("[SECTION] L3 %s loci shared with leiomyoma / T2D / hypertension (min p in lead +-%d kb < %.0e; descriptive)" % (
        PRIMARY, GENE_WIN_BP // 1000, SHARED_P))
    if A is None:
        return
    K = A[A["trait"] == PRIMARY].copy()
    if not len(K):
        return
    win = {}
    for _, r in K.iterrows():
        win.setdefault(r["chr"], []).append((r["bp"] - GENE_WIN_BP, r["bp"] + GENE_WIN_BP, r["lead_snp"]))
    snp_locus = {}
    t0 = time.time()
    with open(G1000_BIM) as fh:
        for line in fh:
            t = line.split()
            c = int(t[0]) if t[0].isdigit() else -1
            if c in win:
                b = int(t[3])
                for lo, hi, ls in win[c]:
                    if lo <= b <= hi:
                        snp_locus[t[1]] = ls
                        break
    say("[L3] reference SNPs in %d windows: %d (%.1f min)" % (len(K), len(snp_locus), (time.time() - t0) / 60))
    for tr in SHARE_TRAITS:
        f = MAGMA_OUT / (tr + ".pval.txt")
        if not f.exists():
            say("[L3] %s pval file missing" % tr)
            continue
        best = {}
        with open(f) as fh:
            hdr = fh.readline().split()
            ip = hdr.index("P") if "P" in hdr else 1
            for line in fh:
                t = line.split()
                ls = snp_locus.get(t[0])
                if ls is None:
                    continue
                try:
                    p = float(t[ip])
                except (ValueError, IndexError):
                    continue
                if p < best.get(ls, 1.0):
                    best[ls] = p
        K["minp_" + tr] = K["lead_snp"].map(best)
        K["shared_" + tr] = K["minp_" + tr] < SHARED_P
    cols = [c for c in K.columns if c.startswith("shared_")]
    if not cols:
        return
    K["shared_any"] = K[cols].any(axis=1)
    K.to_csv(OUT / "l3_shared_loci.csv", index=False)
    for c in cols + ["shared_any"]:
        g = K[K[c]]
        say("[L3] %-20s %2d loci: %s" % (c, len(g), ", ".join("%s(%s)" % (r["prioritized_gene"], r["compartment"][:5]) for _, r in g.iterrows())[:300] or "-"))
    for flag, g in K.groupby("shared_any"):
        a = g[~g["compartment"].isin(["broad", "no expression"])]
        say("[L3] %-9s loci %2d | assigned compartments: %s" % ("shared" if flag else "specific", len(g),
                                                               ", ".join("%s %d" % (k, v) for k, v in a["compartment"].value_counts().items()) or "-"))


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.39 | locus atlas: compartments of risk loci and drivers of cell-type rankings | " + TS)
    say("=" * 72)
    t0 = time.time()
    ent2sym, sym2ent, pos_gene = gene_loc()
    A = None
    for name, fn in (("L1", lambda: l1(ent2sym, sym2ent)), ("L2", None), ("L3", None)):
        try:
            if name == "L1":
                A = fn()[0]
            elif name == "L2":
                l2(A, sym2ent, pos_gene)
            else:
                l3(A)
        except Exception as e:
            import traceback
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v039_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.28 before this run")
    say("[DONE] PAPER9_V039_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
