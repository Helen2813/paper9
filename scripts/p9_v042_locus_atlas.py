#!/usr/bin/env python3
# Paper 9 | v0.0.42 | Locus atlas, final form (rules CHANGELOG v0.2.28 - v0.2.32): v040 + explicit chromosome / mapping
#   checks, L3 sharing against random windows matched on genes and SNPs, L2c (loci removed by compartment, compared
#   with random locus subsets of the same size), K1c (which sample definition gives Koller's 22 / 13).
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
OUT = REPO / "results" / "locus_atlas" / ("p9_v042_" + TS)
WORK = INTERIM / "magma_v042"
LOG = REPO / "logs" / ("p9_v042_locus_atlas_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.28)
GWS = 5e-8
CLUMP_BP = 500_000
GENE_WIN_BP = 250_000
DROP_WIN_BP = 500_000
PCT_ASSIGN = 0.90   # a gene is assigned to the compartment where its specificity share ranks highest across genes, if >= 90th percentile
SHARED_P = 1e-5
N_PERM_L2C = 50
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
SUM = {"task": "PAPER9_V042_LOCUS_ATLAS", "timestamp_utc": TS, "results": {}}


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
    rows = sorted([(p, s) for s, p in sig.items() if s in pos and 1 <= pos[s][0] <= 22])
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
    say("[L1] loci = distance clumping of SNPs with p < 5e-8 (lead = smallest p, +-500 kb; no LD clumping); positions GRCh37 from the")
    say("[L1]   MAGMA 1000G EUR reference via rsID; autosomes 1-22 only (chrX / non-autosomal SNPs and genes excluded explicitly)")
    for t, s in sigs.items():
        nonauto = sum(1 for x in s if x in pos and not (1 <= pos[x][0] <= 22))
        if nonauto:
            say("[L1_CHECK] %s: %d significant SNPs on non-autosomal chromosomes excluded" % (t, nonauto))
        go = MAGMA_OUT / (t + "_35_10.genes.out")
        if not go.exists():
            say("[L1] %s genes.out missing" % t)
            continue
        g = pd.read_csv(go, sep=r"\s+", dtype={"GENE": str, "CHR": str})
        n0 = len(g)
        nonauto_vals = sorted(set(g.loc[~g["CHR"].str.fullmatch(r"\d+"), "CHR"]))
        for c in ("CHR", "START", "STOP", "P"):
            g[c] = pd.to_numeric(g[c], errors="coerce")
        g = g.dropna(subset=["CHR", "START", "STOP", "P"])
        g["CHR"] = g["CHR"].astype(int)
        g = g[(g["CHR"] >= 1) & (g["CHR"] <= 22)]
        say("[L1_CHECK] %s genes.out: %d genes, autosomal %d, excluded %d (chromosome values %s)" % (
            t, n0, len(g), n0 - len(g), ",".join(nonauto_vals) or "none"))
        D = loci_for(t, s, pos, g, ent2sym, ratio, S)
        if not len(D):
            say("[L1] %s: no loci" % t)
            continue
        all_loci.append(D)
        say("[L1_CHECK] %s: loci %d, with >= 1 candidate gene %d, candidate genes total %d, with an expressed (HECA) gene %d" % (
            t, len(D), int((D["n_genes"] > 0).sum()), int(D["n_genes"].sum()), int((D["compartment"] != "no expression").sum())))
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
    # L2b: every genome-wide significant locus removed at once
    drop = set()
    for _, r in K.iterrows():
        drop |= {e for e, (c, s, en) in pos_gene.items() if c == r["chr"] and en >= r["bp"] - DROP_WIN_BP and s <= r["bp"] + DROP_WIN_BP}
    write_sets(drop, sfl)
    res = magma_sets_run(raw, sfl, "drop_all_gws")
    if res:
        keep = [k for k in tracked if np.isfinite(res.get(k, np.nan))]
        say("[L2_ALL_GWS] all %d loci removed (%d genes out of the sets): %s" % (
            len(K), len(drop), "; ".join("%s %s p=%.1e (base %.1e)" % (mode[k], lab[k], res[k], base[k]) for k in keep)))
        for m in ("ABS", "UP", "DOWN"):
            ks = sorted([k for k in res if mode.get(k) == m], key=lambda k: res[k])[:3]
            say("    [L2_ALL_GWS] %-4s new top: %s" % (m, "; ".join("%s p=%.1e" % (lab[k], res[k]) for k in ks)))
        surv = [lab[k] for k in keep if res[k] < 0.05]
        say("[L2_ALL_GWS_VERDICT] tracked sets still p < 0.05 without any genome-wide locus: %d/%d -> %s" % (
            len(surv), len(keep), "signal persists without the genome-wide loci (distributed over sub-threshold regions)" if len(surv) >= len(keep) / 2 else "signal depends on the genome-wide loci jointly"))
        SUM["results"]["l2_all_gws"] = {lab[k]: {"base": base[k], "without_gws": res[k]} for k in keep}
    # L2c: loci removed by compartment; endothelium and stroma compared with random locus subsets of the same size
    say("[L2c] loci removed by assigned compartment (all loci of the compartment at once)")
    rng = np.random.default_rng(20261108)

    def drop_for(lead_ids):
        d_ = set()
        for _, r in K[K["lead_snp"].isin(lead_ids)].iterrows():
            d_ |= {e for e, (c, s, en) in pos_gene.items() if c == r["chr"] and en >= r["bp"] - DROP_WIN_BP and s <= r["bp"] + DROP_WIN_BP}
        return d_

    focus = {"endothelium": next((k for k in tracked if mode[k] == "ABS"), None), "stroma": next((k for k in tracked if mode[k] == "UP"), None)}
    l2c = {}
    for comp_name, g in K.groupby("compartment"):
        ids = list(g["lead_snp"])
        write_sets(drop_for(ids), sfl)
        res = magma_sets_run(raw, sfl, "drop_comp")
        if not res:
            continue
        line = "; ".join("%s %s %.1e" % (mode[k], lab[k], res.get(k, np.nan)) for k in tracked)
        say("    [L2c] %-15s %2d loci (%s) -> %s" % (comp_name, len(ids), ",".join(str(x) for x in g["prioritized_gene"].head(9)), line))
        l2c[comp_name] = {lab[k]: res.get(k, np.nan) for k in tracked}
        if comp_name in focus and focus[comp_name]:
            k = focus[comp_name]
            obs = math.log10(res[k]) - math.log10(base[k])
            others = [x for x in K["lead_snp"] if x not in ids]
            null = []
            for _ in range(N_PERM_L2C):
                pick = list(rng.choice(others, size=min(len(ids), len(others)), replace=False))
                write_sets(drop_for(pick), sfl)
                r2 = magma_sets_run(raw, sfl, "drop_rand")
                if r2 and np.isfinite(r2.get(k, np.nan)):
                    null.append(math.log10(r2[k]) - math.log10(base[k]))
            pperm = (1 + sum(x >= obs for x in null)) / (1 + len(null))
            say("[L2c_VERDICT] removing the %d %s loci raises %s %s p by %.2f log10 units; random %d-locus subsets: median %.2f, 95th pct %.2f (n=%d) -> permutation p %.3f -> %s" % (
                len(ids), comp_name, mode[k], lab[k], obs, len(ids), float(np.median(null)) if null else float("nan"),
                float(np.percentile(null, 95)) if null else float("nan"), len(null), pperm, "SUPPORTED" if pperm < 0.05 else "NOT SUPPORTED"))
            SUM["results"]["l2c_" + comp_name] = {"set": lab[k], "delta_log10": obs, "null_median": float(np.median(null)) if null else None, "p_perm": pperm}
    SUM["results"]["l2c"] = l2c


# ----------------------------------------------------------------------------
# L3 loci shared with other traits (descriptive)
# ----------------------------------------------------------------------------
CHR_LEN_B37 = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516,
               133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]
N_POOL = 800
K_MATCH = 7
RAND_SEED = 20261108


def l3(A, pos_gene):
    from scipy.stats import fisher_exact
    say("")
    say("-" * 72)
    say("[SECTION] L3 %s loci shared with leiomyoma / T2D / hypertension (min p in lead +-%d kb < %.0e) vs random windows matched on genes and SNPs" % (
        PRIMARY, GENE_WIN_BP // 1000, SHARED_P))
    if A is None:
        return
    K = A[A["trait"] == PRIMARY].copy()
    if not len(K):
        return
    rng = np.random.default_rng(RAND_SEED)
    lens = np.array(CHR_LEN_B37, dtype=float)
    taken = {c: [int(b) for b in K.loc[K["chr"] == c, "bp"]] for c in range(1, 23)}
    pool, tries = [], 0
    while len(pool) < N_POOL and tries < 200000:
        tries += 1
        c = int(rng.choice(22, p=lens / lens.sum())) + 1
        b = int(rng.integers(1_000_000, CHR_LEN_B37[c - 1] - 1_000_000))
        lead_d = min([abs(b - x) for x in K.loc[K["chr"] == c, "bp"]] or [10 ** 9])
        if lead_d < 1_000_000 or any(abs(b - x) < 2 * GENE_WIN_BP for x in taken[c]):
            continue
        taken[c].append(b)
        pool.append((c, b, "R%03d" % len(pool)))
    wins = [(int(r["chr"]), int(r["bp"]), r["lead_snp"]) for _, r in K.iterrows()] + pool
    bych = {}
    for c, b, name in wins:
        bych.setdefault(c, []).append((b - GENE_WIN_BP, b + GENE_WIN_BP, name))
    for c in bych:
        bych[c].sort()
    snp_win, nsnp = {}, {}
    t0 = time.time()
    for ch in pd.read_csv(G1000_BIM, sep=r"\s+", header=None, usecols=[0, 1, 3], names=["c", "snp", "bp"], chunksize=2_000_000,
                          dtype={"c": str, "snp": str, "bp": np.int64}):
        ch = ch[ch["c"].str.fullmatch(r"\d+")]
        ch["c"] = ch["c"].astype(int)
        for c, g in ch.groupby("c"):
            if c not in bych:
                continue
            st = np.array([w[0] for w in bych[c]])
            en = np.array([w[1] for w in bych[c]])
            nm = [w[2] for w in bych[c]]
            bp = g["bp"].to_numpy()
            ix = np.searchsorted(st, bp, side="right") - 1
            ok = (ix >= 0) & (bp <= en[np.clip(ix, 0, None)])
            for snp, i in zip(g["snp"].to_numpy()[ok], ix[ok]):
                snp_win[snp] = nm[i]
                nsnp[nm[i]] = nsnp.get(nm[i], 0) + 1
    ngene = {}
    for c, b, name in wins:
        ngene[name] = sum(1 for e, (cc, s_, en_) in pos_gene.items() if cc == c and en_ >= b - GENE_WIN_BP and s_ <= b + GENE_WIN_BP)
    feat = lambda nm: np.array([np.log1p(ngene.get(nm, 0)), np.log1p(nsnp.get(nm, 0))])
    P = np.array([feat(w[2]) for w in pool])
    sd = P.std(axis=0)
    sd[sd == 0] = 1
    used, matched = set(), []
    for _, r in K.sort_values("p").iterrows():
        f0 = feat(r["lead_snp"])
        d = np.sqrt((((P - f0) / sd) ** 2).sum(axis=1))
        order = [i for i in np.argsort(d) if i not in used and nsnp.get(pool[i][2], 0) >= 50][:K_MATCH]
        used |= set(order)
        matched += [pool[i][2] for i in order]
    pool_ok = [w[2] for w in pool if nsnp.get(w[2], 0) >= 50]
    med = lambda names, dct: float(np.median([dct.get(x, 0) for x in names])) if names else float("nan")
    say("[L3] reference SNPs in windows %d (%.1f min); pool %d windows (%d with >= 50 SNPs); matched controls %d" % (
        len(snp_win), (time.time() - t0) / 60, len(pool), len(pool_ok), len(matched)))
    say("[L3_BALANCE] median genes / SNPs per window: endometriosis loci %.0f / %.0f; matched %.0f / %.0f; all pool %.0f / %.0f" % (
        med(list(K["lead_snp"]), ngene), med(list(K["lead_snp"]), nsnp), med(matched, ngene), med(matched, nsnp), med(pool_ok, ngene), med(pool_ok, nsnp)))
    rows = []
    for tr in SHARE_TRAITS:
        f = MAGMA_OUT / (tr + ".pval.txt")
        if not f.exists():
            say("[L3] %s pval file missing" % tr)
            continue
        b_ = {}
        with open(f) as fh:
            hdr = fh.readline().split()
            ip = hdr.index("P") if "P" in hdr else 1
            for line in fh:
                t = line.split()
                w = snp_win.get(t[0])
                if w is None:
                    continue
                try:
                    p = float(t[ip])
                except (ValueError, IndexError):
                    continue
                if p < b_.get(w, 1.0):
                    b_[w] = p
        k_sh = int(sum(b_.get(ls, 1.0) < SHARED_P for ls in K["lead_snp"]))
        m_sh = int(sum(b_.get(w, 1.0) < SHARED_P for w in matched))
        a_sh = int(sum(b_.get(w, 1.0) < SHARED_P for w in pool_ok))
        orm, pm = fisher_exact([[k_sh, len(K) - k_sh], [m_sh, len(matched) - m_sh]], alternative="greater")
        ora, pa = fisher_exact([[k_sh, len(K) - k_sh], [a_sh, len(pool_ok) - a_sh]], alternative="greater")
        K["minp_" + tr] = K["lead_snp"].map(b_)
        K["shared_" + tr] = K["minp_" + tr] < SHARED_P
        rows.append({"trait": tr, "loci_shared": k_sh, "loci": len(K), "matched_shared": m_sh, "matched": len(matched), "or_matched": orm,
                     "p_matched": pm, "pool_shared": a_sh, "pool": len(pool_ok), "or_pool": ora, "p_pool": pa})
        say("[L3] %-12s shared %2d/%d loci (%.0f%%) | matched windows %d/%d (%.0f%%): OR %.2f, p %.2g | all pool %d/%d (%.0f%%): OR %.2f, p %.2g" % (
            tr, k_sh, len(K), 100 * k_sh / len(K), m_sh, len(matched), 100 * m_sh / max(1, len(matched)), orm, pm,
            a_sh, len(pool_ok), 100 * a_sh / max(1, len(pool_ok)), ora, pa))
        g = K[K["shared_" + tr]]
        say("       %s" % ", ".join("%s(%s)" % (r["prioritized_gene"] or r["lead_snp"], str(r["compartment"])[:5]) for _, r in g.iterrows())[:400])
    pd.DataFrame(rows).to_csv(OUT / "l3_sharing_vs_matched.csv", index=False)
    cols = [c for c in K.columns if c.startswith("shared_")]
    if cols:
        K["shared_any"] = K[cols].any(axis=1)
        for flag, g in K.groupby("shared_any"):
            a = g[~g["compartment"].isin(["broad", "no expression"])]
            say("[L3] %-9s loci %2d | assigned compartments: %s" % ("shared" if flag else "specific", len(g),
                                                                   ", ".join("%s %d" % (k, v) for k, v in a["compartment"].value_counts().items()) or "-"))
    K.to_csv(OUT / "l3_shared_loci.csv", index=False)
    SUM["results"]["l3"] = rows


# ----------------------------------------------------------------------------
# K1b Koller cohort audit against the E-MTAB-14039 sample table (descriptive)
# ----------------------------------------------------------------------------
HECA_FILE = REPO / "data" / "raw" / "heca" / "endometriumAtlasV2_cells_with_counts.h5ad"
SDRF = REPO / "data" / "raw" / "arrayexpress" / "E-MTAB-14039.sdrf.txt"
SDRF_URLS = ["https://www.ebi.ac.uk/biostudies/files/E-MTAB-14039/E-MTAB-14039.sdrf.txt",
             "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/039/E-MTAB-14039/Files/E-MTAB-14039.sdrf.txt"]


def obs_cat(f, name):
    import h5py
    g = f["obs"][name]
    dec = lambda x: x.decode() if isinstance(x, bytes) else str(x)
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = np.array([dec(c) for c in g["categories"][:]], dtype=object)
        codes = g["codes"][:]
        out = np.empty(len(codes), dtype=object)
        out[codes >= 0] = cats[codes[codes >= 0]]
        return out
    return np.array([dec(x) for x in g[:]], dtype=object)


def k1b():
    import h5py
    say("")
    say("-" * 72)
    say("[SECTION] K1b Koller cohort audit: E-MTAB-14039 sample table vs HECA libraries (descriptive)")
    SDRF.parent.mkdir(parents=True, exist_ok=True)
    if not SDRF.exists():
        for u in SDRF_URLS:
            rc, o, e = sh(["curl", "-sSL", "--fail", "--max-time", "300", "-A", "paper9", "-o", SDRF, u], timeout=400)
            if rc == 0 and SDRF.exists() and SDRF.stat().st_size > 1000:
                say("[K1b] downloaded %s" % u)
                break
            if SDRF.exists():
                SDRF.unlink()
    if not SDRF.exists():
        say("[K1b] sample table not available (download failed)")
        return
    T = pd.read_csv(SDRF, sep="\t", dtype=str).fillna("")
    say("[K1b] rows %d, columns %d" % (len(T), T.shape[1]))
    keys = re.compile(r"individual|disease|endometriosis|treat|hormon|stage|phase|library|single|nucle|assay name|source name|organism part|sampling", re.I)
    for c in T.columns:
        if keys.search(c):
            vc = T[c].value_counts()
            say("    [SDRF] %-45s n=%-4d %s" % (c[:45], len(vc), " | ".join("%s=%d" % (k[:30], v) for k, v in vc.head(6).items())))
    with h5py.File(HECA_FILE, "r") as f:
        cols = list(f["obs"].keys())
        ds = obs_cat(f, "dataset") if "dataset" in cols else None
        lib = obs_cat(f, "library") if "library" in cols else None
        path = obs_cat(f, "Endometrial_pathology")
        stage = obs_cat(f, "Binary Stage")
        don = obs_cat(f, "genotype")
    if ds is None or lib is None:
        say("[K1b] HECA dataset / library columns missing")
        return
    m = ds == "Mareckova"
    H = pd.DataFrame({"library": lib[m], "pathology": path[m], "stage": stage[m], "donor": don[m]}).drop_duplicates("library")
    text = "\n".join("\t".join(r) for r in T.values.tolist())
    H["in_sdrf"] = H["library"].map(lambda x: bool(x) and x in text)
    say("[K1b] HECA Mareckova libraries %d; found in the sample table %d" % (len(H), int(H["in_sdrf"].sum())))
    for (pa, st), g in H.groupby([H["pathology"].fillna("NA"), H["stage"].fillna("NA")]):
        say("    [HECA_LIB] pathology %-10s stage %-14s libraries %3d donors %3d (in table %d)" % (pa[:10], st[:14], len(g), g["donor"].nunique(), int(g["in_sdrf"].sum())))
    H.to_csv(OUT / "k1b_heca_mareckova_libraries.csv", index=False)
    T.to_csv(OUT / "k1b_sdrf_copy.tsv", sep="\t", index=False)
    # K1c: which definition gives Koller's 22 endometriosis / 13 control samples?
    col = lambda pat: next((c for c in T.columns if re.search(pat, c, re.I)), None)
    c_ind, c_dis, c_lib, c_stage = col(r"Characteristics\[individual\]"), col(r"Characteristics\[disease\]"), col(r"Characteristics\[library_id\]"), col(r"Characteristics\[menstrual_cycle_stage\]$")
    if all([c_ind, c_dis, c_lib, c_stage]):
        U = T[[c_ind, c_dis, c_lib, c_stage]].drop_duplicates()
        U.columns = ["individual", "disease", "library", "stage"]
        U["in_heca_cells"] = U["library"].isin(set(H["library"])) | U["library"].map(lambda x: any(x and x in h for h in H["library"]))
        U["hormones"] = U["stage"] == "Hormones"
        say("[K1c] distinct (library, individual) pairs %d; individuals %d; libraries %d" % (len(U), U["individual"].nunique(), U["library"].nunique()))
        for scope, m_ in (("all libraries", np.ones(len(U), bool)), ("libraries in HECA cells", U["in_heca_cells"].to_numpy()), ("libraries not in HECA cells", ~U["in_heca_cells"].to_numpy())):
            for horm in (False, True):
                V_ = U[m_ & (U["hormones"].to_numpy() if horm else True)] if horm else U[m_]
                if horm:
                    V_ = U[m_]
                else:
                    V_ = U[m_ & ~U["hormones"].to_numpy()]
                e, c_ = V_[V_["disease"] == "endometriosis"], V_[V_["disease"] == "normal"]
                say("    [K1c] %-28s hormones %-8s individuals E %2d C %2d | (library, individual) pairs E %2d C %2d | libraries E %2d C %2d" % (
                    scope, "included" if horm else "excluded", e["individual"].nunique(), c_["individual"].nunique(), len(e), len(c_),
                    e["library"].nunique(), c_["library"].nunique()))
        U.to_csv(OUT / "k1c_sample_pairs.csv", index=False)


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.42 | locus atlas: compartments, drivers by compartment, matched sharing, Koller cohort | " + TS)
    say("=" * 72)
    t0 = time.time()
    ent2sym, sym2ent, pos_gene = gene_loc()
    A = None
    try:
        k1b()
    except Exception as e:
        say("[STEP_ERROR] K1b %s: %s" % (type(e).__name__, str(e)[:300]))
    for name, fn in (("L1", lambda: l1(ent2sym, sym2ent)), ("L2", None), ("L3", None)):
        try:
            if name == "L1":
                A = fn()[0]
            elif name == "L2":
                l2(A, sym2ent, pos_gene)
            else:
                l3(A, pos_gene)
        except Exception as e:
            import traceback
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v042_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.28 - v0.2.32 before this run (L2c formulated after the v040 locus assignments)")
    say("[DONE] PAPER9_V042_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
