#!/usr/bin/env python3
# Paper 9 | v0.0.27 | PLAN v0.3, Option A core (rules CHANGELOG v0.2.17, fixed before this run):
#   R1 covariate structure / VIF   R2 pairwise conditional matrix   R3 delete-one-block jackknife (G4 and G1)
#   R4 simulations with a known source of signal (misassignment)   R5 false-positive rate under a true null
#   B0 fallopian-tube atlas discovery (CELLxGENE listing only)   K0 Koller supplementary tables 9-10 (if present)
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). ~1.5-2 h (about 1,400 short MAGMA runs).

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
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------
# Hardcoded locations
# ----------------------------------------------------------------------------
HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
INTERIM = REPO / "data" / "interim"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
HECA_GS = PROC / "genesets" / "heca_v014"
COV_M1 = INTERIM / "magma_v020" / "C5_M1_specificity.covar.txt"
COVAR_G1 = PROC / "magma_covar" / "g1"
PB213 = PROC / "gse213216_pseudobulk_clusters"
REFS = REPO / "references"
CXG_API = "https://api.cellxgene.cziscience.com/curation/v1/collections"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "resolution" / ("p9_v027_" + TS)
WORK = INTERIM / "magma_v027"
LOG = REPO / "logs" / ("p9_v027_resolution_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.17)
FOCUS = ["dStromal_mid", "ePV_2", "ePV_1b", "Arterial", "uSMCs", "mPV"]
R1_TYPES = ["dStromal_early", "dStromal_mid", "dStromal_late", "ePV_1a", "ePV_1b", "ePV_2", "mPV", "Arterial", "Venous", "uSMCs"]
R2_TRAITS = ["KOLLER_EUR", "GCST90205183", "HYPERTENSION", "T2D", "LEIOMYOMA"]
JK_TRAITS = ["KOLLER_EUR", "GCST90205183", "HYPERTENSION", "T2D"]
JK_TYPES = ["dStromal_mid", "ePV_2", "Arterial", "uSMCs"]
N_BLOCKS = 50
SIM_TARGETS = ["dStromal_mid", "ePV_2", "Arterial", "uSMCs"]
SIM_STRENGTH = [1.0, 2.0]
N_SIM = 30
OBS_BETA_DSMID = 0.947      # KOLLER_EUR dStromal_mid beta (v024 base model)
R5_SETS = [("eutopic", "stromal"), ("eutopic", "epithelial"), ("peritoneal_lesion", "stromal")]
N_SPLIT = 150
SEED = 20261106
DS3 = ["dStromal_early", "dStromal_mid", "dStromal_late"]

for d in (OUT, WORK, LOG.parent):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V027_RESOLUTION", "timestamp_utc": TS, "results": {}}


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


def gene_positions():
    pos = {}
    for line in open(GENE_LOC_B37):
        x = line.split()
        if len(x) >= 6 and x[1] in [str(c) for c in range(1, 23)]:
            pos[x[0]] = (int(x[1]), int(x[2]))
    return pos


def magma_gp(raw, covf, cond, tag):
    raw = Path(raw)
    if not raw.exists():
        return None
    outp = WORK / tag
    gsa = Path(str(outp) + ".gsa.out")
    keys = ["condition-hide=", "condition="] if cond else [None]
    o = e = ""
    for key in keys:
        if gsa.exists():
            gsa.unlink()
        model = ["direction=greater"] if key is None else [key + cond, "direction=greater"]
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--gene-covar", covf, "--model"] + model + ["--out", outp],
                      timeout=3600)
        if rc == 0 and gsa.exists():
            rows = [l for l in open(gsa) if not l.startswith("#") and l.strip()]
            hdr = rows[0].split()
            res = {}
            for l in rows[1:]:
                d = dict(zip(hdr, l.split()))
                name = d.get("FULL_NAME", d.get("VARIABLE"))
                try:
                    res[name] = {"beta": float(d["BETA"]), "beta_std": float(d.get("BETA_STD", "nan")),
                                 "se": float(d["SE"]), "p": float(d["P"])}
                except (KeyError, ValueError):
                    continue
            return res
    say("    [MAGMA_FAIL] %s %s" % (tag, (o + e).strip()[-200:]))
    return None


def raw_of(trait):
    return MAGMA_OUT / (trait + "_35_10.genes.raw")


def heca_tables():
    S = pd.read_csv(HECA_GS / "specificity_by_type.csv.gz", index_col=0)
    M = pd.read_csv(HECA_GS / "mean_cpm_by_type.csv.gz", index_col=0).loc[S.index, S.columns]
    return S, M


# ----------------------------------------------------------------------------
# R1 / R2
# ----------------------------------------------------------------------------
def focus_cov(S, M, sym2):
    keys = {}
    cov = pd.DataFrame(index=S.index)
    for i, t in enumerate([t for t in FOCUS if t in S.columns]):
        k = "F%d" % (i + 1)
        keys[k] = t
        cov[k] = S[t]
    cov["Average"] = np.log2(M.mean(axis=1) + 1.0)
    cov.insert(0, "GENE", [sym2.get(g) for g in cov.index])
    cov = cov.dropna(subset=["GENE"])
    f = WORK / "focus.covar.txt"
    cov.to_csv(f, sep="\t", index=False, float_format="%.6g")
    return f, keys, cov


def vif(y, X):
    X1 = np.column_stack([np.ones(len(y)), X])
    b, *_ = np.linalg.lstsq(X1, y, rcond=None)
    r2 = 1 - np.sum((y - X1 @ b) ** 2) / np.sum((y - y.mean()) ** 2)
    return r2, (1.0 / (1.0 - r2)) if r2 < 1 else float("inf")


def r1_r2(S, M, sym2):
    say("")
    say("-" * 72)
    say("[SECTION] R1 covariate structure (HECA control-donor specificity)")
    types = [t for t in R1_TYPES if t in S.columns]
    C = S[types]
    pr = C.corr(method="pearson")
    sp = C.corr(method="spearman")
    short = {t: t.replace("dStromal_", "dS_") for t in types}
    say("    [R1] Pearson / Spearman with dStromal_mid: " + " | ".join(
        "%s %.2f/%.2f" % (short[t], pr.loc["dStromal_mid", t], sp.loc["dStromal_mid", t]) for t in types if t != "dStromal_mid"))
    for a in ["ePV_2", "Arterial", "uSMCs"]:
        if a in types:
            say("    [R1] Pearson with %-8s: " % a + " ".join("%s=%.2f" % (short[t], pr.loc[a, t]) for t in types if t != a))
    y = S["dStromal_mid"].values
    v3 = [t for t in ["ePV_2", "Arterial", "uSMCs"] if t in S.columns]
    r2a, va = vif(y, S[v3].values)
    others = [t for t in S.columns if t not in DS3]
    r2b, vb = vif(y, S[others].values)
    say("    [R1] dStromal_mid ~ %s: R2=%.3f VIF=%.2f | ~ all %d non-decidual types: R2=%.3f VIF=%.2f" % (
        "+".join(v3), r2a, va, len(others), r2b, vb))
    pr.to_csv(OUT / "R1_pearson.csv")
    SUM["results"]["R1"] = {"vif_vasc3": va, "r2_vasc3": r2a, "vif_all_nondecidual": vb}

    say("")
    say("-" * 72)
    say("[SECTION] R2 pairwise conditional matrix: p(type i | Average + type j)")
    covf, keys, _ = focus_cov(S, M, sym2)
    inv = {v: k for k, v in keys.items()}
    rows = []
    for trait in R2_TRAITS:
        base = magma_gp(raw_of(trait), covf, "Average", "R2_%s_base" % trait)
        if not base:
            continue
        mat = {}
        for tj, kj in inv.items():
            r = magma_gp(raw_of(trait), covf, "Average," + kj, "R2_%s_%s" % (trait, kj))
            if r:
                mat[tj] = {ti: r[ki]["p"] for ti, ki in inv.items() if ki in r}
        say("  [R2] %s  (rows: tested | cols: conditioned on; 'base' = Average only)" % trait)
        say("      %-13s %7s " % ("", "base") + " ".join("%8s" % t[:8] for t in inv))
        for ti, ki in inv.items():
            say("      %-13s %7.3f " % (ti[:13], base[ki]["p"] if ki in base else np.nan) +
                " ".join("%8s" % ("-" if ti == tj else ("%.3f" % mat.get(tj, {}).get(ti, np.nan))) for tj in inv))
            for tj in inv:
                if tj != ti:
                    rows.append({"trait": trait, "tested": ti, "conditioned_on": tj, "p": mat.get(tj, {}).get(ti, np.nan),
                                 "p_base": base[ki]["p"] if ki in base else np.nan})
    pd.DataFrame(rows).to_csv(OUT / "R2_conditional.csv", index=False)
    return covf, keys


# ----------------------------------------------------------------------------
# R3 jackknife
# ----------------------------------------------------------------------------
def blocks_for(genes, pos):
    g = [x for x in genes if x in pos]
    g = sorted(g, key=lambda x: pos[x])
    idx = np.array_split(np.arange(len(g)), N_BLOCKS)
    return [set(np.array(g)[i]) for i in idx]


def jackknife(full, parts):
    parts = np.asarray(parts, dtype=float)
    B = len(parts)
    m = parts.mean()
    se = math.sqrt((B - 1) / B * np.sum((parts - m) ** 2))
    return full, se


def r3(covf_focus, keys):
    say("")
    say("-" * 72)
    say("[SECTION] R3 delete-one-block jackknife (%d contiguous gene blocks)" % N_BLOCKS)
    pos = gene_positions()
    inv = {v: k for k, v in keys.items()}
    jobs = [("G4", covf_focus, "Average", [inv[t] for t in JK_TYPES if t in inv])]
    for lin in ("lymphocyte", "epithelial"):
        f = COVAR_G1 / ("LES_vs_EU__%s.covar.txt" % lin)
        if f.exists():
            jobs.append(("G1_" + lin, f, "AVG", ["D"]))
    store = {}
    t0 = time.time()
    for name, covf, cond, vars_ in jobs:
        cov = pd.read_csv(covf, sep="\t")
        cov["GENE"] = cov["GENE"].astype(str)
        blocks = blocks_for(list(cov["GENE"]), pos)
        for trait in JK_TRAITS:
            full = magma_gp(raw_of(trait), covf, cond, "R3_%s_%s_full" % (name, trait))
            if not full:
                continue
            parts = {v: [] for v in vars_}
            parts_std = {v: [] for v in vars_}
            for bi, bl in enumerate(blocks):
                f = WORK / "R3_block.covar.txt"
                cov[~cov["GENE"].isin(bl)].to_csv(f, sep="\t", index=False, float_format="%.6g")
                r = magma_gp(raw_of(trait), f, cond, "R3_blk")
                for v in vars_:
                    parts[v].append(r[v]["beta"] if r and v in r else np.nan)
                    parts_std[v].append(r[v]["beta_std"] if r and v in r else np.nan)
            store[(name, trait)] = {"full": full, "parts": parts, "parts_std": parts_std}
        say("    [R3] %-14s done (%.0f min)" % (name, (time.time() - t0) / 60))

    def jk_line(tag, full, parts):
        parts = np.array(parts, dtype=float)
        ok = np.isfinite(parts)
        _, se = jackknife(full, parts[ok])
        z = full / se if se > 0 else np.nan
        return "%s %.4f (jk se %.4f, z %.2f)" % (tag, full, se, z), se

    rows = []
    for (name, trait), d in store.items():
        for v, parts in d["parts"].items():
            line, se = jk_line(v if name == "G4" else name, d["full"][v]["beta"], parts)
            label = keys.get(v, v)
            say("    [R3_EST] %-14s %-13s %-14s %s | MAGMA se %.4f" % (name, trait[:13], label[:14], line, d["full"][v]["se"]))
            rows.append({"analysis": name, "trait": trait, "var": label, "beta": d["full"][v]["beta"], "se_magma": d["full"][v]["se"], "se_jk": se})
    # within-trait: dStromal_mid - ePV_2 (standardized)
    if "dStromal_mid" in inv and "ePV_2" in inv:
        a, b = inv["dStromal_mid"], inv["ePV_2"]
        for trait in JK_TRAITS:
            d = store.get(("G4", trait))
            if not d:
                continue
            full = d["full"][a]["beta_std"] - d["full"][b]["beta_std"]
            parts = np.array(d["parts_std"][a]) - np.array(d["parts_std"][b])
            _, se = jackknife(full, parts[np.isfinite(parts)])
            say("    [R3_WITHIN] %-13s BETA_STD(dStromal_mid) - BETA_STD(ePV_2) = %.4f (jk se %.4f, z %.2f)" % (trait[:13], full, se, full / se if se else np.nan))
    # cross-trait differences, block-paired
    pairs = [("KOLLER_EUR", "HYPERTENSION"), ("KOLLER_EUR", "T2D"), ("GCST90205183", "HYPERTENSION"), ("GCST90205183", "T2D")]
    for name in sorted(set(k[0] for k in store)):
        vars_ = list(store[next(k for k in store if k[0] == name)]["parts"].keys())
        for v in vars_:
            for t1, t2 in pairs:
                d1, d2 = store.get((name, t1)), store.get((name, t2))
                if not d1 or not d2:
                    continue
                full = d1["full"][v]["beta"] - d2["full"][v]["beta"]
                parts = np.array(d1["parts"][v]) - np.array(d2["parts"][v])
                _, se = jackknife(full, parts[np.isfinite(parts)])
                z = full / se if se else np.nan
                p = math.erfc(abs(z) / math.sqrt(2)) if np.isfinite(z) else np.nan
                say("    [R3_CROSS] %-14s %-14s %s - %s = %+.4f (jk se %.4f, z %+.2f, p %.3f)" % (
                    name, keys.get(v, v)[:14], t1[:6], t2[:6], full, se, z, p))
                rows.append({"analysis": name, "trait": t1 + "-" + t2, "var": keys.get(v, v), "beta": full, "se_jk": se, "z": z, "p": p})
    pd.DataFrame(rows).to_csv(OUT / "R3_jackknife.csv", index=False)


# ----------------------------------------------------------------------------
# R4 simulations
# ----------------------------------------------------------------------------
def read_raw(trait):
    raw = raw_of(trait)
    lines = raw.read_text().splitlines()
    out = pd.read_csv(MAGMA_OUT / (trait + "_35_10.genes.out"), sep=r"\s+")
    zmap = dict(zip(out["GENE"].astype(str), out["ZSTAT"]))
    data_idx = [i for i, l in enumerate(lines) if l.strip() and not l.startswith("#")]
    votes = {}
    for i in data_idx[:300]:
        tok = lines[i].split()
        z = zmap.get(tok[0])
        if z is None:
            continue
        for j in range(1, min(len(tok), 12)):
            try:
                if abs(float(tok[j]) - z) < 1e-3:
                    votes[j] = votes.get(j, 0) + 1
            except ValueError:
                pass
    zi = max(votes, key=votes.get)
    genes = [lines[i].split()[0] for i in data_idx]
    z = np.array([float(lines[i].split()[zi]) for i in data_idx])
    return lines, data_idx, zi, genes, z


def write_raw(lines, data_idx, zi, znew, path):
    out = list(lines)
    for k, i in enumerate(data_idx):
        t2 = out[i].split()
        t2[zi] = "%.6g" % znew[k]
        out[i] = " ".join(t2)
    Path(path).write_text("\n".join(out) + "\n")


def r4(S, sym2):
    say("")
    say("-" * 72)
    say("[SECTION] R4 simulations with a known source of signal (null = KOLLER gene z circularly shifted)")
    if not COV_M1.exists():
        say("[R4] M1 covariate file missing")
        return
    lines, data_idx, zi, genes, z0 = read_raw("KOLLER_EUR")
    say("    [R4] genes.raw: %d genes; Z column index %d" % (len(genes), zi))
    types = list(S.columns)
    lab = {"T%02d" % (i + 1): t for i, t in enumerate(types)}
    ent = {sym2[g]: g for g in S.index if g in sym2}
    sym_of = np.array([ent.get(g) for g in genes], dtype=object)
    has = np.array([s is not None for s in sym_of])
    sd_ds = float(S["dStromal_mid"].std())
    b_std = OBS_BETA_DSMID * sd_ds
    say("    [R4] observed KOLLER dStromal_mid beta %.3f x sd(S) %.4f -> standardized effect %.4f per SD" % (OBS_BETA_DSMID, sd_ds, b_std))
    rng = np.random.default_rng(SEED)
    scen = [("null", None, 0.0)] + [(t, t, s) for t in SIM_TARGETS if t in S.columns for s in SIM_STRENGTH]
    rows = []
    t0 = time.time()
    simraw = WORK / "SIM_35_10.genes.raw"
    for name, target, strength in scen:
        for rep in range(N_SIM):
            off = int(rng.integers(500, len(z0) - 500))
            zs = np.roll(z0, off).copy()
            if target is not None:
                s = np.zeros(len(genes))
                v = S.loc[[g for g in sym_of[has]], target].to_numpy(dtype=float)
                s[has] = (v - v.mean()) / v.std()
                zs = zs + strength * b_std * s
            write_raw(lines, data_idx, zi, zs, simraw)
            r = magma_gp(simraw, COV_M1, "Average", "R4_sim")
            if not r:
                continue
            df = pd.DataFrame([{"type": lab[k], **r[k]} for k in lab if k in r]).sort_values("p").reset_index(drop=True)
            ds = df[df["type"].isin(DS3)]
            h3 = dict(zip(ds["type"], holm(ds["p"].values)))
            q = dict(zip(df["type"], bh(df["p"].values)))
            row = {"scenario": name, "strength": strength, "rep": rep,
                   "dsmid_top": df["type"].iloc[0] == "dStromal_mid",
                   "dsmid_p": float(df.loc[df["type"] == "dStromal_mid", "p"].iloc[0]),
                   "dsmid_holm3": h3.get("dStromal_mid", np.nan), "dsmid_q32": q.get("dStromal_mid", np.nan)}
            if target is not None:
                row["true_top"] = df["type"].iloc[0] == target
                row["true_p"] = float(df.loc[df["type"] == target, "p"].iloc[0])
                row["true_rank"] = int(df.index[df["type"] == target][0]) + 1
            rows.append(row)
        say("    [R4] scenario %-14s x%.0f done (%.0f min)" % (name, strength, (time.time() - t0) / 60))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "R4_simulations.csv", index=False)
    for (name, strength), g in R.groupby(["scenario", "strength"], sort=False):
        line = "    [R4_SUM] truth=%-13s x%.0f n=%d | dS_mid top %.2f, p<0.05 %.2f, Holm3<0.05 %.2f, q32<0.05 %.2f" % (
            name, strength, len(g), g["dsmid_top"].mean(), (g["dsmid_p"] < 0.05).mean(), (g["dsmid_holm3"] < 0.05).mean(),
            (g["dsmid_q32"] < 0.05).mean())
        if "true_top" in g and g["true_top"].notna().any():
            line += " | true type top %.2f, p<0.05 %.2f, median rank %d" % (g["true_top"].mean(), (g["true_p"] < 0.05).mean(), int(g["true_rank"].median()))
        say(line)


# ----------------------------------------------------------------------------
# R5 null splits
# ----------------------------------------------------------------------------
def r5(sym2):
    say("")
    say("-" * 72)
    say("[SECTION] R5 false-positive rate under a true null (random patient splits within one context; KOLLER_EUR)")
    M = np.load(PB213 / "pseudobulk_counts.npy", mmap_mode="r")
    G = pd.read_csv(PB213 / "pseudobulk_groups.csv").reset_index(drop=True)
    genes = list(pd.read_csv(PB213 / "genes.csv")["symbol"].astype(str))
    mapped = np.array([x in sym2 for x in genes])
    rng = np.random.default_rng(SEED + 1)
    rows = []
    for ctx, lin in R5_SETS:
        sub = G[(G["lineage_ref"] == lin) & (G["context"] == ctx) & (G["agree"])]
        prof = []
        for pid, g in sub.groupby("patient"):
            if int(g["n_cells"].sum()) >= 10:
                v = M[g.index.to_numpy()].sum(axis=0)
                prof.append(np.log2(v / max(v.sum(), 1) * 1e6 + 1))
        n = len(prof)
        if n < 6:
            say("    [R5] %s %s: %d patients -> skip" % (ctx, lin, n))
            continue
        P = np.vstack(prof)
        t0 = time.time()
        for s in range(N_SPLIT):
            perm = rng.permutation(n)
            a, b = perm[: n // 2], perm[n // 2:]
            A, B = P[a].mean(axis=0), P[b].mean(axis=0)
            ex = (np.maximum(2 ** A - 1, 2 ** B - 1) >= 1.0) & mapped
            g = np.asarray(genes)[ex]
            f = WORK / "R5.covar.txt"
            pd.DataFrame({"GENE": [sym2[x] for x in g], "D": (A - B)[ex], "DN": (B - A)[ex], "AVG": ((A + B) / 2)[ex]}).to_csv(
                f, sep="\t", index=False, float_format="%.6g")
            r = magma_gp(raw_of("KOLLER_EUR"), f, "AVG", "R5")
            if r and "D" in r:
                rows.append({"context": ctx, "lineage": lin, "n": n, "split": s, "p_greater": r["D"]["p"], "p_less": r["DN"]["p"],
                             "t": r["D"]["beta"] / r["D"]["se"]})
        g = pd.DataFrame([x for x in rows if x["context"] == ctx and x["lineage"] == lin])
        if len(g):
            say("    [R5] %-17s %-10s patients=%d splits=%d | P(p_greater<0.05)=%.3f P(p_less<0.05)=%.3f | t sd=%.2f (expected 1) | %.0f min" % (
                ctx, lin, n, len(g), (g["p_greater"] < 0.05).mean(), (g["p_less"] < 0.05).mean(), g["t"].std(), (time.time() - t0) / 60))
    pd.DataFrame(rows).to_csv(OUT / "R5_null_splits.csv", index=False)


# ----------------------------------------------------------------------------
# B0 / K0
# ----------------------------------------------------------------------------
def b0():
    say("")
    say("-" * 72)
    say("[SECTION] B0 fallopian-tube single-cell atlases on CELLxGENE (listing only)")
    p = subprocess.run(["curl", "-sSL", "--fail", "--max-time", "120", "-A", "paper9", CXG_API], capture_output=True)
    if p.returncode != 0:
        say("[B0] CELLxGENE API not reachable")
        return
    cols = json.loads(p.stdout.decode("utf-8", "replace"))
    hits = [c for c in cols if re.search(r"fallopian|oviduct|uterine tube", (c.get("name", "") + " " + c.get("description", "")), re.I)]
    say("[B0] collections matching: %d" % len(hits))
    for c in hits[:10]:
        say("  [B0_COLL] %s | %s" % (c.get("collection_id", c.get("id")), c.get("name", "")[:150]))
        q = subprocess.run(["curl", "-sSL", "--fail", "--max-time", "120", "-A", "paper9", CXG_API + "/" + str(c.get("collection_id", c.get("id")))], capture_output=True)
        if q.returncode != 0:
            continue
        det = json.loads(q.stdout.decode("utf-8", "replace"))
        for ds in det.get("datasets", [])[:12]:
            tissues = ", ".join(sorted(set(t.get("label", "") for t in ds.get("tissue", []))))[:80]
            h5 = [a for a in ds.get("assets", []) if a.get("filetype", "").upper() == "H5AD"]
            say("      [B0_DS] %s | cells=%s | %s | tissue: %s | h5ad %.0f MB" % (
                ds.get("dataset_id", ""), ds.get("cell_count", "?"), (ds.get("title", "") or "")[:70], tissues,
                (h5[0].get("filesize", 0) / 1e6) if h5 else 0))


def read_xlsx_sheets(path, pattern):
    z = zipfile.ZipFile(path)
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
          "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root.findall("m:si", ns):
            shared.append("".join(t.text or "" for t in si.iter("{%s}t" % ns["m"])))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target = {r.get("Id"): r.get("Target") for r in rels}
    out = {}
    for sh_ in wb.find("m:sheets", ns):
        name = sh_.get("name")
        if not re.search(pattern, name, re.I):
            continue
        rid = sh_.get("{%s}id" % ns["r"])
        tpath = "xl/" + target[rid].lstrip("/").replace("xl/", "")
        root = ET.fromstring(z.read(tpath))
        rows = []
        for row in root.iter("{%s}row" % ns["m"]):
            vals = []
            for c in row.findall("m:c", ns):
                v = c.find("m:v", ns)
                txt = v.text if v is not None else ""
                if c.get("t") == "s" and txt:
                    txt = shared[int(txt)]
                elif c.get("t") == "inlineStr":
                    txt = "".join(t.text or "" for t in c.iter("{%s}t" % ns["m"]))
                vals.append(txt or "")
            rows.append(vals)
        out[name] = rows
    return out


def k0():
    say("")
    say("-" * 72)
    say("[SECTION] K0 Koller 2026 Supplementary Tables 9-10 (from references/, if present)")
    files = sorted(REFS.rglob("*.xlsx")) if REFS.exists() else []
    if not files:
        say("[K0] no .xlsx under " + str(REFS))
        return
    for f in files:
        try:
            sheets = read_xlsx_sheets(f, r"(^|\D)(9|10)(\D|$)")
        except Exception as e:
            say("[K0] %s unreadable (%s)" % (f.name, type(e).__name__))
            continue
        for name, rows in sheets.items():
            say("  [K0_SHEET] %s :: %s (%d rows)" % (f.name[:40], name, len(rows)))
            for r in rows[:70]:
                say("      " + " | ".join(x[:28] for x in r[:12]))
            pd.DataFrame(rows).to_csv(OUT / ("K0_%s.csv" % re.sub(r"\W+", "_", name)), index=False, header=False)


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.27 | Option A core: resolution, jackknife, simulations, null calibration | " + TS)
    say("=" * 72)
    t0 = time.time()
    sym2 = entrez_unique()
    S, M = heca_tables()
    covf, keys = None, {}
    steps = [("K0", k0), ("B0", b0)]
    for name, fn in steps:
        try:
            fn()
        except Exception as e:
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
    try:
        covf, keys = r1_r2(S, M, sym2)
    except Exception as e:
        import traceback
        say("[STEP_ERROR] R1/R2 %s: %s" % (type(e).__name__, str(e)[:300]))
        say(traceback.format_exc()[-800:])
    for name, fn in (("R4", lambda: r4(S, sym2)), ("R5", lambda: r5(sym2)), ("R3", lambda: r3(covf, keys) if covf else None)):
        try:
            fn()
        except Exception as e:
            import traceback
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-800:])
    with open(OUT / "p9_v027_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.17 before this run (PLAN v0.3)")
    say("[DONE] PAPER9_V027_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
