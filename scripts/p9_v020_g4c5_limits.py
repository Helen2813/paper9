#!/usr/bin/env python3
# Paper 9 | v0.0.20 | rules CHANGELOG v0.2.11 (fixed before this run):
#   (A) Fonseca et al. 2023 (GSE213216) metadata check: authors' GEO-to-patient/class table vs ours,
#       repository file list (per-cell annotations?)
#   (B) G4 check 5: annotation design (specificity | Average vs mean expression | Average vs mean
#       expression without conditioning) on HECA control-donor pseudobulk
#   (C) Detection limits (MDE, 80% power) for G4 and G1
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). ~10-20 min. No large downloads.

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
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------
# Hardcoded locations
# ----------------------------------------------------------------------------
HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
INTERIM = REPO / "data" / "interim"
PROC = REPO / "data" / "processed"
EXT = REPO / "data" / "external"
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
HECA_GS = PROC / "genesets" / "heca_v014"
COVAR_G1 = PROC / "magma_covar" / "g1"
PB_GROUPS = PROC / "gse213216_pseudobulk_clusters" / "pseudobulk_groups.csv"
FONSECA_DIR = EXT / "fonseca2023"
FONSECA_TABLE_URL = "https://raw.githubusercontent.com/lawrenson-lab/AtlasEndometriosis/main/GEO_id_patient_class.csv"
FONSECA_TREE_URL = "https://api.github.com/repos/lawrenson-lab/AtlasEndometriosis/git/trees/main?recursive=1"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "g4c5_limits" / ("p9_v020_" + TS)
WORK = INTERIM / "magma_v020"
LOG = REPO / "logs" / ("p9_v020_g4c5_limits_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Pre-registered (CHANGELOG v0.2.11)
DSTROMAL = ["dStromal_early", "dStromal_mid", "dStromal_late"]
C5_TRAITS = ["GCST90205183", "HYPERTENSION", "LEIOMYOMA", "T2D"]
VASCULAR = ["Arterial", "Venous", "ePV_1a", "ePV_1b", "ePV_2", "mPV", "uSMCs", "Lymphatic"]
LINEAGES = ["stromal", "epithelial", "endothelial", "myeloid", "lymphocyte"]
Z80 = 0.8416212335729143

for d in (OUT, WORK, LOG.parent, FONSECA_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V020_G4C5_LIMITS", "timestamp_utc": TS, "gates": {}, "results": {}}


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


def curl(url, dst=None, max_time=120):
    cmd = ["curl", "-sSL", "--fail", "--max-time", str(max_time), "-A", "paper9"]
    if dst:
        cmd += ["-o", str(dst)]
    cmd.append(url)
    p = subprocess.run(cmd, capture_output=True)
    if p.returncode != 0:
        return None
    return True if dst else p.stdout.decode("utf-8", "replace")


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


def norm_ppf(q):
    # Acklam's rational approximation of the standard normal quantile (no scipy needed)
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02,
         -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01,
         -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00]
    pl = 0.02425
    if q < pl:
        x = np.sqrt(-2 * np.log(q))
        return (((((c[0] * x + c[1]) * x + c[2]) * x + c[3]) * x + c[4]) * x + c[5]) / ((((d[0] * x + d[1]) * x + d[2]) * x + d[3]) * x + 1)
    if q > 1 - pl:
        x = np.sqrt(-2 * np.log(1 - q))
        return -(((((c[0] * x + c[1]) * x + c[2]) * x + c[3]) * x + c[4]) * x + c[5]) / ((((d[0] * x + d[1]) * x + d[2]) * x + d[3]) * x + 1)
    x = q - 0.5
    r = x * x
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * x / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)


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


def magma_gp(trait, covf, cond, tag):
    # Gene-property analysis, direction=greater; cond=None -> no conditioning covariate
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    if not raw.exists():
        return None, "no genes.raw"
    outp = WORK / ("%s__%s" % (trait, tag))
    gsa = Path(str(outp) + ".gsa.out")
    keys = ["condition-hide=", "condition="] if cond else [None]
    err = ""
    for key in keys:
        if gsa.exists():
            gsa.unlink()
        model = ["direction=greater"] if key is None else [key + cond, "direction=greater"]
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--gene-covar", covf, "--model"] + model + ["--out", outp],
                      timeout=3 * 3600)
        if rc == 0 and gsa.exists():
            rows = [l for l in open(gsa) if not l.startswith("#") and l.strip()]
            hdr = rows[0].split()
            res = {}
            for l in rows[1:]:
                d = dict(zip(hdr, l.split()))
                name = d.get("FULL_NAME", d.get("VARIABLE"))
                try:
                    res[name] = {"beta": float(d["BETA"]), "beta_std": float(d.get("BETA_STD", "nan")),
                                 "se": float(d["SE"]), "p": float(d["P"]), "ngenes": int(d["NGENES"])}
                except (KeyError, ValueError):
                    continue
            return res, ""
        err = (o + e).strip()[-300:]
    say("    [MAGMA_FAIL] %s %s %s" % (trait, tag, err))
    return None, err


# ----------------------------------------------------------------------------
# A. Fonseca et al. 2023 metadata
# ----------------------------------------------------------------------------
def fonseca_check():
    say("")
    say("-" * 72)
    say("[SECTION] A. GSE213216 = Fonseca et al. 2023 (Nat Genet): authors' sample table vs ours")
    tree = curl(FONSECA_TREE_URL)
    if tree:
        try:
            items = json.loads(tree).get("tree", [])
            say("[FONSECA_REPO] %d entries" % len(items))
            for it in items[:80]:
                if it.get("type") == "blob":
                    say("    [FILE] %-70s %8.1f KB" % (it["path"][:70], it.get("size", 0) / 1024))
            flag = [it["path"] for it in items if it.get("type") == "blob" and any(
                k in it["path"].lower() for k in ("meta", "annot", "barcode", "cell", "label"))]
            say("[FONSECA_REPO_ANNOT_CANDIDATES] " + (", ".join(flag[:20]) if flag else "none"))
        except Exception as e:
            say("[FONSECA_REPO] unexpected response %s" % type(e).__name__)
    else:
        say("[FONSECA_REPO] GitHub API not reachable")
    dst = FONSECA_DIR / "GEO_id_patient_class.csv"
    if not (dst.exists() and dst.stat().st_size > 0):
        if not curl(FONSECA_TABLE_URL, dst):
            say("[FONSECA_TABLE] download failed")
            return
    T = pd.read_csv(dst, sep=None, engine="python")
    say("[FONSECA_TABLE] rows=%d columns=%s" % (len(T), ", ".join(map(str, T.columns))))
    for line in T.head(6).to_string(index=False).splitlines():
        say("    " + line)
    gcol = next((c for c in T.columns if T[c].astype(str).str.startswith("GSM").mean() > 0.5), None)
    if gcol is None or not PB_GROUPS.exists():
        say("[FONSECA_TABLE] no GSM column or no v016 groups file; comparison skipped")
        return
    G = pd.read_csv(PB_GROUPS)
    ours = G.groupby("gsm").agg(patient=("patient", "first"), context=("context", "first")).reset_index()
    other = [c for c in T.columns if c != gcol]
    M = ours.merge(T.rename(columns={gcol: "gsm"}), on="gsm", how="outer", indicator=True)
    say("[FONSECA_MATCH] GSM in both=%d | only ours=%d | only authors=%d" % (
        (M["_merge"] == "both").sum(), (M["_merge"] == "left_only").sum(), (M["_merge"] == "right_only").sum()))
    B = M[M["_merge"] == "both"]
    for c in other:
        vals = B[c].astype(str)
        if vals.nunique() <= 12:
            tab = pd.crosstab(B["context"], vals)
            say("[FONSECA_XTAB] our context x authors '%s'" % c)
            for line in tab.to_string().splitlines():
                say("    " + line)
        # patient grouping agreement: same partition of GSMs?
        if 3 <= vals.nunique() <= 60 and vals.nunique() != len(B):
            a = B.groupby(c)["patient"].nunique()
            b = B.groupby("patient")[c].nunique()
            say("[FONSECA_PATIENT] column '%s': authors' groups mapping to >1 of our patients=%d; our patients mapping to >1 authors' group=%d" % (
                c, int((a > 1).sum()), int((b > 1).sum())))
    B.drop(columns="_merge").to_csv(OUT / "fonseca_sample_table_vs_ours.csv", index=False)


# ----------------------------------------------------------------------------
# B. G4 check 5: annotation design
# ----------------------------------------------------------------------------
def c5_covars(sym2):
    S = pd.read_csv(HECA_GS / "specificity_by_type.csv.gz", index_col=0)
    M = pd.read_csv(HECA_GS / "mean_cpm_by_type.csv.gz", index_col=0).loc[S.index, S.columns]
    types = list(S.columns)
    keys = ["T%02d" % (i + 1) for i in range(len(types))]
    lab = dict(zip(keys, types))
    gene = [sym2.get(g) for g in S.index]
    avg = np.log2(M.mean(axis=1).values + 1.0)
    f1 = WORK / "C5_M1_specificity.covar.txt"
    d1 = pd.DataFrame(S.values, columns=keys)
    d1["Average"] = avg
    d1.insert(0, "GENE", gene)
    d1.dropna(subset=["GENE"]).to_csv(f1, sep="\t", index=False, float_format="%.6g")
    f2 = WORK / "C5_M23_logmean.covar.txt"
    d2 = pd.DataFrame(np.log2(M.values + 1.0), columns=keys)
    d2["Average"] = avg
    d2.insert(0, "GENE", gene)
    d2.dropna(subset=["GENE"]).to_csv(f2, sep="\t", index=False, float_format="%.6g")
    say("[C5_COVAR] HECA control-donor types=%d genes=%d" % (len(types), int(d1["GENE"].notna().sum())))
    return {"M1": (f1, "Average"), "M2": (f2, "Average"), "M3": (f2, None)}, lab


def c5_table(res, lab):
    df = pd.DataFrame([{"key": k, "type": lab[k], **res[k]} for k in lab if k in res])
    df = df.sort_values("p").reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    m = df["type"].isin(DSTROMAL)
    df["holm_dS"] = np.nan
    df.loc[m, "holm_dS"] = holm(df.loc[m, "p"].values)
    return df


def g4_check5(sym2):
    say("")
    say("-" * 72)
    say("[SECTION] B. G4 CHECK 5 (v0.2.11): annotation design | M1 specificity|Average, M2 logmean|Average, M3 logmean")
    models, lab = c5_covars(sym2)
    tabs = {}
    for trait in C5_TRAITS:
        for mname, (covf, cond) in models.items():
            res, err = magma_gp(trait, covf, cond, "C5_" + mname)
            if not res:
                continue
            df = c5_table(res, lab)
            df["trait"], df["model"] = trait, mname
            tabs[(trait, mname)] = df
            ds = df[df["type"].isin(DSTROMAL)].sort_values("p")
            say("[C5] %-13s %s top3: %s | dStromal: %s" % (
                trait[:13], mname, "; ".join("%s(p=%.1e)" % (t, p) for t, p in zip(df["type"].head(3), df["p"].head(3))),
                "; ".join("%s r%d p=%.2g H=%.2g" % (t.replace("dStromal_", "dS_"), r, p, h)
                          for t, r, p, h in zip(ds["type"], ds["rank"], ds["p"], ds["holm_dS"]))))
    if tabs:
        pd.concat(tabs.values(), ignore_index=True).to_csv(OUT / "G4_check5_results.csv", index=False)
    g = "GCST90205183"
    if (g, "M1") in tabs and (g, "M3") in tabs and ("HYPERTENSION", "M3") in tabs:
        m1 = bool((tabs[(g, "M1")]["holm_dS"] < 0.05).any())
        m3 = bool((tabs[(g, "M3")]["holm_dS"] < 0.05).any())
        shared = len(set(tabs[(g, "M3")]["type"].head(3)) & set(tabs[("HYPERTENSION", "M3")]["type"].head(3)))
        v = "SUPPORTED" if (m3 and not m1 and shared >= 2) else "NOT SUPPORTED"
        say("[C5_VERDICT] annotation design explains discrepancy: %s (dStromal Holm<0.05: M1=%s M3=%s; M3 top-3 shared GCST/HTN=%d)" % (
            v, m1, m3, shared))
        SUM["results"]["G4_check5"] = {"verdict": v, "M1_dS_sig": m1, "M3_dS_sig": m3, "M3_top3_shared_with_HTN": shared}
    return tabs


# ----------------------------------------------------------------------------
# C. Detection limits
# ----------------------------------------------------------------------------
def mde_line(label, r, k):
    zc = norm_ppf(1 - 0.05 / k)
    mde = (zc + Z80) * r["se"]
    ratio = r["beta_std"] / r["beta"] if r["beta"] not in (0.0,) and np.isfinite(r["beta_std"]) else np.nan
    lo, hi = r["beta"] - 1.96 * r["se"], r["beta"] + 1.96 * r["se"]
    say("    [%s] beta=%8.4f (95%% CI %8.4f, %8.4f) se=%.4f | MDE80=%.4f (BETA_STD %.3f; observed BETA_STD %.3f)" % (
        label, r["beta"], lo, hi, r["se"], mde, mde * ratio if np.isfinite(ratio) else np.nan, r["beta_std"]))
    return {"label": label, "beta": r["beta"], "se": r["se"], "ci_lo": lo, "ci_hi": hi, "mde80": mde,
            "mde80_std": mde * ratio if np.isfinite(ratio) else np.nan, "beta_std": r["beta_std"], "k": k}


def limits(tabs):
    say("")
    say("-" * 72)
    say("[SECTION] C. DETECTION LIMITS (80% power, one-sided alpha 0.05/k)")
    rows = []
    t = tabs.get(("GCST90205183", "M1"))
    if t is not None:
        say("[MDE_G4] GCST90205183, HECA M1 (k=3 dStromal types)")
        for _, r in t[t["type"].isin(DSTROMAL)].iterrows():
            rows.append({"analysis": "G4", **mde_line(r["type"], r, 3)})
    pc = tabs.get(("LEIOMYOMA", "M1"))
    if pc is not None and (pc["type"] == "uSMCs").any():
        r = pc[pc["type"] == "uSMCs"].iloc[0]
        say("[MDE_PC] LEIOMYOMA uSMCs (own GWAS): beta=%.4f se=%.4f BETA_STD=%.3f p=%.1e" % (r["beta"], r["se"], r["beta_std"], r["p"]))
        rows.append({"analysis": "PC_leiomyoma_uSMCs", "label": "uSMCs", "beta": r["beta"], "se": r["se"], "beta_std": r["beta_std"]})
    pc = tabs.get(("HYPERTENSION", "M1"))
    if pc is not None:
        v = pc[pc["type"].isin(VASCULAR)].sort_values("p")
        if len(v):
            r = v.iloc[0]
            say("[MDE_PC] HYPERTENSION best vascular %s (own GWAS): beta=%.4f se=%.4f BETA_STD=%.3f p=%.1e" % (
                r["type"], r["beta"], r["se"], r["beta_std"], r["p"]))
            rows.append({"analysis": "PC_hypertension_vascular", "label": r["type"], "beta": r["beta"], "se": r["se"],
                         "beta_std": r["beta_std"]})
    say("[MDE_G1] GCST90205183, lesion vs eutopic (k=5 lineages; MDE applies to either direction)")
    for lin in LINEAGES:
        f = COVAR_G1 / ("LES_vs_EU__%s.covar.txt" % lin)
        if not f.exists():
            continue
        res, _ = magma_gp("GCST90205183", f, "AVG", "MDE_G1_" + lin)
        if res and "D" in res:
            rows.append({"analysis": "G1", **mde_line(lin, res["D"], 5)})
    pd.DataFrame(rows).to_csv(OUT / "detection_limits.csv", index=False)
    SUM["results"]["limits_rows"] = len(rows)


# ----------------------------------------------------------------------------
def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.20 | Fonseca metadata + G4 check 5 (annotation design) + detection limits | " + TS)
    say("=" * 72)
    t0 = time.time()
    sym2 = entrez_unique()
    for name, fn in (("fonseca", fonseca_check),):
        try:
            fn()
        except Exception as e:
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
    tabs = {}
    try:
        tabs = g4_check5(sym2)
    except Exception as e:
        import traceback
        say("[STEP_ERROR] g4_check5 %s: %s" % (type(e).__name__, str(e)[:300]))
        say(traceback.format_exc()[-1200:])
    try:
        limits(tabs)
    except Exception as e:
        import traceback
        say("[STEP_ERROR] limits %s: %s" % (type(e).__name__, str(e)[:300]))
        say(traceback.format_exc()[-1200:])
    with open(OUT / "p9_v020_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.11 before this run")
    say("[DONE] PAPER9_V020_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
