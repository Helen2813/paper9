#!/usr/bin/env python3
# Paper 9 | v0.0.22 | X4 organ-identity hypothesis (exploratory, rules CHANGELOG v0.2.13, fixed before this run).
#   X4a Muellerian-tract expression beyond smooth muscle; X4b by HECA within-endometrium breadth tertile;
#   X4c HECA dStromal specificity conditional on Muellerian expression.
# Runner: if the v021 audit has not finished, it is started first (it runs v019 / v020 itself if needed).
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). X4 itself ~5-15 min.

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
import subprocess
import sys
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
VENV_PY = REPO / ".venv" / "bin" / "python"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
GTEX_FILE = EXT / "gtex" / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
HECA_GS = PROC / "genesets" / "heca_v014"
LOGS = REPO / "logs"
SCRIPTS = REPO / "scripts"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "x4" / ("p9_v022_x4_" + TS)
WORK = INTERIM / "magma_v022"
LOG = LOGS / ("p9_v022_x4_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Pre-registered (CHANGELOG v0.2.13)
MULL = ["Uterus", "Fallopian Tube", "Cervix - Ectocervix", "Cervix - Endocervix", "Vagina"]
SMOOTH = ["Colon - Sigmoid", "Esophagus - Muscularis", "Bladder", "Artery - Tibial"]
TRAITS = ["GCST90205183", "ENDO_ALL", "ENDO_OVARY", "ENDO_PERITONEUM", "ENDO_DEEP", "LEIOMYOMA", "HYPERTENSION", "T2D"]
DSTROMAL = ["dStromal_mid", "dStromal_early"]

_fh = None
SUM = {"task": "PAPER9_V022_X4", "timestamp_utc": TS, "results": {}}


def say(msg=""):
    print(msg, flush=True)
    if _fh:
        _fh.write(msg + "\n")
        _fh.flush()


def sh(cmd, timeout=None):
    try:
        p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout, env=ENV1)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"


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
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    if not raw.exists():
        return None
    outp = WORK / ("%s__%s" % (trait, tag))
    gsa = Path(str(outp) + ".gsa.out")
    o = e = ""
    for key in ("condition-hide=", "condition="):
        if gsa.exists():
            gsa.unlink()
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--gene-covar", covf, "--model", key + cond,
                              "direction=greater", "--out", outp], timeout=3 * 3600)
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
            return res
    say("    [MAGMA_FAIL] %s %s %s" % (trait, tag, (o + e).strip()[-200:]))
    return None


# ----------------------------------------------------------------------------
# Runner
# ----------------------------------------------------------------------------
def run_previous():
    logs = sorted(LOGS.glob("p9_v021_audit_*.log"), key=lambda p: p.stat().st_mtime)
    if any("PAPER9_V021_COMPLETE" in p.read_text(errors="ignore") for p in logs):
        print("[RUNNER] v021 audit already complete", flush=True)
        return
    if logs and time.time() - logs[-1].stat().st_mtime < 3 * 3600:
        print("[RUNNER] v021 has a recent unfinished log (%s): it may still be running; start this again after it ends." % logs[-1].name, flush=True)
        sys.exit(3)
    print("[RUNNER] running v021 audit first (it also runs v019 / v020 if needed) ...", flush=True)
    rc = subprocess.run([str(VENV_PY), str(SCRIPTS / "p9_v021_audit.py")], cwd=str(REPO)).returncode
    print("[RUNNER] v021 exit %d" % rc, flush=True)


# ----------------------------------------------------------------------------
# X4
# ----------------------------------------------------------------------------
def build_tables(sym2):
    df = pd.read_csv(GTEX_FILE, sep="\t", skiprows=2, compression="gzip")
    tissues = [c for c in df.columns[2:] if not c.startswith("Cells")]
    df = df[df["Description"].isin(sym2.keys())].copy()
    df["_m"] = df[tissues].astype(float).mean(axis=1)
    df = df.sort_values("_m", ascending=False).drop_duplicates("Description")
    tpm = df[tissues].astype(float)
    keep = (tpm.values >= 1.0).any(axis=1)
    L = np.log2(tpm[keep] + 1.0)
    L.index = df.loc[keep, "Description"].values
    miss = [t for t in MULL + SMOOTH if t not in tissues]
    if miss:
        say("[X4_WARN] GTEx tissues not found: " + ", ".join(miss))
    G = pd.DataFrame({"MULL": L[[t for t in MULL if t in tissues]].mean(axis=1),
                      "UTERUS": L["Uterus"] if "Uterus" in tissues else np.nan,
                      "SMOOTH": L[[t for t in SMOOTH if t in tissues]].mean(axis=1),
                      "Average": L.mean(axis=1)})
    S = pd.read_csv(HECA_GS / "specificity_by_type.csv.gz", index_col=0)
    M = pd.read_csv(HECA_GS / "mean_cpm_by_type.csv.gz", index_col=0).loc[S.index, S.columns]
    LH = np.log2(M + 1.0)
    mx = LH.max(axis=1)
    tau = ((1.0 - LH.div(mx, axis=0)).sum(axis=1)) / (LH.shape[1] - 1)
    H = pd.DataFrame({"TAU_H": tau, "AVG_H": np.log2(M.mean(axis=1) + 1.0)})
    for t in DSTROMAL:
        if t in S.columns:
            H["S_" + t] = S[t]
    say("[X4_TABLES] GTEx genes=%d | HECA genes=%d | shared=%d | Muellerian tissues=%d smooth=%d" % (
        len(G), len(H), len(G.index.intersection(H.index)), len([t for t in MULL if t in tissues]),
        len([t for t in SMOOTH if t in tissues])))
    return G, H


def write_cov(df, sym2, name):
    d = df.copy()
    d.insert(0, "GENE", [sym2[g] for g in d.index])
    f = WORK / (name + ".covar.txt")
    d.dropna().to_csv(f, sep="\t", index=False, float_format="%.6g")
    return f


def x4(sym2):
    say("")
    say("-" * 72)
    say("[SECTION] X4 organ-identity hypothesis (exploratory, CHANGELOG v0.2.13)")
    G, H = build_tables(sym2)
    rows = []
    # X4a
    fa = write_cov(G[["MULL", "UTERUS", "SMOOTH", "Average"]], sym2, "X4a")
    say("[X4a] MULL and UTERUS | Average + SMOOTH (one-sided greater)")
    for trait in TRAITS:
        r = magma_gp(trait, fa, "Average,SMOOTH", "X4a")
        r0 = magma_gp(trait, fa, "Average", "X4a0")
        if not r or "MULL" not in r:
            continue
        m, u = r["MULL"], r.get("UTERUS", {"beta": np.nan, "p": np.nan})
        m0 = (r0 or {}).get("MULL", {"beta": np.nan, "p": np.nan})
        say("    [X4a] %-15s MULL|Avg p=%.2e -> MULL|Avg+SMOOTH beta=%.4f p=%.2e | UTERUS|Avg+SMOOTH beta=%.4f p=%.2e | genes=%d" % (
            trait[:15], m0["p"], m["beta"], m["p"], u["beta"], u["p"], m["ngenes"]))
        rows.append({"test": "X4a", "trait": trait, "beta": m["beta"], "beta_std": m["beta_std"], "se": m["se"], "p": m["p"],
                     "p_unadj_smooth": m0["p"], "uterus_beta": u["beta"], "uterus_p": u["p"]})
    # X4b
    shared = G.index.intersection(H.index)
    tau = H.loc[shared, "TAU_H"]
    cuts = tau.quantile([1 / 3, 2 / 3]).values
    ter = pd.Series(np.where(tau <= cuts[0], "T1", np.where(tau <= cuts[1], "T2", "T3")), index=shared)
    say("[X4b] HECA tau tertile cuts %.3f / %.3f (T1 broad ... T3 cell-type specific); genes per tertile %s" % (
        cuts[0], cuts[1], ", ".join("%s=%d" % (k, v) for k, v in ter.value_counts().sort_index().items())))
    for trait in TRAITS:
        res = {}
        for k in ("T1", "T2", "T3"):
            idx = ter.index[ter == k]
            f = write_cov(G.loc[idx, ["MULL", "SMOOTH", "Average"]], sym2, "X4b_" + k)
            r = magma_gp(trait, f, "Average,SMOOTH", "X4b_" + k)
            if r and "MULL" in r:
                res[k] = r["MULL"]
                rows.append({"test": "X4b_" + k, "trait": trait, "beta": r["MULL"]["beta"], "beta_std": r["MULL"]["beta_std"],
                             "se": r["MULL"]["se"], "p": r["MULL"]["p"]})
        if "T1" in res and "T3" in res:
            a, b = res["T1"], res["T3"]
            # compare on the semi-standardized scale (covariate SD may differ between tertiles)
            sa = a["se"] * (a["beta_std"] / a["beta"]) if a["beta"] else np.nan
            sb = b["se"] * (b["beta_std"] / b["beta"]) if b["beta"] else np.nan
            z = (a["beta_std"] - b["beta_std"]) / math.sqrt(sa ** 2 + sb ** 2) if np.isfinite(sa) and np.isfinite(sb) else np.nan
            pz = 0.5 * math.erfc(z / math.sqrt(2)) if np.isfinite(z) else np.nan
            say("    [X4b] %-15s BETA_STD T1=%.4f (p=%.1e) T2=%.4f (p=%.1e) T3=%.4f (p=%.1e) | T1>T3 z=%.2f p=%.3f" % (
                trait[:15], a["beta_std"], a["p"], res.get("T2", {}).get("beta_std", np.nan), res.get("T2", {}).get("p", np.nan),
                b["beta_std"], b["p"], z, pz))
            rows.append({"test": "X4b_T1vsT3", "trait": trait, "z": z, "p": pz})
    # X4c
    say("[X4c] HECA dStromal specificity | AVG_H, then | AVG_H + MULL (GTEx)")
    cc = H.loc[shared].copy()
    cc["MULL"] = G.loc[shared, "MULL"]
    cols = [c for c in cc.columns if c.startswith("S_")] + ["AVG_H", "MULL"]
    fc = write_cov(cc[cols], sym2, "X4c")
    for trait in TRAITS:
        r1 = magma_gp(trait, fc, "AVG_H", "X4c1")
        r2 = magma_gp(trait, fc, "AVG_H,MULL", "X4c2")
        if not r1 or not r2:
            continue
        parts = []
        for c in [c for c in cols if c.startswith("S_")]:
            if c in r1 and c in r2:
                parts.append("%s p %.3f -> %.3f (beta %.3f -> %.3f)" % (c[2:], r1[c]["p"], r2[c]["p"], r1[c]["beta"], r2[c]["beta"]))
                rows.append({"test": "X4c_" + c, "trait": trait, "p_before": r1[c]["p"], "p_after": r2[c]["p"],
                             "beta_before": r1[c]["beta"], "beta_after": r2[c]["beta"]})
        say("    [X4c] %-15s %s" % (trait[:15], " | ".join(parts)))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "X4_results.csv", index=False)
    # Fixed reading (v0.2.13)
    def get(test, trait, col="p"):
        s = R[(R["test"] == test) & (R["trait"] == trait)]
        return float(s[col].iloc[0]) if len(s) and col in s else float("nan")
    mull = get("X4a", "GCST90205183") < 0.05
    spec = mull and get("X4a", "HYPERTENSION") >= 0.05 and get("X4a", "T2D") >= 0.05
    organ = get("X4b_T1", "GCST90205183") < 0.05 and get("X4b_T1vsT3", "GCST90205183") < 0.05
    say("[X4_VERDICT] Muellerian beyond smooth muscle: %s | endometriosis-specific (HTN, T2D p>=0.05): %s | organ-identity pattern (T1 > T3): %s" % (
        "YES" if mull else "NO", "YES" if spec else "NO", "YES" if organ else "NO"))
    SUM["results"]["X4"] = {"mullerian": mull, "specific": spec, "organ_identity": organ}


def main():
    global _fh
    run_previous()
    for d in (OUT, WORK, LOGS):
        d.mkdir(parents=True, exist_ok=True)
    _fh = open(LOG, "w", encoding="utf-8")
    say("=" * 72)
    say("PAPER 9 | v0.0.22 | X4 organ-identity hypothesis (exploratory) | " + TS)
    say("=" * 72)
    t0 = time.time()
    try:
        x4(entrez_unique())
    except Exception as e:
        import traceback
        say("[STEP_ERROR] x4 %s: %s" % (type(e).__name__, str(e)[:300]))
        say(traceback.format_exc()[-1200:])
    with open(OUT / "p9_v022_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] X4 rules fixed in CHANGELOG v0.2.13 before this run (exploratory, post hoc)")
    say("[DONE] PAPER9_V022_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
