#!/usr/bin/env python3
# Paper 9 | v0.0.75 | Registered chromatin test: S-LDSC on GSE183771 lineage-specific accessible chromatin
#   Rules: CHANGELOG v0.2.64 (fixed after the GWAS-free v074 audit and before this run).
#   Annotations and LD scores from v074 (data/processed/ldscores/gse183771_v074); baseline v1.2 + ALL control.
#   CH-1  KOLLER_EUR: stroma_S10 coefficient > 0 given baseline + ALL (one-sided p < 0.05)
#   CH-2  KOLLER_EUR: stroma_S10 coefficient > epithelium_S10 coefficient in one model; tested as the coefficient of
#         stroma_S10 in a model with (stroma_S10 + epithelium_S10) and ALL (re-parameterisation, block jackknife SE)
#   CH-3  KOLLER_EUR: stroma_S10 has the largest coefficient z among the evaluable lineages
#   Replication (non-independent): GCST90205183, ENDO_ALL. Controls: HYPERTENSION, T2D (specificity), LEIOMYOMA
#   (descriptive). Sensitivity: top 5% / 20%, mesenchymal_S10. Stopping rule: if CH-1 fails, the chromatin line ends.
# Runs in ~/Desktop/paper9/.venv, low heat (nice 10, 2 BLAS threads). ~30-90 min for 6 traits; reruns use the cache.

import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
try:
    os.nice(10)
except Exception:
    pass

import datetime as dt
import json
import re
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
RAW = REPO / "data" / "raw"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
LDSC_PY = TOOLS / "ldsc_venv" / "bin" / "python"
LDSC = TOOLS / "ldsc_py3" / "ldsc.py"
REF = RAW / "ldsc_ref" / "unpacked"
MUNGED = PROC / "munged"
LD_DIR = PROC / "ldscores" / "gse183771_v074"
SETS = LD_DIR / "sets"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "chromatin" / ("p9_v075_" + TS)
CACHE = REPO / "results" / "chromatin" / "cache_v075"
LOG = REPO / "logs" / ("p9_v075_chromatin_sldsc_" + TS + ".log")
NICE = ["nice", "-n", "10"]
BLAS_THREADS = "2"
L2_ENV = dict(os.environ, OPENBLAS_NUM_THREADS=BLAS_THREADS, OMP_NUM_THREADS=BLAS_THREADS,
              VECLIB_MAXIMUM_THREADS=BLAS_THREADS, MKL_NUM_THREADS=BLAS_THREADS, NUMEXPR_NUM_THREADS=BLAS_THREADS)

# ----------------------------------------------------------------------------
# Fixed rules (CHANGELOG v0.2.64)
# ----------------------------------------------------------------------------
PRIMARY = "KOLLER_EUR"
REPLICATION = ["GCST90205183", "ENDO_ALL"]
CONTROLS = ["HYPERTENSION", "T2D"]
DESCRIPTIVE = ["LEIOMYOMA"]
TRAITS = [PRIMARY] + REPLICATION + CONTROLS + DESCRIPTIVE
LINEAGES = ["stroma", "epithelium", "endothelium", "immune"]
MAIN = ["%s_S10" % x for x in LINEAGES]
SENS = ["stroma_S5", "stroma_S20", "mesenchymal_S10", "epithelium_S5", "epithelium_S20"]
SUMSET = "stroma_plus_epithelium_S10"
CONTRAST = "stroma_minus_epithelium_S10"
ALPHA = 0.05

OUT.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V075_CHROMATIN_SLDSC", "timestamp_utc": TS, "verdicts": {}, "results": {}}


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def sh(cmd, timeout=None, env=None):
    try:
        p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout, env=env)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"
    except FileNotFoundError as e:
        return 127, "", str(e)


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "-"
    return "%.3f" % p if p >= 0.001 else "%.0e" % p


def prefix_of(dirname):
    d = REF / dirname
    for p in sorted(d.rglob("*.l2.ldscore.gz")) if d.exists() else []:
        m = re.match(r"^(.*\.)1\.l2\.ldscore\.gz$", p.name)
        if m:
            return str(p.parent / m.group(1))
    return None


def set_prefix(k):
    return str(SETS / k / (k + "."))


def have_set(k):
    return all((SETS / k / ("%s.%d.l2.ldscore.gz" % (k, c))).exists() for c in range(1, 23))


# ----------------------------------------------------------------------------
# Sum annotation for the contrast (stroma_S10 + epithelium_S10)
# ----------------------------------------------------------------------------
def make_sum_set():
    if have_set(SUMSET):
        return True
    d = SETS / SUMSET
    d.mkdir(parents=True, exist_ok=True)
    for c in range(1, 23):
        a = pd.read_csv(SETS / "stroma_S10" / ("stroma_S10.%d.l2.ldscore.gz" % c), sep="\t")
        b = pd.read_csv(SETS / "epithelium_S10" / ("epithelium_S10.%d.l2.ldscore.gz" % c), sep="\t")
        if len(a) != len(b) or not (a["SNP"].values == b["SNP"].values).all():
            say("[SUM] chr%d: SNP rows differ between stroma_S10 and epithelium_S10 -> stop" % c)
            return False
        out = a[["CHR", "SNP", "BP"]].copy()
        out[SUMSET + "L2"] = a.iloc[:, 3].values + b.iloc[:, 3].values
        out.to_csv(d / ("%s.%d.l2.ldscore.gz" % (SUMSET, c)), sep="\t", index=False, compression="gzip")
        for suf in ("M", "M_5_50"):
            ma = float(open(SETS / "stroma_S10" / ("stroma_S10.%d.l2.%s" % (c, suf))).read().split()[0])
            mb = float(open(SETS / "epithelium_S10" / ("epithelium_S10.%d.l2.%s" % (c, suf))).read().split()[0])
            open(d / ("%s.%d.l2.%s" % (SUMSET, c, suf)), "w").write("%d\n" % int(round(ma + mb)))
    return True


# ----------------------------------------------------------------------------
# h2-cts runs
# ----------------------------------------------------------------------------
def write_ldcts():
    allp = set_prefix("ALL")
    f1 = CACHE / "marginal.ldcts"
    with open(f1, "w") as f:
        for k in MAIN + SENS:
            if have_set(k):
                f.write("%s\t%s,%s\n" % (k, set_prefix(k), allp))
    f2 = CACHE / "contrast.ldcts"
    with open(f2, "w") as f:
        f.write("%s\t%s,%s,%s\n" % (CONTRAST, set_prefix("stroma_S10"), set_prefix(SUMSET), allp))
        f.write("stroma_S10__GIVEN__epithelium_S10\t%s,%s,%s\n" % (set_prefix("stroma_S10"), allp, set_prefix("epithelium_S10")))
        f.write("epithelium_S10__GIVEN__stroma_S10\t%s,%s,%s\n" % (set_prefix("epithelium_S10"), allp, set_prefix("stroma_S10")))
    return {"marginal": f1, "contrast": f2}


def run_cts(trait, tag, ldcts, baseline, weights):
    out = CACHE / ("%s__%s" % (trait, tag))
    res = Path(str(out) + ".cell_type_results.txt")
    if not res.exists():
        cmd = NICE + [str(LDSC_PY), str(LDSC), "--h2-cts", str(MUNGED / (trait + ".sumstats.gz")), "--ref-ld-chr", baseline,
                      "--ref-ld-chr-cts", str(ldcts), "--w-ld-chr", weights, "--out", str(out)]
        t0 = time.time()
        rc, o, e = sh(cmd, timeout=8 * 3600, env=L2_ENV)
        if rc != 0 or not res.exists():
            say("    [CTS] %-14s %-9s rc=%d FAILED" % (trait, tag, rc))
            for line in (o + e).strip().splitlines()[-6:]:
                say("       " + line[:200])
            return None
        say("    [CTS] %-14s %-9s done (%.0f min)" % (trait, tag, (time.time() - t0) / 60))
    df = pd.read_csv(res, sep="\t")
    df["trait"] = trait
    df["z"] = df["Coefficient"] / df["Coefficient_std_error"]
    return df


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
def figure(M):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    sub = M[M["Name"].isin(MAIN)]
    Z = sub.pivot(index="Name", columns="trait", values="z").reindex(index=MAIN, columns=[t for t in TRAITS if t in set(sub["trait"])])
    P = sub.pivot(index="Name", columns="trait", values="Coefficient_P_value").reindex(index=Z.index, columns=Z.columns)
    fig, ax = plt.subplots(figsize=(1.2 * Z.shape[1] + 2.4, 0.55 * Z.shape[0] + 1.6))
    v = float(np.nanmax(np.abs(Z.values))) if np.isfinite(Z.values).any() else 1.0
    im = ax.imshow(Z.values, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            p = P.values[i, j]
            ax.text(j, i, "%.2f%s" % (Z.values[i, j], "*" if p < ALPHA else ""), ha="center", va="center", fontsize=8,
                    color="white" if abs(Z.values[i, j]) > 0.6 * v else "black")
    ax.set_xticks(range(Z.shape[1]))
    ax.set_xticklabels(Z.columns, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(Z.shape[0]))
    ax.set_yticklabels([x.replace("_S10", "") for x in Z.index], fontsize=9)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.06, pad=0.28)
    cb.set_label("S-LDSC coefficient z (baseline v1.2 + all accessible peaks); * one-sided p < 0.05", fontsize=7)
    ax.set_title("GSE183771 lineage-specific accessible chromatin (top 10%)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "S7_chromatin_sldsc.png", dpi=300)
    fig.savefig(OUT / "S7_chromatin_sldsc.pdf")
    plt.close(fig)


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.75 | registered chromatin S-LDSC (CH-1..CH-3) on GSE183771 annotations | " + TS)
    say("=" * 72)
    t0 = time.time()
    baseline, weights = prefix_of("baseline_v1.2"), prefix_of("1000G_Phase3_weights_hm3_no_MHC")
    missing = [k for k in ["ALL"] + MAIN if not have_set(k)]
    say("[SETUP] baseline %s | weights %s | LDSC %s" % ("OK" if baseline else "MISSING", "OK" if weights else "MISSING",
                                                     "OK" if LDSC_PY.exists() and LDSC.exists() else "MISSING"))
    say("[SETUP] v074 sets: %s" % ("all present" if not missing else "missing " + ", ".join(missing)))
    tr = [t for t in TRAITS if (MUNGED / (t + ".sumstats.gz")).exists()]
    say("[SETUP] traits: %s%s" % (", ".join(tr), "" if len(tr) == len(TRAITS) else " (missing: %s)" % ", ".join(sorted(set(TRAITS) - set(tr)))))
    if not baseline or not weights or missing or PRIMARY not in tr or not LDSC_PY.exists():
        say("[STOP] inputs incomplete")
        return finish(t0)
    if not make_sum_set():
        return finish(t0)
    ld = write_ldcts()
    say("")
    say("-" * 72)
    say("[SECTION] S-LDSC runs (h2-cts; each annotation with ALL and baseline v1.2)")
    rows = []
    for t in tr:
        for tag in ("marginal", "contrast"):
            df = run_cts(t, tag, ld[tag], baseline, weights)
            if df is not None:
                df["model"] = tag
                rows.append(df)
    if not rows:
        say("[STOP] no S-LDSC results")
        return finish(t0)
    M = pd.concat(rows, ignore_index=True)
    M.to_csv(OUT / "chromatin_sldsc_all.csv", index=False)

    say("")
    say("-" * 72)
    say("[SECTION] Results (coefficient z; one-sided p)")
    mg = M[M["model"] == "marginal"]
    say("  %-22s %s" % ("annotation", " ".join("%16s" % t[:16] for t in tr)))
    for k in MAIN + SENS:
        cells = []
        for t in tr:
            r = mg[(mg["trait"] == t) & (mg["Name"] == k)]
            cells.append("%16s" % ("%.2f (%s)" % (r["z"].iloc[0], fmt_p(r["Coefficient_P_value"].iloc[0])) if len(r) else "-"))
        say("  %-22s %s" % (k, " ".join(cells)))
    ct = M[M["model"] == "contrast"]
    for k in [CONTRAST, "stroma_S10__GIVEN__epithelium_S10", "epithelium_S10__GIVEN__stroma_S10"]:
        cells = []
        for t in tr:
            r = ct[(ct["trait"] == t) & (ct["Name"] == k)]
            cells.append("%16s" % ("%.2f (%s)" % (r["z"].iloc[0], fmt_p(r["Coefficient_P_value"].iloc[0])) if len(r) else "-"))
        say("  %-22s %s" % (k.replace("__GIVEN__", " | ")[:22], " ".join(cells)))

    say("")
    say("-" * 72)
    say("[SECTION] Verdicts (rules v0.2.64)")

    def get(t, k, model="marginal"):
        r = M[(M["trait"] == t) & (M["Name"] == k) & (M["model"] == model)]
        return (float(r["z"].iloc[0]), float(r["Coefficient_P_value"].iloc[0])) if len(r) else (np.nan, np.nan)

    z1, p1 = get(PRIMARY, "stroma_S10")
    ch1 = bool(p1 < ALPHA)
    say("[CH-1] %s stroma_S10: z %.2f, p %s -> %s" % (PRIMARY, z1, fmt_p(p1), "SUPPORTED" if ch1 else "NOT SUPPORTED"))
    z2, p2 = get(PRIMARY, CONTRAST, "contrast")
    ch2 = bool(p2 < ALPHA)
    say("[CH-2] %s stroma_S10 - epithelium_S10 (joint model): z %.2f, p %s -> %s" % (PRIMARY, z2, fmt_p(p2), "SUPPORTED" if ch2 else "NOT SUPPORTED"))
    zs = {k: get(PRIMARY, k)[0] for k in MAIN}
    top = max(zs, key=lambda k: zs[k] if np.isfinite(zs[k]) else -1e9)
    ch3 = top == "stroma_S10"
    say("[CH-3] %s rank by z: %s -> %s" % (PRIMARY, ", ".join("%s %.2f" % (k.replace("_S10", ""), zs[k])
                                                            for k in sorted(zs, key=lambda k: -zs[k])), "SUPPORTED" if ch3 else "NOT SUPPORTED"))
    rep = {t: get(t, "stroma_S10") for t in REPLICATION if t in tr}
    say("[REPLICATION] stroma_S10: %s" % "; ".join("%s z %.2f p %s" % (t, z, fmt_p(p)) for t, (z, p) in rep.items()))
    ctl = {t: get(t, "stroma_S10") for t in CONTROLS if t in tr}
    spec = all(p >= ALPHA for z, p in ctl.values())
    say("[CONTROLS] stroma_S10: %s -> %s" % ("; ".join("%s z %.2f p %s" % (t, z, fmt_p(p)) for t, (z, p) in ctl.items()),
                                            "trait-specific" if spec else "NOT trait-specific (stroma also enriched for a control trait)"))
    for t in DESCRIPTIVE:
        if t in tr:
            zz = {k: get(t, k)[0] for k in MAIN}
            say("[DESCRIPTIVE] %s rank by z: %s" % (t, ", ".join("%s %.2f" % (k.replace("_S10", ""), zz[k]) for k in sorted(zz, key=lambda k: -zz[k]))))
    sens = {k: get(PRIMARY, k) for k in ["stroma_S5", "stroma_S20", "mesenchymal_S10"]}
    say("[SENSITIVITY] %s: %s" % (PRIMARY, "; ".join("%s z %.2f p %s" % (k, z, fmt_p(p)) for k, (z, p) in sens.items())))
    if not ch1:
        say("[STOP RULE] CH-1 not supported -> the chromatin line ends here; reported as a negative result")
    SUM["verdicts"] = {"CH-1": ch1, "CH-2": ch2, "CH-3": ch3, "controls_trait_specific": spec}
    SUM["results"] = {"CH-1": [z1, p1], "CH-2": [z2, p2], "CH-3_z": zs, "replication": rep, "controls": ctl, "sensitivity": sens}
    figure(mg)
    finish(t0)


def finish(t0):
    with open(OUT / "p9_v075_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules CHANGELOG v0.2.64 fixed before this run")
    say("[DONE] PAPER9_V075_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
