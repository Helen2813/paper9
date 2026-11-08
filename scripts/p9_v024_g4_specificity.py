#!/usr/bin/env python3
# Paper 9 | v0.0.24 | G4-S specificity of the decidual signal (HECA M1) + G1 epithelial permutation.
#   Rules: CHANGELOG v0.2.15 (fixed before this run; post hoc, descriptive).
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). ~15-25 min.

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
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
HECA_GS = PROC / "genesets" / "heca_v014"
COV_X4A = INTERIM / "magma_v022" / "X4a.covar.txt"
PB213 = PROC / "gse213216_pseudobulk_clusters"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "g4s" / ("p9_v024_g4s_" + TS)
WORK = INTERIM / "magma_v024"
LOG = REPO / "logs" / ("p9_v024_g4s_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.15)
FOCUS = ["KOLLER_EUR", "GCST90205183"]
PANEL = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL", "ENDO_OVARY", "ENDO_PERITONEUM", "ENDO_DEEP",
         "LEIOMYOMA", "HYPERTENSION", "T2D"]
CONTROLS_Z = ["T2D", "HYPERTENSION", "LEIOMYOMA"]
DS = {"dStromal_mid": "DSMID", "dStromal_early": "DSEARLY", "dStromal_late": "DSLATE"}
VASC = {"ePV_2": "EPV2", "Arterial": "ART", "uSMCs": "USMC"}
N_PERM = 100
PERM_SEED = 20261106

for d in (OUT, WORK, LOG.parent):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V024_G4S", "timestamp_utc": TS, "results": {}}


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
                    res[name] = {"beta": float(d["BETA"]), "se": float(d["SE"]), "p": float(d["P"])}
                except (KeyError, ValueError):
                    continue
            return res
    say("    [MAGMA_FAIL] %s %s %s" % (trait, tag, (o + e).strip()[-200:]))
    return None


def gene_z(trait):
    f = MAGMA_OUT / (trait + "_35_10.genes.out")
    if not f.exists():
        return None
    g = pd.read_csv(f, sep=r"\s+")
    return pd.Series(g["ZSTAT"].values, index=g["GENE"].astype(str))


# ----------------------------------------------------------------------------
# G4-S
# ----------------------------------------------------------------------------
def build_cov(sym2):
    S = pd.read_csv(HECA_GS / "specificity_by_type.csv.gz", index_col=0)
    M = pd.read_csv(HECA_GS / "mean_cpm_by_type.csv.gz", index_col=0).loc[S.index, S.columns]
    cov = pd.DataFrame(index=S.index)
    for lab, k in list(DS.items()) + list(VASC.items()):
        if lab in S.columns:
            cov[k] = S[lab]
        else:
            say("[G4S_WARN] type missing in HECA table: " + lab)
    cov["Average"] = np.log2(M.mean(axis=1) + 1.0)
    cov["GENE"] = [sym2.get(g) for g in cov.index]
    cov = cov.dropna(subset=["GENE"])
    cov["GENE"] = cov["GENE"].astype(str)
    x4 = pd.read_csv(COV_X4A, sep="\t")
    x4["GENE"] = x4["GENE"].astype(str)
    cov = cov.merge(x4[["GENE", "MULL"]], on="GENE", how="left")
    for t in CONTROLS_Z:
        z = gene_z(t)
        cov["Z_" + t] = cov["GENE"].map(z) if z is not None else np.nan
    full = cov.dropna()
    say("[G4S_COV] HECA genes=%d; with MULL and all control gene z=%d" % (len(cov), len(full)))
    cols = ["GENE"] + [c for c in full.columns if c != "GENE"]
    f = WORK / "G4S.covar.txt"
    full[cols].to_csv(f, sep="\t", index=False, float_format="%.6g")
    return f


def g4s(sym2):
    say("")
    say("-" * 72)
    say("[SECTION] G4-S specificity of the decidual signal (HECA control donors, M1; post hoc, v0.2.15)")
    covf = build_cov(sym2)
    ds_keys = [k for k in DS.values()]
    models = [("base", "Average"),
              ("+zT2D", "Average,Z_T2D"), ("+zHTN", "Average,Z_HYPERTENSION"), ("+zLEIO", "Average,Z_LEIOMYOMA"),
              ("+vasc", "Average,EPV2,ART,USMC"), ("+MULL", "Average,MULL")]
    have = set(open(covf).readline().rstrip("\n").split("\t"))
    models = [(m, ",".join(v for v in c.split(",") if v in have)) for m, c in models]
    rows = []
    for trait in FOCUS:
        for mname, cond in models:
            r = magma_gp(trait, covf, cond, "G4S_" + mname.replace("+", "p"))
            if not r:
                continue
            for lab, k in DS.items():
                if k in r:
                    rows.append({"trait": trait, "model": mname, "type": lab, **r[k]})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "G4S_conditional.csv", index=False)
    for trait in FOCUS:
        sub = df[df["trait"] == trait]
        if not len(sub):
            continue
        say("  [S5] %s (dStromal_mid beta / p; same genes in every model)" % trait)
        for mname, _ in models:
            s = sub[sub["model"] == mname].set_index("type")
            if "dStromal_mid" in s.index:
                say("    [S5] %-7s dS_mid beta=%.3f p=%.4f | dS_early p=%.3f | dS_late p=%.3f" % (
                    mname, s.loc["dStromal_mid", "beta"], s.loc["dStromal_mid", "p"],
                    s.loc["dStromal_early", "p"] if "dStromal_early" in s.index else np.nan,
                    s.loc["dStromal_late", "p"] if "dStromal_late" in s.index else np.nan))
        ok_z = all(((sub["model"] == m) & (sub["type"] == "dStromal_mid") & (sub["p"] < 0.05)).any() for m in ("+zT2D", "+zHTN", "+zLEIO"))
        ok_v = ((sub["model"] == "+vasc") & (sub["type"] == "dStromal_mid") & (sub["p"] < 0.05)).any()
        ok_m = ((sub["model"] == "+MULL") & (sub["type"] == "dStromal_mid") & (sub["p"] < 0.05)).any()
        say("  [S5_READING] %s: not explained by generic gene-level signal: %s | separable from perivascular/SMC: %s | beyond MULL: %s" % (
            trait, "YES" if ok_z else "NO", "YES" if ok_v else "NO", "YES" if ok_m else "NO"))
        SUM["results"]["S5_" + trait] = {"generic": ok_z, "vascular": ok_v, "mull": ok_m}

    # S5d panel with the original M1 file (all 32 types)
    say("  [S5d] M1 panel (all HECA types; dStromal Holm over 3; rank among types)")
    types = list(pd.read_csv(HECA_GS / "specificity_by_type.csv.gz", index_col=0, nrows=1).columns)
    lab = {"T%02d" % (i + 1): t for i, t in enumerate(types)}
    m1 = INTERIM / "magma_v020" / "C5_M1_specificity.covar.txt"
    for trait in PANEL:
        r = magma_gp(trait, m1, "Average", "S5d")
        if not r:
            continue
        t = pd.DataFrame([{"type": lab[k], **r[k]} for k in lab if k in r]).sort_values("p").reset_index(drop=True)
        ds = t[t["type"].isin(DS.keys())].copy()
        ds["holm"] = holm(ds["p"].values)
        say("    [S5d] %-15s top=%s(p=%.1e) | %s" % (trait[:15], t["type"].iloc[0], t["p"].iloc[0], "; ".join(
            "%s r%d p=%.3f H=%.3f" % (x.replace("dStromal_", "dS_"), i + 1, p, h) for i, x, p, h in zip(ds.index, ds["type"], ds["p"], ds["holm"]))))


# ----------------------------------------------------------------------------
# G1-S epithelial permutation
# ----------------------------------------------------------------------------
def unit_profiles(M, G, lineage, context):
    sub = G[(G["lineage_ref"] == lineage) & (G["context"] == context) & (G["agree"])]
    prof = {}
    for pid, g in sub.groupby("patient"):
        if int(g["n_cells"].sum()) < 10:
            continue
        v = M[g.index.to_numpy()].sum(axis=0)
        prof[(pid, context)] = np.log2(v / max(v.sum(), 1) * 1e6 + 1)
    return prof


def write_contrast(A_list, B_list, genes, mapped, sym2, tag):
    A = np.vstack(A_list).mean(axis=0)
    B = np.vstack(B_list).mean(axis=0)
    ex = (np.maximum(2 ** A - 1, 2 ** B - 1) >= 1.0) & mapped
    g = np.asarray(genes)[ex]
    covf = WORK / ("G1_%s.covar.txt" % tag)
    pd.DataFrame({"GENE": [sym2[x] for x in g], "D": (A - B)[ex], "DN": (B - A)[ex], "AVG": ((A + B) / 2)[ex]}).to_csv(
        covf, sep="\t", index=False, float_format="%.6g")
    return covf


def g1_epithelial(sym2):
    say("")
    say("-" * 72)
    say("[SECTION] G1-S epithelial contrast: patient-label permutation (B=%d)" % N_PERM)
    M = np.load(PB213 / "pseudobulk_counts.npy", mmap_mode="r")
    G = pd.read_csv(PB213 / "pseudobulk_groups.csv").reset_index(drop=True)
    genes = list(pd.read_csv(PB213 / "genes.csv")["symbol"].astype(str))
    mapped = np.array([x in sym2 for x in genes])
    pa = unit_profiles(M, G, "epithelial", "peritoneal_lesion")
    pb = unit_profiles(M, G, "epithelial", "eutopic")
    say("    [G1E] units lesion=%d eutopic=%d" % (len(pa), len(pb)))
    if len(pa) < 3 or len(pb) < 3:
        return
    units = list(pa.values()) + list(pb.values())
    lab = np.array([1] * len(pa) + [0] * len(pb))
    obs = write_contrast(list(pa.values()), list(pb.values()), genes, mapped, sym2, "obs")
    traits = [t for t in ["KOLLER_EUR", "GCST90205183", "T2D"] if (MAGMA_OUT / (t + "_35_10.genes.raw")).exists()]
    tobs = {}
    for t in traits:
        r = magma_gp(t, obs, "AVG", "G1E_obs")
        if r and "D" in r:
            tobs[t] = (r["D"]["beta"] / r["D"]["se"], r["DN"]["p"])
    rng = np.random.default_rng(PERM_SEED)
    perm = {t: [] for t in tobs}
    t0 = time.time()
    for b in range(N_PERM):
        pl = rng.permutation(lab)
        covf = write_contrast([units[j] for j in range(len(units)) if pl[j] == 1],
                              [units[j] for j in range(len(units)) if pl[j] == 0], genes, mapped, sym2, "perm")
        for t in tobs:
            r = magma_gp(t, covf, "AVG", "G1E_perm")
            if r and "D" in r:
                perm[t].append(r["D"]["beta"] / r["D"]["se"])
    for t, (to, pp) in tobs.items():
        ts = np.array(perm[t])
        pe = (1 + np.sum(ts <= to)) / (1 + len(ts)) if len(ts) else float("nan")
        say("    [G1E_PERM] %-13s t_obs=%.2f permutation p_less=%.3f (B=%d) vs parametric %.2e | perm t sd=%.2f" % (
            t[:13], to, pe, len(ts), pp, ts.std() if len(ts) else float("nan")))
        SUM["results"]["G1E_perm_" + t] = {"t_obs": to, "p_perm": pe}
    say("    [G1E_PERM] %.0f min" % ((time.time() - t0) / 60))


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.24 | G4-S decidual specificity + G1 epithelial permutation | " + TS)
    say("=" * 72)
    t0 = time.time()
    sym2 = entrez_unique()
    for name, fn in (("g4s", lambda: g4s(sym2)), ("g1e", lambda: g1_epithelial(sym2))):
        try:
            fn()
        except Exception as e:
            import traceback
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:300]))
            say(traceback.format_exc()[-1000:])
    with open(OUT / "p9_v024_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.15 before this run (post hoc, descriptive)")
    say("[DONE] PAPER9_V024_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
