#!/usr/bin/env python3
# Paper 9 | v0.0.50 | independent endometrial atlas and candidate cells of origin (rules CHANGELOG v0.2.40)
#   R2 GSE266265 (60 donors, not part of HECA): donor-paired pseudobulk DE per refined type vs the rest of the donor's
#      cells; seven marker families; six GWAS; rule flip and hypertension control
#   O1 candidate cells of origin inside the same biopsies: endometrial stroma, perivascular (MSC-like), mesothelium,
#      epithelium, endothelium; rank of each group under up-type rules; trait contrast vs T2D / hypertension
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). ~20-40 min; no downloads.

import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
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
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
PB266 = PROC / "gse266265_pseudobulk"
DE_DIR = PROC / "genesets" / "gse266265_v050_de"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "review" / ("p9_v050_" + TS)
WORK = INTERIM / "magma_v050"
LOG = REPO / "logs" / ("p9_v050_atlas2_origin_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.40)
MIN_CELLS_DONOR = 10
MIN_DONORS = 10
MIN_CPM = 1.0
LFC, PADJ, MIN_SET, N_MATCH = 1.2, 0.01, 20, 400
FAMILIES = ["ABS", "UP", "DOWN", "ABSN", "UPN", "DOWNN", "SPEC"]
TRAITS = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL", "LEIOMYOMA", "HYPERTENSION", "T2D"]
GROUPS = [("endometrial stroma", r"^(eF\d|dS\d|fib C7)"), ("perivascular / MSC-like", r"^(Prv|VSMC)"), ("mesothelium", r"^mesothelial"),
          ("epithelium", r"^(glandular|lumenal|ciliated|mid-secretory|MUC5B|TP63)"), ("endothelium", r"^(EC-|LEC)")]
VASC_RE = r"^(EC-artery|EC-|VSMC|Prv)"
MESO_MARKERS = ["MSLN", "UPK3B", "CALB2", "PRG4", "ITLN1"]
CONTRAST_TRAITS = [("KOLLER_EUR", "T2D"), ("KOLLER_EUR", "HYPERTENSION"), ("KOLLER_EUR", "LEIOMYOMA")]

for d in (OUT, WORK, LOG.parent, DE_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V050_ATLAS2_ORIGIN", "timestamp_utc": TS, "results": {}}
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


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def parse_full(gsa):
    rows = [l for l in open(gsa) if not l.startswith("#") and l.strip()]
    hdr = rows[0].split()
    out = {}
    for l in rows[1:]:
        d = dict(zip(hdr, l.split()))
        k = d.get("FULL_NAME", d.get("VARIABLE"))
        try:
            out[k] = {"p": float(d["P"]), "beta": float(d.get("BETA", "nan")), "beta_std": float(d.get("BETA_STD", "nan")), "se": float(d.get("SE", "nan"))}
        except ValueError:
            continue
    return out


def magma_covar(raw, cf, tag, cond):
    outp = WORK / tag
    gsa = Path(str(outp) + ".gsa.out")
    model = (["condition-hide=" + ",".join(cond)] if cond else []) + ["direction=greater"]
    if gsa.exists():
        gsa.unlink()
    rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--gene-covar", cf, "--model"] + model + ["--out", outp], timeout=7200)
    if gsa.exists():
        return parse_full(gsa)
    say("    [MAGMA_FAIL] %s %s" % (tag, (o + e).strip()[-300:]))
    return None


def tdist_sf(t, df):
    # two-sided p for Student t via the regularized incomplete beta function (scipy if available)
    try:
        from scipy.stats import t as T
        return 2.0 * T.sf(np.abs(t), df)
    except Exception:
        z = np.abs(t)
        return np.array([math.erfc(v / math.sqrt(2)) for v in z])


def clean(label):
    return re.sub(r"\$|\\|\{|\}|_\{|\^", "", str(label)).replace("Phi", "Phi")


def group_of(label):
    for g, rx in GROUPS:
        if re.search(rx, str(label)):
            return g
    return "immune / other"


# ----------------------------------------------------------------------------
# DE from donor pseudobulk
# ----------------------------------------------------------------------------
def de_tables(sym2ent):
    f_de, f_cpm = DE_DIR / "de_all.csv.gz", DE_DIR / "mean_cpm_by_type.csv.gz"
    if f_de.exists() and f_cpm.exists():
        say("[R2] DE tables cached (%s)" % DE_DIR.name)
        return pd.read_csv(f_de), pd.read_csv(f_cpm, index_col=0)
    G = pd.read_csv(PB266 / "pseudobulk_groups.csv")
    genes = np.array(pd.read_csv(PB266 / "genes.csv")["symbol"].astype(str))
    pb = np.load(PB266 / "pseudobulk_counts.npy", mmap_mode="r")
    say("[R2] GSE266265 pseudobulk: %d sample x type rows, %d genes, %d donors (%s)" % (
        len(G), len(genes), G["sample"].nunique(), ", ".join("%s %d" % (k, v) for k, v in G.groupby("status")["sample"].nunique().items())))
    t0 = time.time()
    tot = {}
    for s, idx in G.groupby("sample").groups.items():
        tot[s] = np.asarray(pb[np.sort(np.asarray(idx)), :], dtype=np.float64).sum(axis=0)
    ok = G[G["n_cells"] >= MIN_CELLS_DONOR]
    cnt = ok.groupby("refined")["sample"].nunique()
    types = sorted(cnt[cnt >= MIN_DONORS].index)
    say("[R2] eligible types (>= %d donors with >= %d cells): %d of %d" % (MIN_DONORS, MIN_CELLS_DONOR, len(types), G["refined"].nunique()))
    cpm_t, parts = {}, []
    for ty in types:
        rows = ok[ok["refined"] == ty]
        X = np.asarray(pb[rows.index.to_numpy(), :], dtype=np.float64)
        Rr = np.vstack([tot[s] for s in rows["sample"]]) - X
        Ct = X / np.maximum(X.sum(axis=1, keepdims=True), 1) * 1e6
        Cr = Rr / np.maximum(Rr.sum(axis=1, keepdims=True), 1) * 1e6
        cpm_t[ty] = Ct.mean(axis=0)
        expr = (Ct.mean(axis=0) >= MIN_CPM) | (Cr.mean(axis=0) >= MIN_CPM)
        D = np.log2(Ct[:, expr] + 1) - np.log2(Cr[:, expr] + 1)
        n = D.shape[0]
        m, sd = D.mean(axis=0), D.std(axis=0, ddof=1)
        tt = m / np.maximum(sd / math.sqrt(n), 1e-9)
        p = tdist_sf(tt, n - 1)
        parts.append(pd.DataFrame({"group": clean(ty), "gene": genes[expr], "lfc": m, "score": tt, "p": p, "padj": bh(p), "n_donors": n}))
    DE = pd.concat(parts, ignore_index=True)
    CPM = pd.DataFrame({clean(k): v for k, v in cpm_t.items()}, index=genes)
    CPM = CPM[(CPM >= MIN_CPM).any(axis=1)]
    DE.to_csv(f_de, index=False)
    CPM.to_csv(f_cpm)
    say("[R2] DE done: %d types, median %d donors per type (%.0f min)" % (len(types), int(DE.groupby("group")["n_donors"].first().median()), (time.time() - t0) / 60))
    return DE, CPM


def build(DE, CPM, sym2ent):
    meta, cols = {}, {}
    uni = sorted({g for g in DE["gene"].astype(str) if g in sym2ent})
    uset = set(uni)
    for g, d in DE.groupby("group", sort=False):
        d = d.assign(gene=d["gene"].astype(str))
        sig = d["padj"] < PADJ
        sel = {"ABS": d.loc[(d["lfc"].abs() > LFC) & sig, "gene"], "UP": d.loc[(d["lfc"] > LFC) & sig, "gene"],
               "DOWN": d.loc[(d["lfc"] < -LFC) & sig, "gene"],
               "ABSN": d.reindex(d["score"].abs().sort_values(ascending=False).index)["gene"].head(N_MATCH),
               "UPN": d.sort_values("score", ascending=False)["gene"].head(N_MATCH),
               "DOWNN": d.sort_values("score", ascending=True)["gene"].head(N_MATCH)}
        for fam, genes in sel.items():
            gl = [x for x in genes if x in uset]
            if len(gl) >= MIN_SET:
                k = "V%04d" % (len(meta) + 1)
                meta[k], cols[k] = (fam, g), gl
    T = pd.DataFrame({k: np.isin(uni, gl).astype(int) for k, gl in cols.items()}, index=[sym2ent[g] for g in uni])
    T = T[~T.index.duplicated()]
    C = CPM.copy()
    C.index = C.index.astype(str)
    C = C.groupby(level=0).mean()
    C = C.loc[[g for g in C.index if g in sym2ent]]
    S = C.div(C.sum(axis=1), axis=0)
    spec, sd = {}, {}
    for i, c in enumerate(S.columns):
        k = "S%03d" % (i + 1)
        spec[k], sd[k] = c, S[c].values
    sd["Average"] = np.log2(C.mean(axis=1).values + 1.0)
    Ts = pd.DataFrame(sd, index=[sym2ent[g] for g in C.index])
    Ts = Ts[~Ts.index.duplicated()]
    return T, Ts, meta, spec


def ranked(res, meta, spec, fam):
    if fam == "SPEC":
        it = [(spec[k], v) for k, v in res.items() if k in spec]
    else:
        it = [(meta[k][1], v) for k, v in res.items() if k in meta and meta[k][0] == fam]
    it.sort(key=lambda x: x[1]["p"])
    q = bh([v["p"] for _, v in it]) if it else []
    return [(l, v, qq) for (l, v), qq in zip(it, q)]


def main_run(sym2ent):
    DE, CPM = de_tables(sym2ent)
    say("")
    say("-" * 72)
    say("[SECTION] mesothelial label check (pre-specified: MSLN + UPK3B mean CPM >= 5x the median over other types)")
    meso = [c for c in CPM.columns if re.search(GROUPS[2][1], c)]
    meso_ok = False
    if meso:
        mk = [g for g in MESO_MARKERS if g in CPM.index]
        say("  [MESO] markers (mean CPM): " + "; ".join("%s %s" % (g, ", ".join("%s %.0f" % (c[:12], CPM.loc[g, c]) for c in meso)) for g in mk))
        core = [g for g in ("MSLN", "UPK3B") if g in CPM.index]
        if core:
            v = CPM.loc[core, meso[0]].sum()
            oth = CPM.loc[core, [c for c in CPM.columns if c not in meso]].sum(axis=0).median()
            meso_ok = bool(v >= 5 * max(oth, 1e-9))
            say("  [MESO] MSLN + UPK3B in %s: %.0f vs median of other types %.1f -> %s" % (meso[0], v, oth, "mesothelial identity CONFIRMED" if meso_ok else "NOT confirmed (group reported but not interpreted)"))
    else:
        say("  [MESO] no mesothelial type eligible")
    SUM["results"]["mesothelium_confirmed"] = meso_ok
    T, Ts, meta, spec = build(DE, CPM, sym2ent)
    cf, cfs = WORK / "covar_sets.txt", WORK / "covar_spec.txt"
    T.to_csv(cf, sep="\t", index_label="GENE")
    Ts.to_csv(cfs, sep="\t", float_format="%.6g", index_label="GENE")
    say("[R2] %d genes; %d set indicators; %d specificity columns" % (len(T), len(meta), len(spec)))
    res = {}
    t0 = time.time()
    for t in TRAITS:
        raw = MAGMA_OUT / (t + "_35_10.genes.raw")
        if not raw.exists():
            say("[R2] %s genes.raw missing" % t)
            continue
        r0 = magma_covar(raw, cf, "sets_%s" % t, [])
        r1 = magma_covar(raw, cfs, "spec_%s" % t, ["Average"])
        if r0 and r1:
            res[t] = {fam: ranked(r1 if fam == "SPEC" else r0, meta, spec, fam) for fam in FAMILIES}
    say("[R2] MAGMA runs done (%.0f min)" % ((time.time() - t0) / 60))
    rows = []
    for t, F in res.items():
        for fam, it in F.items():
            for i, (l, v, q) in enumerate(it):
                rows.append({"trait": t, "family": fam, "rank": i + 1, "label": l, "group": group_of(l), "p": v["p"], "q": q,
                             "beta_std": v["beta_std"], "beta": v["beta"], "se": v["se"]})
    A = pd.DataFrame(rows)
    A.to_csv(OUT / "r2_gse266265_family_rankings.csv", index=False)
    # R2 report
    say("")
    say("-" * 72)
    say("[SECTION] R2 seven families on an independent endometrial atlas (GSE266265; top 3 per family: p / BH q)")
    for t in TRAITS:
        if t not in res:
            continue
        say("  %s" % t)
        for fam in FAMILIES:
            it = res[t][fam]
            say("    %-6s %s" % (fam, "; ".join("%s %.0e/%.2f" % (l[:20], v["p"], q) for l, v, q in it[:3])))
    first = {(t, fam): res[t][fam][0][0] for t in res for fam in FAMILIES if res[t][fam]}
    k_up, k_abs = first.get(("KOLLER_EUR", "UP"), ""), first.get(("KOLLER_EUR", "ABS"), "")
    flip = group_of(k_up) == "endometrial stroma" and group_of(k_abs) != "endometrial stroma"
    say("[H_REPL_FLIP_VERDICT] KOLLER_EUR: up-only first %s (%s) | |log2FC| first %s (%s) -> %s" % (
        k_up, group_of(k_up), k_abs, group_of(k_abs), "SUPPORTED" if flip else "NOT SUPPORTED"))
    hv = {fam: bool(re.search(VASC_RE, first.get(("HYPERTENSION", fam), "^$"))) for fam in FAMILIES if ("HYPERTENSION", fam) in first}
    ctrl = all(hv.get(f, False) for f in ("UP", "UPN", "SPEC"))
    say("[H_REPL_CTRL_VERDICT] HYPERTENSION vascular / perivascular first: %s -> %s (up-type families UP, UPN, SPEC required)" % (
        ", ".join("%s %s" % (f, "YES" if v else "no") for f, v in hv.items()), "SUPPORTED" if ctrl else "NOT SUPPORTED"))
    lv = {fam: first.get(("LEIOMYOMA", fam), "") for fam in FAMILIES if ("LEIOMYOMA", fam) in first}
    say("  [R2] LEIOMYOMA first (no myometrium in these biopsies; descriptive): " + "; ".join("%s %s" % (f, l[:16]) for f, l in lv.items()))
    SUM["results"]["r2"] = {"flip": flip, "ctrl": ctrl, "first": {"%s_%s" % k: v for k, v in first.items()}}
    # O1 report
    say("")
    say("-" * 72)
    say("[SECTION] O1 candidate cells of origin: best rank of each group (KOLLER_EUR, GCST90205183, ENDO_ALL and controls)")
    best = {}
    for t in TRAITS:
        if t not in res:
            continue
        for fam in ("UPN", "UP", "SPEC"):
            it = res[t][fam]
            br = {}
            for i, (l, v, q) in enumerate(it):
                g = group_of(l)
                if g not in br:
                    br[g] = (i + 1, l, v["p"], q)
            best[(t, fam)] = br
            say("  %-12s %-4s %s" % (t[:12], fam, " | ".join("%s %d (%s %.0e)" % (g.split(" ")[0], br[g][0], br[g][1][:12], br[g][2]) for g, _ in GROUPS if g in br)))
    b = best.get(("KOLLER_EUR", "UPN"), {})
    rs = b.get("endometrial stroma", (99,))[0]
    others = [b[g][0] for g in ("perivascular / MSC-like", "mesothelium", "epithelium") if g in b]
    origin = bool(others) and all(rs < r for r in others)
    say("[H_ORIGIN_VERDICT] KOLLER_EUR UPN: best endometrial-stroma type rank %d vs perivascular / mesothelium / epithelium best ranks %s -> %s" % (
        rs, ", ".join(str(r) for r in others), "SUPPORTED (stroma ranks above every other candidate group)" if origin else "NOT SUPPORTED"))
    meso_rank = b.get("mesothelium", (None,))[0]
    say("  [O1] mesothelium (coelomic-metaplasia candidate) best rank under UPN: %s; label identity %s" % (
        meso_rank, "confirmed by markers" if meso_ok else "not confirmed -> not interpreted"))
    SUM["results"]["o1"] = {"supported": origin, "stroma_rank": rs, "other_ranks": others, "mesothelium_rank": meso_rank}
    # trait contrast of the best stromal type (UPN)
    if "endometrial stroma" in b:
        lab = b["endometrial stroma"][1]
        est = {}
        for t in TRAITS:
            d = A[(A["trait"] == t) & (A["family"] == "UPN") & (A["label"] == lab)]
            if len(d):
                r = d.iloc[0]
                ses = abs(r["se"] * r["beta_std"] / r["beta"]) if (np.isfinite(r["beta"]) and r["beta"] != 0) else np.nan
                est[t] = (r["beta_std"], ses)
        parts, pc = [], {}
        for a, c in CONTRAST_TRAITS:
            if a in est and c in est:
                z = (est[a][0] - est[c][0]) / math.sqrt(est[a][1] ** 2 + est[c][1] ** 2)
                p1 = 0.5 * math.erfc(z / math.sqrt(2))
                pc[c] = p1
                parts.append("%s>%s z %.2f p %.3f" % (a[:6], c[:6], z, p1))
        ok = pc.get("T2D", 1) < 0.05 and pc.get("HYPERTENSION", 1) < 0.05
        say("  [O1] trait contrast, UPN %s: %s" % (lab, " | ".join(parts)))
        say("[H_ORIGIN_SPECIFIC_VERDICT] best stromal type more enriched for KOLLER_EUR than T2D and hypertension -> %s" % ("SUPPORTED" if ok else "NOT SUPPORTED"))
        SUM["results"]["o1_contrast"] = {"label": lab, "p": pc, "supported": ok}


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.50 | independent endometrial atlas (GSE266265) and candidate cells of origin | " + TS)
    say("=" * 72)
    t0 = time.time()
    try:
        main_run(gene_loc())
    except Exception as e:
        import traceback
        ERR.append("main")
        say("[STEP_ERROR] %s: %s" % (type(e).__name__, str(e)[:300]))
        say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v050_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.40 before this run")
    if ERR:
        say("[INCOMPLETE] steps with errors: " + ", ".join(ERR))
    else:
        say("[DONE] PAPER9_V050_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
