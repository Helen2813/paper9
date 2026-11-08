#!/usr/bin/env python3
# Paper 9 | v0.0.46 | (rules CHANGELOG v0.2.36)
#   D1  composition of |log2FC| sets vs their significance (depleted-gene share vs -log10 p, per trait, both atlases)
#   D2  sign decomposition of |log2FC| hits with MAGMA (K7 rule: UP part | DOWN part and DOWN part | UP part)
#       for the top-5 |log2FC| hits per trait, the KOLLER key sets, and positive-control rankings
#   D3  H-MASK: why the |log2FC| rule misses the leiomyoma control (sets ranked above the expected group)
# Runs in ~/Desktop/paper9/.venv, low-heat (MAGMA single thread, nice 10). ~15-30 min. Uses cached v041 / v044 DE tables
# and the gnomAD constraint file downloaded by v045.

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
RAW = REPO / "data" / "raw"
INTERIM = REPO / "data" / "interim"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
GTEX_DE = PROC / "genesets" / "gtex_v044_de"
HECA_DE = PROC / "genesets" / "heca_koller_v041_de"
GNOMAD_FILE = RAW / "gnomad" / "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "specificity" / ("p9_v046_" + TS)
WORK = INTERIM / "magma_v046"
LOG = REPO / "logs" / ("p9_v046_sign_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.36)
LFC, PADJ, MIN_SET = 1.2, 0.01, 20
TOP_K = 5
MAX_ABOVE = 10
TRAITS = ["KOLLER_EUR", "GCST90205183", "ENDO_ALL", "LEIOMYOMA", "HYPERTENSION", "T2D"]
EXPECT = {"GTEx": {"LEIOMYOMA": r"smooth muscle", "HYPERTENSION": r"smooth muscle|endothel|pericyte|vascular|^contractile cell$"},
          "HECA": {"LEIOMYOMA": r"^(uSMCs|mPV|ePV)", "HYPERTENSION": r"^(Arterial|Venous|ePV|mPV|uSMCs)"}}
KEY_HECA = ["Venous/CTL", "Venous/ENDO", "Arterial/CTL", "Arterial/ENDO", "dStromal_early/ENDO", "dStromal_mid/ENDO", "dStromal_mid/CTL",
            "dStromal_early/CTL", "Luminal/ENDO", "Ciliated/ENDO"]

for d in (OUT, WORK, LOG.parent):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V046_SIGN", "timestamp_utc": TS, "results": {}}
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


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return np.nan
    return float(pd.Series(x[m]).rank().corr(pd.Series(y[m]).rank()))


def short(label, n=26):
    s = str(label).replace("SOX9_functionalis_", "SOX9f_").replace("dStromal_", "dS_")
    return s[:n]


def fp(p):
    return "%.1e" % p if np.isfinite(p) else "-"


def loeuf_table():
    if not GNOMAD_FILE.exists():
        say("[D] gnomAD constraint file missing (run v045 first); universe without the LOEUF filter")
        return None
    G = pd.read_csv(GNOMAD_FILE, sep="\t", compression="gzip")
    sc = next(c for c in ("gene", "gene_symbol") if c in G.columns)
    lc = next(c for c in ("oe_lof_upper", "LOEUF") if c in G.columns)
    return pd.to_numeric(G.set_index(sc)[lc], errors="coerce").dropna().groupby(level=0).min()


def atlas_tables(DE, CPM, sym2ent, loeuf):
    # Universe as v045 C4: CPM >= 1 in >= 1 group, mapped, with LOEUF
    C = CPM.copy()
    C.index = C.index.astype(str)
    C = C.groupby(level=0).mean()
    C = C[(C >= 1.0).any(axis=1)]
    uni = [g for g in C.index if g in sym2ent and (loeuf is None or g in loeuf.index)]
    uset = set(uni)
    parts = {}
    for g, d in DE.groupby("group", sort=False):
        s = d[(d["lfc"].abs() > LFC) & (d["padj"] < PADJ)]
        s = s[s["gene"].astype(str).isin(uset)]
        up = sorted(set(s.loc[s["lfc"] > 0, "gene"].astype(str)))
        dn = sorted(set(s.loc[s["lfc"] < 0, "gene"].astype(str)))
        if len(up) + len(dn) >= MIN_SET:
            parts[g] = (up, dn)
    ent = [sym2ent[g] for g in uni]
    return uni, ent, parts


def write_covar(path, uni, ent, cols):
    T = pd.DataFrame({k: np.isin(uni, v).astype(int) for k, v in cols.items()}, index=ent)
    T = T[~T.index.duplicated()]
    T.to_csv(path, sep="\t", index_label="GENE")


def ranked(res, names):
    items = sorted([(names[k], p) for k, p in (res or {}).items() if k in names], key=lambda x: x[1])
    return items


def classify(pu, pd_, nu, nd):
    if nu < MIN_SET or nd < MIN_SET:
        return "single-direction (%s)" % ("UP" if nu >= MIN_SET else "DOWN")
    return {(True, False): "UP-carried", (False, True): "DOWN-carried", (True, True): "both", (False, False): "neither"}[(pu < 0.05, pd_ < 0.05)]


def run_atlas(atlas, DE, CPM, sym2ent, loeuf):
    say("")
    say("-" * 72)
    say("[SECTION] %s" % atlas)
    uni, ent, parts = atlas_tables(DE, CPM, sym2ent, loeuf)
    labs = list(parts)
    names = {}
    cols = {}
    for i, g in enumerate(labs):
        up, dn = parts[g]
        cols["A%03d" % i], names["A%03d" % i] = up + dn, ("ABS", g)
        if len(up) >= MIN_SET:
            cols["U%03d" % i], names["U%03d" % i] = up, ("UP", g)
        if len(dn) >= MIN_SET:
            cols["D%03d" % i], names["D%03d" % i] = dn, ("DOWN", g)
    cf = WORK / ("%s_marginal.txt" % atlas)
    write_covar(cf, uni, ent, cols)
    fd = {g: len(parts[g][1]) / float(len(parts[g][0]) + len(parts[g][1])) for g in labs}
    say("[D] %s: %d genes; %d |log2FC| sets (median depleted share %.2f)" % (atlas, len(set(ent)), len(labs), float(np.median(list(fd.values())))))
    rows, cache = [], {}
    t0 = time.time()
    for t in TRAITS:
        raw = MAGMA_OUT / (t + "_35_10.genes.raw")
        if not raw.exists():
            say("[D] %s genes.raw missing" % t)
            continue
        r = magma_covar(raw, cf, "%s_%s_marg" % (atlas, t), [])
        if not r:
            continue
        marg = {}
        for k, p in r.items():
            if k in names:
                marg[names[k]] = p
        absr = sorted([(g, marg.get(("ABS", g), np.nan)) for g in labs], key=lambda x: (np.nan_to_num(x[1], nan=2.0)))
        rank = {g: i + 1 for i, (g, p) in enumerate(absr)}
        # D1
        rho = spearman([fd[g] for g, p in absr], [-np.log10(max(p, 1e-300)) if np.isfinite(p) else np.nan for g, p in absr])
        say("  [D1] %-12s rho(depleted share, -log10 p ABS) across %d sets: %+.2f" % (t, len(absr), rho))
        SUM["results"].setdefault("d1_" + atlas, {})[t] = rho
        # analysed labels
        todo = [g for g, p in absr[:TOP_K]]
        why = {g: "top%d" % rank[g] for g in todo}
        if t in EXPECT[atlas]:
            hit = [g for g, p in absr if re.search(EXPECT[atlas][t], g, re.I)]
            if hit:
                e = hit[0]
                for g, p in absr[:min(rank[e], MAX_ABOVE + 1)]:
                    if g not in todo:
                        todo.append(g)
                        why[g] = "above-expected"
                why[e] = "expected (ABS rank %d)" % rank[e]
        if atlas == "HECA" and t == "KOLLER_EUR":
            for g in KEY_HECA:
                if g in parts and g not in todo:
                    todo.append(g)
                    why[g] = "key"
        for g in todo:
            up, dn = parts[g]
            pu = pd_ = np.nan
            if len(up) >= MIN_SET and len(dn) >= MIN_SET:
                ck = (t, g)
                if ck not in cache:
                    jf = WORK / ("%s_joint.txt" % atlas)
                    write_covar(jf, uni, ent, {"U": up, "D": dn})
                    ru = magma_covar(raw, jf, "%s_%s_j_u" % (atlas, t), ["D"])
                    rd = magma_covar(raw, jf, "%s_%s_j_d" % (atlas, t), ["U"])
                    cache[ck] = ((ru or {}).get("U", np.nan), (rd or {}).get("D", np.nan))
                pu, pd_ = cache[ck]
            pa = marg.get(("ABS", g), np.nan)
            cls = classify(pu, pd_, len(up), len(dn)) if pa < 0.05 else "ABS not nominal"
            rows.append({"atlas": atlas, "trait": t, "label": g, "why": why[g], "is_top": rank[g] <= TOP_K, "abs_rank": rank[g], "p_abs": pa,
                         "p_up": marg.get(("UP", g), np.nan), "p_down": marg.get(("DOWN", g), np.nan), "p_up_given_down": pu,
                         "p_down_given_up": pd_, "n_up": len(up), "n_down": len(dn), "depleted_share": fd[g], "class": cls})
    say("[D] %s MAGMA runs done (%.0f min)" % (atlas, (time.time() - t0) / 60))
    return pd.DataFrame(rows)


def report(R):
    say("")
    say("-" * 72)
    say("[SECTION] D2  sign decomposition of |log2FC| hits (K7 rule, MAGMA)")
    res = {}
    for atlas, A in R.groupby("atlas"):
        for t in TRAITS:
            d = A[(A["trait"] == t) & A["is_top"]].sort_values("abs_rank")
            if len(d):
                say("  [D2] %s %-12s top %d: %s" % (atlas, t, len(d), "; ".join("%s %s (%.0f%% down) %s" % (
                    short(r["label"]), fp(r["p_abs"]), 100 * r["depleted_share"], r["class"]) for _, r in d.iterrows())))
        top = A[A["is_top"] & (A["p_abs"] < 0.05)]
        cnt = top["class"].value_counts().to_dict()
        share = (cnt.get("DOWN-carried", 0) + cnt.get("single-direction (DOWN)", 0)) / max(1, len(top))
        say("  [D2] %s top-%d nominal |log2FC| hits pooled over traits (n %d): %s -> DOWN-carried share %.2f" % (
            atlas, TOP_K, len(top), ", ".join("%s %d" % (k, v) for k, v in sorted(cnt.items())), share))
        res[atlas] = {"share_down": share, "counts": cnt, "n": int(len(top))}
        SUM["results"]["d2_" + atlas] = res[atlas]
        if atlas == "HECA":
            k = A[(A["trait"] == "KOLLER_EUR") & A["label"].isin(KEY_HECA)].drop_duplicates("label")
            for _, r in k.sort_values("abs_rank").iterrows():
                say("    [D2_KEY] KOLLER_EUR %-20s ABS %s (rank %d, %.0f%% down) | UP %s (|DOWN %s) | DOWN %s (|UP %s) -> %s" % (
                    short(r["label"], 20), fp(r["p_abs"]), r["abs_rank"], 100 * r["depleted_share"], fp(r["p_up"]), fp(r["p_up_given_down"]),
                    fp(r["p_down"]), fp(r["p_down_given_up"]), r["class"]))
    ok = all(v["share_down"] >= 0.5 for v in res.values()) and len(res) == 2
    say("[H_K7REP_VERDICT] top |log2FC| hits mostly DOWN-carried in both atlases -> %s (%s)" % (
        "SUPPORTED" if ok else "NOT SUPPORTED", ", ".join("%s %.2f" % (a, v["share_down"]) for a, v in res.items())))
    SUM["results"]["h_k7rep"] = ok
    say("")
    say("-" * 72)
    say("[SECTION] D3  positive controls under the |log2FC| rule: masking of the expected group and DOWN-carried sets ranked above it")
    mask = {}
    for atlas, A in R.groupby("atlas"):
        for t in ("LEIOMYOMA", "HYPERTENSION"):
            d = A[A["trait"] == t]
            e = d[d["why"].str.startswith("expected")]
            if not len(e):
                continue
            e = e.iloc[0]
            above = d[(d["abs_rank"] < e["abs_rank"])].sort_values("abs_rank")
            ac = above["class"].value_counts().to_dict()
            nd = int(ac.get("DOWN-carried", 0) + ac.get("single-direction (DOWN)", 0))
            say("  [D3] %s %-12s expected %s: ABS %s (rank %d, %.0f%% down) | UP part %s (|DOWN %s) | DOWN part %s (|UP %s)" % (
                atlas, t, short(e["label"]), fp(e["p_abs"]), e["abs_rank"], 100 * e["depleted_share"], fp(e["p_up"]), fp(e["p_up_given_down"]),
                fp(e["p_down"]), fp(e["p_down_given_up"])))
            for _, r in above.iterrows():
                say("       above: %-26s ABS %s (%.0f%% down) | UP|DOWN %s | DOWN|UP %s -> %s" % (
                    short(r["label"]), fp(r["p_abs"]), 100 * r["depleted_share"], fp(r["p_up_given_down"]), fp(r["p_down_given_up"]), r["class"]))
            if t == "LEIOMYOMA":
                mask[atlas] = bool(e["abs_rank"] > 1 and len(above) and nd >= len(above) / 2.0)
                say("  [D3] %s LEIOMYOMA: expected group not first under |log2FC| and >= half of the sets above it DOWN-carried (%d/%d) -> %s" % (
                    atlas, nd, len(above), "YES" if mask[atlas] else "NO"))
    ok = len(mask) == 2 and all(mask.values())
    say("[H_MASK_VERDICT] leiomyoma control missed by |log2FC| because DOWN-carried sets outrank the expected group, in both atlases -> %s" % (
        "SUPPORTED" if ok else "NOT SUPPORTED"))
    SUM["results"]["h_mask"] = {"supported": ok, **mask}


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.46 | composition and sign decomposition of |log2FC| hits; positive-control masking | " + TS)
    say("=" * 72)
    t0 = time.time()
    sym2ent = gene_loc()
    R = []
    try:
        loeuf = loeuf_table()
        for atlas, ddir, cpmf in (("HECA", HECA_DE, "mean_cpm_by_group.csv.gz"), ("GTEx", GTEX_DE, "mean_cpm_by_type.csv.gz")):
            if not ((ddir / "de_all.csv.gz").exists() and (ddir / cpmf).exists()):
                say("[D] %s DE tables missing" % atlas)
                continue
            try:
                R.append(run_atlas(atlas, pd.read_csv(ddir / "de_all.csv.gz"), pd.read_csv(ddir / cpmf, index_col=0), sym2ent, loeuf))
            except Exception as e:
                import traceback
                ERR.append(atlas)
                say("[STEP_ERROR] %s %s: %s" % (atlas, type(e).__name__, str(e)[:300]))
                say(traceback.format_exc()[-1500:])
        if R:
            RR = pd.concat(R, ignore_index=True)
            RR.to_csv(OUT / "d_sign_decomposition.csv", index=False)
            report(RR)
    except Exception as e:
        import traceback
        ERR.append("report")
        say("[STEP_ERROR] %s: %s" % (type(e).__name__, str(e)[:300]))
        say(traceback.format_exc()[-1500:])
    with open(OUT / "p9_v046_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    say("[TIME] %.0f min" % ((time.time() - t0) / 60))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] rules fixed in CHANGELOG v0.2.36 before this run (D1-D3 formulated after v045; K7 rule from v0.2.23)")
    if ERR:
        say("[INCOMPLETE] steps with errors: " + ", ".join(ERR))
    else:
        say("[DONE] PAPER9_V046_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
