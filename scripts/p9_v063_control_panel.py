#!/usr/bin/env python3
# Paper 9 | v0.0.63 | known-source control panel in GTEx v9, completed (rules CHANGELOG v0.2.52; deviations v0.2.53)
#   Same design as v062. Fixes: (1) p-values below 1e-300 are written as 1e-300 (MAGMA cannot read subnormal values;
#   AF and Dupuytren failed in v062); (2) two label rules corrected after the v062 label list was printed, reported next
#   to the registered rules: lymphocytes use word boundaries ("T cell" also matched "mast cell"); prostate epithelium also
#   matches "prostatic duct" ("basal epithelial cell of prostatic duct" was missed).
#   S6 figure: rank of the expected group (corrected rules), trait x family.
# Re-uses downloaded FinnGen files and finished gene analyses; runs MAGMA gene analysis for AF and Dupuytren
# (30-40 min each, single thread, nice 10), then the set analyses (~10 min).

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
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RAW = REPO / "data" / "raw"
INTERIM = REPO / "data" / "interim"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
MAGMA = TOOLS / "magma" / "magma_bin"
G1000 = TOOLS / "magma" / "ref_eur" / "g1000_eur"
GENE_LOC_B37 = TOOLS / "magma" / "genes_b37" / "NCBI37.3.gene.loc"
MAGMA_OUT = INTERIM / "magma"
ANNOT = MAGMA_OUT / "annot_b37_35_10.genes.annot"
FINNGEN_DIR = RAW / "finngen"
FINNGEN_MANIFEST = FINNGEN_DIR / "finngen_R13_manifest.tsv"
GTEX_DE = PROC / "genesets" / "gtex_v044_de"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "external_atlas" / ("p9_v063_" + TS)
WORK = INTERIM / "magma_v063"
LOG = REPO / "logs" / ("p9_v063_control_panel_" + TS + ".log")
NICE = ["nice", "-n", "10"]
ENV1 = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# Fixed (CHANGELOG v0.2.52)
LFC, PADJ, MIN_SET, N_MATCH = 1.2, 0.01, 20, 400            # as v044
FAMILIES = ["ABS", "UP", "DOWN", "ABSN", "UPN", "DOWNN", "SPEC"]
FAM_NAME = {"ABS": "|log2FC|", "UP": "up-only", "DOWN": "down-only", "ABSN": "|log2FC|, N=400", "UPN": "up-only, N=400",
            "DOWNN": "down-only, N=400", "SPEC": "specificity"}
NEW = {  # trait: (candidate FinnGen codes in order, name regex fallback)
    "AF": (["I9_AF", "I9_AF_REIMB"], r"atrial fibrillation"),
    "DUPUYTREN": (["M13_DUPUTRYEN", "M13_DUPUTRENS", "M13_DUPUYTREN"], r"dupuytren"),
    "PROSTATE_CA": (["C3_PROSTATE_EXALLC", "C3_PROSTATE"], r"malignant neoplasm of prostate"),
    "AIHYPO": (["E4_HYTHY_AI_STRICT", "E4_HYTHY_AI_DX", "E4_HYTHY_AI_EXCL"], r"hypothyroidism, autoimmune|autoimmune hypothyroid"),
}
EXPECT_REG = {  # registered in v0.2.52: (regex on GTEx v9 cell-type label, description)
    "AF": (r"cardiac muscle|cardiomyocyte", "cardiomyocytes"),
    "DUPUYTREN": (r"fibroblast", "fibroblasts"),
    "PROSTATE_CA": (r"prostate epithelium|epithelial cell of prostate|prostate gland epithel", "prostate epithelium"),
    "AIHYPO": (r"T cell|natural killer|B cell|plasma cell|lymphocyte", "lymphocytes (T, NK, B, plasma)"),
    "LEIOMYOMA": (r"smooth muscle", "smooth muscle"),
    "HYPERTENSION": (r"smooth muscle|endothel|pericyte|vascular|^contractile cell$", "vascular cells (v045 label rule)"),
}
EXPECT = dict(EXPECT_REG)  # corrected in v0.2.53 (primary from v063; registered rules reported alongside)
EXPECT["PROSTATE_CA"] = (r"prostate epithelium|epithelial cell of prostate|prostatic duct", "prostate epithelium")
EXPECT["AIHYPO"] = (r"\bT cell|natural killer|\bNK\b|\bB cell|plasma cell|lymphocyte", "lymphocytes (T, NK, B, plasma)")
P_FLOOR = 1e-300
TRAIT_NAME = {"AF": "Atrial fibrillation", "DUPUYTREN": "Dupuytren contracture", "PROSTATE_CA": "Prostate cancer",
              "AIHYPO": "Autoimmune hypothyroidism", "LEIOMYOMA": "Uterine leiomyoma", "HYPERTENSION": "Hypertension"}
MIN_FREE_GB = 5.0
MIN_TRAITS = 5             # verdicts need at least 5 evaluable traits

for d in (OUT, WORK, LOG.parent, FINNGEN_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V063_CONTROL_PANEL", "timestamp_utc": TS, "results": {}}


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


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


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


def magma_sets(raw, sf, tag):
    outp = WORK / tag
    gsa = Path(str(outp) + ".gsa.out")
    o = e = ""
    for model in (["--model", "condition-hide=ALLGENES"], []):
        if gsa.exists():
            gsa.unlink()
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--set-annot", sf] + model + ["--out", outp], timeout=3600)
        if gsa.exists():
            return parse_gsa(gsa)
    say("    [MAGMA_FAIL] sets %s" % (o + e).strip()[-200:])
    return None


def magma_covar(raw, cf, tag):
    outp = WORK / tag
    gsa = Path(str(outp) + ".gsa.out")
    o = e = ""
    for key in ("condition-hide=", "condition="):
        if gsa.exists():
            gsa.unlink()
        rc, o, e = sh(NICE + [MAGMA, "--gene-results", raw, "--gene-covar", cf, "--model", key + "Average", "direction=greater",
                              "--out", outp], timeout=3600)
        if gsa.exists():
            return parse_gsa(gsa)
    say("    [MAGMA_FAIL] covar %s" % (o + e).strip()[-200:])
    return None


# ----------------------------------------------------------------------------
# FinnGen endpoints and MAGMA gene analysis
# ----------------------------------------------------------------------------
def manifest():
    lines = FINNGEN_MANIFEST.read_text().splitlines()
    hdr = lines[0].split("\t")
    rows = [dict(zip(hdr, l.split("\t"))) for l in lines[1:] if l.strip()]
    code_col = hdr[0]
    return rows, code_col


def resolve(trait, rows, code_col):
    cands, rx = NEW[trait]
    by_code = {r[code_col]: r for r in rows}
    for c in cands:
        if c in by_code:
            return by_code[c], "code"
    hits = [r for r in rows if re.search(rx, " ".join(str(v) for v in r.values()), re.I)]
    say("  [FG] %s: none of %s in the manifest; name matches: %s" % (
        trait, ",".join(cands), ", ".join("%s (%s cases)" % (r[code_col], r.get("num_cases", "?")) for r in hits[:8]) or "none"))
    if not hits:
        return None, "none"
    best = max(hits, key=lambda r: int(r.get("num_cases", "0") or 0))
    return best, "name"


def download(url, dst):
    if dst.exists() and dst.stat().st_size > 1_000_000:
        say("  [DL] present %s (%.0f MB)" % (dst.name, dst.stat().st_size / 1e6))
        return True
    part = dst.with_name(dst.name + ".part")
    t0 = time.time()
    rc, o, e = sh(["curl", "-L", "-sS", "--fail", "--retry", "3", "--connect-timeout", "30", "-C", "-", "-A", "paper9", "-o", part, url],
                  timeout=4 * 3600)
    if rc == 0 and part.exists() and part.stat().st_size > 1_000_000:
        part.rename(dst)
        say("  [DL] OK %s %.0f MB in %.0f min" % (dst.name, dst.stat().st_size / 1e6, (time.time() - t0) / 60))
        return True
    say("  [DL] FAILED rc=%s %s" % (rc, (e or o).strip()[:200]))
    return False


def gene_analysis(trait, row, code_col):
    raw = MAGMA_OUT / (trait + "_35_10.genes.raw")
    if raw.exists():
        say("  [MAGMA_GENE] %s present" % trait)
        return raw
    url = row.get("path_https", "")
    if not url:
        say("  [FG] %s: no path_https in manifest" % trait)
        return None
    src = FINNGEN_DIR / url.split("/")[-1]
    if not download(url, src):
        return None
    rc, _, _ = sh(["gzip", "-t", src], timeout=1800)
    if rc != 0:
        say("  [DL] %s gzip test failed; file removed, re-run to download again" % src.name)
        src.unlink()
        return None
    pfile = MAGMA_OUT / (trait + ".pval.txt")
    n = 0
    tmp = pfile.with_name(pfile.name + ".part")
    with gzip.open(src, "rt") as f, open(tmp, "w") as o:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {h.lstrip("#"): i for i, h in enumerate(hdr)}
        o.write("SNP\tP\n")
        for line in f:
            t = line.rstrip("\n").split("\t")
            rs = t[ix["rsids"]].split(",")[0]
            if rs.startswith("rs"):
                try:
                    pv = float(t[ix["pval"]])
                except ValueError:
                    continue
                if not (pv == pv):
                    continue
                o.write(rs + "\t" + (repr(P_FLOOR) if pv < P_FLOOR else t[ix["pval"]]) + "\n")
                n += 1
    tmp.rename(pfile)
    ncas, ncon = int(row["num_cases"]), int(row["num_controls"])
    say("  [MAGMA_GENE] %s: %d SNPs, N = %d (%d cases); running (20-40 min, quiet)" % (trait, n, ncas + ncon, ncas))
    t0 = time.time()
    rc, o, e = sh(NICE + [MAGMA, "--bfile", G1000, "--pval", pfile, "N=%d" % (ncas + ncon), "--gene-annot", ANNOT,
                          "--out", MAGMA_OUT / (trait + "_35_10")], timeout=8 * 3600)
    say("  [MAGMA_GENE] %s rc=%d %.0f min" % (trait, rc, (time.time() - t0) / 60))
    if raw.exists():
        pfile.unlink()
        return raw
    say("  [MAGMA_GENE] %s failed: %s" % (trait, (o + e).strip()[-300:]))
    return None


# ----------------------------------------------------------------------------
# Marker families from the cached v044 GTEx DE
# ----------------------------------------------------------------------------
def families(DE, CPM, sym2ent):
    fams, meta = {}, {}
    for g, d in DE.groupby("group", sort=False):
        sig = d["padj"] < PADJ
        sel = {"ABS": d.loc[(d["lfc"].abs() > LFC) & sig, "gene"], "UP": d.loc[(d["lfc"] > LFC) & sig, "gene"],
               "DOWN": d.loc[(d["lfc"] < -LFC) & sig, "gene"],
               "ABSN": d.reindex(d["score"].abs().sort_values(ascending=False).index)["gene"].head(N_MATCH),
               "UPN": d.sort_values("score", ascending=False)["gene"].head(N_MATCH),
               "DOWNN": d.sort_values("score", ascending=True)["gene"].head(N_MATCH)}
        for fam, genes in sel.items():
            gl = sorted({sym2ent[x] for x in genes.astype(str) if x in sym2ent})
            if len(gl) >= MIN_SET:
                key = ("%s_%s" % (fam, re.sub(r"[^A-Za-z0-9]+", "_", g)))[:60]
                fams[key], meta[key] = gl, (fam, g)
    fams["ALLGENES"] = sorted({sym2ent[x] for x in DE["gene"].astype(str) if x in sym2ent})
    C = CPM[(CPM >= 1.0).any(axis=1)]
    C = C.loc[[g for g in C.index if g in sym2ent]]
    S = C.div(C.sum(axis=1), axis=0)
    names = {"C%03d" % (i + 1): c for i, c in enumerate(S.columns)}
    S.columns = list(names)
    S["Average"] = np.log2(C.mean(axis=1) + 1.0)
    S.index = [sym2ent[g] for g in S.index]
    S = S[~S.index.duplicated()]
    return fams, meta, S, names


def depleted_share(DE, label):
    d = DE[DE["group"] == label]
    a = d[(d["lfc"].abs() > LFC) & (d["padj"] < PADJ)]
    return (float((a["lfc"] < 0).mean()) if len(a) else float("nan")), len(a)


# ----------------------------------------------------------------------------
# S6 figure
# ----------------------------------------------------------------------------
def s6(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    INK, INK2, MUTED, AXIS = "#0b0b0b", "#52514e", "#898781", "#c3c2b7"
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7,
                         "axes.titlesize": 7.5, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.edgecolor": AXIS,
                         "text.color": INK, "axes.titlecolor": INK, "xtick.color": MUTED, "ytick.color": INK2,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "figure.facecolor": "white", "savefig.facecolor": "white"})
    traits = [t for t in EXPECT if t in R["trait"].unique()]
    M = R.pivot(index="trait", columns="family", values="rank").reindex(index=traits, columns=FAMILIES)
    Q = R.pivot(index="trait", columns="family", values="q").reindex(index=traits, columns=FAMILIES)
    cls = M.copy()
    cls[:] = 3
    cls[M <= 5] = 2
    cls[M <= 2] = 1
    cls[M == 1] = 0
    cls[M.isna()] = 3
    cmap = ListedColormap(["#eb6834", "#f5b48f", "#fbe3d4", "#f0efec"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], 4)
    fig, ax = plt.subplots(figsize=(5.6, 0.36 * len(traits) + 1.3))
    ax.imshow(cls.values.astype(float), cmap=cmap, norm=norm, aspect="auto")
    for i in range(len(traits)):
        for j in range(len(FAMILIES)):
            r, q = M.iat[i, j], Q.iat[i, j]
            txt = "-" if not np.isfinite(r) else ("%d" % r) + ("*" if (r == 1 and np.isfinite(q) and q < 0.05) else "")
            ax.text(j, i, txt, ha="center", va="center", fontsize=6.6, color=INK if (np.isfinite(r) and r <= 2) else INK2)
    ax.set_xticks(range(len(FAMILIES)))
    ax.set_xticklabels([FAM_NAME[f] for f in FAMILIES], rotation=35, ha="right")
    ax.set_yticks(range(len(traits)))
    ax.set_yticklabels(["%s → %s" % (TRAIT_NAME[t], EXPECT[t][1].split(" (")[0]) for t in traits])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    ax.set_xticks(np.arange(-0.5, len(FAMILIES), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(traits), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", length=0)
    ax.set_title("Rank of the expected cell group (GTEx v9; 1 = first)", loc="left", fontweight="bold")
    for ext in ("png", "pdf"):
        fig.savefig(OUT / ("S6_control_panel.%s" % ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    M.to_csv(OUT / "S6_control_panel_data.csv")
    leg = ("Known-source control panel in the GTEx v9 snRNA-seq atlas (Eraslan et al. 2022; 45 cell types). Each row is a FinnGen "
           "R13 GWAS with a cell group fixed before the run; each cell gives the rank of the best label of that group under one "
           "marker family (MAGMA; * first with BH q < 0.05; - set not formed or absent). Label rules for prostate epithelium and lymphocytes as corrected after the label list was printed (registered-rule ranks are given in the source data). Thresholded |log2FC| > 1.2 and adjusted "
           "p < 0.01 as in Koller et al.; N = 400 families are size-matched; specificity is the continuous share of expression.")
    (OUT / "S6_legend.md").write_text("## S6_control_panel\n\n" + leg + "\n", encoding="utf-8")
    say("[FIG] S6_control_panel written")


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.63 | known-source control panel, completed (GTEx v9; six FinnGen traits) | " + TS)
    say("=" * 72)
    free = shutil.disk_usage(str(REPO)).free / 1e9
    say("[DISK] free %.1f GB" % free)
    need = [t for t in NEW if not (MAGMA_OUT / (t + "_35_10.genes.raw")).exists()]
    if need and free < MIN_FREE_GB:
        say("[FATAL] less than %.0f GB free; free some space and re-run" % MIN_FREE_GB)
        _fh.close()
        return
    for f in (GTEX_DE / "de_all.csv.gz", GTEX_DE / "mean_cpm_by_type.csv.gz", FINNGEN_MANIFEST, ANNOT):
        if not f.exists():
            say("[FATAL] missing input " + str(f))
            _fh.close()
            return
    sym2ent = gene_loc()
    DE = pd.read_csv(GTEX_DE / "de_all.csv.gz")
    CPM = pd.read_csv(GTEX_DE / "mean_cpm_by_type.csv.gz", index_col=0)
    labels = sorted(DE["group"].astype(str).unique())
    say("[GTEX] %d cell types: %s" % (len(labels), "; ".join(labels)))
    evaluable = {}
    for t, (rx, desc) in EXPECT.items():
        m = [l for l in labels if re.search(rx, l, re.I)]
        m0 = [l for l in labels if re.search(EXPECT_REG[t][0], l, re.I)]
        evaluable[t] = bool(m)
        say("[EXPECT] %-13s %-32s labels: %s%s" % (t, desc, ", ".join(m) or "NONE (trait not evaluable)",
                                                   "" if m == m0 else " | registered rule: " + ", ".join(m0)))

    say("")
    say("-" * 72)
    say("[SECTION] FinnGen R13 endpoints and MAGMA gene analysis")
    rows, code_col = manifest()
    chosen = {}
    for t in NEW:
        r, how = resolve(t, rows, code_col)
        if r is None:
            say("  [FG] %s: no endpoint found; trait skipped" % t)
            continue
        chosen[t] = r[code_col]
        say("  [FG] %-12s %-24s cases %s controls %s (%s)" % (t, r[code_col], r.get("num_cases"), r.get("num_controls"),
                                                              "pre-registered code" if how == "code" else "name fallback"))
        gene_analysis(t, r, code_col)
    SUM["results"]["endpoints"] = chosen

    say("")
    say("-" * 72)
    say("[SECTION] Marker families and rankings")
    fams, meta, S, names = families(DE, CPM, sym2ent)
    say("[GTEX] sets: %s; continuous covariates %d" % (", ".join("%s %d" % (f, sum(1 for k in meta if meta[k][0] == f))
                                                               for f in FAMILIES if f != "SPEC"), len(names)))
    sf, cf = WORK / "gtex_sets.txt", WORK / "gtex_covar.txt"
    with open(sf, "w") as o:
        for k, g in fams.items():
            o.write(k + " " + " ".join(g) + "\n")
    S.to_csv(cf, sep="\t", index_label="GENE")
    RANKED, top_rows = {}, []
    for t in EXPECT:
        raw = MAGMA_OUT / (t + "_35_10.genes.raw")
        if not raw.exists() or not evaluable[t]:
            say("  [RANK] %s not evaluable (%s)" % (t, "genes.raw missing" if not raw.exists() else "no expected label"))
            continue
        rs, rc = magma_sets(raw, sf, "sets_" + t), magma_covar(raw, cf, "cov_" + t)
        say("  %s (%s)" % (TRAIT_NAME[t], EXPECT[t][1]))
        for fam in FAMILIES:
            if fam == "SPEC":
                items = [(names[k], p) for k, p in (rc or {}).items() if k in names]
            else:
                items = [(meta[k][1], p) for k, p in (rs or {}).items() if k in meta and meta[k][0] == fam]
            if not items:
                continue
            items.sort(key=lambda x: x[1])
            q = bh(np.array([p for _, p in items]))
            ranked = [(l, p, qq) for (l, p), qq in zip(items, q)]
            RANKED[(t, fam)] = ranked
            for i, (l, p, qq) in enumerate(ranked[:5]):
                top_rows.append({"trait": t, "family": fam, "rank": i + 1, "label": l, "p": p, "q": qq})
            hit = [(i + 1, l, p, qq) for i, (l, p, qq) in enumerate(ranked) if re.search(EXPECT[t][0], str(l), re.I)]
            best = hit[0] if hit else (np.nan, "-", np.nan, np.nan)
            say("    %-6s expected rank %-3s (%s, p %.1e, q %.2f) | first: %s %.1e" % (
                fam, "-" if not np.isfinite(best[0]) else int(best[0]), str(best[1])[:30], best[2] if np.isfinite(best[2]) else np.nan,
                best[3] if np.isfinite(best[3]) else np.nan, str(ranked[0][0])[:30], ranked[0][1]))
    pd.DataFrame(top_rows).to_csv(OUT / "panel_top5.csv", index=False)
    if not RANKED:
        say("[FATAL] no trait evaluable")
        _fh.close()
        return

    def ranks(rules):
        rows = []
        for (t, fam), ranked in RANKED.items():
            hit = [(i + 1, l, p, qq) for i, (l, p, qq) in enumerate(ranked) if re.search(rules[t][0], str(l), re.I)]
            best = hit[0] if hit else (np.nan, "-", np.nan, np.nan)
            rows.append({"trait": t, "family": fam, "rank": best[0], "label": best[1], "p": best[2], "q": best[3],
                         "n_labels": len(ranked), "top_label": ranked[0][0], "top_p": ranked[0][1], "top_q": ranked[0][2]})
        return pd.DataFrame(rows)

    def verdicts(R, tag):
        traits = list(dict.fromkeys(R["trait"]))
        hits = {f: int(((R["family"] == f) & (R["rank"] == 1)).sum()) for f in FAMILIES}
        n_t = len(traits)
        need_up = int(np.ceil(2 * n_t / 3))
        v1 = ("NOT EVALUABLE" if n_t < MIN_TRAITS else
              "SUPPORTED" if (hits["UP"] >= need_up and hits["UP"] > hits["ABS"]) else "NOT SUPPORTED")
        v2 = ("NOT EVALUABLE" if n_t < MIN_TRAITS else "SUPPORTED" if hits["UPN"] > hits["ABSN"] else "NOT SUPPORTED")
        rk = lambda t, f: np.nan_to_num(R[(R.trait == t) & (R.family == f)]["rank"].min(), nan=99)
        better = sum(1 for t in traits if rk(t, "UP") < rk(t, "ABS"))
        worse = sum(1 for t in traits if rk(t, "UP") > rk(t, "ABS"))
        say("[PANEL%s] traits %d; expected group first per family: %s" % (tag, n_t, ", ".join("%s %d" % (f, hits[f]) for f in FAMILIES)))
        say("[H_PANEL_VERDICT%s] up-only first in %d of %d traits (%d required), |log2FC| first in %d -> %s" % (
            tag, hits["UP"], n_t, need_up, hits["ABS"], v1))
        say("[H_PANEL_N_VERDICT%s] size-matched up-only first in %d, size-matched |log2FC| first in %d -> %s" % (
            tag, hits["UPN"], hits["ABSN"], v2))
        say("[PANEL%s] expected group ranked higher with up-only than with |log2FC| in %d traits, lower in %d" % (tag, better, worse))
        return {"hits": hits, "H_PANEL": v1, "H_PANEL_N": v2, "up_better": better, "up_worse": worse}

    say("")
    say("-" * 72)
    say("[SECTION] Verdicts")
    R0, R = ranks(EXPECT_REG), ranks(EXPECT)
    R0.to_csv(OUT / "panel_expected_ranks_registered_rules.csv", index=False)
    R.to_csv(OUT / "panel_expected_ranks.csv", index=False)
    for t in dict.fromkeys(R["trait"]):
        r = R[R["trait"] == t].set_index("family")["rank"]
        r0 = R0[R0["trait"] == t].set_index("family")["rank"]
        lab = [l for f in ("UP", "ABS") for l in R[(R["trait"] == t) & (R["family"] == f)]["label"] if l != "-"]
        dep, n = depleted_share(DE, lab[0]) if lab else (np.nan, 0)
        diff = [f for f in FAMILIES if np.nan_to_num(r.get(f, np.nan), nan=99) != np.nan_to_num(r0.get(f, np.nan), nan=99)]
        say("  [PANEL] %-13s ranks %s | depleted share of the |log2FC| set of the expected label %.2f (%d genes)%s" % (
            t, " ".join("%s:%s" % (f, "-" if not np.isfinite(r.get(f, np.nan)) else int(r.get(f))) for f in FAMILIES), dep, n,
            ("" if not diff else " | registered rule differs: " + ", ".join(
                "%s %s" % (f, "-" if not np.isfinite(r0.get(f, np.nan)) else int(r0.get(f))) for f in diff))))
    SUM["results"]["registered_rules"] = verdicts(R0, "_REGISTERED")
    SUM["results"]["corrected_rules"] = verdicts(R, "")
    try:
        s6(R)
    except Exception as e:
        say("[S6] figure failed: %s: %s" % (type(e).__name__, e))
    (OUT / "summary.json").write_text(json.dumps(SUM, indent=1, default=str), encoding="utf-8")
    say("")
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V063_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
