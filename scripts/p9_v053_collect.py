#!/usr/bin/env python3
# Paper 9 | v0.0.53 | collect the manuscript figure set (no new analyses)
#   Copies F1-F6, S1 (v048), F7 (v049), F8 (v051) into one folder with combined legends and an input manifest.
# Runs in seconds; no network.

import datetime as dt
import shutil
from pathlib import Path

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RES = REPO / "results"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = RES / "paper_figures" / ("final_v053_" + TS)
LOG = REPO / "logs" / ("p9_v053_collect_" + TS + ".log")
OUT.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")

EXTRA_LEGENDS = {
    "F7_locus_genes_heca": ("Expression of the gene used for each KOLLER_EUR genome-wide locus (v042 assignment) across HECA cell states "
                            "(row z-score of log2(CPM + 1), control and endometriosis donors averaged), grouped by assigned compartment. "
                            "Genes at the nine stromal loci (WNT4, GREB1, IGF1 region, TBX3, WT1, TSHZ3, HOXA9-11, OSR2 region, SFRP1) peak "
                            "in proliferative / hormone-responsive and decidual endometrial stroma."),
    "F8_candidate_origin": ("Best rank of each candidate-origin group per GWAS in an endometrial atlas independent of HECA (GSE266265, 60 donors; "
                            "donor-paired pseudobulk; size-matched up-only marker sets, N = 400; MAGMA). Endometrial stroma (eF3; dS3 for FinnGen) "
                            "ranks first for all three endometriosis GWAS; perivascular cells are second; the 'mesothelial' label of this atlas is "
                            "not mesothelium by marker expression (see Results) and is not interpreted."),
}


def say(m=""):
    print(m, flush=True)
    _fh.write(m + "\n")
    _fh.flush()


def latest(base, pattern):
    fs = sorted(base.glob(pattern))
    return fs[-1] if fs else None


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.53 | collect manuscript figures | " + TS)
    say("=" * 72)
    manifest, legends = [], []
    d48 = latest(RES / "paper_figures", "p9_v048_*/figure_legends.md")
    if d48 is not None:
        d48 = d48.parent
        for f in sorted(d48.glob("*")):
            if f.suffix in (".png", ".pdf", ".csv"):
                shutil.copy2(f, OUT / f.name)
        legends.append(d48.joinpath("figure_legends.md").read_text().split("\n", 2)[-1])
        manifest.append("F1-F6, S1: %s" % d48.relative_to(RES))
        say("[COPY] F1-F6, S1 from %s" % d48.name)
    else:
        say("[MISS] v048 figures not found")
    for name, base, pat in (("F7_locus_genes_heca", RES / "review", "p9_v049_*/F7_locus_genes_heca.png"),
                            ("F8_candidate_origin", RES / "review", "p9_v051_*/F8_candidate_origin.png")):
        f = latest(base, pat)
        if f is None:
            say("[MISS] %s not found" % name)
            continue
        for ext in ("png", "pdf"):
            src = f.with_suffix("." + ext)
            if src.exists():
                shutil.copy2(src, OUT / src.name)
        for extra in f.parent.glob("*.csv"):
            if name.startswith("F7") and extra.name.startswith("o2_"):
                shutil.copy2(extra, OUT / ("F7_locus_genes_heca_data.csv"))
            if name.startswith("F8") and extra.name.startswith("F8_"):
                shutil.copy2(extra, OUT / extra.name)
        legends.append("## %s\n\n%s\n" % (name, EXTRA_LEGENDS[name]))
        manifest.append("%s: %s" % (name, f.parent.relative_to(RES)))
        say("[COPY] %s from %s" % (name, f.parent.name))
    (OUT / "figure_legends.md").write_text("# Paper 9 figure legends (v053 collection)\n\n" + "\n".join(legends))
    (OUT / "MANIFEST.txt").write_text("\n".join(manifest) + "\n")
    say("")
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V053_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
