#!/usr/bin/env python3
# Paper 9 | v0.0.69 | sample sizes of the Koller ancestry GWAS used in v068 (no analyses)
#   Prints the README of the Koller Zenodo record in full, the maximum and median per-SNP N of each summary-statistics
#   file, and the per-gene N used by MAGMA (genes.out). Runs in 1-3 minutes; no downloads.

import datetime as dt
import gzip
from pathlib import Path

import pandas as pd

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
KDIR = REPO / "data" / "raw" / "gwas" / "Koller2026"
MAGMA_OUT = REPO / "data" / "interim" / "magma"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
LOG = REPO / "logs" / ("p9_v069_ancestry_n_" + TS + ".log")
LOG.parent.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
FILES = {"EUR": "endometriosis_EUR_wo_23andMe.txt.gz", "EAS": "endometriosis_EAS.txt.gz", "AFR": "endometriosis_AFR.txt.gz",
         "AMR": "endometriosis_AMR.txt.gz"}


def say(m=""):
    print(m, flush=True)
    _fh.write(m + "\n")
    _fh.flush()


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.69 | Koller ancestry GWAS sample sizes | " + TS)
    say("=" * 72)
    rd = KDIR / "README_GWAS.txt"
    if rd.exists():
        say("[README]")
        for line in rd.read_text(errors="ignore").splitlines():
            if line.strip():
                say("  " + line.rstrip()[:220])
    for anc, name in FILES.items():
        f = KDIR / name
        if not f.exists():
            say("[N] %s file missing" % anc)
            continue
        cols = pd.read_csv(f, sep=None, engine="python", nrows=5).columns
        ncols = [c for c in cols if c.lower() in ("n", "neff", "n_eff", "ncase", "n_case", "ncases", "ncontrol", "n_control", "ncontrols", "weight")]
        if not ncols:
            say("[N] %s no N column (%s)" % (anc, ", ".join(cols)))
            continue
        sep = "\t"
        with gzip.open(f, "rt") as fh:
            sep = "\t" if "\t" in fh.readline() else r"\s+"
        mx, med, n = {}, {}, 0
        parts = {c: [] for c in ncols}
        for ch in pd.read_csv(f, sep=sep, usecols=ncols, chunksize=1_000_000):
            n += len(ch)
            for c in ncols:
                parts[c].append(pd.to_numeric(ch[c], errors="coerce").describe(percentiles=[0.5])[["max", "50%"]])
        for c in ncols:
            d = pd.concat(parts[c], axis=1)
            say("[N] %-3s %-8s rows %d | max %.0f | median of chunk medians %.0f" % (anc, c, n, d.loc["max"].max(), d.loc["50%"].median()))
        g = MAGMA_OUT / ("KOLLER_%s_35_10.genes.out" % anc)
        if g.exists():
            G = pd.read_csv(g, sep=r"\s+")
            say("[N] %-3s MAGMA genes %d | per-gene N median %.0f, max %.0f | genes with P < 2.6e-6: %d" % (
                anc, len(G), G["N"].median(), G["N"].max(), int((G["P"] < 2.6e-6).sum())))
    say("")
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V069_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
