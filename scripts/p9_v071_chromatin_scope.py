#!/usr/bin/env python3
# Paper 9 | v0.0.71 | scope of a chromatin-based (marker-free) test of the stromal signal (no analyses, no large downloads)
#   1. Local S-LDSC resources (plink, frq, HapMap3 list, baseline, weights, munged GWAS, LDSC-SEG LD scores if present)
#   2. GEO supplementary files and sample titles of endometrial stromal ATAC-seq series (GSE104720, GSE246585, GSE205473)
#   3. Remote availability of the LDSC-SEG multi-tissue chromatin LD scores and a hg38->hg19 chain file (size only)
#   4. Python packages for coordinate lift-over
# Runs in 1-3 minutes.

import datetime as dt
import gzip
import importlib.util
import re
import subprocess
from pathlib import Path

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RAW = REPO / "data" / "raw"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
REF = RAW / "ldsc_ref" / "unpacked"
SCOPE = RAW / "chromatin_scope"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
LOG = REPO / "logs" / ("p9_v071_chromatin_scope_" + TS + ".log")
SERIES = ["GSE104720", "GSE246585", "GSE205473"]
REMOTE = {
    "LDSC-SEG chromatin LD scores": "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/LDSC_SEG_ldscores/Multi_tissue_chromatin_1000Gv3_ldscores.tgz",
    "LDSC-SEG chromatin LD scores (alt)": "https://storage.googleapis.com/broad-alkesgroup-public-requester-pays/LDSCORE/LDSC_SEG_ldscores/Multi_tissue_chromatin_1000Gv3_ldscores.tgz",
    "UCSC hg38ToHg19 chain": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
}
LOG.parent.mkdir(parents=True, exist_ok=True)
SCOPE.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")


def say(m=""):
    print(m, flush=True)
    _fh.write(m + "\n")
    _fh.flush()


def sh(cmd, timeout=120):
    try:
        p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"


def geo_base(acc):
    return "https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/" % (acc[:-3], acc)


def listing(url):
    rc, o, e = sh(["curl", "-sSL", "--fail", "--max-time", "60", "-A", "paper9", url])
    if rc != 0:
        return None
    rows = []
    for m in re.finditer(r'<a href="([^"?/][^"]*)">[^<]*</a>\s*([0-9]{4}-[0-9]{2}-[0-9]{2} [0-9:]+)?\s*([0-9.]+[KMG]?)?', o):
        rows.append((m.group(1), m.group(3) or ""))
    return rows


def head_size(url):
    rc, o, e = sh(["curl", "-sSIL", "--max-time", "60", "-A", "paper9", url])
    if rc != 0:
        return "not reachable"
    codes = re.findall(r"HTTP/[0-9.]+ (\d+)", o)
    size = re.findall(r"(?i)content-length:\s*(\d+)", o)
    return "HTTP %s, %s" % (codes[-1] if codes else "?", ("%.1f GB" % (int(size[-1]) / 1e9)) if size else "size unknown")


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.71 | chromatin-test scope (local LDSC resources, stromal ATAC series) | " + TS)
    say("=" * 72)
    say("[LOCAL] LDSC python %s | ldsc.py %s" % ((TOOLS / "ldsc_venv" / "bin" / "python").exists(), (TOOLS / "ldsc_py3" / "ldsc.py").exists()))
    if REF.exists():
        for d in sorted(p for p in REF.iterdir() if p.is_dir()):
            n = sum(1 for _ in d.rglob("*") if _.is_file())
            say("  [REF] %-50s %5d files" % (d.name[:50], n))
    for name in ("hm3_no_MHC.list.txt", "w_hm3.snplist"):
        f = RAW / "ldsc_ref" / name
        say("  [REF] %-30s %s" % (name, f.exists()))
    seg = [p for p in RAW.rglob("*") if p.is_dir() and re.search(r"Multi_tissue|LDSC_SEG|Roadmap|ENTEx|chromatin", p.name, re.I)]
    say("  [REF] LDSC-SEG / chromatin folders: %s" % (", ".join(str(p.relative_to(REPO)) for p in seg[:10]) or "none"))
    mg = sorted((PROC / "munged").glob("*.sumstats.gz"))
    say("  [MUNGED] %s" % (", ".join(p.name.replace(".sumstats.gz", "") for p in mg) or "none"))
    say("")
    for acc in SERIES:
        base = geo_base(acc)
        sup = listing(base + "suppl/")
        say("[GEO] %s supplementary: %s" % (acc, "listing failed" if sup is None else (", ".join("%s (%s)" % (a, b) for a, b in sup if not a.startswith("http")) or "none")))
        mats = listing(base + "matrix/") or []
        for m, _ in mats:
            if not m.endswith("_series_matrix.txt.gz"):
                continue
            f = SCOPE / m
            if not f.exists():
                sh(["curl", "-sSL", "--fail", "--max-time", "120", "-A", "paper9", "-o", f, base + "matrix/" + m], timeout=180)
            if not f.exists():
                continue
            txt = gzip.open(f, "rt", errors="ignore").read().split("!series_matrix_table_begin")[0]
            keep = {}
            for line in txt.splitlines():
                k = line.split("\t")[0]
                if k in ("!Series_title", "!Series_platform_id", "!Sample_title", "!Sample_library_strategy", "!Sample_source_name_ch1",
                         "!Sample_characteristics_ch1", "!Sample_supplementary_file_1"):
                    keep.setdefault(k, []).append([v.strip().strip('"') for v in line.split("\t")[1:]])
            say("  [%s] %s | platform %s" % (acc, " ".join(keep.get("!Series_title", [[""]])[0])[:120], " ".join(keep.get("!Series_platform_id", [[""]])[0])))
            titles = keep.get("!Sample_title", [[]])[0]
            strat = keep.get("!Sample_library_strategy", [[]])[0]
            srcs = keep.get("!Sample_source_name_ch1", [[]])[0]
            chars = keep.get("!Sample_characteristics_ch1", [])
            sfile = keep.get("!Sample_supplementary_file_1", [[]])[0]
            for i, t in enumerate(titles):
                ch = " | ".join(c[i] for c in chars if i < len(c))[:110]
                say("    %-34s %-10s %-26s %s | %s" % (t[:34], strat[i] if i < len(strat) else "", (srcs[i] if i < len(srcs) else "")[:26], ch,
                                                      (sfile[i].split("/")[-1] if i < len(sfile) else "")[:60]))
        say("")
    for k, u in REMOTE.items():
        say("[REMOTE] %-36s %s" % (k, head_size(u)))
    for pkg in ("pyliftover", "liftover", "pybedtools"):
        say("[PY] %-12s %s" % (pkg, "installed" if importlib.util.find_spec(pkg) else "not installed"))
    rc, o, e = sh(["bash", "-lc", "command -v bedtools liftOver"])
    say("[BIN] %s" % (o.strip().replace("\n", ", ") or "bedtools / liftOver not on PATH"))
    say("")
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V071_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
