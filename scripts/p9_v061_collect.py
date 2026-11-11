#!/usr/bin/env python3
# Paper 9 | v0.0.61 | final figure set v2 and Overleaf package of manuscript v2 (no analyses; CHANGELOG v0.2.51)
#   Collect: every file of the newest results/paper_figures/final_v056_* folder, F2 and F9 replaced by the relabelled
#   versions (x labels no longer clipped; F2c title without "true source"; same data), plus S5 from the newest
#   results/review/p9_v060_* run -> results/paper_figures/final_v061_<TS> with combined legends.
#   Package: manuscript/overleaf/main.tex (v2: lesion localization, somatic-mutation note) + references.bib + included
#   figures -> ~/Desktop/paper9_overleaf_v2_<TS>.zip (outside the repository). Compiles once if latexmk exists.
# Runs in seconds; no network.

import base64
import datetime as dt
import json
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RES = REPO / "results"
SRC = REPO / "manuscript" / "overleaf"
ASSETS = REPO / "manuscript" / "assets" / "figures_fixed_v061.json"
FIGROOT = RES / "paper_figures"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
FINAL = FIGROOT / ("final_v061_" + TS)
BUILD = RES / "overleaf" / ("overleaf_v2_build_" + TS)
ZIP = HOME / "Desktop" / ("paper9_overleaf_v2_" + TS + ".zip")
LOG = REPO / "logs" / ("p9_v061_collect_" + TS + ".log")
for d in (FINAL, BUILD / "figures", LOG.parent):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")


def say(m=""):
    print(m, flush=True)
    _fh.write(m + "\n")
    _fh.flush()


def latest(base, pattern):
    fs = sorted(base.glob(pattern))
    return fs[-1] if fs else None


def collect():
    src = latest(FIGROOT, "final_v056_*/figure_legends.md")
    if src is None:
        say("[FATAL] no results/paper_figures/final_v056_* folder")
        return False
    n = 0
    for f in sorted(src.parent.glob("*")):
        if f.is_file() and f.suffix in (".png", ".pdf", ".csv"):
            shutil.copy2(f, FINAL / f.name)
            n += 1
    say("[COPY] %d files from %s" % (n, src.parent.name))
    fixed = json.loads(ASSETS.read_text()) if ASSETS.exists() else {}
    for name, b64 in fixed.items():
        (FINAL / name).write_bytes(base64.b64decode(b64))
    say("[FIXED] relabelled F2 / F9 written: %s" % (", ".join(sorted(fixed)) or "none (assets file missing)"))
    run = latest(RES / "review", "p9_v060_*/S5_lesion_locus_genes.pdf")
    s5_leg = ""
    if run is None:
        say("[MISS] S5 not found (run v060 first)")
    else:
        for f in sorted(run.parent.glob("S5*")):
            if f.suffix in (".png", ".pdf", ".csv"):
                shutil.copy2(f, FINAL / f.name)
        lf = run.parent / "S5_legend.md"
        s5_leg = lf.read_text() if lf.exists() else ""
        say("[COPY] S5 from %s" % run.parent.name)
    legends = src.read_text().split("\n", 2)[-1]
    note = ("\n## Note (v061)\n\nF2 and F9 were redrawn from the same data with shorter axis labels; the F2c title no longer "
            "uses the words \"true source\".\n")
    (FINAL / "figure_legends.md").write_text("# Paper 9 figure legends (v061 collection)\n\n" + legends + "\n" + s5_leg + note)
    (FINAL / "MANIFEST.txt").write_text("F1-F9, S1-S4: %s (F2, F9 relabelled in v061)\nS5: %s\n" % (
        src.parent.name, run.parent.name if run else "-"))
    return True


def package():
    tex, bib = SRC / "main.tex", SRC / "references.bib"
    if not tex.exists() or not bib.exists():
        say("[FATAL] manuscript/overleaf/main.tex or references.bib missing")
        return []
    shutil.copy2(tex, BUILD / "main.tex")
    shutil.copy2(bib, BUILD / "references.bib")
    text = tex.read_text(encoding="utf-8")
    wanted = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{figures/([^}]+)\}", text)
    missing = []
    for name in wanted:
        f = FINAL / name
        if f.exists():
            shutil.copy2(f, BUILD / "figures" / name)
        else:
            missing.append(name)
    say("[FIGURES] %d of %d included figures copied%s" % (len(wanted) - len(missing), len(wanted),
                                                          ("; MISSING: " + ", ".join(missing)) if missing else ""))
    (BUILD / "README_OVERLEAF.txt").write_text(
        "Overleaf: New Project -> Upload Project -> choose the zip.\n"
        "Compiler: pdfLaTeX (Menu -> Settings). Main document: main.tex. Bibliography: BibTeX (unsrtnat).\n", encoding="utf-8")
    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(BUILD.rglob("*")):
            if f.is_file():
                z.write(f, f.relative_to(BUILD))
    say("[ZIP] " + str(ZIP) + " (%.1f MB)" % (ZIP.stat().st_size / 1e6))
    if shutil.which("latexmk"):
        rc = subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "main.tex"], cwd=str(BUILD),
                            capture_output=True, text=True).returncode
        say("[COMPILE] latexmk rc=%d" % rc)
    else:
        say("[COMPILE] no local LaTeX; compile in Overleaf")
    return missing


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.61 | final figure set v2 + Overleaf package v2 | " + TS)
    say("=" * 72)
    if not collect():
        _fh.close()
        return
    missing = package()
    say("")
    say("[FIGURE SET] " + str(FINAL))
    say("[OUTPUT] " + str(ZIP))
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V061_COMPLETE" if not missing else "[INCOMPLETE] missing figures: " + ", ".join(missing))
    _fh.close()


if __name__ == "__main__":
    main()
