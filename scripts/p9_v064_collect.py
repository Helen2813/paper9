#!/usr/bin/env python3
# Paper 9 | v0.0.64 | final figure set and Overleaf package of manuscript v2 with the control panel (no analyses; v0.2.54)
#   Collect: every file of the newest results/paper_figures/final_v061_* folder (F1-F9 with relabelled F2 / F9, S1-S5)
#   plus S6 from the newest results/external_atlas/p9_v063_* run -> results/paper_figures/final_v064_<TS>.
#   Package: manuscript/overleaf/main.tex (v2 final: lesion localization, known-source panel) + references.bib +
#   included figures -> ~/Desktop/paper9_overleaf_v2_final_<TS>.zip (outside the repository).
# Runs in seconds; no network.

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
FIGROOT = RES / "paper_figures"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
FINAL = FIGROOT / ("final_v064_" + TS)
BUILD = RES / "overleaf" / ("overleaf_v2_final_build_" + TS)
ZIP = HOME / "Desktop" / ("paper9_overleaf_v2_final_" + TS + ".zip")
LOG = REPO / "logs" / ("p9_v064_collect_" + TS + ".log")
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
    src = latest(FIGROOT, "final_v061_*/figure_legends.md")
    if src is None:
        say("[FATAL] no results/paper_figures/final_v061_* folder (run v061 first)")
        return False
    n = 0
    for f in sorted(src.parent.glob("*")):
        if f.is_file() and f.suffix in (".png", ".pdf", ".csv"):
            shutil.copy2(f, FINAL / f.name)
            n += 1
    say("[COPY] %d files from %s" % (n, src.parent.name))
    run = latest(RES / "external_atlas", "p9_v063_*/S6_control_panel.pdf")
    s6_leg = ""
    if run is None:
        say("[MISS] S6 not found (run v063 first)")
    else:
        for f in sorted(run.parent.glob("S6*")):
            if f.suffix in (".png", ".pdf", ".csv"):
                shutil.copy2(f, FINAL / f.name)
        lf = run.parent / "S6_legend.md"
        s6_leg = lf.read_text() if lf.exists() else ""
        say("[COPY] S6 from %s" % run.parent.name)
    legends = src.read_text().split("\n", 2)[-1]
    (FINAL / "figure_legends.md").write_text("# Paper 9 figure legends (v064 collection)\n\n" + legends + "\n" + s6_leg)
    (FINAL / "MANIFEST.txt").write_text("F1-F9, S1-S5: %s\nS6: %s\n" % (src.parent.name, run.parent.name if run else "-"))
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
    say("PAPER 9 | v0.0.64 | final figure set + Overleaf package v2 final | " + TS)
    say("=" * 72)
    if not collect():
        _fh.close()
        return
    missing = package()
    say("")
    say("[FIGURE SET] " + str(FINAL))
    say("[OUTPUT] " + str(ZIP))
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V064_COMPLETE" if not missing else "[INCOMPLETE] missing figures: " + ", ".join(missing))
    _fh.close()


if __name__ == "__main__":
    main()
