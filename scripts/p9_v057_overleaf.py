#!/usr/bin/env python3
# Paper 9 | v0.0.57 | Overleaf package of the manuscript (no analyses)
#   Takes manuscript/overleaf/main.tex and references.bib (written by this installer), copies every figure that
#   main.tex includes from the newest results/paper_figures/final_v056_* folder into figures/, and zips the folder
#   into results/overleaf/paper9_overleaf_<TS>.zip (git-ignored) for Overleaf (New Project -> Upload Project).
#   If latexmk is installed locally, it also compiles once as a check. Runs in seconds; no network.

import datetime as dt
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
SRC = REPO / "manuscript" / "overleaf"
FIGROOT = REPO / "results" / "paper_figures"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
BUILD = REPO / "results" / "overleaf" / ("overleaf_build_" + TS)
ZIP = REPO / "results" / "overleaf" / ("paper9_overleaf_" + TS + ".zip")
LOG = REPO / "logs" / ("p9_v057_overleaf_" + TS + ".log")
LOG.parent.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")


def say(m=""):
    print(m, flush=True)
    _fh.write(m + "\n")
    _fh.flush()


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.57 | Overleaf package | " + TS)
    say("=" * 72)
    tex, bib = SRC / "main.tex", SRC / "references.bib"
    if not tex.exists() or not bib.exists():
        say("[FATAL] manuscript/overleaf/main.tex or references.bib missing")
        return
    finals = sorted(FIGROOT.glob("final_v056_*"))
    if not finals:
        say("[FATAL] no results/paper_figures/final_v056_* folder (run v056 first)")
        return
    fig_dir = finals[-1]
    say("[FIGURES] source " + str(fig_dir.relative_to(REPO)))
    (BUILD / "figures").mkdir(parents=True, exist_ok=True)
    shutil.copy2(tex, BUILD / "main.tex")
    shutil.copy2(bib, BUILD / "references.bib")
    text = tex.read_text(encoding="utf-8")
    wanted = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{figures/([^}]+)\}", text)
    missing = []
    for name in wanted:
        src = fig_dir / name
        if src.exists():
            shutil.copy2(src, BUILD / "figures" / name)
        else:
            missing.append(name)
    say("[FIGURES] %d of %d included figures copied%s" % (len(wanted) - len(missing), len(wanted),
                                                          ("; MISSING: " + ", ".join(missing)) if missing else ""))
    legends = fig_dir / "figure_legends.md"
    if legends.exists():
        shutil.copy2(legends, BUILD / "figure_legends_auto.md")
    (BUILD / "README_OVERLEAF.txt").write_text(
        "Overleaf: New Project -> Upload Project -> choose the zip.\n"
        "Compiler: pdfLaTeX (Menu -> Settings). Main document: main.tex. Bibliography: BibTeX (unsrtnat).\n"
        "Red [TODO: ...] markers in main.tex need a decision before submission.\n", encoding="utf-8")
    todos = sorted(set(re.findall(r"\\todo\{([^}]*)\}", text)))
    say("[TODO] %d open items: %s" % (len(todos), "; ".join(todos)))
    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(BUILD.rglob("*")):
            if f.is_file():
                z.write(f, f.relative_to(BUILD))
    say("[ZIP] " + str(ZIP.relative_to(REPO)) + " (%.1f MB)" % (ZIP.stat().st_size / 1e6))
    if shutil.which("latexmk"):
        rc = subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "main.tex"], cwd=str(BUILD),
                            capture_output=True, text=True).returncode
        say("[COMPILE] latexmk rc=%d%s" % (rc, " -> results/overleaf/" + BUILD.name + "/main.pdf" if rc == 0 else " (see main.log in the build folder)"))
    else:
        say("[COMPILE] no local LaTeX; compile in Overleaf")
    say("")
    say("[OUTPUT] " + str(ZIP))
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V057_COMPLETE" if not missing else "[INCOMPLETE] missing figures: " + ", ".join(missing))
    _fh.close()


if __name__ == "__main__":
    main()
