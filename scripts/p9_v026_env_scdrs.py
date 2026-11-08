#!/usr/bin/env python3
# Paper 9 | v0.0.26 | technical (no rule change): the scdrs install in v025 left an extension compiled for a different
#   numpy ABI. (1) probe the main .venv and repair only modules that now fail to import; (2) build a separate
#   .venv_scdrs with a pinned, mutually compatible set (numpy 1.26); (3) re-run the v025 analysis unchanged in it.
# Standard library only (the main venv may be partly broken). ~10 min setup + v025 runtime (~60-90 min).

import datetime as dt
import json
import subprocess
import sys
from pathlib import Path

REPO = Path.home() / "Desktop" / "paper9"
VENV_PY = REPO / ".venv" / "bin" / "python"
SC_VENV = REPO / ".venv_scdrs"
SC_PY = SC_VENV / "bin" / "python"
SCRIPT_V025 = REPO / "scripts" / "p9_v025_scdrs_ochoa.py"
TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
LOG = REPO / "logs" / ("p9_v026_env_scdrs_" + TS + ".log")
LOG.parent.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")

MAIN_MODULES = ["numpy", "scipy", "pandas", "sklearn", "h5py", "anndata", "scanpy", "statsmodels", "numba", "leidenalg", "skmisc", "scdrs"]
PIP_NAME = {"sklearn": "scikit-learn"}
# Pinned set for Python 3.9 / macOS arm64 (compiled packages from wheels; fire is pure Python, sdist only)
# Cause of the v025 failure: scikit-misc (scdrs dependency) wheel for cp39 is built against numpy 1.x,
# while the main venv runs numpy 2.x -> "numpy.dtype size changed" on import.
SC_PINS = ["numpy==1.26.4", "scipy==1.11.4", "pandas==2.1.4", "scikit-learn==1.3.2", "statsmodels==0.14.1",
           "numba==0.58.1", "llvmlite==0.41.1", "h5py==3.11.0", "anndata==0.10.9", "scanpy==1.10.3",
           "igraph==0.11.6", "leidenalg==0.10.2", "scikit-misc==0.3.1", "fire==0.7.0", "scdrs==1.0.2"]
PROBE = r"""
import importlib, json, sys
out = {}
for m in sys.argv[1:]:
    try:
        mod = importlib.import_module(m)
        out[m] = "OK " + str(getattr(mod, "__version__", "?"))
    except Exception as e:
        out[m] = "FAIL %s: %s" % (type(e).__name__, str(e)[:150])
print("PROBE_JSON " + json.dumps(out))
"""


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def run(cmd, timeout=3600):
    try:
        p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"


def probe(py, mods):
    rc, o, e = run([py, "-c", PROBE] + mods, timeout=900)
    for line in o.splitlines():
        if line.startswith("PROBE_JSON "):
            return json.loads(line[len("PROBE_JSON "):])
    return {m: "FAIL probe rc=%d %s" % (rc, (o + e).strip()[-150:]) for m in mods}


def show(tag, res):
    for m, v in res.items():
        say("    [%s] %-12s %s" % (tag, m, v))


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.26 | environment repair + scDRS venv + v025 re-run | " + TS)
    say("=" * 72)
    rc, o, e = run([VENV_PY, "-m", "pip", "freeze"], timeout=300)
    (REPO / "logs" / ("pip_freeze_main_before_v026_" + TS + ".txt")).write_text(o)
    res = probe(VENV_PY, MAIN_MODULES)
    say("[MAIN_VENV] before repair")
    show("MAIN", res)
    core = [m for m in MAIN_MODULES if m not in ("scdrs", "statsmodels", "skmisc") and res.get(m, "").startswith("FAIL")]
    for m in core:
        pkg = PIP_NAME.get(m, m)
        rc, o, e = run([VENV_PY, "-m", "pip", "install", "-q", "--upgrade", "--prefer-binary", pkg], timeout=1800)
        say("[MAIN_FIX] pip upgrade %s rc=%d %s" % (pkg, rc, (e or o).strip()[-200:]))
    if core:
        res = probe(VENV_PY, MAIN_MODULES)
        say("[MAIN_VENV] after repair")
        show("MAIN", res)
    bad = [m for m in MAIN_MODULES if m not in ("scdrs", "statsmodels", "skmisc") and res.get(m, "").startswith("FAIL")]
    say("[MAIN_VENV_STATUS] %s" % ("OK (pipeline modules import)" if not bad else "STILL FAILING: " + ", ".join(bad)))

    # Dedicated scDRS environment
    if not SC_PY.exists():
        rc, o, e = run(["/usr/bin/python3", "-m", "venv", SC_VENV], timeout=600)
        say("[SC_VENV] create rc=%d %s" % (rc, (e or o).strip()[-200:]))
        run([SC_PY, "-m", "pip", "install", "-q", "--upgrade", "pip"], timeout=900)
    sres = probe(SC_PY, ["numpy", "scanpy", "anndata", "scdrs", "leidenalg"])
    if any(v.startswith("FAIL") for v in sres.values()):
        say("[SC_VENV] installing pinned set (prefer wheels, ~300 MB, a few minutes)")
        rc, o, e = run([SC_PY, "-m", "pip", "install", "-q", "--prefer-binary"] + SC_PINS, timeout=3600)
        say("[SC_VENV] pip rc=%d %s" % (rc, (e or o).strip()[-400:]))
        sres = probe(SC_PY, ["numpy", "scipy", "pandas", "sklearn", "statsmodels", "anndata", "scanpy", "scdrs", "leidenalg"])
    say("[SC_VENV] modules")
    show("SC", sres)
    if any(v.startswith("FAIL") for v in sres.values()):
        say("[STOP] scDRS environment not usable; send this log")
    else:
        say("")
        say("[RUN] v025 analysis (unchanged rules, CHANGELOG v0.2.16) in .venv_scdrs ...")
        _fh.flush()
        rc = subprocess.run([str(SC_PY), str(SCRIPT_V025)], cwd=str(REPO)).returncode
        say("[RUN] v025 exit %d (its own log: logs/p9_v025_scdrs_*.log, newest)" % rc)
    say("[LOG] " + str(LOG))
    say("[DONE] PAPER9_V026_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
