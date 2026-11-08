#!/usr/bin/env python3
# Paper 9 | v0.0.28 | runner (standard library only): completes v027 if needed, then runs the Ochoa-faithful
#   scDRS core (scripts/p9_v028_scdrs_core.py) in .venv_scdrs. Rules: CHANGELOG v0.2.17 / v0.2.18.

import subprocess
import sys
import time
from pathlib import Path

REPO = Path.home() / "Desktop" / "paper9"
LOGS = REPO / "logs"
MAIN_PY = REPO / ".venv" / "bin" / "python"
SC_PY = REPO / ".venv_scdrs" / "bin" / "python"


def done(prefix, marker):
    logs = sorted(LOGS.glob(prefix + "*.log"), key=lambda p: p.stat().st_mtime)
    if any(marker in p.read_text(errors="ignore") for p in logs):
        return "done"
    if logs and time.time() - logs[-1].stat().st_mtime < 2 * 3600:
        return "running?"
    return "missing"


def main():
    st = done("p9_v027_resolution_", "PAPER9_V027_COMPLETE")
    if st == "running?":
        print("[RUNNER] v027 has a recent unfinished log: it may still be running; start this again after it ends.", flush=True)
        sys.exit(3)
    if st == "missing":
        print("[RUNNER] running v027 first (~1.5-2 h) ...", flush=True)
        rc = subprocess.run([str(MAIN_PY), str(REPO / "scripts" / "p9_v027_resolution.py")], cwd=str(REPO)).returncode
        print("[RUNNER] v027 exit %d" % rc, flush=True)
    else:
        print("[RUNNER] v027 already complete", flush=True)
    if not SC_PY.exists():
        print("[RUNNER] .venv_scdrs missing: run paper9_v026.py first", flush=True)
        sys.exit(4)
    print("[RUNNER] running Ochoa-faithful scDRS core in .venv_scdrs (~20-40 min) ...", flush=True)
    rc = subprocess.run([str(SC_PY), str(REPO / "scripts" / "p9_v028_scdrs_core.py")], cwd=str(REPO)).returncode
    print("[RUNNER] v028 core exit %d" % rc, flush=True)
    print("[DONE] PAPER9_V028_RUNNER_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
