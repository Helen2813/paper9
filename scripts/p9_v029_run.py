#!/usr/bin/env python3
# Paper 9 | v0.0.29 | runner (standard library only): completes v028 (Ochoa-faithful scDRS) if needed, then runs B1.

import subprocess
import sys
import time
from pathlib import Path

REPO = Path.home() / "Desktop" / "paper9"
LOGS = REPO / "logs"
MAIN_PY = REPO / ".venv" / "bin" / "python"


def state(prefix, marker):
    logs = sorted(LOGS.glob(prefix + "*.log"), key=lambda p: p.stat().st_mtime)
    if any(marker in p.read_text(errors="ignore") for p in logs):
        return "done"
    if logs and time.time() - logs[-1].stat().st_mtime < 2 * 3600:
        return "running?"
    return "missing"


def main():
    st = state("p9_v028_scdrs_core_", "PAPER9_V028_COMPLETE")
    if st == "running?":
        print("[RUNNER] v028 scDRS has a recent unfinished log: it may still be running; start this again after it ends.", flush=True)
        sys.exit(3)
    if st == "missing":
        print("[RUNNER] running v028 first (Ochoa-faithful scDRS, ~20-40 min) ...", flush=True)
        rc = subprocess.run([str(MAIN_PY), str(REPO / "scripts" / "p9_v028_run.py")], cwd=str(REPO)).returncode
        print("[RUNNER] v028 exit %d" % rc, flush=True)
    else:
        print("[RUNNER] v028 already complete", flush=True)
    print("[RUNNER] running B1 (fallopian tube; ~1.6 GB download, ~30-50 min) ...", flush=True)
    rc = subprocess.run([str(MAIN_PY), str(REPO / "scripts" / "p9_v029_b1_ft.py")], cwd=str(REPO)).returncode
    print("[RUNNER] B1 exit %d" % rc, flush=True)
    print("[DONE] PAPER9_V029_RUNNER_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
