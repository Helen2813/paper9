#!/usr/bin/env python3
# Paper 9 | v0.0.30 | runner (standard library only): B1 (v029) first if not complete, then the Koller 2026
#   method-matched reproduction (long, resumable: re-running continues the LD scores where they stopped).

import subprocess
import sys
from pathlib import Path

REPO = Path.home() / "Desktop" / "paper9"
LOGS = REPO / "logs"
MAIN_PY = REPO / ".venv" / "bin" / "python"


def state(prefix, marker, script_name):
    logs = sorted(LOGS.glob(prefix + "*.log"), key=lambda p: p.stat().st_mtime)
    if any(marker in p.read_text(errors="ignore") for p in logs):
        return "done"
    try:
        if subprocess.run(["pgrep", "-f", script_name], capture_output=True).returncode == 0:
            return "running"
    except FileNotFoundError:
        pass
    return "missing"


def main():
    st = state("p9_v029_b1_ft_", "PAPER9_V029_COMPLETE", "p9_v029_b1_ft.py")
    if st == "running":
        print("[RUNNER] B1 (v029) is still running in another window; start this again after it ends.", flush=True)
        sys.exit(3)
    if st == "missing":
        print("[RUNNER] running B1 first (fallopian tube; ~1.6 GB download, ~30-60 min) ...", flush=True)
        rc = subprocess.run([str(MAIN_PY), str(REPO / "scripts" / "p9_v029_b1_ft.py")], cwd=str(REPO)).returncode
        print("[RUNNER] B1 exit %d" % rc, flush=True)
    else:
        print("[RUNNER] B1 already complete", flush=True)
    print("[RUNNER] running the Koller 2026 reproduction (markers ~20-40 min, LD scores ~4-8 h, S-LDSC ~1 h) ...", flush=True)
    rc = subprocess.run([str(MAIN_PY), str(REPO / "scripts" / "p9_v030_koller_repro.py")], cwd=str(REPO)).returncode
    print("[RUNNER] Koller reproduction exit %d" % rc, flush=True)
    print("[DONE] PAPER9_V030_RUNNER_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
