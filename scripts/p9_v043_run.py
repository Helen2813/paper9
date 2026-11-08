#!/usr/bin/env python3
# Paper 9 | v0.0.43 | runner (standard library only): locus atlas v042 first if not complete, then v043.

import subprocess
from pathlib import Path

REPO = Path.home() / "Desktop" / "paper9"
PY = REPO / ".venv" / "bin" / "python"


def complete(prefix, marker):
    return any(marker in p.read_text(errors="ignore") for p in (REPO / "logs").glob(prefix + "*.log"))


def main():
    if not complete("p9_v042_locus_atlas_", "PAPER9_V042_COMPLETE"):
        print("[RUNNER] running the locus atlas v042 first (~20 min) ...", flush=True)
        rc = subprocess.run([str(PY), str(REPO / "scripts" / "p9_v042_locus_atlas.py")], cwd=str(REPO)).returncode
        print("[RUNNER] v042 exit %d" % rc, flush=True)
    else:
        print("[RUNNER] v042 already complete", flush=True)
    print("[RUNNER] running v043 (~5-10 min) ...", flush=True)
    rc = subprocess.run([str(PY), str(REPO / "scripts" / "p9_v043_real_families.py")], cwd=str(REPO)).returncode
    print("[RUNNER] v043 exit %d" % rc, flush=True)
    print("[DONE] PAPER9_V043_RUNNER_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
