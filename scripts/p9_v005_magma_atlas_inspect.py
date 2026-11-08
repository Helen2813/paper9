#!/usr/bin/env python3
# Paper 9 | v0.0.5 | MAGMA install (surfsara links) + MAGMA gene analysis of GCST90205183 (sanity check)
#                   + structure of the Paper 8 GSE266265 export + HECA download and metadata inspection
# No disease effects are computed. Standard library only (h5ad inspection runs in .venv).

import datetime as dt
import gzip
import json
import re
import shutil
import subprocess
import time
import zipfile
from pathlib import Path

# ----------------------------------------------------------------------------
# Hardcoded locations
# ----------------------------------------------------------------------------
HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RAW = REPO / "data" / "raw"
INTERIM = REPO / "data" / "interim"
EXT = REPO / "data" / "external"
TOOLS = REPO / "tools"
VENV_PY = REPO / ".venv" / "bin" / "python"
MAGMA_DIR = TOOLS / "magma"
GWAS_H = RAW / "gwas" / "GCST90205183" / "harmonised" / "GCST90205183.h.tsv.gz"
GWAS_N = 474160
P8_EXPORT = EXT / "p8_GSE266265_v232c_export"
HECA_DIR = RAW / "heca"
HECA_URL = "https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/endometriumAtlasV2_cells_with_counts.h5ad"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "magma" / ("p9_v005_magma_atlas_inspect_" + TS)
LOG = REPO / "logs" / ("p9_v005_magma_atlas_inspect_" + TS + ".log")
MAGMA_OUT = INTERIM / "magma"

SURF = "https://vu.data.surfsara.nl/index.php/s/"
MAGMA_SHARES = {
    "mac": "1M1d9vHtVidEwvU",
    "source": "1OOi7bxLWef0GwY",
    "genes_b37": "Pj2orwuF2JYyKxq",
    "genes_b38": "yj952iHqy5anYhH",
    "ref_eur": "VZNByNwpD8qqINe",
}
KNOWN_ENDO_GENES = ["WNT4", "CDC42", "GREB1", "ESR1", "FN1", "VEZT", "ID4", "KDR", "IL1A", "SYNE1",
                    "CDKN2B", "FSHB", "ETAA1", "RSPO3", "SKAP1", "BSN", "MLLT10", "HOXC10"]

for d in (OUT, LOG.parent, MAGMA_DIR, MAGMA_OUT, HECA_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V005_MAGMA_ATLAS_INSPECT", "timestamp_utc": TS, "gates": {}}


def say(msg=""):
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def human(n):
    if n is None:
        return "NA"
    n = float(n)
    for u in ["B", "KB", "MB", "GB", "TB"]:
        if n < 1024:
            return "%.1f %s" % (n, u)
        n /= 1024.0
    return "%.1f PB" % n


def sh(cmd, cwd=None, timeout=None):
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout, p.stderr
    except FileNotFoundError as e:
        return 127, "", str(e)
    except OSError as e:
        return 126, "", str(e)
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"


def remote_size(url):
    rc, out, err = sh(["curl", "-sSIL", "--max-time", "60", "-A", "paper9", url])
    blocks = [b for b in re.split(r"\r?\n\r?\n", out) if b.strip().startswith("HTTP/")]
    if not blocks:
        return None
    last = blocks[-1].strip()
    code = re.match(r"HTTP/\S+\s+(\d+)", last)
    if not code or not code.group(1).startswith("2"):
        return None
    m = re.search(r"^content-length:\s*(\d+)", last, flags=re.M | re.I)
    size = int(m.group(1)) if m else None
    return size if size else None


def download(url, dst, label):
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    expect = remote_size(url)
    if dst.exists() and dst.stat().st_size > 0 and (expect is None or dst.stat().st_size == expect):
        say("[DL] %-48s %10s SKIP_PRESENT" % (label[:48], human(dst.stat().st_size)))
        return True
    part = dst.with_name(dst.name + ".part")
    base = ["curl", "-L", "-sS", "--fail", "--retry", "4", "--retry-delay", "5",
            "--connect-timeout", "30", "-A", "paper9", "-o", str(part), url]
    t0 = time.time()
    say("[DL] %-48s %10s START" % (label[:48], human(expect)))
    rc, out, err = sh(base[:1] + ["-C", "-"] + base[1:])
    if rc == 33:
        if part.exists():
            part.unlink()
        rc, out, err = sh(base)
    if rc != 0 or not part.exists() or part.stat().st_size == 0:
        say("[DL] %-48s FAILED rc=%s %s" % (label[:48], rc, err.strip()[:140]))
        return False
    if expect and part.stat().st_size != expect:
        say("[DL] %-48s SIZE_MISMATCH (re-run to resume)" % label[:48])
        return False
    part.rename(dst)
    say("[DL] %-48s %10s OK %.0fs" % (label[:48], human(dst.stat().st_size), time.time() - t0))
    return True


# ----------------------------------------------------------------------------
# A. MAGMA
# ----------------------------------------------------------------------------
def magma_version(binpath):
    txt = ""
    for args in (["--version"], []):
        rc, out, err = sh([str(binpath)] + args, timeout=60)
        txt = out + err
        m = re.search(r"v1\.\d+\w*", txt)
        if m:
            return m.group(0), txt
    return None, txt


def setup_magma():
    say("")
    say("-" * 72)
    say("[SECTION] MAGMA (surfsara share links + /download)")
    for f in list(MAGMA_DIR.glob("*.zip")):
        if not zipfile.is_zipfile(f):
            f.unlink()
    for key, share in MAGMA_SHARES.items():
        z = MAGMA_DIR / (key + ".zip")
        if not (z.exists() and zipfile.is_zipfile(z)):
            download(SURF + share + "/download", z, "MAGMA " + key)
        if z.exists() and not zipfile.is_zipfile(z):
            head = z.read_bytes()[:200]
            say("[MAGMA_DL] %s not a zip: %r" % (key, head[:120]))
            z.unlink()
            continue
        sub = MAGMA_DIR / key
        if z.exists() and not sub.exists():
            with zipfile.ZipFile(z) as zf:
                zf.extractall(sub)
        if sub.exists():
            files = sorted(str(p.relative_to(sub)) for p in sub.rglob("*") if p.is_file())
            say("[MAGMA_FILES] %-10s %s" % (key, ", ".join(files[:8])))

    status, version, binpath = "NOT_INSTALLED", None, None
    for cand in sorted((MAGMA_DIR / "mac").rglob("magma")) if (MAGMA_DIR / "mac").exists() else []:
        if cand.is_file():
            cand.chmod(0o755)
            rc, out, err = sh(["file", str(cand)])
            say("[MAGMA_MAC_FILE] " + out.strip()[-120:])
            version, txt = magma_version(cand)
            if version:
                status, binpath = "MAC_BINARY_OK", cand
            else:
                say("[MAGMA_MAC_RUN] " + (txt.strip() or "not runnable (likely x86_64 without Rosetta)")[:200])
            break
    if not binpath and (MAGMA_DIR / "source").exists():
        built = next((p for p in (MAGMA_DIR / "source").rglob("magma") if p.is_file() and p.stat().st_size > 100000), None)
        if built is None:
            mk = next(iter((MAGMA_DIR / "source").rglob("Makefile")), None) or next(iter((MAGMA_DIR / "source").rglob("makefile")), None)
            if mk:
                text = mk.read_text(errors="ignore")
                patched = text
                for flag in ["-static-libstdc++", "-static-libgcc", "-static", "-fopenmp"]:
                    patched = re.sub(r"(?<=\s)" + re.escape(flag) + r"(?=\s|$)", "", patched, flags=re.M)
                if patched != text:
                    mk.write_text(patched)
                    say("[MAGMA_COMPILE] Makefile patched (static/OpenMP flags removed for macOS clang)")
                say("[MAGMA_COMPILE] make in " + str(mk.parent))
                rc, out, err = sh(["make"], cwd=str(mk.parent), timeout=3600)
                lines = (out + err).strip().splitlines()
                say("[MAGMA_COMPILE] rc=%d" % rc)
                if rc != 0:
                    for l in lines[-12:]:
                        say("       " + l[:200])
                built = mk.parent / "magma" if (mk.parent / "magma").exists() else None
            else:
                say("[MAGMA_COMPILE] no Makefile; source files: " +
                    ", ".join(sorted(p.name for p in (MAGMA_DIR / "source").rglob("*") if p.is_file())[:20]))
        if built:
            built.chmod(0o755)
            version, txt = magma_version(built)
            if version:
                status, binpath = "COMPILED_FROM_SOURCE_OK", built
            else:
                say("[MAGMA_BUILT_RUN] " + txt.strip()[:200])
    say("[MAGMA] %s %s %s" % (status, version or "", binpath or ""))
    SUM["gates"]["MAGMA"] = status
    return binpath


def magma_sanity(binpath):
    say("")
    say("-" * 72)
    say("[SECTION] MAGMA GENE ANALYSIS: GCST90205183 (window 35kb/10kb, 1000G EUR)")
    bfile = next((p.with_suffix("") for p in (MAGMA_DIR / "ref_eur").rglob("*.bed")), None)
    gloc = next((p for p in (MAGMA_DIR / "genes_b37").rglob("*.gene.loc")), None)
    if not bfile or not gloc:
        say("[MAGMA_SANITY] SKIPPED reference=%s geneloc=%s" % (bfile, gloc))
        return
    annot = MAGMA_OUT / "annot_b37_35_10"
    if not Path(str(annot) + ".genes.annot").exists():
        rc, out, err = sh([str(binpath), "--annotate", "window=35,10", "--snp-loc", str(bfile) + ".bim",
                           "--gene-loc", str(gloc), "--out", str(annot)], timeout=3600)
        say("[MAGMA_ANNOTATE] rc=%d" % rc)
    pfile = MAGMA_OUT / "GCST90205183.pval.txt"
    if not pfile.exists():
        n = 0
        with gzip.open(GWAS_H, "rt") as f, open(pfile, "w") as o:
            hdr = f.readline().rstrip("\n").split("\t")
            ix = {h: i for i, h in enumerate(hdr)}
            o.write("SNP\tP\n")
            for line in f:
                t = line.rstrip("\n").split("\t")
                rs, p = t[ix["rsid"]], t[ix["p_value"]]
                if rs.startswith("rs") and p not in ("", "NA"):
                    o.write(rs + "\t" + p + "\n")
                    n += 1
        say("[MAGMA_PFILE] %d SNPs" % n)
    outp = MAGMA_OUT / "GCST90205183_35_10"
    if not Path(str(outp) + ".genes.out").exists():
        t0 = time.time()
        rc, out, err = sh([str(binpath), "--bfile", str(bfile), "--pval", str(pfile), "N=%d" % GWAS_N,
                           "--gene-annot", str(annot) + ".genes.annot", "--out", str(outp)], timeout=4 * 3600)
        say("[MAGMA_GENE] rc=%d %.0f min" % (rc, (time.time() - t0) / 60))
        if rc != 0:
            for l in (out + err).strip().splitlines()[-8:]:
                say("       " + l[:200])
    gout = Path(str(outp) + ".genes.out")
    if not gout.exists():
        return
    # Map Entrez IDs to symbols from the gene location file
    sym = {}
    for line in open(gloc):
        t = line.split()
        if len(t) >= 6:
            sym[t[0]] = t[5]
    rows = []
    with open(gout) as f:
        hdr = f.readline().split()
        ix = {h: i for i, h in enumerate(hdr)}
        for line in f:
            t = line.split()
            rows.append((float(t[ix["P"]]), sym.get(t[ix["GENE"]], t[ix["GENE"]]), float(t[ix["ZSTAT"]])))
    rows.sort()
    bonf = 0.05 / max(1, len(rows))
    nsig = sum(1 for r in rows if r[0] < bonf)
    say("[MAGMA_GENES] tested=%d bonferroni(p<%.1e)=%d" % (len(rows), bonf, nsig))
    for p, s, z in rows[:15]:
        say("    [TOP] %-12s z=%6.2f p=%.2e" % (s, z, p))
    pos = {s: i for i, (p, s, z) in enumerate(rows)}
    hits = [(g, pos[g] + 1, rows[pos[g]][0]) for g in KNOWN_ENDO_GENES if g in pos]
    say("[MAGMA_KNOWN_LOCI] " + "; ".join("%s rank=%d p=%.1e" % h for h in hits))
    n_known_sig = sum(1 for h in hits if h[2] < bonf)
    SUM["magma_known_bonferroni"] = n_known_sig
    SUM["gates"]["MAGMA_SANITY"] = "PASS" if n_known_sig >= 5 else "CHECK"
    say("[GATE] MAGMA_SANITY %s (%d known endometriosis genes Bonferroni-significant)" % (SUM["gates"]["MAGMA_SANITY"], n_known_sig))


# ----------------------------------------------------------------------------
# B. Paper 8 GSE266265 export: structure only
# ----------------------------------------------------------------------------
def inspect_p8_export():
    say("")
    say("-" * 72)
    say("[SECTION] GSE266265 EXPORT FROM PAPER 8 (read-only structure)")
    if not P8_EXPORT.exists():
        say("[P8_EXPORT] MISSING " + str(P8_EXPORT))
        return
    root = P8_EXPORT.resolve()
    files = sorted(p for p in root.rglob("*") if p.is_file())
    say("[P8_EXPORT] %s files=%d total=%s" % (root, len(files), human(sum(p.stat().st_size for p in files))))
    names = None
    for p in files:
        if p.name.endswith("meta.data_names.json"):
            names = json.loads(p.read_text())
    if names:
        say("[META_COLUMNS] %d" % len(names))
        for i, n in enumerate(names):
            say("    [COL] %2d %s" % (i, n))
    for p in files:
        rel = str(p.relative_to(root))
        size = p.stat().st_size
        if size > 50_000_000:
            say("    [BIG] %-70s %10s" % (rel[:70], human(size)))
    for key in ["sparse_vectors.tsv", "manifest.json", "run_manifest.json", "metadata_manifest_excerpt.json"]:
        for p in files:
            if p.name == key:
                txt = p.read_text(errors="ignore")
                say("[FILE] %s (%s)" % (key, human(p.stat().st_size)))
                for line in txt.splitlines()[:45]:
                    say("    " + line[:220])
    for p in files:
        m = re.match(r"^\d+_meta\.data_(\d+)__levels\.json$", p.name)
        if m:
            idx = int(m.group(1))
            try:
                lv = json.loads(p.read_text())
            except Exception:
                lv = []
            name = names[idx] if names and idx < len(names) else "?"
            say("[LEVELS] col %d %-28s n=%d: %s" % (idx, name, len(lv), " | ".join(str(x) for x in lv[:40])[:600]))
    other = [p for p in files if re.match(r"^\d+_meta\.data_\d+", p.name) and "levels" not in p.name]
    say("[META_FILES] %d per-column files, e.g.: %s" % (len(other), ", ".join(p.name for p in other[:12])))


# ----------------------------------------------------------------------------
# C. HECA
# ----------------------------------------------------------------------------
H5AD_INSPECT = r"""
import sys, h5py, numpy as np
p = sys.argv[1]
f = h5py.File(p, "r")
print("[H5AD] keys:", list(f.keys()))
X = f["X"]
if isinstance(X, h5py.Group):
    print("[H5AD] X:", dict(X.attrs).get("encoding-type"), "shape", X.attrs.get("shape"), "nnz", X["data"].shape[0], "dtype", X["data"].dtype)
    d = X["data"][:2000000]
    print("[H5AD] X sample integer-valued:", bool(np.all(np.mod(d, 1) == 0)), "max", float(d.max()))
else:
    print("[H5AD] X dense", X.shape, X.dtype)
if "layers" in f:
    print("[H5AD] layers:", list(f["layers"].keys()))
if "raw" in f:
    print("[H5AD] raw keys:", list(f["raw"].keys()))
obs = f["obs"]
cols = [k for k in obs.keys() if k not in ("_index", "__categories")]
print("[H5AD] n_obs", obs[obs.attrs.get("_index", "_index")].shape[0], "obs columns:", len(cols))
oldcats = obs["__categories"] if "__categories" in obs else None
for k in cols:
    g = obs[k]
    cats = None
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in g["categories"][:]]
        codes = g["codes"][:]
    elif oldcats is not None and k in oldcats:
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in oldcats[k][:]]
        codes = g[:]
    if cats is not None:
        cnt = np.bincount(codes[codes >= 0], minlength=len(cats))
        pairs = sorted(zip(cats, cnt), key=lambda x: -x[1])
        shown = " | ".join("%s=%d" % (c, n) for c, n in pairs[:45])
        print("[OBS] %-28s categorical n=%d :: %s" % (k, len(cats), shown[:900]))
    else:
        ds = g if not isinstance(g, h5py.Group) else None
        print("[OBS] %-28s %s" % (k, (ds.dtype if ds is not None else "group")))
var = f["var"]
print("[H5AD] var columns:", [k for k in var.keys()][:20])
idx = var[var.attrs.get("_index", "_index")][:5]
print("[H5AD] var index sample:", [x.decode() if isinstance(x, bytes) else x for x in idx])
"""


def heca():
    say("")
    say("-" * 72)
    say("[SECTION] HECA (endometriumAtlasV2_cells_with_counts.h5ad)")
    dst = HECA_DIR / "endometriumAtlasV2_cells_with_counts.h5ad"
    ok = download(HECA_URL, dst, "HECA cells_with_counts")
    SUM["gates"]["HECA_DOWNLOAD"] = ok
    if not ok:
        return
    rc, out, err = sh([str(VENV_PY), "-c", H5AD_INSPECT, str(dst)], timeout=3600)
    for line in (out + err).strip().splitlines():
        say(line[:1000])
    SUM["gates"]["HECA_READ"] = rc == 0


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.5 | MAGMA + ATLAS INSPECTION | " + TS)
    say("=" * 72)
    steps = [("magma", None), ("p8_export", inspect_p8_export), ("heca", heca)]
    binpath = None
    try:
        binpath = setup_magma()
    except Exception as e:
        say("[STEP_ERROR] magma %s: %s" % (type(e).__name__, str(e)[:200]))
    if binpath:
        try:
            magma_sanity(binpath)
        except Exception as e:
            say("[STEP_ERROR] magma_sanity %s: %s" % (type(e).__name__, str(e)[:200]))
    for name, fn in steps[1:]:
        try:
            fn()
        except Exception as e:
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:200]))
    with open(OUT / "p9_v005_magma_atlas_inspect_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    for k, v in SUM["gates"].items():
        say("[GATE] %-14s %s" % (k, v))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] No disease effects computed; Paper 8 export read-only")
    say("[DONE] PAPER9_V005_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
