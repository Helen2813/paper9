#!/usr/bin/env python3
# Paper 9 | v0.0.2 | TOOLS + CORE DATA (no analysis, no disease effects)
# 1) Python environments (main + pinned LDSC env), LDSC, S-PrediXcan, MAGMA
# 2) Downloads: GWAS Catalog, FinnGen R13 endpoints, LDSC references, PredictDB, GSE213216 RAW
# 3) Read-only links to Paper 7 / Paper 8 data in data/external
# 4) Annotation-source probes (HECA, Tan 2022), integrity checks, gates
# Re-run safe: finished downloads are skipped, interrupted downloads resume.

import datetime as dt
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time
import zipfile
from pathlib import Path

# ----------------------------------------------------------------------------
# Hardcoded locations
# ----------------------------------------------------------------------------
HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
RAW = REPO / "data" / "raw"
EXT = REPO / "data" / "external"
TOOLS = REPO / "tools"
VENV = REPO / ".venv"
LDSC_VENV = TOOLS / "ldsc_venv"
P8 = HOME / "Desktop" / "paper8" / "endo-research"
P7 = HOME / "Desktop" / "paper7"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "setup" / ("p9_v002_tools_and_data_" + TS)
LOG = REPO / "logs" / ("p9_v002_tools_and_data_" + TS + ".log")

BASE_PY = "/usr/bin/python3"

MAIN_PKGS = ["numpy", "pandas", "scipy", "statsmodels", "pyarrow", "openpyxl", "xlrd",
             "h5py", "bitarray", "patsy"]
SC_PKGS = ["anndata", "scanpy"]
LDSC_PKGS = ["numpy==1.23.5", "pandas==1.5.3", "scipy==1.10.1", "bitarray==2.8.2", "nose"]

LDSC_TARBALL = "https://github.com/bulik/ldsc/archive/refs/heads/master.tar.gz"
METAXCAN_TARBALL = "https://github.com/hakyimlab/MetaXcan/archive/refs/heads/master.tar.gz"
MAGMA_PAGE = "https://ctg.thebluebus.nl/software/magma"

GWASCAT = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
           "GCST90205001-GCST90206000/GCST90205183/")
GWASCAT_FILES = ["GCST90205183_buildGRCh37.tsv.gz", "GCST90205183_buildGRCh37.tsv.gz-meta.yaml",
                 "README.txt", "md5sum.txt",
                 "harmonised/GCST90205183.h.tsv.gz", "harmonised/GCST90205183.h.tsv.gz-meta.yaml",
                 "harmonised/md5sum.txt"]
GWASCAT_API = "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90205183"

FINNGEN_MANIFEST = "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv"
FINNGEN_ENDPOINTS = [
    # code, role
    ("N14_ENDOMETRIOSIS", "all"),
    ("N14_ENDOMETRIOSIS_OVARY", "G3 subtype"),
    ("N14_ENDOMETRIOSIS_PELVICPERITONEUM", "G3 subtype"),
    ("N14_ENDOMETRIOSIS_DEEP", "G3 subtype"),
    ("N14_ENDOMETRIOSIS_RECTPVAGSEPT_VAGINA", "G3 subtype"),
    ("N14_ENDOMETRIOSIS_ASRM_STAGE3_4", "G3 severity"),
    ("N14_ENDOMETRIOSIS_UTERUS", "G5 adenomyosis"),
    ("CD2_BENIGN_LEIOMYOMA_UTERI", "G5 fibroid"),
    ("I9_HYPTENS", "G5 non-reproductive control"),
    ("T2D", "G5 non-reproductive control"),
]

LDSC_ZENODO = "10515792"
LDSC_FILES = ["readme_baseline_versions.txt", "1000G_Phase3_baseline_v1.2_ldscores.tgz",
              "1000G_Phase3_baselineLD_v2.2_ldscores.tgz", "1000G_Phase3_weights_hm3_no_MHC.tgz",
              "1000G_Phase3_plinkfiles.tgz", "1000G_Phase3_frq.tgz", "hm3_no_MHC.list.txt",
              "1000G_Phase3_ldscores.tgz"]
PREDICTDB_ZENODO = "3518299"
PREDICTDB_FILES = ["mashr_eqtl.tar", "gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz",
                   "eur_ld.bed.gz", "README.md"]

GSE213216_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213216/suppl/"
GSE141549_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141549/suppl/"
HECA_PAGE = "https://www.reproductivecellatlas.org/endometrium_reference.html"
HECA_GITHUB_API = "https://api.github.com/repos/ventolab/HECA-Human-Endometrial-Cell-Atlas/contents"

LINKS = {
    "p8_GSE141549": P8 / "data" / "p8_multicohort_full_sources" / "GSE141549",
    "p8_GSE51981_series_matrix.txt.gz": P8 / "data" / "p8_multicohort_full_sources" / "GSE51981" / "GSE51981_series_matrix.txt.gz",
    "p8_GSE266265_v232c_export": P8 / "results" / "paper8_v232c_full_rna_disk_export_20261103T230602Z",
    "p8_GSE213216_sample_metadata.tsv": P8 / "research" / "p8_01_sources" / "gse213216_sample_metadata.tsv",
    "p8_GSE213216_family.soft": P8 / "research" / "p8_01_sources" / "gse213216_family.soft",
    "p7_GSE246591_verified_sample_manifest.csv": P7 / "publication" / "recomb2027_paper7" / "evidence" / "GSE246591_verified_sample_manifest.csv",
}

# ----------------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------------
for d in (OUT, LOG.parent, RAW, EXT, TOOLS):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V002_TOOLS_AND_DATA", "timestamp_utc": TS, "downloads": [], "tools": {}, "links": {}}


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


def sh(cmd, cwd=None, timeout=None, env=None):
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout, env=env)
        return p.returncode, p.stdout, p.stderr
    except FileNotFoundError as e:
        return 127, "", str(e)
    except OSError as e:
        return 126, "", str(e)
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"


def curl_text(url, max_time=60):
    rc, out, err = sh(["curl", "-sSL", "--fail", "--max-time", str(max_time), "-A", "paper9-setup", url])
    return out if rc == 0 else None


def remote_size(url):
    # Size from the final response only (redirect hops carry their own content-length)
    rc, out, err = sh(["curl", "-sSIL", "--max-time", "60", "-A", "paper9-setup", url])
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


def download(url, dst, label, expect=None):
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if expect is None:
        expect = remote_size(url)
    if dst.exists() and (expect is None or dst.stat().st_size == expect):
        say("[DL] %-58s %10s SKIP_PRESENT" % (label[:58], human(dst.stat().st_size)))
        SUM["downloads"].append({"label": label, "path": str(dst), "bytes": dst.stat().st_size, "status": "PRESENT"})
        return True
    part = dst.with_name(dst.name + ".part")
    t0 = time.time()
    say("[DL] %-58s %10s START" % (label[:58], human(expect)))
    cmd = ["curl", "-L", "-sS", "--fail", "--retry", "5", "--retry-delay", "5",
           "-A", "paper9-setup", "-o", str(part), url]
    rc, out, err = sh(cmd[:1] + ["-C", "-"] + cmd[1:])
    if rc == 33:  # server does not support resume: restart from zero
        if part.exists():
            part.unlink()
        rc, out, err = sh(cmd)
    if rc not in (0,) and part.exists() and expect and part.stat().st_size == expect:
        rc = 0
    if rc != 0:
        say("[DL] %-58s FAILED rc=%s %s" % (label[:58], rc, err.strip()[:160]))
        SUM["downloads"].append({"label": label, "path": str(dst), "status": "FAILED", "error": err.strip()[:200]})
        return False
    size = part.stat().st_size
    if expect and size != expect:
        say("[DL] %-58s SIZE_MISMATCH got=%s expect=%s (re-run to resume)" % (label[:58], size, expect))
        SUM["downloads"].append({"label": label, "path": str(dst), "status": "SIZE_MISMATCH"})
        return False
    part.rename(dst)
    dt_s = time.time() - t0
    say("[DL] %-58s %10s OK %.0fs" % (label[:58], human(size), dt_s))
    SUM["downloads"].append({"label": label, "path": str(dst), "bytes": size, "status": "OK", "seconds": round(dt_s)})
    return True


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_ok(path):
    rc, out, err = sh(["gzip", "-t", str(path)], timeout=1800)
    return rc == 0


# ----------------------------------------------------------------------------
# A. Python environments
# ----------------------------------------------------------------------------
def make_venv(path, label):
    py = path / "bin" / "python"
    if not py.exists():
        rc, out, err = sh([BASE_PY, "-m", "venv", str(path)])
        if rc != 0:
            say("[VENV] %s CREATE_FAILED %s" % (label, err.strip()[:200]))
            return None
    sh([str(py), "-m", "pip", "install", "-q", "--upgrade", "pip", "setuptools", "wheel"], timeout=900)
    return py


def pip_install(py, pkgs, label):
    rc, out, err = sh([str(py), "-m", "pip", "install", "-q"] + pkgs, timeout=3600)
    status = "OK" if rc == 0 else "FAILED"
    say("[PIP] %-10s %-6s %s" % (label, status, " ".join(pkgs)))
    if rc != 0:
        say("      " + err.strip().splitlines()[-1][:200] if err.strip() else "      (no stderr)")
    return rc == 0


def versions(py, mods):
    code = ("import importlib\n"
            "for m in %r:\n"
            "    try:\n"
            "        x = importlib.import_module(m); print(m, getattr(x, '__version__', '?'))\n"
            "    except Exception as e:\n"
            "        print(m, 'IMPORT_FAILED', type(e).__name__)\n") % (mods,)
    rc, out, err = sh([str(py), "-c", code], timeout=300)
    for line in out.strip().splitlines():
        say("      [VER] " + line)


def setup_python():
    say("")
    say("-" * 72)
    say("[SECTION] PYTHON ENVIRONMENTS")
    rc, out, err = sh([BASE_PY, "--version"])
    say("[BASE_PYTHON] " + BASE_PY + " " + (out or err).strip())
    py = make_venv(VENV, "main")
    if py:
        say("[VENV] main " + str(py))
        pip_install(py, MAIN_PKGS, "main")
        pip_install(py, SC_PKGS, "singlecell")
        versions(py, ["numpy", "pandas", "scipy", "statsmodels", "pyarrow", "openpyxl", "h5py", "anndata", "scanpy"])
    SUM["tools"]["venv_main"] = str(py) if py else None
    lpy = make_venv(LDSC_VENV, "ldsc")
    if lpy:
        say("[VENV] ldsc " + str(lpy))
        if not pip_install(lpy, LDSC_PKGS, "ldsc"):
            pip_install(lpy, ["numpy<2", "pandas<2", "scipy", "bitarray", "nose"], "ldsc_alt")
        versions(lpy, ["numpy", "pandas", "scipy", "bitarray"])
    SUM["tools"]["venv_ldsc"] = str(lpy) if lpy else None
    return py, lpy


# ----------------------------------------------------------------------------
# B. LDSC, S-PrediXcan, MAGMA
# ----------------------------------------------------------------------------
def fetch_github_tarball(url, dest_dir, label):
    if dest_dir.exists() and any(dest_dir.iterdir()):
        say("[TOOL_SRC] %s present %s" % (label, dest_dir))
        return True
    tgz = TOOLS / (label + "_master.tar.gz")
    if not download(url, tgz, label + " source (GitHub tarball)"):
        return False
    tmp = TOOLS / (label + "_extract_tmp")
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir()
    with tarfile.open(tgz, "r:gz") as tf:
        for m in tf.getmembers():
            if m.name.startswith("/") or ".." in Path(m.name).parts:
                continue
            tf.extract(m, tmp)
    inner = [p for p in tmp.iterdir() if p.is_dir()]
    if not inner:
        return False
    if dest_dir.exists():
        shutil.rmtree(dest_dir)
    inner[0].rename(dest_dir)
    shutil.rmtree(tmp)
    say("[TOOL_SRC] %s extracted to %s" % (label, dest_dir))
    return True


def setup_ldsc(lpy):
    ok = False
    d = TOOLS / "ldsc"
    if lpy and fetch_github_tarball(LDSC_TARBALL, d, "ldsc"):
        for script in ["ldsc.py", "munge_sumstats.py"]:
            rc, out, err = sh([str(lpy), str(d / script), "-h"], timeout=300)
            txt = (out + err).strip().splitlines()
            say("[LDSC] %-18s rc=%d %s" % (script, rc, (txt[-1] if txt else "")[:150]))
            if script == "ldsc.py":
                ok = rc == 0
            if rc != 0:
                for line in (out + err).strip().splitlines()[-6:]:
                    say("       " + line[:200])
    SUM["tools"]["ldsc_ok"] = ok
    return ok


def setup_metaxcan(py):
    ok = False
    d = TOOLS / "MetaXcan"
    if py and fetch_github_tarball(METAXCAN_TARBALL, d, "MetaXcan"):
        script = d / "software" / "SPrediXcan.py"
        rc, out, err = sh([str(py), str(script), "--help"], timeout=300)
        txt = (out + err).strip().splitlines()
        say("[SPREDIXCAN] rc=%d %s" % (rc, (txt[0] if txt else "")[:150]))
        if rc != 0:
            for line in (out + err).strip().splitlines()[-6:]:
                say("       " + line[:200])
        ok = rc == 0
    SUM["tools"]["spredixcan_ok"] = ok
    return ok


def magma_run(binpath):
    rc, out, err = sh([str(binpath), "--version"], timeout=60)
    txt = (out + err)
    m = re.search(r"MAGMA\s+version[:\s]*([\w.]+)", txt, flags=re.I) or re.search(r"v1\.\d+\w*", txt)
    return rc, txt, (m.group(0) if m else None)


def setup_magma():
    say("")
    say("-" * 72)
    say("[SECTION] MAGMA")
    mdir = TOOLS / "magma"
    mdir.mkdir(exist_ok=True)
    html = curl_text(MAGMA_PAGE)
    links = []  # (url, context text preceding the link)
    if html:
        for m in re.finditer(r'href="([^"]+)"', html):
            href = m.group(1)
            ctx = re.sub(r"<[^>]+>", " ", html[max(0, m.start() - 160):m.start()])
            ctx = re.sub(r"\s+", " ", ctx).strip()[-90:]
            if re.search(r"\.zip($|\?)|download|surf|/s/", href, flags=re.I):
                url = href if href.startswith("http") else "https://ctg.thebluebus.nl" + (href if href.startswith("/") else "/" + href)
                links.append((url, ctx))
    say("[MAGMA_PAGE] %s candidate links=%d" % ("OK" if html else "FAILED", len(links)))
    for u, c in links:
        say("    [LINK] %s  <<%s>>" % (u, c))

    def pick(name_pat, ctx_pat):
        for u, c in links:
            if re.search(name_pat, u, flags=re.I):
                return u
        for u, c in links:
            if re.search(ctx_pat, c, flags=re.I):
                return u
        return None

    wanted = {
        "mac": pick(r"magma_v1\.10[^/]*mac[^/]*\.zip", r"mac\s*os[^|]*$"),
        "source": pick(r"magma_v1\.10[^/]*(source|src)[^/]*\.zip", r"source\s*code[^|]*$"),
        "genes_b37": pick(r"NCBI37[^/]*\.zip", r"build\s*37[^|]*$"),
        "genes_b38": pick(r"NCBI38[^/]*\.zip", r"build\s*38[^|]*$"),
        "ref_eur": pick(r"g1000_eur[^/]*\.zip", r"european[^|]*$"),
    }
    for k, u in wanted.items():
        say("[MAGMA_WANT] %-10s %s" % (k, u or "NOT_FOUND_ON_PAGE"))
    for k, u in wanted.items():
        if not u:
            continue
        z = mdir / (k + ".zip")
        if download(u, z, "MAGMA " + k):
            sub = mdir / k
            if not sub.exists():
                try:
                    with zipfile.ZipFile(z) as zf:
                        zf.extractall(sub)
                except Exception as e:
                    say("[MAGMA_UNZIP] %s FAILED %s" % (k, e))
    # Locate and test binary
    status = "NOT_INSTALLED"
    version = None
    binpath = None
    for cand in sorted((mdir / "mac").rglob("magma")) if (mdir / "mac").exists() else []:
        if cand.is_file():
            cand.chmod(0o755)
            binpath = cand
            break
    if binpath:
        rc, txt, version = magma_run(binpath)
        if version:
            status = "MAC_BINARY_OK"
        else:
            say("[MAGMA_MAC_BINARY] not runnable: " + txt.strip()[:200])
    if status != "MAC_BINARY_OK" and (mdir / "source").exists():
        mk = None
        for cand in (mdir / "source").rglob("Makefile"):
            mk = cand
            break
        if mk:
            rc, out, err = sh(["make"], cwd=str(mk.parent), timeout=1800)
            say("[MAGMA_COMPILE] make rc=%d %s" % (rc, (err or out).strip().splitlines()[-1][:200] if (err or out).strip() else ""))
            built = mk.parent / "magma"
            if built.exists():
                binpath = built
                rc, txt, version = magma_run(built)
                if version:
                    status = "COMPILED_FROM_SOURCE_OK"
        else:
            say("[MAGMA_COMPILE] no Makefile in source archive")
    if binpath and version:
        link = mdir / "magma_bin"
        if link.exists() or link.is_symlink():
            link.unlink()
        link.symlink_to(binpath)
    say("[MAGMA] %s version=%s binary=%s" % (status, version, binpath))
    SUM["tools"]["magma_status"] = status
    SUM["tools"]["magma_binary"] = str(binpath) if binpath else None
    return status.endswith("_OK")


# ----------------------------------------------------------------------------
# C. Downloads
# ----------------------------------------------------------------------------
def zenodo_files(rec):
    txt = curl_text("https://zenodo.org/api/records/" + rec)
    if not txt:
        return {}
    js = json.loads(txt)
    out = {}
    for f in js.get("files", []):
        key = f.get("key")
        link = (f.get("links") or {}).get("self") or ("https://zenodo.org/records/%s/files/%s?download=1" % (rec, key))
        out[key] = (link, f.get("size"), (f.get("checksum") or "").replace("md5:", ""))
    return out


def download_gwas_catalog():
    say("")
    say("-" * 72)
    say("[SECTION] GWAS CATALOG GCST90205183")
    d = RAW / "gwas" / "GCST90205183"
    ok = True
    for f in GWASCAT_FILES:
        ok &= download(GWASCAT + f, d / f, "GCST90205183/" + f)
    for md5file, folder in [(d / "md5sum.txt", d), (d / "harmonised" / "md5sum.txt", d / "harmonised")]:
        if md5file.exists():
            for line in md5file.read_text().splitlines():
                parts = line.split()
                if len(parts) == 2 and (folder / parts[1]).exists() and parts[1].endswith(".gz"):
                    good = md5(folder / parts[1]) == parts[0]
                    say("[MD5] %-50s %s" % (parts[1], "OK" if good else "MISMATCH"))
                    ok &= good
    api = curl_text(GWASCAT_API)
    if api:
        try:
            js = json.loads(api)
            say("[GWASCAT_API] initialSampleSize: " + str(js.get("initialSampleSize"))[:300])
            say("[GWASCAT_API] replicationSampleSize: " + str(js.get("replicationSampleSize"))[:300])
            say("[GWASCAT_API] trait: " + str((js.get("diseaseTrait") or {}).get("trait")))
            SUM["gwascat_initial_sample"] = js.get("initialSampleSize")
        except Exception as e:
            say("[GWASCAT_API] parse failed " + str(e))
    else:
        say("[GWASCAT_API] not reachable")
    SUM["gwascat_ok"] = ok
    return ok


def download_finngen():
    say("")
    say("-" * 72)
    say("[SECTION] FINNGEN R13 (about 765 MB per endpoint)")
    txt = curl_text(FINNGEN_MANIFEST)
    if not txt:
        say("[FINNGEN] manifest not reachable")
        SUM["finngen_ok"] = 0
        return 0
    lines = txt.splitlines()
    hdr = lines[0].split("\t")
    idx = {h: i for i, h in enumerate(hdr)}
    rows = {}
    for line in lines[1:]:
        r = line.split("\t")
        rows[r[0]] = r
    (RAW / "finngen").mkdir(parents=True, exist_ok=True)
    (RAW / "finngen" / "finngen_R13_manifest.tsv").write_text(txt)
    n_ok = 0
    chosen = []
    for code, role in FINNGEN_ENDPOINTS:
        r = rows.get(code)
        if r is None:
            alt = [c for c in rows if code.split("_")[-1] in c and c.startswith(code.split("_")[0])]
            say("[FG] %-40s NOT_IN_MANIFEST; similar: %s" % (code, ", ".join(alt[:6]) or "none"))
            continue
        cases = r[idx["num_cases"]]
        ctrls = r[idx["num_controls"]]
        url = r[idx["path_https"]]
        say("[FG] %-40s cases=%-7s controls=%-7s role=%s" % (code, cases, ctrls, role))
        chosen.append({"code": code, "role": role, "cases": cases, "controls": ctrls, "url": url})
        dst = RAW / "finngen" / url.split("/")[-1]
        if download(url, dst, "FinnGen " + code):
            n_ok += 1
    for c in chosen:
        dst = RAW / "finngen" / c["url"].split("/")[-1]
        c["gzip_ok"] = gzip_ok(dst) if dst.exists() else False
        say("[GZIP_TEST] %-40s %s" % (c["code"], "OK" if c["gzip_ok"] else "FAILED"))
    SUM["finngen_endpoints"] = chosen
    SUM["finngen_ok"] = n_ok
    return n_ok


def download_zenodo_set(rec, files, subdir, label):
    say("")
    say("-" * 72)
    say("[SECTION] " + label + " (Zenodo " + rec + ")")
    meta = zenodo_files(rec)
    if not meta:
        say("[ZENODO] %s metadata not reachable" % rec)
        return False
    ok = True
    d = RAW / subdir
    for f in files:
        if f not in meta:
            say("[ZENODO] %-50s NOT_IN_RECORD" % f)
            ok = False
            continue
        link, size, checksum = meta[f]
        got = download(link, d / f, label + "/" + f, expect=size)
        if got and checksum and (size or 0) < 3 * 1024 ** 3:
            good = md5(d / f) == checksum
            say("[MD5] %-50s %s" % (f, "OK" if good else "MISMATCH"))
            got = got and good
        ok &= got
    return ok


def unpack_ldsc_refs():
    d = RAW / "ldsc_ref"
    u = d / "unpacked"
    u.mkdir(parents=True, exist_ok=True)
    for f in sorted(d.glob("*.tgz")):
        marker = u / (f.name + ".done")
        if marker.exists():
            continue
        try:
            with tarfile.open(f, "r:gz") as tf:
                for m in tf.getmembers():
                    if m.name.startswith("/") or ".." in Path(m.name).parts:
                        continue
                    tf.extract(m, u)
            marker.write_text("ok")
            say("[UNPACK] " + f.name)
        except Exception as e:
            say("[UNPACK] %s FAILED %s" % (f.name, e))
    tops = sorted(p.name for p in u.iterdir() if p.is_dir())
    say("[LDSC_REF_DIRS] " + ", ".join(tops))
    SUM["ldsc_ref_dirs"] = tops


def unpack_predictdb():
    d = RAW / "predictdb"
    t = d / "mashr_eqtl.tar"
    u = d / "unpacked"
    if t.exists() and not (u / "mashr_eqtl.tar.done").exists():
        u.mkdir(parents=True, exist_ok=True)
        with tarfile.open(t, "r:") as tf:
            for m in tf.getmembers():
                if m.name.startswith("/") or ".." in Path(m.name).parts:
                    continue
                tf.extract(m, u)
        (u / "mashr_eqtl.tar.done").write_text("ok")
    dbs = sorted(p.name for p in u.rglob("*.db")) if u.exists() else []
    keep = [x for x in dbs if re.search(r"Uterus|Ovary|Vagina|Fallopian|Whole_Blood|Spleen|Adipose|Colon|Small_Intestine", x)]
    say("[PREDICTDB] models=%d; reproductive/related: %s" % (len(dbs), ", ".join(keep)))
    SUM["predictdb_models"] = len(dbs)
    SUM["predictdb_selected"] = keep


def download_gse213216():
    say("")
    say("-" * 72)
    say("[SECTION] GSE213216 (Tan 2022) supplementary")
    d = RAW / "geo" / "GSE213216"
    download(GSE213216_SUPPL + "filelist.txt", d / "filelist.txt", "GSE213216/filelist.txt")
    fl = d / "filelist.txt"
    if fl.exists():
        lines = fl.read_text(errors="ignore").splitlines()
        say("[GSE213216_FILELIST] %d lines" % len(lines))
        kinds = {}
        for line in lines[1:]:
            parts = line.split("\t")
            name = parts[1] if len(parts) > 1 else parts[0]
            ext = re.sub(r"^GSM\d+_", "", name)
            ext = re.sub(r"^[^_]*_", "", ext) if "_" in ext else ext
            kinds[ext] = kinds.get(ext, 0) + 1
        for line in lines[:12]:
            say("    " + line[:160])
        for k, v in sorted(kinds.items(), key=lambda kv: -kv[1])[:15]:
            say("    [KIND] %4d  %s" % (v, k[:100]))
        annot = [l for l in lines if re.search(r"annot|cell_?type|metadata|cluster|label", l, re.I)]
        say("[GSE213216_ANNOTATION_FILES_IN_RAW] %d" % len(annot))
        for l in annot[:10]:
            say("    [ANNOT] " + l[:160])
        SUM["gse213216_annotation_files"] = len(annot)
    say("[INFO] GSE213216_RAW.tar is 15.7 GB; this is the longest step and can be resumed by re-running")
    ok = download(GSE213216_SUPPL + "GSE213216_RAW.tar", d / "GSE213216_RAW.tar", "GSE213216_RAW.tar")
    if ok:
        rc, out, err = sh(["tar", "-tf", str(d / "GSE213216_RAW.tar")], timeout=3600)
        names = out.strip().splitlines()
        say("[TAR_LIST] entries=%d rc=%d" % (len(names), rc))
        for n in names[:8]:
            say("    " + n)
        SUM["gse213216_tar_entries"] = len(names)
    SUM["gse213216_raw_ok"] = ok
    return ok


# ----------------------------------------------------------------------------
# D. Links to other papers' data (read-only)
# ----------------------------------------------------------------------------
def make_links():
    say("")
    say("-" * 72)
    say("[SECTION] READ-ONLY LINKS TO PAPER 7 / PAPER 8 DATA")
    for name, target in LINKS.items():
        link = EXT / name
        if not target.exists():
            say("[LINK] %-45s MISSING_TARGET %s" % (name, target))
            SUM["links"][name] = None
            continue
        if link.is_symlink() or link.exists():
            if link.is_symlink() and os.readlink(link) == str(target):
                say("[LINK] %-45s present" % name)
                SUM["links"][name] = str(target)
                continue
            if link.is_symlink():
                link.unlink()
        link.symlink_to(target)
        say("[LINK] %-45s -> %s" % (name, target))
        SUM["links"][name] = str(target)
    # Inspect GSE141549 folder from Paper 8
    g = LINKS["p8_GSE141549"]
    if g.exists():
        for f in sorted(g.iterdir()):
            if f.is_file():
                say("    [P8_GSE141549] %-75s %10s" % (f.name[:75], human(f.stat().st_size)))
        if not (g / "GSE141549_Sample_link.xlsx").exists():
            download(GSE141549_SUPPL + "GSE141549_Sample_link.xlsx",
                     RAW / "geo" / "GSE141549" / "GSE141549_Sample_link.xlsx", "GSE141549_Sample_link.xlsx")
        for mtx in ["GSE141549-GPL13376_series_matrix.txt.gz", "GSE141549-GPL10558_series_matrix.txt.gz"]:
            if not (g / mtx).exists():
                download("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141549/matrix/" + mtx,
                         RAW / "geo" / "GSE141549" / mtx, "GSE141549/" + mtx)
    # Inspect Paper 8 GSE266265 export (structure only)
    e = LINKS["p8_GSE266265_v232c_export"]
    if e.exists():
        files = sorted(p for p in e.rglob("*") if p.is_file())
        say("[P8_V232C_EXPORT] files=%d" % len(files))
        for f in files[:40]:
            say("    %-80s %10s" % (str(f.relative_to(e))[:80], human(f.stat().st_size)))
        for f in files:
            if f.suffix.lower() in (".json", ".csv", ".tsv", ".txt") and f.stat().st_size < 300_000:
                try:
                    head = f.read_text(errors="ignore").splitlines()[:3]
                except Exception:
                    head = []
                if head:
                    say("    [HEAD] %s :: %s" % (f.name[:50], head[0][:220]))


# ----------------------------------------------------------------------------
# E. Annotation-source probes
# ----------------------------------------------------------------------------
def probe_annotations():
    say("")
    say("-" * 72)
    say("[SECTION] CELL-ANNOTATION SOURCES")
    html = curl_text(HECA_PAGE)
    if html:
        hrefs = sorted(set(h for h in re.findall(r'href="([^"]+)"', html)
                           if re.search(r"\.(h5ad|rds|csv|tsv|zip|gz|h5)(\?|$)", h, re.I)))
        say("[HECA_PAGE] OK data links=%d" % len(hrefs))
        for h in hrefs[:40]:
            say("    [HECA_LINK] " + h[:200])
        SUM["heca_links"] = hrefs
    else:
        say("[HECA_PAGE] not reachable")
    gh = curl_text(HECA_GITHUB_API)
    if gh:
        try:
            names = [x.get("name") for x in json.loads(gh)]
            say("[HECA_GITHUB] " + ", ".join(names)[:400])
        except Exception:
            say("[HECA_GITHUB] unexpected response")
    soft = LINKS["p8_GSE213216_family.soft"]
    if soft.exists():
        txt = soft.read_text(errors="ignore")
        sup = sorted(set(re.findall(r"!Series_supplementary_file\s*=\s*(\S+)", txt)))
        say("[GSE213216_SERIES_SUPPL] " + " | ".join(sup)[:400])
        hits = sorted(set(re.findall(r"(https?://\S*(?:zenodo|cellxgene|figshare|github|synapse)\S*)", txt, re.I)))
        for h in hits[:10]:
            say("    [GSE213216_EXTERNAL_LINK] " + h[:200])


# ----------------------------------------------------------------------------
# F. Gates
# ----------------------------------------------------------------------------
def gates():
    say("")
    say("-" * 72)
    say("[SECTION] GATES")
    t = SUM["tools"]
    fg_total = len(SUM.get("finngen_endpoints", []))
    fg_gz = sum(1 for c in SUM.get("finngen_endpoints", []) if c.get("gzip_ok"))
    g = {
        "TOOLS_LDSC": bool(t.get("ldsc_ok")),
        "TOOLS_MAGMA": str(t.get("magma_status", "")).endswith("_OK"),
        "TOOLS_SPREDIXCAN": bool(t.get("spredixcan_ok")),
        "DATA_GWASCAT": bool(SUM.get("gwascat_ok")),
        "DATA_FINNGEN": fg_total > 0 and fg_gz == fg_total,
        "DATA_LDSC_REF": len(SUM.get("ldsc_ref_dirs", [])) > 0,
        "DATA_PREDICTDB": SUM.get("predictdb_models", 0) > 0,
        "DATA_GSE213216_RAW": bool(SUM.get("gse213216_raw_ok")),
        "LINKS_P7_P8": all(v for v in SUM["links"].values()) if SUM["links"] else False,
    }
    for k, v in g.items():
        say("[GATE] %-22s %s" % (k, "PASS" if v else "FAIL"))
    say("[FINNGEN_FILES_VERIFIED] %d/%d" % (fg_gz, fg_total))
    SUM["gates"] = g
    ready_v003 = g["TOOLS_LDSC"] and g["DATA_FINNGEN"] and g["DATA_LDSC_REF"] and g["DATA_GWASCAT"]
    say("[NEXT] v003 (munge + h2 + subtype rg, G3/G5 first) " + ("READY" if ready_v003 else "BLOCKED - send this log"))


def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.2 | TOOLS + CORE DATA | " + TS)
    say("=" * 72)
    du = shutil.disk_usage(str(HOME))
    say("[DISK_FREE_START] " + human(du.free) + " (about 30 GB will be downloaded)")
    steps = [
        ("python", lambda: SUM.__setitem__("_py", setup_python())),
        ("ldsc", lambda: setup_ldsc(Path(SUM["tools"]["venv_ldsc"]) if SUM["tools"].get("venv_ldsc") else None)),
        ("metaxcan", lambda: setup_metaxcan(Path(SUM["tools"]["venv_main"]) if SUM["tools"].get("venv_main") else None)),
        ("magma", setup_magma),
        ("links", make_links),
        ("annotations", probe_annotations),
        ("gwascat", download_gwas_catalog),
        ("ldsc_ref", lambda: download_zenodo_set(LDSC_ZENODO, LDSC_FILES, "ldsc_ref", "LDSC_REF")),
        ("ldsc_unpack", unpack_ldsc_refs),
        ("predictdb", lambda: download_zenodo_set(PREDICTDB_ZENODO, PREDICTDB_FILES, "predictdb", "PREDICTDB")),
        ("predictdb_unpack", unpack_predictdb),
        ("finngen", download_finngen),
        ("gse213216", download_gse213216),
    ]
    for name, fn in steps:
        try:
            fn()
        except Exception as e:
            say("[STEP_ERROR] %s %s: %s" % (name, type(e).__name__, str(e)[:240]))
            SUM.setdefault("errors", []).append({"step": name, "error": repr(e)[:300]})
    SUM.pop("_py", None)
    gates()
    du = shutil.disk_usage(str(HOME))
    say("[DISK_FREE_END] " + human(du.free))
    with open(OUT / "p9_v002_tools_and_data_summary.json", "w", encoding="utf-8") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] No analysis; no disease effects; Paper 7/8 data linked read-only, never modified")
    say("[DONE] PAPER9_V002_TOOLS_AND_DATA_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
