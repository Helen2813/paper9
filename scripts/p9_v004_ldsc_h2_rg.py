#!/usr/bin/env python3
# Paper 9 | v0.0.4 | = v003 with two fixes: plain-text input to munge_sumstats (gzip header bug
#                   in the Python 3 port) and MAGMA link discovery from the page the old URLs return.
#                   Steps: LDSC validation, MAGMA, munge, SNP-heritability, genetic correlations (G3, G5)
# Standard library only; LDSC runs in tools/ldsc_venv. Re-run safe (finished steps are reused).

import csv
import datetime as dt
import gzip
import json
import math
import os
import re
import shutil
import subprocess
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
INTERIM = REPO / "data" / "interim"
PROC = REPO / "data" / "processed"
TOOLS = REPO / "tools"
LDSC_PY = TOOLS / "ldsc_venv" / "bin" / "python"
LDSC_DIR = TOOLS / "ldsc_py3"
LDSC_TARBALL = "https://github.com/belowlab/ldsc/archive/refs/heads/2-to-3.tar.gz"
LDSC_PYPI = "ldsc==2.0.1"

REF_UNPACKED = RAW / "ldsc_ref" / "unpacked"
HM3_LIST = RAW / "ldsc_ref" / "hm3_no_MHC.list.txt"
FINNGEN_DIR = RAW / "finngen"
FINNGEN_MANIFEST = FINNGEN_DIR / "finngen_R13_manifest.tsv"
GWAS_DIR = RAW / "gwas" / "GCST90205183"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "ldsc" / ("p9_v004_ldsc_h2_rg_" + TS)
LOG = REPO / "logs" / ("p9_v004_ldsc_h2_rg_" + TS + ".log")
HM3_DIR = INTERIM / "sumstats_hm3"
MUNGED = PROC / "munged"
MERGE_ALLELES = INTERIM / "ldsc" / "hm3_merge_alleles.tsv"

# Traits: code, label, role, population prevalence for liability scale (None = observed scale only)
TRAITS = [
    ("N14_ENDOMETRIOSIS", "ENDO_ALL", "all", 0.10),
    ("N14_ENDOMETRIOSIS_OVARY", "ENDO_OVARY", "G3", None),
    ("N14_ENDOMETRIOSIS_PELVICPERITONEUM", "ENDO_PERITONEUM", "G3", None),
    ("N14_ENDOMETRIOSIS_DEEP", "ENDO_DEEP", "G3", None),
    ("N14_ENDOMETRIOSIS_RECTPVAGSEPT_VAGINA", "ENDO_RECTOVAG", "G3", None),
    ("N14_ENDOMETRIOSIS_ASRM_STAGE3_4", "ENDO_STAGE34", "G3_secondary", None),
    ("N14_ENDOMETRIOSIS_UTERUS", "ADENOMYOSIS", "G5", None),
    ("CD2_BENIGN_LEIOMYOMA_UTERI", "LEIOMYOMA", "G5", None),
    ("I9_HYPTENS", "HYPERTENSION", "G5_control", None),
    ("T2D", "T2D", "G5_control", None),
]
G3_PRIMARY = ["ENDO_OVARY", "ENDO_PERITONEUM", "ENDO_DEEP", "ENDO_RECTOVAG"]
H2_Z_GATE = 4.0

GWASCAT_PATHS = ["harmonised/GCST90205183.h.tsv.gz", "harmonised/md5sum.txt",
                 "GCST90205183_buildGRCh37.tsv.gz-meta.yaml", "README.txt"]
GWASCAT_MIRRORS = [
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90205001-GCST90206000/GCST90205183/",
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90205001-GCST90206000/GCST90205183/",
    "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90205001-GCST90206000/GCST90205183/",
]
GWASCAT_API = "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90205183"
GWASCAT_TOTAL_N = 474160

MAGMA_DIR = TOOLS / "magma"
MAGMA_URLS = {
    "mac": ["https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10_mac.zip"],
    "source": ["https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10_source.zip",
               "https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10_src.zip"],
    "genes_b37": ["https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip"],
    "genes_b38": ["https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI38.zip"],
    "ref_eur": ["https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip"],
}

# ----------------------------------------------------------------------------
# Logging and helpers
# ----------------------------------------------------------------------------
for d in (OUT, LOG.parent, HM3_DIR, MUNGED, MERGE_ALLELES.parent, MAGMA_DIR):
    d.mkdir(parents=True, exist_ok=True)
_fh = open(LOG, "w", encoding="utf-8")
SUM = {"task": "PAPER9_V004_LDSC_H2_RG", "timestamp_utc": TS, "gates": {}}


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
    expect = remote_size(url) if url.startswith("http") else None
    if dst.exists() and (expect is None or dst.stat().st_size == expect) and dst.stat().st_size > 0:
        say("[DL] %-52s %10s SKIP_PRESENT" % (label[:52], human(dst.stat().st_size)))
        return True
    part = dst.with_name(dst.name + ".part")
    base = ["curl", "-L", "-sS", "--fail", "--retry", "4", "--retry-delay", "5",
            "--connect-timeout", "30", "-A", "paper9", "-o", str(part), url]
    rc, out, err = sh(base[:1] + ["-C", "-"] + base[1:])
    if rc == 33:
        if part.exists():
            part.unlink()
        rc, out, err = sh(base)
    if rc != 0 or not part.exists() or part.stat().st_size == 0:
        say("[DL] %-52s FAILED rc=%s %s" % (label[:52], rc, err.strip()[:140]))
        return False
    if expect and part.stat().st_size != expect:
        say("[DL] %-52s SIZE_MISMATCH (re-run to resume)" % label[:52])
        return False
    part.rename(dst)
    say("[DL] %-52s %10s OK" % (label[:52], human(dst.stat().st_size)))
    return True


def norm_sf(z):
    # Two-sided p-value for a standard normal statistic
    return math.erfc(abs(z) / math.sqrt(2.0))


def holm(pvals):
    order = sorted(range(len(pvals)), key=lambda i: pvals[i])
    adj = [None] * len(pvals)
    running = 0.0
    m = len(pvals)
    for rank, i in enumerate(order):
        val = min(1.0, (m - rank) * pvals[i])
        running = max(running, val)
        adj[i] = running
    return adj


def safe_extract_tar(tgz, dest, mode="r:gz"):
    with tarfile.open(tgz, mode) as tf:
        for m in tf.getmembers():
            if m.name.startswith("/") or ".." in Path(m.name).parts:
                continue
            tf.extract(m, dest)


# ----------------------------------------------------------------------------
# A. LDSC Python 3 port + validation by the bundled test suite
# ----------------------------------------------------------------------------
def setup_ldsc():
    say("")
    say("-" * 72)
    say("[SECTION] LDSC PYTHON 3 PORT (belowlab 2-to-3, = PyPI ldsc 2.0.1)")
    if not LDSC_PY.exists():
        say("[FATAL_STEP] ldsc venv missing: " + str(LDSC_PY))
        return None, None
    if not (LDSC_DIR / "ldsc.py").exists():
        tgz = TOOLS / "ldsc_py3_2-to-3.tar.gz"
        if download(LDSC_TARBALL, tgz, "ldsc 2-to-3 source"):
            tmp = TOOLS / "ldsc_py3_tmp"
            if tmp.exists():
                shutil.rmtree(tmp)
            tmp.mkdir()
            safe_extract_tar(tgz, tmp)
            inner = [p for p in tmp.iterdir() if p.is_dir()]
            if inner:
                if LDSC_DIR.exists():
                    shutil.rmtree(LDSC_DIR)
                inner[0].rename(LDSC_DIR)
            shutil.rmtree(tmp, ignore_errors=True)
    ldsc_py = LDSC_DIR / "ldsc.py"
    munge_py = LDSC_DIR / "munge_sumstats.py"
    if not ldsc_py.exists():
        say("[LDSC] tarball unavailable -> pip " + LDSC_PYPI)
        sh([str(LDSC_PY), "-m", "pip", "install", "-q", LDSC_PYPI], timeout=1800)
        cands = list((TOOLS / "ldsc_venv").rglob("ldsc.py"))
        mcands = list((TOOLS / "ldsc_venv").rglob("munge_sumstats.py"))
        ldsc_py = cands[0] if cands else None
        munge_py = mcands[0] if mcands else None
    if not ldsc_py or not Path(ldsc_py).exists():
        say("[LDSC] NOT_AVAILABLE")
        return None, None
    rc, out, err = sh([str(LDSC_PY), "-m", "pip", "install", "-q", "nose", "bitarray==2.8.2"], timeout=900)
    for script in [ldsc_py, munge_py]:
        rc, out, err = sh([str(LDSC_PY), str(script), "-h"], timeout=300)
        tail = (out + err).strip().splitlines()
        say("[LDSC_HELP] %-18s rc=%d %s" % (Path(script).name, rc, (tail[0] if tail else "")[:100]))
        if rc != 0:
            for line in tail[-5:]:
                say("       " + line[:200])
    # Bundled nose test suite (upstream test data)
    tests_status = "NOT_RUN"
    if (LDSC_DIR / "test").exists():
        t0 = time.time()
        rc, out, err = sh([str(LDSC_PY), "-m", "nose", "-v", "test"], cwd=str(LDSC_DIR), timeout=3600)
        txt = (out + err)
        ran = re.search(r"Ran (\d+) tests? in", txt)
        failed = re.search(r"FAILED \(([^)]*)\)", txt)
        tests_status = "PASS" if rc == 0 else "FAIL"
        say("[LDSC_TESTS] %s rc=%d ran=%s %s (%.0fs)" % (tests_status, rc, ran.group(1) if ran else "?",
                                                        failed.group(0) if failed else "", time.time() - t0))
        if rc != 0:
            fails = [l for l in txt.splitlines() if l.startswith(("FAIL:", "ERROR:"))]
            for l in fails[:25]:
                say("       " + l[:200])
            (OUT / "ldsc_tests_output.txt").write_text(txt)
    SUM["ldsc_tests"] = tests_status
    SUM["gates"]["LDSC_TESTS"] = tests_status
    return Path(ldsc_py), (Path(munge_py) if munge_py else None)


# ----------------------------------------------------------------------------
# B. MAGMA (needed later for G1/G4; installed now)
# ----------------------------------------------------------------------------
def magma_links_from_html(html, base_host):
    found = []
    for m in re.finditer(r'href="([^"]+)"', html):
        href = m.group(1)
        ctx = re.sub(r"<[^>]+>", " ", html[max(0, m.start() - 200):m.start()])
        ctx = re.sub(r"\s+", " ", ctx).strip()[-100:]
        if href.startswith("//"):
            url = "https:" + href
        elif href.startswith("http"):
            url = href
        elif href.startswith("/"):
            url = base_host + href
        else:
            continue
        if re.search(r"\.zip|download|surf|/s/|magma", url, flags=re.I):
            found.append((url, ctx))
    return found


def setup_magma():
    say("")
    say("-" * 72)
    say("[SECTION] MAGMA")
    # Remove non-zip files saved by v003 (HTML pages) but keep one copy for link discovery
    pages = []
    for z in sorted(MAGMA_DIR.glob("*.zip")):
        if not zipfile.is_zipfile(z):
            html = z.read_text(errors="ignore")
            pages.append(html)
            z.rename(z.with_suffix(".html"))
    for h in MAGMA_DIR.glob("*.html"):
        pages.append(h.read_text(errors="ignore"))
    for page_url in ["https://ctg.thebluebus.nl/software/magma", "https://cncr.nl/research/magma/"]:
        p = MAGMA_DIR / ("page_" + re.sub(r"\W+", "_", page_url) + ".html")
        if not p.exists():
            download(page_url, p, "MAGMA page " + page_url.split("/")[2])
        if p.exists():
            pages.append(p.read_text(errors="ignore"))
    title = ""
    for html in pages:
        m = re.search(r"<title>(.*?)</title>", html, flags=re.I | re.S)
        if m:
            title = m.group(1).strip()[:100]
            break
    say("[MAGMA_PAGE_TITLE] " + (title or "none"))
    links = []
    for html in pages:
        for host in ["https://ctg.thebluebus.nl", "https://cncr.nl"]:
            for u, c in magma_links_from_html(html, host):
                if (u, c) not in links:
                    links.append((u, c))
    say("[MAGMA_LINKS] %d candidates" % len(links))
    for u, c in links[:60]:
        say("    [LINK] %s  <<%s>>" % (u[:160], c))

    def pick(name_pat, ctx_pat):
        for u, c in links:
            if re.search(name_pat, u, flags=re.I):
                return u
        for u, c in links:
            if re.search(ctx_pat, c, flags=re.I):
                return u
        return None

    wanted = {
        "mac": pick(r"magma_v1\.10[^/]*mac[^/]*\.zip", r"mac\s*os[^|]{0,40}$"),
        "source": pick(r"magma_v1\.10[^/]*(source|src)[^/]*\.zip", r"source\s*code[^|]{0,40}$"),
        "genes_b37": pick(r"NCBI37[^/]*\.zip", r"build\s*37[^|]{0,40}$"),
        "genes_b38": pick(r"NCBI38[^/]*\.zip", r"build\s*38[^|]{0,40}$"),
        "ref_eur": pick(r"g1000_eur[^/]*\.zip", r"european[^|]{0,40}$"),
    }
    for key, u in wanted.items():
        say("[MAGMA_WANT] %-10s %s" % (key, u or "NOT_FOUND"))
        z = MAGMA_DIR / (key + ".zip")
        if u and not (z.exists() and zipfile.is_zipfile(z)):
            if z.exists():
                z.unlink()
            download(u, z, "MAGMA " + key)
        if z.exists() and not zipfile.is_zipfile(z):
            say("[MAGMA_DL] %s is not a zip (%s)" % (key, human(z.stat().st_size)))
            z.rename(z.with_suffix(".notzip.html"))
        sub = MAGMA_DIR / key
        if z.exists() and zipfile.is_zipfile(z) and not sub.exists():
            with zipfile.ZipFile(z) as zf:
                zf.extractall(sub)
    rc, out, err = sh(["arch", "-x86_64", "/usr/bin/true"])
    rosetta = rc == 0
    say("[ROSETTA] " + ("INSTALLED" if rosetta else "NOT_INSTALLED"))

    def try_bin(p):
        p.chmod(0o755)
        rc, out, err = sh([str(p), "--version"], timeout=60)
        txt = out + err
        m = re.search(r"MAGMA\s+version[:\s]*v?([\d.]+\w*)", txt, flags=re.I) or re.search(r"v1\.\d+\w*", txt)
        return (m.group(0) if m else None), txt

    status, version, binpath = "NOT_INSTALLED", None, None
    for cand in sorted((MAGMA_DIR / "mac").rglob("magma")) if (MAGMA_DIR / "mac").exists() else []:
        if cand.is_file():
            rc, out, err = sh(["file", str(cand)])
            say("[MAGMA_MAC_FILE] " + out.strip()[:200])
            version, txt = try_bin(cand)
            if version:
                status, binpath = "MAC_BINARY_OK", cand
            else:
                say("[MAGMA_MAC_RUN] " + txt.strip()[:200])
            break
    if not version and (MAGMA_DIR / "source").exists():
        mk = next(iter((MAGMA_DIR / "source").rglob("Makefile")), None)
        if mk:
            text = mk.read_text(errors="ignore")
            patched = re.sub(r"\s-static\b", "", text)
            if patched != text:
                mk.write_text(patched)
                say("[MAGMA_COMPILE] removed -static from Makefile (not supported on macOS)")
            rc, out, err = sh(["make"], cwd=str(mk.parent), timeout=3600)
            lines = (out + err).strip().splitlines()
            say("[MAGMA_COMPILE] make rc=%d %s" % (rc, (lines[-1] if lines else "")[:200]))
            if rc != 0:
                for l in lines[-8:]:
                    say("       " + l[:200])
            built = mk.parent / "magma"
            if built.exists():
                version, txt = try_bin(built)
                if version:
                    status, binpath = "COMPILED_FROM_SOURCE_OK", built
        else:
            say("[MAGMA_COMPILE] no Makefile found in source archive")
    if binpath:
        link = MAGMA_DIR / "magma_bin"
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(binpath)
    for key in ("genes_b37", "genes_b38", "ref_eur"):
        files = sorted(p.name for p in (MAGMA_DIR / key).rglob("*") if p.is_file()) if (MAGMA_DIR / key).exists() else []
        say("[MAGMA_DATA] %-10s %s" % (key, ", ".join(files[:6]) if files else "MISSING"))
    say("[MAGMA] %s version=%s" % (status, version))
    SUM["gates"]["MAGMA"] = status


# ----------------------------------------------------------------------------
# C. GWAS Catalog retry
# ----------------------------------------------------------------------------
def gwascat_retry():
    say("")
    say("-" * 72)
    say("[SECTION] GWAS CATALOG GCST90205183 (retry with mirrors)")
    ok_main = False
    for rel in GWASCAT_PATHS:
        got = False
        for base in GWASCAT_MIRRORS:
            got = download(base + rel, GWAS_DIR / rel, "GCST " + rel + " via " + base.split(":")[0])
            if got:
                break
        if rel.endswith(".h.tsv.gz"):
            ok_main = got
    md5f = GWAS_DIR / "harmonised" / "md5sum.txt"
    if ok_main and md5f.exists():
        import hashlib
        for line in md5f.read_text().splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[1].endswith(".h.tsv.gz"):
                h = hashlib.md5()
                with open(GWAS_DIR / "harmonised" / parts[1], "rb") as f:
                    for chunk in iter(lambda: f.read(8 << 20), b""):
                        h.update(chunk)
                good = h.hexdigest() == parts[0]
                say("[MD5] %s %s" % (parts[1], "OK" if good else "MISMATCH"))
                ok_main = ok_main and good
    rc, out, err = sh(["curl", "-sSL", "--fail", "--max-time", "60", GWASCAT_API])
    if rc == 0:
        try:
            js = json.loads(out)
            say("[GWASCAT_API] initialSampleSize: " + str(js.get("initialSampleSize"))[:300])
            SUM["gwascat_initial_sample"] = js.get("initialSampleSize")
        except Exception:
            say("[GWASCAT_API] unparsable")
    else:
        say("[GWASCAT_API] not reachable")
    SUM["gates"]["GWASCAT"] = ok_main
    return ok_main


# ----------------------------------------------------------------------------
# D. HapMap3 merge-alleles file and pre-filtered sumstats
# ----------------------------------------------------------------------------
def load_hm3():
    s = set()
    with open(HM3_LIST, encoding="utf-8", errors="ignore") as f:
        for line in f:
            tok = line.split()
            if tok and tok[0].startswith("rs"):
                s.add(tok[0])
    return s


def build_merge_alleles(hm3):
    if MERGE_ALLELES.exists() and MERGE_ALLELES.stat().st_size > 1000:
        n = sum(1 for _ in open(MERGE_ALLELES)) - 1
        say("[MERGE_ALLELES] present (%d SNPs)" % n)
        return n
    bims = sorted(REF_UNPACKED.rglob("*.bim"))
    n = 0
    with open(MERGE_ALLELES, "w") as out:
        out.write("SNP\tA1\tA2\n")
        for b in bims:
            with open(b) as f:
                for line in f:
                    t = line.split()
                    if len(t) >= 6 and t[1] in hm3:
                        out.write("%s\t%s\t%s\n" % (t[1], t[4], t[5]))
                        n += 1
    say("[MERGE_ALLELES] built from %d bim files: %d SNPs" % (len(bims), n))
    return n


def prefilter_finngen(code, label, hm3):
    src = None
    for p in FINNGEN_DIR.glob("*" + code + ".gz"):
        if p.name.endswith("_" + code + ".gz"):
            src = p
    dst = HM3_DIR / (label + ".tsv")
    if dst.exists() and dst.stat().st_size > 1000:
        return dst
    if src is None:
        say("[PREFILTER] %s SOURCE_MISSING" % code)
        return None
    t0 = time.time()
    n_in = n_out = 0
    with gzip.open(src, "rt") as f, open(dst, "w") as o:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        o.write("SNP\tA1\tA2\tBETA\tSE\tP\tFRQ\n")
        for line in f:
            n_in += 1
            t = line.rstrip("\n").split("\t")
            rs = t[ix["rsids"]]
            if not rs:
                continue
            hit = None
            for r in rs.split(","):
                if r in hm3:
                    hit = r
                    break
            if hit is None:
                continue
            o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (hit, t[ix["alt"]], t[ix["ref"]], t[ix["beta"]],
                                                     t[ix["sebeta"]], t[ix["pval"]], t[ix["af_alt"]]))
            n_out += 1
    say("[PREFILTER] %-16s in=%d hm3=%d %.0fs" % (label, n_in, n_out, time.time() - t0))
    return dst


def prefilter_gwascat(hm3):
    src = GWAS_DIR / "harmonised" / "GCST90205183.h.tsv.gz"
    dst = HM3_DIR / "GCST90205183.tsv"
    if dst.exists() and dst.stat().st_size > 1000:
        return dst
    if not src.exists():
        return None
    n_in = n_out = 0
    with gzip.open(src, "rt") as f, open(dst, "w") as o:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        o.write("SNP\tA1\tA2\tBETA\tSE\tP\tFRQ\n")
        for line in f:
            n_in += 1
            t = line.rstrip("\n").split("\t")
            rs = t[ix["rsid"]]
            if rs not in hm3:
                continue
            o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (rs, t[ix["effect_allele"]], t[ix["other_allele"]],
                                                     t[ix["beta"]], t[ix["standard_error"]], t[ix["p_value"]],
                                                     t[ix["effect_allele_frequency"]]))
            n_out += 1
    say("[PREFILTER] %-16s in=%d hm3=%d" % ("GCST90205183", n_in, n_out))
    return dst


def manifest_counts():
    out = {}
    if FINNGEN_MANIFEST.exists():
        lines = FINNGEN_MANIFEST.read_text().splitlines()
        hdr = lines[0].split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        for line in lines[1:]:
            t = line.split("\t")
            out[t[0]] = (int(t[ix["num_cases"]]), int(t[ix["num_controls"]]))
    return out


def munge(label, src, ncas=None, ncon=None, n=None, munge_py=None):
    outp = MUNGED / label
    if Path(str(outp) + ".sumstats.gz").exists():
        return Path(str(outp) + ".sumstats.gz")
    cmd = [str(LDSC_PY), str(munge_py), "--sumstats", str(src), "--snp", "SNP", "--a1", "A1", "--a2", "A2",
           "--p", "P", "--signed-sumstats", "BETA,0", "--frq", "FRQ",
           "--merge-alleles", str(MERGE_ALLELES), "--chunksize", "500000", "--out", str(outp)]
    if ncas is not None:
        cmd += ["--N-cas", str(ncas), "--N-con", str(ncon)]
    else:
        cmd += ["--N", str(n)]
    rc, out, err = sh(cmd, timeout=3600)
    logf = Path(str(outp) + ".log")
    txt = logf.read_text(errors="ignore") if logf.exists() else (out + err)
    kept = re.search(r"Writing summary statistics for (\d+) SNPs", txt)
    chi = re.search(r"Mean chi\^2 = ([\d.]+)", txt)
    say("[MUNGE] %-16s rc=%d snps=%s mean_chi2=%s" % (label, rc, kept.group(1) if kept else "?",
                                                    chi.group(1) if chi else "?"))
    if rc != 0:
        for l in (out + err).strip().splitlines()[-6:]:
            say("       " + l[:200])
        return None
    return Path(str(outp) + ".sumstats.gz")


# ----------------------------------------------------------------------------
# E. LDSC h2 and rg
# ----------------------------------------------------------------------------
def ld_prefixes():
    def find(dirname, must=None):
        d = REF_UNPACKED / dirname
        for p in sorted(d.rglob("*.l2.ldscore.gz")) if d.exists() else []:
            m = re.match(r"^(.*\.)1\.l2\.ldscore\.gz$", p.name)
            if m and (must is None or must in p.name):
                return str(p.parent / m.group(1))
        return None
    ref = find("LDscore")
    w = find("1000G_Phase3_weights_hm3_no_MHC")
    say("[LD_REF] ref=%s" % ref)
    say("[LD_WEIGHTS] w=%s" % w)
    return ref, w


def parse_h2(logtxt):
    r = {}
    for key, pat in [("h2_obs", r"Total Observed scale h2:\s*([-\d.e]+)\s*\(([-\d.e]+)\)"),
                     ("h2_liab", r"Total Liability scale h2:\s*([-\d.e]+)\s*\(([-\d.e]+)\)"),
                     ("intercept", r"Intercept:\s*([-\d.e]+)\s*\(([-\d.e]+)\)")]:
        m = re.search(pat, logtxt)
        if m:
            r[key] = float(m.group(1))
            r[key + "_se"] = float(m.group(2))
    for key, pat in [("lambda_gc", r"Lambda GC:\s*([-\d.e]+)"), ("mean_chi2", r"Mean Chi\^2:\s*([-\d.e]+)"),
                     ("ratio", r"Ratio:\s*([-\d.e]+)")]:
        m = re.search(pat, logtxt)
        if m:
            r[key] = float(m.group(1))
    return r


def run_h2(label, ss, ref, w, prev=None, ncas=None, ncon=None, ldsc_py=None):
    outp = OUT / ("h2_" + label)
    cmd = [str(LDSC_PY), str(ldsc_py), "--h2", str(ss), "--ref-ld-chr", ref, "--w-ld-chr", w, "--out", str(outp)]
    if prev is not None and ncas:
        cmd += ["--samp-prev", "%.6f" % (ncas / float(ncas + ncon)), "--pop-prev", str(prev)]
    rc, out, err = sh(cmd, timeout=3600)
    logf = Path(str(outp) + ".log")
    txt = logf.read_text(errors="ignore") if logf.exists() else out + err
    r = parse_h2(txt)
    r["rc"] = rc
    if rc != 0:
        for l in (out + err).strip().splitlines()[-6:]:
            say("       " + l[:200])
    return r


def run_rg(labels, files, ref, w, ldsc_py):
    rows = []
    for i in range(len(labels) - 1):
        anchor = labels[i]
        others = labels[i + 1:]
        outp = OUT / ("rg_" + anchor)
        cmd = [str(LDSC_PY), str(ldsc_py), "--rg", ",".join(str(files[x]) for x in [anchor] + others),
               "--ref-ld-chr", ref, "--w-ld-chr", w, "--out", str(outp)]
        t0 = time.time()
        rc, out, err = sh(cmd, timeout=7200)
        logf = Path(str(outp) + ".log")
        txt = logf.read_text(errors="ignore") if logf.exists() else out + err
        if "Summary of Genetic Correlation Results" in txt:
            block = txt.split("Summary of Genetic Correlation Results")[1].strip().splitlines()
            hdr = block[0].split()
            for line in block[1:]:
                t = line.split()
                if len(t) != len(hdr):
                    break
                d = dict(zip(hdr, t))
                p1 = Path(d["p1"]).name.replace(".sumstats.gz", "")
                p2 = Path(d["p2"]).name.replace(".sumstats.gz", "")
                row = {"trait1": p1, "trait2": p2}
                for k in ["rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se"]:
                    try:
                        row[k] = float(d.get(k, "nan"))
                    except ValueError:
                        row[k] = float("nan")
                rows.append(row)
        say("[RG] anchor=%-16s partners=%d rc=%d %.0fs" % (anchor, len(others), rc, time.time() - t0))
        if rc != 0:
            for l in (out + err).strip().splitlines()[-5:]:
                say("       " + l[:200])
    return rows


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    return ("%." + str(nd) + "f") % x


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    say("=" * 72)
    say("PAPER 9 | v0.0.4 | LDSC h2 + rg (G3 primary, G5 descriptive) | " + TS)
    say("=" * 72)
    ldsc_py, munge_py = setup_ldsc()
    setup_magma()
    gwas_ok = gwascat_retry()
    if ldsc_py is None or munge_py is None:
        say("[STOP] LDSC unavailable; analysis not run")
        finish()
        return

    say("")
    say("-" * 72)
    say("[SECTION] PREPARE SUMMARY STATISTICS (HapMap3, no MHC)")
    hm3 = load_hm3()
    say("[HM3] %d rsIDs" % len(hm3))
    build_merge_alleles(hm3)
    counts = manifest_counts()
    files = {}
    meta = {}
    for code, label, role, prev in TRAITS:
        src = prefilter_finngen(code, label, hm3)
        if not src:
            continue
        ncas, ncon = counts.get(code, (None, None))
        ss = munge(label, src, ncas=ncas, ncon=ncon, munge_py=munge_py)
        if ss:
            files[label] = ss
            meta[label] = {"code": code, "role": role, "prev": prev, "ncas": ncas, "ncon": ncon}
    if gwas_ok:
        src = prefilter_gwascat(hm3)
        ss = munge("GCST90205183", src, n=GWASCAT_TOTAL_N, munge_py=munge_py) if src else None
        if ss:
            files["GCST90205183"] = ss
            meta["GCST90205183"] = {"code": "GCST90205183", "role": "primary_G1_G2_G4", "prev": None,
                                    "ncas": None, "ncon": None}

    ref, w = ld_prefixes()
    if not ref or not w:
        say("[STOP] LD reference prefixes not found")
        finish()
        return

    say("")
    say("-" * 72)
    say("[SECTION] SNP-HERITABILITY")
    h2rows = []
    for label in files:
        m = meta[label]
        r = run_h2(label, files[label], ref, w, prev=m["prev"], ncas=m["ncas"], ncon=m["ncon"], ldsc_py=ldsc_py)
        z = r["h2_obs"] / r["h2_obs_se"] if r.get("h2_obs_se") else float("nan")
        r.update({"trait": label, "h2_z": z, "gate_h2z_ge4": bool(z >= H2_Z_GATE), **m})
        h2rows.append(r)
        say("[H2] %-16s cases=%-6s h2_obs=%s (%s) z=%s int=%s (%s) ratio=%s liab=%s gate=%s" % (
            label, m["ncas"], fmt(r.get("h2_obs"), 4), fmt(r.get("h2_obs_se"), 4), fmt(z, 1),
            fmt(r.get("intercept"), 3), fmt(r.get("intercept_se"), 3), fmt(r.get("ratio"), 3),
            fmt(r.get("h2_liab"), 3), "PASS" if r["gate_h2z_ge4"] else "FAIL"))
    with open(OUT / "h2_table.csv", "w", newline="") as f:
        cols = ["trait", "code", "role", "ncas", "ncon", "h2_obs", "h2_obs_se", "h2_z", "h2_liab", "h2_liab_se",
                "intercept", "intercept_se", "ratio", "lambda_gc", "mean_chi2", "gate_h2z_ge4"]
        wr = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        wr.writeheader()
        for r in h2rows:
            wr.writerow(r)
    gate = {r["trait"]: r["gate_h2z_ge4"] for r in h2rows}

    say("")
    say("-" * 72)
    say("[SECTION] GENETIC CORRELATIONS (all pairs)")
    labels = [x for x in files]
    rg_rows = run_rg(labels, files, ref, w, ldsc_py)
    with open(OUT / "rg_all_pairs.csv", "w", newline="") as f:
        cols = ["trait1", "trait2", "rg", "se", "z", "p", "gcov_int", "gcov_int_se", "h2_obs", "h2_obs_se"]
        wr = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        wr.writeheader()
        for r in rg_rows:
            wr.writerow(r)
    for r in rg_rows:
        say("[RG_PAIR] %-16s %-16s rg=%s (%s) p=%s gcov_int=%s" % (r["trait1"], r["trait2"], fmt(r["rg"]),
                                                                    fmt(r["se"]), "%.2e" % r["p"] if not math.isnan(r["p"]) else "NA",
                                                                    fmt(r["gcov_int"])))

    say("")
    say("-" * 72)
    say("[SECTION] G3 PRIMARY TEST: rg = 1 between anatomical subtypes (Holm)")
    eligible = [t for t in G3_PRIMARY if gate.get(t)]
    say("[G3_ELIGIBLE] " + ", ".join(eligible) + " | excluded by h2 z<4: " +
        (", ".join(t for t in G3_PRIMARY if t in gate and not gate[t]) or "none"))
    prim = []
    for r in rg_rows:
        if r["trait1"] in eligible and r["trait2"] in eligible and not math.isnan(r["se"]) and r["se"] > 0:
            z1 = (1.0 - r["rg"]) / r["se"]
            prim.append({"pair": r["trait1"] + " ~ " + r["trait2"], "rg": r["rg"], "se": r["se"],
                         "lo95": r["rg"] - 1.96 * r["se"], "hi95": r["rg"] + 1.96 * r["se"],
                         "z_vs_1": z1, "p_vs_1": norm_sf(z1), "gcov_int": r["gcov_int"]})
    if prim:
        adj = holm([x["p_vs_1"] for x in prim])
        for x, a in zip(prim, adj):
            x["p_vs_1_holm"] = a
            x["verdict"] = "rg<1 SUPPORTED" if (a < 0.05 and x["rg"] < 1) else "NOT REJECTED"
            say("[G3] %-34s rg=%s [%s, %s] z(1-rg)=%s p=%.2e Holm=%.2e %s" % (
                x["pair"], fmt(x["rg"]), fmt(x["lo95"]), fmt(x["hi95"]), fmt(x["z_vs_1"], 2),
                x["p_vs_1"], a, x["verdict"]))
        with open(OUT / "G3_primary_rg_vs_1.csv", "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=list(prim[0].keys()))
            wr.writeheader()
            for x in prim:
                wr.writerow(x)
    else:
        say("[G3] no eligible pairs")
    SUM["G3_primary"] = prim

    say("")
    say("-" * 72)
    say("[SECTION] G5 DESCRIPTIVE: endometriosis vs comparators (cell-type tests come later)")
    for r in rg_rows:
        pair = {r["trait1"], r["trait2"]}
        if "ENDO_ALL" in pair and pair & {"ADENOMYOSIS", "LEIOMYOMA", "HYPERTENSION", "T2D", "GCST90205183"}:
            other = (pair - {"ENDO_ALL"}).pop()
            say("[G5] ENDO_ALL ~ %-14s rg=%s (%s) gcov_int=%s" % (other, fmt(r["rg"]), fmt(r["se"]), fmt(r["gcov_int"])))
    SUM["h2"] = h2rows
    SUM["rg_pairs"] = len(rg_rows)
    finish()


def finish():
    with open(OUT / "p9_v004_ldsc_h2_rg_summary.json", "w") as f:
        json.dump(SUM, f, indent=2, default=str)
    say("")
    for k, v in SUM["gates"].items():
        say("[GATE] %-12s %s" % (k, v))
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] G3 primary family and h2 gate fixed in PLAN v0.2 / CHANGELOG v0.2.1 before this run")
    say("[DONE] PAPER9_V004_LDSC_H2_RG_COMPLETE")
    _fh.close()


if __name__ == "__main__":
    main()
