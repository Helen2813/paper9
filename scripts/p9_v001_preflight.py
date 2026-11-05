#!/usr/bin/env python3
# Paper 9 | v0.0.1 | PREFLIGHT (data + tools feasibility, no analysis, no large downloads)
# Self-contained, standard library only. Uses system curl for all network access.
# Reads only small byte ranges of remote files. Writes nothing outside the repo.

import csv
import datetime as dt
import gzip
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import zlib
from pathlib import Path

# ----------------------------------------------------------------------------
# Hardcoded locations
# ----------------------------------------------------------------------------
HOME = Path.home()
REPO = HOME / "Desktop" / "paper9"
PAPER8_ROOT = HOME / "Desktop" / "paper8" / "endo-research"
DESKTOP = HOME / "Desktop"

TS = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT = REPO / "results" / "preflight" / ("p9_v001_preflight_" + TS)
LOG = REPO / "logs" / ("p9_v001_preflight_" + TS + ".log")

# GWAS sources
GWAS_CATALOG_ID = "GCST90205183"
GWAS_CATALOG_DIR = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
                    "GCST90205001-GCST90206000/GCST90205183/")
FINNGEN_MANIFESTS = [
    ("R13", "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv"),
    ("R12", "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv"),
    ("R12", "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/R12_manifest.tsv"),
]
FINNGEN_KEYWORDS = ["ENDOMET", "LEIOMYOMA", "ADENOMYO"]
FINNGEN_PROBE_HEADERS = ["N14_ENDOMETRIOSIS", "N14_ENDOMETRIOSIS_OVARY",
                         "N14_ENDOMETRIOSIS_PELVICPERITONEUM", "N14_ENDOMETRIOSIS_DEEP"]

# Expression sources (GEO)
GEO_SERIES = {
    "GSE213216": "G1 lesion/eutopic/peritoneum single-cell atlas (Tan 2022)",
    "GSE141549": "G2 eutopic + lesion bulk cohort (Gabriel 2020)",
    "GSE266265": "G4 phase-annotated eutopic single-cell atlas (already exported in Paper 8)",
    "GSE51981": "G2 replication bulk cohort (already used in Paper 7/8)",
}

# Model / reference resources
ZENODO_RECORDS = {
    "3518299": "PredictDB GTEx v8 MASHR eQTL models (S-PrediXcan, G2)",
    "3519321": "PredictDB GTEx v8 Elastic Net models (sensitivity, G2)",
    "3569954": "GTEx v8 harmonization support data (G2)",
    "10515792": "LDSC reference resources (S-LDSC, G1/G3/G4)",
}
WEB_PAGES = {
    "MAGMA": "https://cncr.nl/research/magma/",
    "HECA_portal": "https://www.reproductivecellatlas.org/",
    "Mortlock_endometrial_eQTL": "http://reproductivegenomics.com.au/shiny/endo_eqtl_rna/",
}

TOOLS = ["git", "curl", "magma", "ldsc.py", "munge_sumstats.py", "plink", "plink2",
         "R", "Rscript", "conda", "mamba", "micromamba", "brew"]

CURL_TIMEOUT = "60"
RANGE_BYTES = 4_000_000

# ----------------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------------
OUT.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)
_log_fh = open(LOG, "w", encoding="utf-8")


def say(msg=""):
    print(msg, flush=True)
    _log_fh.write(msg + "\n")
    _log_fh.flush()


SUMMARY = {"task": "PAPER9_V001_PREFLIGHT", "timestamp_utc": TS, "sources": [], "gates": {}}
SOURCE_ROWS = []


def record(source_id, kind, url, status, detail):
    row = {"source_id": source_id, "kind": kind, "url": url, "status": status, "detail": detail}
    SOURCE_ROWS.append(row)
    SUMMARY["sources"].append(row)


# ----------------------------------------------------------------------------
# Network helpers (curl)
# ----------------------------------------------------------------------------
def curl_bytes(url, max_bytes=None):
    cmd = ["curl", "-sSL", "--fail", "--max-time", CURL_TIMEOUT, "-A", "paper9-preflight"]
    if max_bytes:
        cmd += ["-r", "0-" + str(max_bytes - 1)]
    cmd.append(url)
    try:
        p = subprocess.run(cmd, capture_output=True, timeout=int(CURL_TIMEOUT) + 15)
    except Exception as e:
        return None, "EXC " + type(e).__name__
    if p.returncode != 0:
        return None, "CURL_RC_" + str(p.returncode) + " " + p.stderr.decode("utf-8", "ignore").strip()[:160]
    data = p.stdout
    if max_bytes and len(data) > max_bytes:
        data = data[:max_bytes]
    return data, "OK"


def curl_head_size(url):
    cmd = ["curl", "-sSIL", "--max-time", CURL_TIMEOUT, "-A", "paper9-preflight", url]
    try:
        p = subprocess.run(cmd, capture_output=True, timeout=int(CURL_TIMEOUT) + 15)
    except Exception as e:
        return None, "EXC " + type(e).__name__
    txt = p.stdout.decode("utf-8", "ignore")
    codes = re.findall(r"^HTTP/\S+\s+(\d+)", txt, flags=re.M)
    sizes = re.findall(r"^content-length:\s*(\d+)", txt, flags=re.M | re.I)
    code = codes[-1] if codes else "NA"
    size = int(sizes[-1]) if sizes else None
    return size, "HTTP_" + code


def gunzip_prefix(data):
    # Decompress as much as possible from a truncated gzip stream
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    try:
        return d.decompress(data)
    except zlib.error:
        try:
            return gzip.decompress(data)
        except Exception:
            return b""


def human(n):
    if n is None:
        return "NA"
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if n < 1024:
            return "%.1f %s" % (n, unit)
        n /= 1024.0
    return "%.1f PB" % n


def parse_dir_listing(html, base):
    # Returns list of (name, url) from an Apache/nginx style listing
    out = []
    for href in re.findall(r'href="([^"?#]+)"', html):
        if href.startswith("/") or href.startswith("http") or href in ("../", "./"):
            if not href.startswith(base):
                continue
        name = href.rstrip("/").split("/")[-1]
        if not name or name in (".", ".."):
            continue
        url = href if href.startswith("http") else base + href
        out.append((name + ("/" if href.endswith("/") else ""), url))
    seen = set()
    uniq = []
    for n, u in out:
        if u not in seen:
            seen.add(u)
            uniq.append((n, u))
    return uniq


# ----------------------------------------------------------------------------
# 1. Environment
# ----------------------------------------------------------------------------
def check_environment():
    say("=" * 72)
    say("PAPER 9 | v0.0.1 | PREFLIGHT | " + TS)
    say("=" * 72)
    say("[REPO] " + str(REPO))
    say("[PYTHON] " + sys.version.split()[0] + " " + sys.executable)
    say("[PLATFORM] " + platform.platform() + " | machine=" + platform.machine())
    mem = None
    try:
        p = subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True)
        mem = int(p.stdout.strip())
    except Exception:
        pass
    say("[RAM] " + human(mem))
    du = shutil.disk_usage(str(HOME))
    say("[DISK_FREE] " + human(du.free) + " of " + human(du.total))
    SUMMARY["env"] = {"python": sys.version.split()[0], "machine": platform.machine(),
                      "ram_bytes": mem, "disk_free_bytes": du.free}
    tools = {}
    for t in TOOLS:
        path = shutil.which(t)
        tools[t] = path
        say("[TOOL] %-18s %s" % (t, path if path else "MISSING"))
    SUMMARY["tools"] = tools
    # Git state of repo
    if (REPO / ".git").exists() and shutil.which("git"):
        p = subprocess.run(["git", "-C", str(REPO), "remote", "-v"], capture_output=True, text=True)
        say("[GIT_REMOTE] " + (p.stdout.strip().splitlines()[0] if p.stdout.strip() else "NONE"))
        p = subprocess.run(["git", "-C", str(REPO), "log", "--oneline", "-n", "3"], capture_output=True, text=True)
        for line in p.stdout.strip().splitlines():
            say("[GIT_LOG] " + line)
    # Python packages useful later (not required now)
    for mod in ["numpy", "pandas", "scipy", "statsmodels", "anndata", "scanpy", "h5py", "pyarrow"]:
        try:
            __import__(mod)
            say("[PYMOD] %-12s present (this interpreter)" % mod)
        except Exception:
            say("[PYMOD] %-12s absent (this interpreter)" % mod)
    venv8 = PAPER8_ROOT / ".venv" / "bin" / "python"
    say("[PAPER8_VENV] " + (str(venv8) if venv8.exists() else "NOT FOUND"))


# ----------------------------------------------------------------------------
# 2. GWAS Catalog (Rahmioglu 2023, excl. 23andMe)
# ----------------------------------------------------------------------------
def check_gwas_catalog():
    say("")
    say("-" * 72)
    say("[SECTION] GWAS CATALOG " + GWAS_CATALOG_ID)
    data, st = curl_bytes(GWAS_CATALOG_DIR, 200_000)
    if data is None:
        say("[GWASCAT] LISTING_FAILED " + st)
        record(GWAS_CATALOG_ID, "gwas", GWAS_CATALOG_DIR, "FAIL", st)
        SUMMARY["gwas_catalog_ok"] = False
        return
    entries = parse_dir_listing(data.decode("utf-8", "ignore"), GWAS_CATALOG_DIR)
    files = []
    for name, url in entries:
        if name.endswith("/"):
            sub, st2 = curl_bytes(url, 200_000)
            if sub is not None:
                for n2, u2 in parse_dir_listing(sub.decode("utf-8", "ignore"), url):
                    if not n2.endswith("/"):
                        files.append((name + n2, u2))
        else:
            files.append((name, url))
    ok_any = False
    for name, url in files:
        size, hst = curl_head_size(url)
        say("[GWASCAT_FILE] %-60s %10s %s" % (name, human(size), hst))
        record(GWAS_CATALOG_ID, "gwas_file", url, hst, human(size))
        low = name.lower()
        if low.endswith((".yaml", ".yml", ".txt", ".md5")) and (size or 0) < 200_000:
            meta, _ = curl_bytes(url, 200_000)
            if meta:
                for line in meta.decode("utf-8", "ignore").splitlines()[:40]:
                    if re.search(r"build|assembly|sample|case|control|ancestry|effect|coordinate|file_type", line, re.I):
                        say("    [META] " + line.strip()[:150])
        if low.endswith((".tsv.gz", ".txt.gz", ".gz")) and "md5" not in low:
            raw, st3 = curl_bytes(url, 262_144)
            if raw:
                txt = gunzip_prefix(raw).decode("utf-8", "ignore").splitlines()
                if txt:
                    ok_any = True
                    say("    [HEADER] " + txt[0][:400])
                    for line in txt[1:3]:
                        say("    [ROW] " + line[:400])
        elif low.endswith(".tsv") and "md5" not in low:
            raw, st3 = curl_bytes(url, 65_536)
            if raw:
                txt = raw.decode("utf-8", "ignore").splitlines()
                if txt:
                    ok_any = True
                    say("    [HEADER] " + txt[0][:400])
                    for line in txt[1:3]:
                        say("    [ROW] " + line[:400])
    SUMMARY["gwas_catalog_ok"] = ok_any
    say("[GWASCAT] STATUS " + ("READABLE" if ok_any else "NO_READABLE_TABLE_FOUND"))


# ----------------------------------------------------------------------------
# 3. FinnGen manifest and endpoints
# ----------------------------------------------------------------------------
def check_finngen():
    say("")
    say("-" * 72)
    say("[SECTION] FINNGEN PUBLIC SUMMARY STATISTICS")
    manifest = None
    used = None
    for rel, url in FINNGEN_MANIFESTS:
        data, st = curl_bytes(url)
        say("[FINNGEN_MANIFEST] %s %s %s" % (rel, st, url))
        record("FinnGen_" + rel + "_manifest", "manifest", url, st, "")
        if data:
            manifest = data.decode("utf-8", "ignore")
            used = (rel, url)
            break
    if manifest is None:
        SUMMARY["finngen_ok"] = False
        say("[FINNGEN] MANIFEST_NOT_REACHABLE (download form may be required; see finngen.gitbook.io/data-download)")
        return
    lines = manifest.splitlines()
    header = lines[0].split("\t")
    say("[FINNGEN_RELEASE] " + used[0])
    say("[FINNGEN_COLUMNS] " + " | ".join(header))
    say("[FINNGEN_ENDPOINTS_TOTAL] " + str(len(lines) - 1))
    hidx = {h.lower(): i for i, h in enumerate(header)}

    def col(row, *names):
        for n in names:
            for k, i in hidx.items():
                if n in k and i < len(row):
                    return row[i]
        return ""

    hits = []
    for line in lines[1:]:
        row = line.split("\t")
        code = row[0]
        text = line.upper()
        if any(k in text for k in FINNGEN_KEYWORDS):
            hits.append({
                "phenocode": code,
                "name": col(row, "phenotype", "name"),
                "n_cases": col(row, "num_cases", "cases"),
                "n_controls": col(row, "num_controls", "controls"),
                "url": next((c for c in row if c.startswith("https://") and c.endswith(".gz")), ""),
            })
    with open(OUT / "finngen_endpoints.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["phenocode", "name", "n_cases", "n_controls", "url"])
        w.writeheader()
        for h in hits:
            w.writerow(h)
    say("[FINNGEN_MATCHED] " + str(len(hits)) + " endpoints for keywords " + ",".join(FINNGEN_KEYWORDS))
    for h in hits:
        say("[FG] %-45s cases=%-7s controls=%-8s %s" % (h["phenocode"], h["n_cases"], h["n_controls"], h["name"][:60]))
    SUMMARY["finngen_endpoints"] = hits
    SUMMARY["finngen_release"] = used[0]
    # Header and size of a few key files
    by_code = {h["phenocode"]: h for h in hits}
    for code in FINNGEN_PROBE_HEADERS:
        h = by_code.get(code)
        if not h or not h["url"]:
            say("[FG_FILE] %-40s NOT_IN_MANIFEST_OR_NO_URL" % code)
            continue
        size, hst = curl_head_size(h["url"])
        say("[FG_FILE] %-40s %10s %s" % (code, human(size), hst))
        raw, st = curl_bytes(h["url"], 131_072)
        if raw:
            txt = gunzip_prefix(raw).decode("utf-8", "ignore").splitlines()
            if txt:
                say("    [HEADER] " + txt[0][:400])
                if len(txt) > 1:
                    say("    [ROW] " + txt[1][:400])
        record("FinnGen_" + code, "gwas_file", h["url"], hst, human(size))
    SUMMARY["finngen_ok"] = len(hits) > 0


# ----------------------------------------------------------------------------
# 4. GEO series: supplementary files and sample characteristics
# ----------------------------------------------------------------------------
def geo_base(gse):
    num = gse[3:]
    stub = "GSE" + (num[:-3] if len(num) > 3 else "") + "nnn"
    return "https://ftp.ncbi.nlm.nih.gov/geo/series/" + stub + "/" + gse + "/"


def summarize_series_matrix(gse, url, fh):
    raw, st = curl_bytes(url, RANGE_BYTES)
    if raw is None:
        say("    [MATRIX_READ] FAIL " + st)
        return None
    text = gunzip_prefix(raw).decode("utf-8", "ignore")
    head = text.split("!series_matrix_table_begin")[0]
    n_samples = None
    fields = {}
    for line in head.splitlines():
        if line.startswith("!Sample_geo_accession"):
            n_samples = len(line.split("\t")) - 1
        if line.startswith(("!Sample_characteristics_ch1", "!Sample_source_name_ch1",
                            "!Sample_molecule_ch1", "!Sample_library_strategy",
                            "!Sample_platform_id", "!Sample_description")):
            parts = [p.strip().strip('"') for p in line.split("\t")[1:]]
            tag = line.split("\t")[0]
            for p in parts:
                if ":" in p and tag == "!Sample_characteristics_ch1":
                    k, v = p.split(":", 1)
                    key = k.strip()
                    val = v.strip()
                else:
                    key = tag.replace("!Sample_", "")
                    val = p
                fields.setdefault(key, {})
                fields[key][val] = fields[key].get(val, 0) + 1
    say("    [SAMPLES] " + str(n_samples))
    fh.write("### " + gse + " " + url + "\n")
    fh.write("samples=" + str(n_samples) + "\n")
    for key, vals in fields.items():
        top = sorted(vals.items(), key=lambda kv: -kv[1])
        fh.write("  " + key + " (" + str(len(vals)) + " distinct)\n")
        for v, c in top[:60]:
            fh.write("     " + str(c).rjust(5) + "  " + v[:160] + "\n")
        shown = "; ".join("%s=%d" % (v[:40], c) for v, c in top[:8])
        say("    [FIELD] %s (%d distinct): %s" % (key[:40], len(vals), shown))
    return n_samples


def check_geo():
    say("")
    say("-" * 72)
    say("[SECTION] GEO SERIES")
    fh = open(OUT / "geo_characteristics.txt", "w", encoding="utf-8")
    geo_status = {}
    for gse, why in GEO_SERIES.items():
        say("")
        say("[GEO] " + gse + " | " + why)
        base = geo_base(gse)
        info = {"suppl": [], "matrix": []}
        for sub in ["suppl/", "matrix/"]:
            data, st = curl_bytes(base + sub, 300_000)
            if data is None:
                say("  [%s] LISTING_FAILED %s" % (sub.rstrip("/").upper(), st))
                record(gse + "_" + sub.rstrip("/"), "geo_listing", base + sub, "FAIL", st)
                continue
            entries = [(n, u) for n, u in parse_dir_listing(data.decode("utf-8", "ignore"), base + sub)
                       if not n.endswith("/")]
            for name, url in entries:
                size, hst = curl_head_size(url)
                say("  [%s] %-70s %10s" % (sub.rstrip("/").upper(), name[:70], human(size)))
                record(gse + "_" + sub.rstrip("/"), "geo_file", url, hst, human(size))
                info[sub.rstrip("/")].append({"name": name, "url": url, "bytes": size})
        for m in info["matrix"]:
            if m["name"].endswith("_series_matrix.txt.gz"):
                say("  [MATRIX_HEADER] " + m["name"])
                summarize_series_matrix(gse, m["url"], fh)
        geo_status[gse] = info
    fh.close()
    SUMMARY["geo"] = geo_status


# ----------------------------------------------------------------------------
# 5. Zenodo records and web resources
# ----------------------------------------------------------------------------
def check_resources():
    say("")
    say("-" * 72)
    say("[SECTION] MODEL AND REFERENCE RESOURCES")
    zen = {}
    for rec, why in ZENODO_RECORDS.items():
        url = "https://zenodo.org/api/records/" + rec
        data, st = curl_bytes(url)
        if data is None:
            say("[ZENODO] %s FAIL %s | %s" % (rec, st, why))
            record("zenodo_" + rec, "zenodo", url, "FAIL", st)
            zen[rec] = None
            continue
        try:
            js = json.loads(data.decode("utf-8", "ignore"))
        except Exception:
            say("[ZENODO] %s BAD_JSON | %s" % (rec, why))
            zen[rec] = None
            continue
        title = js.get("metadata", {}).get("title", "")
        say("[ZENODO] %s OK | %s | title: %s" % (rec, why, title[:90]))
        files = js.get("files", [])
        zen[rec] = []
        for f in files:
            key = f.get("key") or f.get("filename")
            size = f.get("size") or f.get("filesize")
            say("    [ZFILE] %-60s %s" % (str(key)[:60], human(size)))
            zen[rec].append({"key": key, "bytes": size})
        record("zenodo_" + rec, "zenodo", url, "OK", title[:80])
    SUMMARY["zenodo"] = zen
    for name, url in WEB_PAGES.items():
        size, hst = curl_head_size(url)
        say("[WEB] %-28s %s %s" % (name, hst, url))
        record(name, "web", url, hst, "")


# ----------------------------------------------------------------------------
# 6. Local assets from Paper 7 / Paper 8 (read-only discovery, nothing copied)
# ----------------------------------------------------------------------------
def check_local_assets():
    say("")
    say("-" * 72)
    say("[SECTION] LOCAL ASSETS FROM OTHER PAPERS (read-only)")
    roots = []
    if PAPER8_ROOT.exists():
        roots.append(PAPER8_ROOT)
    for p in sorted(DESKTOP.iterdir()) if DESKTOP.exists() else []:
        if p.is_dir() and re.search(r"paper\s*_?7", p.name, re.I):
            roots.append(p)
    say("[LOCAL_ROOTS] " + (" | ".join(str(r) for r in roots) if roots else "NONE FOUND"))
    pat = re.compile(r"(GSE\d{5,6}|v232|pseudobulk|sparse|series_matrix)", re.I)
    skip = {".venv", "venv", "__pycache__", ".git", "node_modules", "site-packages"}
    found = []
    for root in roots:
        base_depth = len(root.parts)
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = [d for d in dirnames if d not in skip]
            if len(Path(dirpath).parts) - base_depth > 6:
                dirnames[:] = []
                continue
            for d in dirnames:
                if pat.search(d):
                    found.append({"path": str(Path(dirpath) / d), "type": "dir", "bytes": None})
            for f in filenames:
                if pat.search(f):
                    fp = Path(dirpath) / f
                    try:
                        sz = fp.stat().st_size
                    except Exception:
                        sz = None
                    found.append({"path": str(fp), "type": "file", "bytes": sz})
    with open(OUT / "local_assets.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["path", "type", "bytes"])
        w.writeheader()
        for r in found:
            w.writerow(r)
    say("[LOCAL_MATCHES] " + str(len(found)) + " (full list: local_assets.csv)")
    # Show the most relevant: GSE-named dirs and large files
    key_ids = ["GSE266265", "GSE51981", "GSE246591", "GSE141549", "GSE213216", "GSE6364"]
    for gid in key_ids:
        rows = [r for r in found if gid.lower() in r["path"].lower()]
        big = sorted([r for r in rows if r["bytes"]], key=lambda r: -r["bytes"])[:4]
        dirs = [r for r in rows if r["type"] == "dir"][:3]
        say("[LOCAL] %-10s matches=%d" % (gid, len(rows)))
        for r in dirs:
            say("    [DIR]  " + r["path"])
        for r in big:
            say("    [FILE] %-10s %s" % (human(r["bytes"]), r["path"]))
    v232 = [r for r in found if "v232" in r["path"].lower() and r["type"] == "dir"]
    for r in v232[:6]:
        say("[LOCAL_P8_V232] " + r["path"])
    SUMMARY["local_assets_count"] = len(found)
    SUMMARY["local_v232_dirs"] = [r["path"] for r in v232]
    SUMMARY["local_by_id"] = {gid: len([r for r in found if gid.lower() in r["path"].lower()]) for gid in key_ids}


# ----------------------------------------------------------------------------
# 7. Gates per hypothesis
# ----------------------------------------------------------------------------
def geo_has(gse, pattern):
    info = SUMMARY.get("geo", {}).get(gse, {})
    return any(re.search(pattern, f["name"], re.I) for f in info.get("suppl", []) + info.get("matrix", []))


def gates():
    say("")
    say("-" * 72)
    say("[SECTION] FEASIBILITY GATES")
    tools = SUMMARY.get("tools", {})
    gwas_ok = SUMMARY.get("gwas_catalog_ok", False)
    fg_ok = SUMMARY.get("finngen_ok", False)
    codes = [h["phenocode"] for h in SUMMARY.get("finngen_endpoints", [])]
    subtype = [c for c in codes if re.search(r"OVARY|PELVICPERITONEUM|DEEP|RECT|INTESTINE", c)]
    fibroid = [c for c in codes if "LEIOMYOMA" in c]
    uterus = [c for c in codes if c.endswith("_UTERUS") or "ADENOMYO" in c]
    zen = SUMMARY.get("zenodo", {})
    g = {}

    def gate(name, ok, partial, notes):
        status = "READY_DATA" if ok else ("PARTIAL" if partial else "BLOCKED")
        g[name] = {"status": status, "notes": notes}
        say("[GATE] %-4s %-11s %s" % (name, status, notes))

    gate("G1", gwas_ok and geo_has("GSE213216", r"(count|matrix|h5|rds|tar|mtx)"),
         gwas_ok or fg_ok,
         "GWAS=%s GSE213216_matrix=%s" % (gwas_ok, geo_has("GSE213216", r"(count|matrix|h5|rds|tar|mtx)")))
    gate("G2", (gwas_ok or fg_ok) and bool(zen.get("3518299")) and geo_has("GSE141549", "series_matrix"),
         bool(zen.get("3518299")) or geo_has("GSE141549", "series_matrix"),
         "PredictDB=%s GSE141549_matrix=%s" % (bool(zen.get("3518299")), geo_has("GSE141549", "series_matrix")))
    gate("G3", fg_ok and len(subtype) >= 2, fg_ok, "subtype_endpoints=" + (",".join(subtype) or "NONE"))
    gate("G4", (gwas_ok or fg_ok) and (bool(SUMMARY.get("local_v232_dirs")) or geo_has("GSE266265", r"rds|h5|mtx")),
         gwas_ok or fg_ok,
         "P8_v232_export=%s GSE266265_remote=%s" % (bool(SUMMARY.get("local_v232_dirs")), geo_has("GSE266265", r"rds|h5|mtx")))
    gate("G5", fg_ok and bool(fibroid), fg_ok,
         "fibroid=%s uterus/adenomyosis=%s" % (",".join(fibroid) or "NONE", ",".join(uterus) or "NONE"))
    missing = [t for t in ["magma", "ldsc.py", "munge_sumstats.py", "plink"] if not tools.get(t)]
    say("[TOOLS_TO_INSTALL_NEXT] " + (", ".join(missing) if missing else "NONE"))
    SUMMARY["gates"] = g
    SUMMARY["tools_missing"] = missing


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    if not shutil.which("curl"):
        say("[FATAL] curl not found")
        sys.exit(2)
    sections = [check_environment, check_gwas_catalog, check_finngen, check_geo,
                check_resources, check_local_assets, gates]
    for fn in sections:
        try:
            fn()
        except Exception as e:
            say("[SECTION_ERROR] %s %s: %s" % (fn.__name__, type(e).__name__, str(e)[:200]))
            SUMMARY.setdefault("errors", []).append({"section": fn.__name__, "error": repr(e)[:300]})
    with open(OUT / "sources.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["source_id", "kind", "url", "status", "detail"])
        w.writeheader()
        for r in SOURCE_ROWS:
            w.writerow(r)
    with open(OUT / "p9_v001_preflight_summary.json", "w", encoding="utf-8") as fh:
        json.dump(SUMMARY, fh, indent=2, default=str)
    say("")
    say("[OUTPUT] " + str(OUT))
    say("[LOG] " + str(LOG))
    say("[POLICY] No analysis run; no disease effects computed; only metadata and file headers read")
    say("[DONE] PAPER9_V001_PREFLIGHT_COMPLETE")
    _log_fh.close()


if __name__ == "__main__":
    main()
