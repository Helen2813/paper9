"""Research-only P9-DARC profile functions. NO inferential GWAS p-values.

Pure standard library; designed for small donor-pseudobulk TSV proof-of-concept files.
Not a replacement for MAGMA, LDSC, CELLECT or a calibrated method.
"""
from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean
from typing import Dict, Iterable, List, Optional, Tuple

MAX_ROWS = 300_000  # guardrail: prototype is not a large-matrix loader


def _rows(path: Path, fields: Iterable[str]):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(fields) - set(reader.fieldnames or ())
        if missing:
            raise ValueError("Missing TSV columns in {}: {}".format(path, sorted(missing)))
        for i, record in enumerate(reader, 1):
            if i > MAX_ROWS:
                raise ValueError("Too many rows for prototype; use a streaming data pipeline")
            yield record


def read_pseudobulk(path: Path) -> List[dict]:
    """One already-normalized donor x cell type x gene CPM per row."""
    seen = set()
    out = []
    for r in _rows(path, ("gene", "donor", "cell_type", "cpm")):
        gene, donor, ct = (r[k].strip() for k in ("gene", "donor", "cell_type"))
        if not (gene and donor and ct):
            raise ValueError("gene, donor, cell_type must be nonempty")
        key = (gene, donor, ct)
        n_cells = None
        if "n_cells" in r and r["n_cells"].strip():
            n_cells = float(r["n_cells"])
            if not math.isfinite(n_cells) or n_cells <= 0 or n_cells != int(n_cells):
                raise ValueError("n_cells must be a positive integer for {}".format(key))
            n_cells = int(n_cells)
        if key in seen:
            raise ValueError("Duplicate donor-state-gene key: {}".format(key))
        seen.add(key)
        value = float(r["cpm"])
        if not math.isfinite(value) or value < 0:
            raise ValueError("CPM must be finite and nonnegative: {}".format(key))
        out.append({"gene": gene, "donor": donor, "cell_type": ct, "cpm": value, "n_cells": n_cells})
    if not out:
        raise ValueError("Empty pseudobulk file")
    return out


def read_gene_z(path: Path) -> Dict[str, float]:
    genes = {}
    for r in _rows(path, ("gene", "z")):
        gene = r["gene"].strip()
        z = float(r["z"])
        if not gene or not math.isfinite(z) or gene in genes:
            raise ValueError("Missing/duplicate gene or non-finite z in gene_z.tsv")
        genes[gene] = z
    if not genes:
        raise ValueError("Empty gene z-score file")
    return genes


def read_hierarchy(path: Optional[Path]) -> Dict[str, str]:
    if path is None:
        return {}
    out = {}
    for r in _rows(path, ("cell_type", "group")):
        ct, group = r["cell_type"].strip(), r["group"].strip()
        if not ct or not group or ct in out:
            raise ValueError("Missing/duplicate hierarchy mapping for {}".format(ct))
        out[ct] = group
    return out


def donor_paired_profiles(rows: List[dict], min_donors: int = 5) -> Dict[str, Dict[str, dict]]:
    """Compute per-donor contrasts with same-donor comparators.

    No phase, batch, diagnosis or anatomical context adjustment. Users must
    provide a single compatible biological stratum. Profiles with insufficient
    matched donors are withheld.
    """
    if min_donors < 3:
        raise ValueError("min_donors must be >= 3")
    by_donor_gene = defaultdict(dict)
    for r in rows:
        by_donor_gene[(r["donor"], r["gene"])][r["cell_type"]] = r["cpm"]
    diffs = defaultdict(list)
    specs = defaultdict(list)
    for (donor, gene), cs in by_donor_gene.items():
        if len(cs) < 2:
            continue
        expression_logs = {c: math.log2(1.0 + x) for c, x in cs.items()}
        sum_exp = sum(cs.values())
        for c, val in expression_logs.items():
            others = [v for other_c, v in expression_logs.items() if other_c != c]
            if not others:
                continue
            diffs[(c, gene)].append(val - mean(others))
            specs[(c, gene)].append(cs[c] / sum_exp if sum_exp > 0 else 0.0)
    result = defaultdict(dict)
    for (ct, gene), values in diffs.items():
        if len(values) < min_donors:
            continue
        avg = mean(values)
        result[ct][gene] = {"delta": avg, "spec": mean(specs[(ct, gene)]), "n_donors": len(values)}
    return dict(result)


def select_markers(profiles: Dict[str, Dict[str, dict]], k: int = 40) -> Dict[str, Dict[str, List[str]]]:
    if k < 1:
        raise ValueError("k must be >=1")
    out = {}
    for ct, per_gene in sorted(profiles.items()):
        up = sorted(((g, p["delta"]) for g, p in per_gene.items() if p["delta"] > 0), key=lambda x: (-x[1], x[0]))
        down = sorted(((g, p["delta"]) for g, p in per_gene.items() if p["delta"] < 0), key=lambda x: (x[1], x[0]))
        abs_ = sorted(((g, abs(p["delta"])) for g, p in per_gene.items()), key=lambda x: (-x[1], x[0]))
        spec = sorted(((g, p["spec"]) for g, p in per_gene.items()), key=lambda x: (-x[1], x[0]))
        out[ct] = {"UP": [x[0] for x in up[:k]],
                   "DOWN": [x[0] for x in down[:k]],
                   "ABS": [x[0] for x in abs_[:k]],
                   "SPEC": [x[0] for x in spec[:k]]}
    return out


def diagnostic_ranks(markers: Dict[str, Dict[str, List[str]]], gene_z: Dict[str, float],
                     min_genes: int = 10) -> dict:
    """Mean gene z by marker set: TOY DIAGNOSTIC ONLY, no LD correction or p values."""
    result = {}
    families = ("UP", "DOWN", "ABS", "SPEC")
    for family in families:
        ranked = []
        for cell, choices in markers.items():
            present = [g for g in choices[family] if g in gene_z]
            if len(present) < min_genes:
                continue
            ranked.append({"cell_type": cell, "mean_gene_z": round(mean(gene_z[g] for g in present), 6),
                           "n_genes": len(present)})
        result[family] = sorted(ranked, key=lambda x: (-x["mean_gene_z"], x["cell_type"]))
    return result


def leave_one_donor_out(rows: List[dict], gene_z: Dict[str, float], k: int = 40,
                        min_donors: int = 5, hierarchy: Optional[Dict[str, str]] = None) -> dict:
    donors = sorted({r["donor"] for r in rows})
    hierarchy = hierarchy or {}
    by_family = {f: Counter() for f in ("UP", "DOWN", "ABS", "SPEC")}
    by_group = {f: Counter() for f in by_family}
    skipped = Counter()
    for donor in donors:
        loo_rows = [r for r in rows if r["donor"] != donor]
        rankings = diagnostic_ranks(select_markers(donor_paired_profiles(loo_rows, min_donors), k), gene_z)
        for family in by_family:
            if rankings[family]:
                winner = rankings[family][0]["cell_type"]
                by_family[family][winner] += 1
                by_group[family][hierarchy.get(winner, "UNMAPPED")] += 1
            else:
                skipped[family] += 1
    return {"donors_left_out": len(donors),
            "top_cell_counts": {family: dict(by_family[family]) for family in by_family},
            "top_group_counts": {family: dict(by_group[family]) for family in by_group},
            "skipped": dict(skipped),
            "warning": "Selection frequencies are descriptive; not probabilities of a causal cell."}


def donor_balance(rows: List[dict]) -> dict:
    """Descriptive cell sampling balance (not an inferential donor correction).

    n_eff = (sum n_cells)^2 / sum(n_cells^2) across independent donors.
    If no n_cells column is present, no result is invented.
    """
    counts = {}
    for r in rows:
        n = r.get("n_cells")
        if n is None:
            continue
        key = (r["cell_type"], r["donor"])
        if key in counts and counts[key] != n:
            raise ValueError("Different n_cells values for one donor/cell_type: {}".format(key))
        counts[key] = n
    cells = sorted({r["cell_type"] for r in rows})
    result = {}
    for ct in cells:
        vals = [v for (cell, _), v in counts.items() if cell == ct]
        if not vals:
            result[ct] = {"status": "UNKNOWN_NO_CELL_COUNTS"}
            continue
        total = sum(vals)
        eff = total**2 / sum(v**2 for v in vals)
        result[ct] = {"status": "DESCRIPTIVE_ONLY", "n_donors": len(vals),
                      "largest_donor_share": round(max(vals)/total, 5),
                      "effective_donor_count_from_cell_weights": round(eff, 4),
                      "donor_imbalance_warning": bool(eff < 5 or max(vals)/total > 0.5)}
    return result


def make_report(rows: List[dict], gene_z: Dict[str, float], hierarchy: Optional[Dict[str, str]] = None,
                k: int = 40, min_donors: int = 5, lodo: bool = True) -> dict:
    hierarchy = hierarchy or {}
    cells = sorted({r["cell_type"] for r in rows})
    if hierarchy and set(cells) - set(hierarchy):
        raise ValueError("Hierarchy missing cell types: {}".format(sorted(set(cells)-set(hierarchy))))
    profiles = donor_paired_profiles(rows, min_donors)
    markers = select_markers(profiles, k)
    rankings = diagnostic_ranks(markers, gene_z)
    return {
        "version": "P9-DARC 0.1.0-dev", "inferential_status": "NOT_CALIBRATED",
        "has_formal_genetic_inference": False, "has_validated_abstention": False,
        "n_donors": len(set(r["donor"] for r in rows)), "n_cells": len(cells),
        "gene_universe_in_pseudobulk": len(set(r["gene"] for r in rows)),
        "n_gene_z": len(gene_z), "marker_top_k": k, "minimum_paired_donors": min_donors,
        "hierarchy": hierarchy, "rankings": rankings,
        "donor_sampling_balance": donor_balance(rows),
        "eligible_profile_genes": {ct: len(p) for ct, p in profiles.items()},
        "marker_counts": {ct: {fam: len(gs) for fam, gs in fams.items()} for ct, fams in markers.items()},
        "lodo": leave_one_donor_out(rows, gene_z, k, min_donors, hierarchy) if lodo else None,
        "cautions": [
            "Mean gene z ranking has NO LD correction and NO calibrated p-values; diagnostic only.",
            "Donor-paired expression contrasts require matched biological strata and correct normalization.",
            "This prototype has not executed MAGMA, S-LDSC, CELLECT, scDRS or CATCH.",
            "Do not infer an etiologic cell, statistical significance or posterior probability from this report."
        ],
    }
