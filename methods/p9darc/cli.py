"""P9-DARC prototype CLI. Python 3.9+; no optional dependencies."""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from pathlib import Path
from .core import read_pseudobulk, read_gene_z, read_hierarchy, make_report


def _write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def generate_demo(output: Path):
    rng = random.Random(20261111)
    cells = [("eF3", "stroma"), ("dS3", "stroma"), ("VSMC", "perivascular"), ("immune", "immune")]
    donors = ["D{:02d}".format(i) for i in range(1, 11)]
    rows, z_rows = [], []
    for j in range(240):
        g = "TOY_G{:04d}".format(j)
        target = "eF3" if j < 65 else ("VSMC" if j < 120 else ("dS3" if j < 170 else None))
        z = (3.0 if j < 65 else (1.7 if j < 120 else (1.3 if j < 170 else 0.4))) + rng.gauss(0, 0.22)
        z_rows.append({"gene": g, "z": "{:.6f}".format(z)})
        for donor in donors:
            donor_effect = 1 + 0.14 * (int(donor[1:]) - 5.5) / 5.5
            for cell, _ in cells:
                base = 4.0 + (j % 11) * 0.7
                extra = 24.0 if cell == target else 0.0
                cpm = max(0, (base + extra) * donor_effect * (1 + rng.uniform(-0.09, 0.09)))
                rows.append({"gene": g, "donor": donor, "cell_type": cell, "cpm": "{:.6f}".format(cpm), "n_cells": 110 + int(donor[1:]) * 4})
    output.mkdir(parents=True, exist_ok=True)
    _write_tsv(output/"toy_pseudobulk.tsv", ["gene", "donor", "cell_type", "cpm", "n_cells"], rows)
    _write_tsv(output/"toy_gene_z.tsv", ["gene", "z"], z_rows)
    _write_tsv(output/"toy_hierarchy.tsv", ["cell_type", "group"], [{"cell_type": c, "group": group} for c, group in cells])


def main(argv=None):
    parser = argparse.ArgumentParser(description="P9-DARC research-only donor/direction audit (NOT CALIBRATED)")
    cmd = parser.add_subparsers(dest="command", required=True)
    demo = cmd.add_parser("demo", help="Generate and audit SYNTHETIC toy pseudobulk")
    demo.add_argument("--out", default="results/methods/p9darc_v001_demo")
    audit = cmd.add_parser("audit", help="Audit small precomputed donor pseudobulk; NO p values")
    audit.add_argument("--pseudobulk", required=True)
    audit.add_argument("--gene-z", required=True)
    audit.add_argument("--hierarchy")
    audit.add_argument("--out", required=True)
    audit.add_argument("--min-donors", type=int, default=5)
    audit.add_argument("--k", type=int, default=40)
    args = parser.parse_args(argv)
    output = Path(args.out)
    if args.command == "demo":
        # Avoid overwriting existing results; a demo is deterministic and separate.
        if output.exists():
            raise SystemExit("Output already exists; refusing to overwrite: {}".format(output))
        generate_demo(output)
        p, z, h, k, min_d = output/"toy_pseudobulk.tsv", output/"toy_gene_z.tsv", output/"toy_hierarchy.tsv", 40, 5
    else:
        if output.exists():
            raise SystemExit("Output already exists; refusing to overwrite: {}".format(output))
        p, z, h, k, min_d = Path(args.pseudobulk), Path(args.gene_z), Path(args.hierarchy) if args.hierarchy else None, args.k, args.min_donors
    data = read_pseudobulk(p)
    zscores = read_gene_z(z)
    hierarchy = read_hierarchy(h)
    report = make_report(data, zscores, hierarchy, k=k, min_donors=min_d)
    output.mkdir(parents=True, exist_ok=True)
    with (output/"audit_report.json").open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    print("[P9-DARC] PROTOTYPE ONLY / NO GWAS INFERENCE")
    print("[DONORS] {}  [CELLS] {}  [GENES] {}".format(report["n_donors"], report["n_cells"], report["gene_universe_in_pseudobulk"]))
    for family, ranking in report["rankings"].items():
        print("[DIAGNOSTIC] {} top: {}".format(family, ", ".join("{}={:.3f}".format(r["cell_type"], r["mean_gene_z"]) for r in ranking[:3]) or "NONE"))
    print("[STATUS] {} [OUTPUT] {}".format(report["inferential_status"], output.resolve()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
