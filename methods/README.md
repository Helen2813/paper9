# Paper 9 / P9-DARC methods workspace

**Research-only prototype.** Created as an isolated method-development namespace after P9 v071; all previous analyses and logs stay untouched.

## Quick start on the Mac (no package installation)

From `~/Desktop/paper9`:

```bash
PYTHONPATH="$PWD/methods" python3 -m unittest discover -s methods/p9darc/tests -v
PYTHONPATH="$PWD/methods" python3 -m p9darc.cli demo --out results/methods/p9darc_v001_demo
```

The demo generates **synthetic, non-biological data** and saves `toy_pseudobulk.tsv`, `toy_gene_z.tsv`, and `audit_report.json`. It will never report formal GWAS significance, an etiologic cell, or calibrated causal probability.

## Files

- `METHOD_SPEC_v0.1.md` — scientific targets, equations, data constraints, uncertainty and refusal rules.
- `BENCHMARK_PLAN_v0.1.md` — calibration plan and clean external benchmark design.
- `NOVELTY_AND_BOUNDARIES.md` — literature overlap and Paper 7/8 separation.
- `p9darc/core.py` — donor-aware directional profile and diagnostic ranking prototype.
- `p9darc/cli.py` — self-contained toy/demo and small TSV audit interface.
- `p9darc/tests/` — standard-library unit tests.

## Input TSV schema for the audit-only prototype

`pseudobulk.tsv`: `gene`, `donor`, `cell_type`, `cpm`; optional `n_cells` for a descriptive donor-imbalance warning. Each `(donor,gene,cell_type)` key must be unique. The script expects **already-normalized donor-level pseudobulk**, not raw scRNA cell rows. Each file must represent one compatible biological context; no implicit phase/batch/disease adjustment is implemented.

`gene_z.tsv`: `gene`, `z` (for example output MAGMA gene z, for ranking **only**).

`hierarchy.tsv` (optional): `cell_type`, `group`.

Example invocation (NO genetic inference or calibration):

```bash
PYTHONPATH="$PWD/methods" python3 -m p9darc.cli audit \
  --pseudobulk /path/to/pseudobulk.tsv \
  --gene-z /path/to/gene_z.tsv \
  --hierarchy /path/to/hierarchy.tsv \
  --out results/methods/audit_demo --min-donors 6 --k 40
```

The report includes a diagnostic donor effective-count warning when `n_cells` is supplied. The report status is always `NOT_CALIBRATED` until official MAGMA/S-LDSC integration and held-out validation are implemented.
