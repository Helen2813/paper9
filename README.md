# Paper 9 — Where does inherited endometriosis risk act?

Genetic localisation of endometriosis risk across tissue context (eutopic endometrium,
lesion niche, outside the uterus), menstrual-cycle cell states and disease subtypes.

Plan: `docs/plan/PLAN_v0.1.md` · Hypothesis registry: `docs/registry/P9_HYPOTHESES.csv`

## Layout

```
scripts/        versioned, self-contained run scripts (p9_vNNN_*.py)
docs/plan/      research plan and change log (frozen versions)
docs/registry/  hypothesis registry
data/raw/       downloaded source files (git-ignored)
data/external/  references and links to other projects' data (git-ignored)
data/interim/   intermediate tables (git-ignored)
data/processed/ analysis-ready tables (git-ignored)
results/        per-run outputs; only *_summary.json is committed
logs/           console logs per run (git-ignored)
figures/        final figures
manuscript/     LaTeX sources
tools/          MAGMA, LDSC, reference panels (git-ignored)
```

## Running

Each script is run without arguments from Terminal, e.g.

```
python3 ~/Desktop/paper9/scripts/p9_v001_preflight.py
```

Every run writes `logs/<script>_<timestamp>.log` and `results/<step>/<script>_<timestamp>/`.

## Related projects

- Paper 7 (preprint): identifiability of delay vs attenuation in decidualization.
- Paper 8: specificity and transportability of endometriosis markers.

Datasets shared with those projects are referenced by path, not copied.
