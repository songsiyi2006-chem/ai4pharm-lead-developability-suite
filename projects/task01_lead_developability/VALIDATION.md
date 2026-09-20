# Delivery validation

## Automated tests

The 16 tests in [tests/test_task1.py](../../tests/test_task1.py) passed locally on both environments during the original scientific delivery:

| Environment | Scientific stack | Outcome |
|---|---|---|
| Windows, CPython 3.14.6 | RDKit 2026.03.5; NumPy 2.5.2; pandas 3.0.5; matplotlib 3.11.1; seaborn 0.13.2 | 16 passed |
| Windows, CPython 3.10.21 | RDKit 2026.03.6; NumPy 2.2.6; pandas 2.3.3; matplotlib 3.10.9; seaborn 0.13.2 | 16 passed |

```bash
python -m unittest discover -s tests -p test_task1.py -v
```

Coverage includes score endpoints, interior interpolation and continuity, custom/original MPO separation, ESOL coefficients and unit conversion, chemical identity, malformed/duplicate inputs, ionization direction, amide exclusion, override validation, label independence, reproducible small-molecule geometry, isolated-worker timeout, RDKit's negative timeout sentinel, invalid CLI options, a full offline descriptor-only run, all three PNGs at 300 DPI, and artifact hashes.

These tests do not validate biological predictivity. GitHub Actions is configured for Python 3.10 and 3.12; local pass results above do not imply that a remote workflow has already passed.

## Delivered full panel

The complete per-compound results, precise execution settings, geometry status counts and checksums are in [run_manifest.json](results_task1/run_manifest.json). The manifest is the authoritative record of the delivered run. Reproduce the extended macrocycle budget with:

```bash
python projects/task01_lead_developability/run_task1_mpo_admet_developability.py --conformer-timeout 120 --strict-3d --output-dir work/task01_reproduction
```

The delivered run completed for all 30 compounds: 26 ensembles fully converged and four partially converged (cyclosporine A, paclitaxel, sirolimus and daclatasvir). Every compound has at least one converged conformer; no volume fallback was needed. Cyclosporine A contributed three converged conformers. All artifact hashes and local document links passed the final integrity check. PNG dimensions are 3000×1980, 4500×1680 and 3300×2310 pixels, each with 300 DPI metadata.

The default 45-second worker budget is intentionally bounded and may yield missing cyclosporine geometry on this machine. A longer budget gives difficult cyclic peptides more embedding time. `--strict-3d` requires at least one converged conformer per compound; partial convergence is reported and is not equivalent to complete sampling. The retry permits cis amides and does not establish equilibrium cis/trans populations.

The three figures were visually inspected for label clipping, visible boxplot medians, units, legends and interpretability. A final integrity check verifies the script checksum, all result and image hashes, the 30-record count and local document links. JSON retains full numeric precision; CSV values are rounded to eight significant digits.

## Known scientific limitations

- The historical panel is purposive and small; no independent assay test set or clinical classifier is provided.
- pKa/logD assumptions and the hERG, Caco-2, HIA and oral-score formulas are uncalibrated.
- Successful geometry does not prove solvent exposure or chameleonicity; `absolute_chameleonic_hb_index` remains null.
- ARV-110's selected source structure has unspecified stereochemistry; its 3D realization is not a verified clinical stereoisomer ensemble.
- RDKit/version/platform differences and time limits can change geometry and dependent proxy scores. A random seed is not a cross-platform guarantee.
- No finite collection of tests can guarantee absence of every software defect.
