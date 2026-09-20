# Task 10 — RNA ensembles and splicing thermodynamics

This project executes real ViennaRNA partition functions and analyzes 40 public NMR conformers, then explores four **hypothetical** ligand models with a thermodynamically consistent finite-pool RNA–ligand–U1 equilibrium. Binding energies and inclusion proxies are not clinical risdiplam predictions.

| Start here | Contents |
|---|---|
| [Chinese report](RNA_TARGETED_CADD_REPORT_ZH.md) / [English report](RNA_TARGETED_CADD_REPORT_EN.md) | Equations, results, limitations, source-backed interpretation |
| [Driver](run_task10_rna_targeted_small_molecule_dynamics.py) | All four calculation and plotting modules in one Python artifact |
| [Original request](TASK.md) | Preserved original bytes; aspirational claims are not evidence |
| [Sources](SOURCES.md) / [example configuration](inputs/example_config.json) | Public structures and explicitly assumed parameters |
| [Results summary](results/results_summary.json) / [verification](results/verification.json) | Actual executed results and numerical checks |
| [Figures](figures_task10/) | Exactly four 300-DPI PNGs |
| [Manifest](manifest.json) | SHA256 hashes of project files |

## Reproduce from repository root

Use Python with the repository requirements, including `ViennaRNA==2.7.2`. No network access is required during computation because both PDB files are archived.

```powershell
python projects/task10_rna_splicing/run_task10_rna_targeted_small_molecule_dynamics.py --out work/task10_reproduce
python projects/task10_rna_splicing/run_task10_rna_targeted_small_molecule_dynamics.py --write-example work/task10_config.json
python projects/task10_rna_splicing/run_task10_rna_targeted_small_molecule_dynamics.py --config work/task10_config.json --out work/task10_custom
python projects/task10_rna_splicing/run_task10_rna_targeted_small_molecule_dynamics.py --self-test
python -m unittest discover -s tests -p test_task10.py -v
```

This Windows run reused the chemistry Python 3.12 environment and a separately installed ViennaRNA wheel. On that machine add `--rna-library work/optional_rna_py312`; the library may also be provided through `PYTHONPATH`. This local directory is ignored and is not a portable dependency location. `--write-example` does not need ViennaRNA.

Output directories must be new or empty. Existing outputs require explicit `--overwrite`; normal reproduction should use a new directory. A fresh run puts figures and a manifest directly under `--out`, with numerical tables alongside them. In the curated repository, tables are under `results/` and figures under `figures_task10/`; the project manifest records this arrangement. Reports describe the frozen default run and are not automatically reinterpreted when parameters change.

## Evidence and interfaces

- Secondary structure: a 26-nt **engineered tethered hairpin**, not full SMN2 or a native two-strand spliceosome. Pseudouridine is approximated as U here, but retained in 3D.
- Geometry: 6HMI apo and 6HMO SMN-C5 bound; 20 NMR conformers each. Ensemble spread is not a thermodynamic population or an MD trajectory. SMN-C5 is a risdiplam analogue, not risdiplam.
- Eight equilibrium states per RNA: closed/open, each with no ligand, ligand, U1, or both. Free ligand and U1 are jointly solved with mass conservation in target and one off-target pool.
- U1 occupancy is an **inclusion proxy assumption**. There is no clinical PK, global transcriptome screen, validated dosing window, or fitted efficacy.
- Module 10A, experimental geometry, and ligand-specific competent substates are complementary analyses. Their energies are not claimed to be inferred from one another.

The default calculation performs 14 partition-function evaluations, produces 560 groove-facing descriptor rows, solves 976 main dose-response conditions and 120 sensitivity conditions, and retains negative results: no position exceeds 1.2 bits at 37 C; only one hypothetical scaffold exceeds 85% asymptotic inclusion proxy. Fourteen dedicated tests cover physical limits and regressions.

[Back to projects](../README.md) · [Repository home](../../README.md)
