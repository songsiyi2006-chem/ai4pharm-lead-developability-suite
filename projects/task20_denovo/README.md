# Task 20 · Three-dimensional pocket-conditioned molecular generation

The local CPU run evolves **100 candidates over 20 generations** against the real unliganded SARS-CoV-2 Mpro pocket in **PDB 6YB7**. Fragment attachment constructs molecular graphs, ETKDG/MMFF constructs conformers, and a dimensionless steric/contact score ranks actual three-dimensional poses. NSGA-II survival balances geometry, a CNS-MPO proxy, QED and RDKit SA score. An equal-budget random-search control tests whether the evolutionary run offers an advantage in this restricted space.

AutoDock Vina separately evaluates 12 generated candidates. Two independently seeded X77 redocking controls use the experimental 6W63 complex, preserve its bond graph and compare heavy-atom coordinates without alignment. The Vina score is **not converted to Kd**. A single rigid monomer, one protonation assignment and a finite fragment vocabulary do not establish biological potency, novelty or synthetic feasibility.

- [中文技术报告](REPORT_ZH.md) · [English technical report](REPORT_EN.md)
- [Source provenance](sources.md) · [Executed summary](outputs/summary.json)
- [Population history](outputs/population_history.csv) · [3D final structures](outputs/final_candidates.sdf)
- [Equal-budget random control](outputs/random_search_summary.json)
- [Independent Vina scores](outputs/docking/candidate_vina_scores.csv)
- [300 DPI figure](outputs/figures/fig20_denovo_pareto_lead_optimization.png)

```bash
# Reuse an existing NumPy/SciPy/RDKit (with Contrib SA_Score) environment.
python projects/task20_denovo/driver.py --out work/task20_new

# Actual docking additionally uses a Vina executable and a Meeko Python runtime.
python projects/task20_denovo/driver.py --out work/task20_docking \
  --vina /path/to/vina --preparation-python /path/to/meeko/python
```

No installation or download occurs during execution. Inputs are archived. Default seed: 20260921; numerical threads: one. The output directory must not already exist. The run without `--vina` explicitly reports docking as not executed. A JSON `--config` can reduce dimensions for software smoke checks, which must not be presented as the full study.
