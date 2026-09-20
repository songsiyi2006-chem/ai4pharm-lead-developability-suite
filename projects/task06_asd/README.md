# Task 6 — ASD thermodynamics and supersaturation

[返回仓库首页](../../README.md) · [项目列表](../README.md)

**Evidence level: executed hypothetical mechanistic scenarios; no compound-specific experimental fit or clinical prediction.**

[TASK.md](TASK.md) preserves the original user request. Its proposed equations, targets and expected outcomes are historical requirements; the reports below document the implemented corrections and actual results.

[中文技术报告](ASD_FORMULATION_KINETICS_REPORT_ZH.md) · [English report](ASD_FORMULATION_KINETICS_REPORT_EN.md) · [Computed summary](results/summary.json) · [Assumptions and units](inputs/parameters.json) · [Primary-source registry](inputs/sources.json)

The standalone driver implements four-carrier Flory–Huggins binodal/spinodal calculations, dry/humid Gordon–Taylor scenarios, an explicitly uncalibrated VFT domain check, and finite-dose dissolution–precipitation–absorption ODEs. All four requested figures are in [figures_task6](figures_task6/).

From the repository root, in an environment with NumPy, SciPy and Matplotlib:

```bash
python projects/task06_asd/run_task6_asd_formulation_supersaturation_kinetics.py --self-test
python projects/task06_asd/run_task6_asd_formulation_supersaturation_kinetics.py --out work/task6_reproduction
python projects/task06_asd/run_task6_asd_formulation_supersaturation_kinetics.py --write-example work/task6_parameters.json
python projects/task06_asd/run_task6_asd_formulation_supersaturation_kinetics.py --config work/task6_parameters.json --out work/task6_custom
```

Use a **new** output directory. Computation is offline and deterministic; there is no stochastic fitting. Inputs are embedded in `DEFAULT` and exported verbatim for provenance. The checked-in `inputs/`, `results/`, `figures_task6/` and `manifest.json` are one completed run. Re-running writes the same structure beneath `--out`.

`--config` reads the complete exported schema and rejects nonphysical values or ignored design changes. The four carriers, drug reference segment Nd=1, prescribed loading grid and endpoint thresholds define this experiment. Users may edit the physical parameters, horizon and two storage conditions. Such changes define new hypothetical calculations unless supplied with new calibration evidence. The script's summary is regenerated; the checked-in narrative reports describe the default run only.

The run contains 324 phase-grid rows, 32 storage scenarios, 69 ODE integrations and four 300-dpi figures. The 20% HPMC-AS candidate gives a hypothetical closed-vessel AUC ratio of 8.00 and an absorption-sink proxy ratio of 8.80. The 30% candidate has a larger closed-vessel AUC under these assumptions, so **20% is not established as optimal**. Twelve storage rows are outside the chosen mathematical VFT domain; none of the storage calculations establishes shelf life.

The 900 mL, pH 6.8 vessel is an assumed medium. Its bile-salt composition, micellar binding and permeability are not calibrated, so it is not a validated FaSSIF test or a human exposure model.
