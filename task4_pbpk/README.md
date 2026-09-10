<!-- TASK4:BEGIN -->
# Developability Suite — Task 4

本目录是任务四的独立模块；默认结果均使用明确标注的示例输入。通过 JSON 模板接入指定先导化合物的参数，不自动把上游基准面板当作实测 PBPK 数据。

| Module | Input / output |
|---|---|
| 4A: IVIVE | Microsomal clearance and binding to whole-liver clearance. |
| 4B: Distribution | Seven amount states with a GI depot and serial pulmonary circulation. |
| 4C: Exposure | IV/oral profiles, AUC to infinity, multiple doses and periodic steady state. |
| 4D: Dose screening | QD/BID free-plasma coverage subject to a total-plasma peak constraint. |

## Run

Requires Python 3.10+, NumPy, SciPy and Matplotlib.

```sh
python -m pip install -r requirements.txt
python run_task4_pbpk_pharmacokinetics_dose_prediction.py --self-test --out .
python run_task4_pbpk_pharmacokinetics_dose_prediction.py --write-example compound.json
python run_task4_pbpk_pharmacokinetics_dose_prediction.py --config compound.json --out results
```

The JSON schema is the exported template: unknown keys and invalid values are rejected. Partial configurations override explicit defaults; omitted fields remain assumptions and are recorded in the resolved configuration. Set `km_unbound_mg_l` to a positive value for optional saturable hepatic metabolism. Nonlinear runs can take substantially longer and may have no finite steady state. Dose bounds are **per administration**; the default lattice is 1 mg. `target_tissue` defaults to `rest`, an explicit proxy. The toxicity threshold is **total plasma mg/L**; the efficacy threshold is **free plasma nM**.

## Deliverables

- [Self-contained script](run_task4_pbpk_pharmacokinetics_dose_prediction.py)
- [中文报告](PBPK_DOSE_PREDICTION_REPORT_ZH.md) / [English report](PBPK_DOSE_PREDICTION_REPORT_EN.md)
- [Resolved inputs](parameters_used.json) / [Editable example](compound_example.json)
- [PK and dose results](results_summary.json) / [Numerical checks](verification.json)
- `single_iv.csv`, `single_oral.csv`, `multidose_bid_7days.csv`, `steady_state_bid.csv`
- `dose_titration.csv`, `sensitivity_analysis.csv`, `manifest.json`
- `figures_task4/`: four requested publication-format PNGs at 300 DPI

The original brief mixes blood flow with plasma concentration and treats lung as a parallel organ. This implementation uses Rb consistently and places lung in series. The seven states include a GI **depot**, not a separate gut-wall tissue. Kp is a documented Poulin–Theil composition approximation with a logD ionization adaptation, not full Rodgers–Rowland. AUC includes its infinite tail; day 7 and true steady state are distinct. Clinical recommendation remains unset because the required validation and toxicology evidence are absent.

## Git delivery

To stage, commit and push generated outputs in an existing repository on local `main` with `origin` configured, run the script with `--git-sync`. It stages only this output directory and the script, refuses pre-existing staged changes, uses the requested commit message, and pushes `origin main` without force. Failures are saved in `git_delivery_status.json` and return a nonzero exit status. The ordinary run has no Git side effects.
<!-- TASK4:END -->
