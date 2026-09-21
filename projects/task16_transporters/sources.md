# Task 16 source and assumption register

Sources checked on 2026-09-21.

| Source | Supported use | Limitation |
|---|---|---|
| [FDA/ICH M12, final August 2024](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/m12-drug-interaction-studies) | Transporter-mediated DDI belongs in drug-development assessment | Does not validate this compound or its parameter values |
| [FDA substrate/inhibitor tables](https://www.fda.gov/drugs/drug-interactions-labeling/drug-development-and-drug-interactions-table-substrates-inhibitors-and-inducers) | OATP1B1/1B3, P-gp and BCRP identities and assay context | No specific clinical inhibitor is represented by imposed I/Ki |
| [Existing Task 4 driver](../task04_pbpk/run_task4_pbpk_pharmacokinetics_dose_prediction.py) | Actual imported circulation, partition factors, amounts and glomerular filtration | Its physiology and partition approximation remain uncalibrated scenario inputs |

Every value in [config.json](outputs/config.json) is an assumption unless explicitly inherited from the Task 4 scenario. There are no measured transporter kinetic inputs.

| Parameter family | Values with units | Status |
|---|---|---|
| OATP1B1/1B3 Vmax | 12 / 4 mg h^-1 | Assumed |
| OATP Km; competitive inhibition | 0.1 mg L^-1; I/Ki 0–100 | Assumed |
| Hepatic passive exchange | 2 L h^-1 | Assumed bidirectional |
| Hepatocyte / renal-cell / enterocyte volumes | 1.3 / 0.15 / 0.25 L | Assumed effective cellular pools |
| Intracellular unbound fraction | 0.2 | Assumed |
| Hepatic intracellular metabolic clearance | 5 L h^-1 | Assumed; replaces original hepatic loss |
| Biliary efflux Vmax / Km | 3 mg h^-1 / 0.1 mg L^-1 | Assumed aggregate carrier |
| Renal uptake Vmax / Km / passive exchange | 3 mg h^-1 / 0.2 mg L^-1 / 0.5 L h^-1 | Assumed; not assigned to hepatic OATP |
| Renal P-gp / BCRP Vmax; Km | 0.8 / 0.6 mg h^-1; 0.15 mg L^-1 | Assumed |
| Gut P-gp / BCRP Vmax; Km | 2 / 1 mg h^-1; 0.2 mg L^-1 | Assumed |
| Enterocyte basolateral transfer / fecal transit | 0.8 / 0.1 h^-1 | Assumed |
| Dose, observation window | 100 mg IV or oral; 96 h | Scenario design, not clinical prescription |

Task 4 organ amounts act as effective perfused exchange spaces; adding cell volumes does not establish a validated human organ-volume decomposition. The finite-time AUC cannot be called total exposure when much drug remains at 96 h.
