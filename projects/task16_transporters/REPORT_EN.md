# Task 16 — Transporter-limited cellular access in the existing PBPK circulation

This extension imports and calls the existing Task 4 `PBPK.derivative`. It preserves that model's series lung/systemic circulation and filtration while adding finite hepatocyte, renal-cell and enterocyte amounts. [Sources](sources.md), [summary](outputs/summary.json), [verification](outputs/verification.json), and [figure](outputs/figures/fig16_transporter_oatp_pgp_kinetics.png) document the calculation. All transporter kinetics are assumptions.

## Balance derivation and coupling

For each ordinary perfused tissue, Task 4 supplies Q(Carterial−Cvenous), with Cvenous determined by its existing partition approximation. To avoid two parallel hepatic elimination mechanisms, the extension calculates the exact old hepatic flux H=CLint·Cu,liver and reverses its two bookkeeping contributions: add H back to the liver derivative and subtract H from the hepatic-loss ledger. Only then is intracellular uptake/metabolism introduced.

For unbound donor concentration Cu, competitive uptake and passive exchange are:

$$J_{in}=\sum_{j\in\{1B1,1B3\}}\frac{V_{max,j}C_{u,ext}}{K_{m,j}(1+I/K_i)+C_{u,ext}}+PS(C_{u,ext}-C_{u,cell}).$$

The passive term changes sign; it is not an irreversible additional sink. Cu,cell=fu,cell·Acell/Vcell. Hepatocyte conservation gives:

$$\dot A_H=J_{in}-CL_{met}C_{u,H}-\frac{V_{max,bile}C_{u,H}}{K_{m,bile}+C_{u,H}}.$$

Uptake is subtracted from the donor liver space. Metabolic and biliary outputs enter separate loss ledgers. Renal cells similarly receive finite uptake and passive exchange, and their saturable P-gp/BCRP efflux enters the urinary ledger. Existing glomerular filtration remains exactly once in the Task 4 balance; renal uptake is not mislabeled hepatic OATP.

The old gut-depot→liver absorption shortcut is removed from the liver derivative. Its flux is instead transferred to enterocytes. Enterocyte basolateral export enters the liver, while P-gp/BCRP export returns to the lumen. Luminal transit produces fecal loss. This permits recycling without creating or destroying drug.

Summing all seven existing amount states, three cellular states and three loss ledgers cancels exchange and elimination transfers:

$$\sum A_i+L_{met}+L_{urine}+L_{feces}=Dose.$$

The AUC state is an integral with different units and is excluded from this sum. These balances make the hybrid model auditable despite its uncalibrated physiology.

## Executed scenarios

Eight IV/oral scenarios use I/Ki=0,1,10,100 at a 100 mg hypothetical input. Another 24 dose/inhibition combinations test saturation at 10,100,300 mg. IV AUC0–96 values are **7.266, 10.250, 16.399, 18.687 mg·h/L**; the corresponding AUCR values are **1, 1.411, 2.257, 2.572**. Oral AUCR values are **1, 1.495, 2.533, 2.919**.

The apparently larger exposure is specifically **96-hour truncated AUC**. At IV I/Ki=100, **55.72 mg** remains in modeled amounts at 96 h. Reporting this ratio as AUC0–infinity would be incorrect. The same imposed inhibitor concentration acts on both OATP pathways and has no perpetrator PK model.

Independent checks include summed derivative conservation, nonnegative boundary flux, absence of carrier transport from an empty donor, competitive Km increase without changing limiting Vmax, and tighter-tolerance integration. Eight detailed trajectories contain 7,688 rows.

## Translational boundary

This calculation can help prioritize uptake/efflux and intracellular concentration assays. It needs measured Vmax/Km, intracellular binding, transporter abundance, species scaling and external PK evaluation. Adding effective cell volumes to inherited exchange spaces does not establish anatomically correct organ decomposition. The old partition approximation also remains an approximation. The output is a conditional transport experiment in silico; clinical recommendations are null.
