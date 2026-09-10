# Clinical interpretation and evidence

This is a selected historical reference panel, not a current prescribing catalogue, a complete withdrawal registry, or a collection of assay labels. Structure identity is sourced separately in `reference_panel.json`. Clinical outcome is metadata only and never enters the equations.

## Toxicity stratum

| Compound | Interpretation | Primary evidence / regulator |
|---|---|---|
| Terfenadine | Parent-drug hERG block; acquired QT/arrhythmia risk | [Roy et al., 1996](https://pubmed.ncbi.nlm.nih.gov/8772706/) |
| Astemizole | Parent and desmethyl metabolite block hERG | [Zhou et al., 1999](https://pubmed.ncbi.nlm.nih.gov/10376921/) |
| Cisapride | Potent hERG blocker; distinct from hepatic or thrombotic liabilities | [Comparative patch-clamp study](https://pubmed.ncbi.nlm.nih.gov/11714889/) |
| Cerivastatin | Withdrawal associated with rhabdomyolysis, including interaction risk | [FDA review](https://www.accessdata.fda.gov/drugsatfda_docs/nda/2003/21-366_Crestor_Admindocs_P1.pdf) |
| Troglitazone | Severe hepatic injury; not a hERG-positive label | [FDA medical review](https://www.accessdata.fda.gov/drugsatfda_docs/nda/99/20719S12_Prelay_medr_P2.pdf) |
| Rofecoxib | Cardiovascular thrombotic risk, not equivalent to QT prolongation | [APPROVe trial](https://pubmed.ncbi.nlm.nih.gov/15713943/) |
| Nefazodone | Severe hepatocellular injury; withdrawals/restrictions are jurisdiction-dependent | [Primary case series](https://pubmed.ncbi.nlm.nih.gov/10068386/) |
| Benoxaprofen | Fatal cholestatic injury reported in elderly patients | [Primary clinical report](https://www.bmj.com/content/284/6326/1372) |
| Ximelagatran | Global withdrawal after severe liver injury report | [EMA announcement](https://www.ema.europa.eu/en/documents/press-release/astrazeneca-withdraws-its-application-ximelagatran-36-mg-film-coated-tablets_en.pdf) |
| Fialuridine | Investigational trial stopped after severe hepatic/mitochondrial toxicity | [Primary clinical trial report](https://pubmed.ncbi.nlm.nih.gov/7565947/) |

Low solubility or high cLogP can coexist with these liabilities, but these sources do not establish solubility as the cause of every withdrawal. The ten compounds do not form a uniform hERG-positive test set. A low hERG heuristic cannot clear rofecoxib's thrombotic risk or fialuridine's mitochondrial toxicity.

## bRo5 examples and interpretation

Venetoclax is a nonmacrocyclic BCL-2 inhibitor. Its oral product employs an amorphous solid dispersion; a human formulation study documents the delivery strategy. [Primary formulation/bioavailability study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9338003/).

ARV-110 is bavdegalutamide, represented here by PubChem CID 134414307. The selected record does not specify complete stereochemistry; the database's depositor notes explicitly acknowledge this. Geometry results therefore refer to a sampled stereochemical realization, not a verified clinical stereoisomer ensemble. [PubChem compound](https://pubchem.ncbi.nlm.nih.gov/compound/134414307), [IUPHAR depositor record](https://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?sid=441604884).

Paclitaxel is a taxane used here as a challenging conventional IV comparator, not as proof of oral macrocycle success. Cyclosporine A, tacrolimus, sirolimus, erythromycin, azithromycin and rifampicin provide macrocyclic references; daclatasvir and venetoclax broaden the panel to nonmacrocyclic bRo5 chemistry. Route/modality annotations are contextual; measured oral bioavailability is not assigned to these rows.

## Method references

- [Wager et al. CNS-MPO, 2010](https://doi.org/10.1021/cn100008c).
- [CNS-MPO implementation comparison, 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8260158/).
- [Delaney ESOL, 2004](https://doi.org/10.1021/ci034243x).
- [hERG binding-site mutagenesis](https://pubmed.ncbi.nlm.nih.gov/18987434/).
- [RDKit conformer and descriptor documentation](https://www.rdkit.org/docs/RDKit_Book.html).

The references support method definitions and clinical interpretation, not validation of the repository's new heuristic coefficients.
