# Task 20 · Input and method provenance

Retrieved and inspected on 2026-09-21. Public structural observations, algorithm definitions and analyst assumptions are separate below.

| Source | Direct URL | Use and scope |
|---|---|---|
| RCSB PDB 6YB7 | https://www.rcsb.org/structure/6YB7 | Experimental unliganded Mpro, 1.25 Å; chain A has 305 modeled residues. Used for generation pocket coordinates. Biological assembly is dimeric; this pilot scores a rigid single chain. |
| Raw 6YB7 input | https://files.rcsb.org/download/6YB7.pdb | Archived byte-for-byte in `inputs/6YB7.pdb`. |
| RCSB PDB 6W63 | https://www.rcsb.org/structure/6W63 | Experimental Mpro–X77 complex, 2.10 Å. Used only for redocking control. RCSB lists its primary citation as to be published; no journal paper is invented. |
| Raw 6W63 input | https://files.rcsb.org/download/6W63.pdb | Archived byte-for-byte in `inputs/6W63.pdb`. |
| X77 chemical component | https://www.rcsb.org/ligand/X77 | Bond-order template and stereochemical identity. `https://files.rcsb.org/ligands/download/X77_ideal.sdf` is archived; ideal coordinates are not the experimental pose. |
| Douangamath et al., 2020 | https://doi.org/10.1038/s41467-020-18709-w | Associated Mpro fragment-screening study; not evidence for generated molecules. |
| RDKit documentation | https://www.rdkit.org/docs/GettingStartedInPython.html | Sanitization, conformers, descriptors, chemical features and atom mapping. |
| Ertl & Schuffenhauer, 2009 | https://doi.org/10.1186/1758-2946-1-8 | SA fragment/complexity score. Low SA is not a demonstrated synthesis route. |
| AutoDock Vina official tutorial | https://autodock-vina.readthedocs.io/en/latest/docking_basic.html | Meeko preparation, rigid-receptor docking and structure export. |
| AutoDock Vina 1.2.7 official release | https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7 | Executed portable Windows binary; original license Apache-2.0. Executable is not redistributed in this repository. |

Vina executable URL: `https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe`. SHA256: `e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5`. Input hashes are recorded in `outputs/summary.json`; all generated data hashes are recorded in `outputs/manifest.json`.

Analyst choices: eight cores, twenty left and sixteen right caps; 150–550 g/mol admission domain; 40 rigid orientations per unique graph; a 13 Å local receptor sphere; a 24 Å Vina box; pKa values 0/5/9 from a simple structural heuristic; contact and feature weights; one evolutionary seed and one random-order control. No trained affinity model, measured pKa, proprietary assay or patient data is used. Protein alternate locations blank/A are selected explicitly. Meeko residue templates determine protonation; crystal waters, biological dimer, receptor flexibility and microstate ensembles are absent.
