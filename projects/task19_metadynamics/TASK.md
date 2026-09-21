# Current Task 19 instruction

Source attachment SHA256: `cf0eab249b9496c13b141684eea48995108fcbe6a840bbf6abf6320be668bd58`

#### TASK 19: All-Atom Ternary Complex Metadynamics & Interface Cooperativity
- Construct an atomistic interface model of an E3-Ligand-Target ternary degradation complex:
  - Map the non-covalent protein-protein interaction (PPI) contacts at the induced composite binding interface.
- Implement a 2D Well-Tempered Metadynamics simulation engine:
  - Collective Variable 1 ($CV_1$): Protein-Protein Center-of-Mass distance ($d_{\text{E3-Target}}$).
  - Collective Variable 2 ($CV_2$): Relative inter-protein orientation / dihedral angle ($\phi$).
- Reconstruct the 2D Potential of Mean Force (PMF) landscape: Quantify the thermodynamic free energy of cooperativity ($\Delta\Delta G_{\text{coop}} = -R T \ln \alpha$) and map energetic basins corresponding to productive degradation geometries.

Requested figure: `fig19_ternary_metadynamics_pmf_landscape.png` (300 DPI).
