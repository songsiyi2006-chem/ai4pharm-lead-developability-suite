# Task 19: explicit atomistic ternary-complex pilot and free-energy acceptance gates

## Outcome and evidence

The code works from the first biological assembly of experimental 5T35: BRD4 BD2 (A), Elongin B (B), Elongin C (C), VHL (D), and MZ1 (CCD 759, author chain D). The crystal heavy-atom audit finds 112 interprotein atom pairs at or below 4.5 angstrom, involving 14 VHL–BRD4 residue pairs; the closest pair is 2.79836 angstrom. A distance cutoff establishes geometric proximity, not an interaction free energy or a hydrogen bond.

The executable atomistic calculation and its actual status are recorded in [summary.json](outputs/summary.json). All-atom molecular mechanics is used: Amber14 protein parameters, TIP3P explicit water, Sage 2.3 ligand parameters, and the Ash/NAGL 1.0 charge model declared by Sage 2.3. The charges are learned approximations to AM1-BCC ELF10 charges. Crystal and CCD stereochemical identities are checked before adding ligand hydrogens. Ligand net charge is zero within numerical precision.

The sampling schedule is deliberately bounded to the local CPU: 300 K, 1 fs timestep, 200 preliminary integration steps, 1,000 biased steps, and a hill every 50 biased steps. This is **0.2 ps initial thermalization plus 1.0 ps biased sampling**, not a thermally equilibrated binding ensemble. There are no claimed interbasin round trips, stationary state weights, productive-degradation basins or calibrated binding free energies. Successful execution of this pilot meets an engineering criterion only.

The actual run completed on 2026-09-21. Its measured results are:

| Quantity | Executed result |
| --- | --- |
| Explicit particles / deposited hills | 150,838 / 20 |
| CPU threads / wall time, including output generation | 1 / 1,013.86 s (16.90 min) |
| Distance at the recorded biased frames | 3.948636–3.972654 nm |
| Orientation at the recorded biased frames | −1.060035 to −0.989082 rad |
| Initial / minimized potential energy | −1,406,165.703 / −2,010,957.599 kJ/mol |
| Minimization callbacks | 300; limit 50 iterations per constraint-penalty stage |
| Finite atom positions and recorded energies | Passed |
| Serialized hill parameters checked against the live history | 60 / 60 |
| Converged PMF / cooperativity / productive basins | Unavailable / not computed / not identified |

The large absolute potential energy includes the solvent and is not a binding energy. The energy decrease during minimization is a preparation diagnostic. The moving distance and narrow angular excursion demonstrate local dynamics only.

## Structure preparation and interactions

SEQRES-aware PDBFixer preparation builds the internal unresolved Elongin C segment and missing atoms; unresolved terminal tails remain omitted. These built coordinates are modeled, not experimentally resolved. Protein hydrogens use pH 7 default variants. No protonation-state ensemble is sampled. Salt is placed at 0.15 M in a dodecahedral solvent box, with 1.0 nm nominal padding and 0.9 nm nonbonded cutoff. Waters that overlap the ligand by less than 0.30 nm are removed after protein solvation. An overlapping ion causes failure rather than changing net charge silently.

The merged system retains all protein/solvent terms, every ligand bond/angle/torsion/constraint, and all ligand nonbonded particles and exceptions. The combined NonbondedForce generates cross protein–ligand interactions; there is no artificial attractive PPI potential. Unlike a protein-only demonstration, MZ1 has a complete bonded and nonbonded parameter set. PME treats electrostatics; Lorentz–Berthelot mixing treats ordinary cross Lennard–Jones terms. The serialized system fixes the actual numerical force parameters.

## The biased coordinates

For a protein atom group with masses \(m_i\),

\[
\mathbf R=\frac{\sum_i m_i\mathbf r_i}{\sum_i m_i},\qquad
d=|\mathbf R_{\mathrm{VHL}}-\mathbf R_{\mathrm{BRD4}}|.
\]

Here the distance uses the heavy-atom mass centers of the VHL and BRD4 chains and has units nm. The orientation is the signed four-centroid dihedral formed by C-alpha anchor groups D:65–80, D:120–135, A:350–365 and A:430–445, in radians. A single dihedral is a reduced orientation coordinate; it does not specify all three rotational degrees of freedom, nor linker conformations.

The interaction system is periodic, but CV positions are whole-solute continuous Cartesian positions. The pilot never wraps the solute and aborts if any solute atom moves more than one-quarter of the shortest box-vector length after initial thermalization. This is a bounded-run validity check, **not a production solution to dissociation across periodic boundaries**. Production unbinding work must define molecule-whole imaging, avoid periodic copies and track orientation singularities. Distances outside the configured sampling interval also stop execution.

## Actual well-tempered hills

Let \(\gamma=(T+\Delta T)/T>1\), \(\delta\phi=\operatorname{atan2}(\sin(\phi-\phi_n),\cos(\phi-\phi_n))\). The deposited bias is

\[
V(d,\phi,t)=\sum_{t_n<t}w_n
\exp\left[-\frac{(d-d_n)^2}{2\sigma_d^2}-\frac{\delta\phi^2}{2\sigma_\phi^2}\right],\qquad
w_n=w_0\exp\left[-\frac{V(d_n,\phi_n,t_n^-)}{RT(\gamma-1)}\right].
\]

The analytic centroid force applies this potential directly to the real atoms. The force engine differentiates the mass centers and dihedral; it is not a Monte Carlo walker on a invented surface. The shortest wrapped angular displacement is continuous in bias value at the periodic boundary; its antipodal slope discontinuity is suppressed by \(\exp[-\pi^2/(2\sigma_\phi^2)]\), far below machine precision for the chosen 0.15 rad width. Parameters are \(w_0=0.25\) kJ/mol, \(\gamma=10\), \(\sigma_d=0.05\) nm and \(\sigma_\phi=0.15\) rad.

Each deposition records its prior bias and decreased height. An independent NumPy sum generates the diagnostic grid. In the long-time, adequate-sampling limit, well-tempered theory gives

\[
F(d,\phi)=-\frac{\gamma}{\gamma-1}V(d,\phi)+C.
\]

Applying that algebra to 20 local hills does **not** meet the limit. The figure labels it an **unconverged bias diagnostic, not a PMF**. Its untouched grid is unsampled, not a flat physical landscape. Grid basins cannot be called cooperativity basins. Full trajectories and hill logs retain what actually occurred.

## Why cooperativity remains null

For a degrader P, ligase E and target T, a conventional equilibrium definition is

\[
\alpha=\frac{K_D(E+P\rightleftharpoons EP)}{K_D(E+PT\rightleftharpoons EPT)},\qquad
\Delta\Delta G_{\mathrm{coop}}=-RT\ln\alpha.
\]

Since \(\Delta G_{\mathrm{bind}}^\circ=RT\ln(K_D/C^\circ)\),

\[
\Delta\Delta G_{\mathrm{coop}}=
\Delta G^\circ_{E|PT}-\Delta G^\circ_{E|P}.
\]

Both legs must share standard-state, restraint, solvation, protonation and volume conventions. A second, independent cycle through the target-binding legs should close within uncertainty. A ternary distance/orientation bias alone lacks the binary reference and standard-state corrections. Even a converged ternary PMF would therefore be insufficient by itself.

Production work requires independent replicas, reversible exploration of relevant basins, time-block stability, sampling of hidden slow variables, angular/translational Jacobians, restraint-release and standard-concentration corrections, appropriate protonation alternatives, and sensitivity to the rebuilt loop and force field. A structural basin cannot be labeled productive degradation from this VCB–BRD4 domain construct alone: the complete catalytic assembly, target lysine accessibility, ubiquitination and cellular degradation evidence are missing. Biophysical binding and cellular degradation assays are separate validation gates.

## Reproduction and audit

See [README](README.md), [sources](SOURCES.md), [configuration](outputs/config.json), [CV definitions](outputs/cv_definition.json), [preparation](outputs/preparation.json), [verification](outputs/verification.json) and [output hashes](outputs/manifest.json). The full atomic parameter cache is checked before reuse. Cached coordinates are restarted at PDB precision; velocity and Langevin seeds are declared. Numerical files and a verified checkpoint are written after each hill, with configuration/step/hash checks on restart. The CPU calculation requires no supercomputer. It supplies reusable atomistic infrastructure and a measurable local pilot while leaving publication-level free-energy objectives explicitly unfinished.
