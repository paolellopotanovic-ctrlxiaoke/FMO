# Default scientific workflow

1. Normalize PDB/mmCIF, select an explicit model for multi-model files, resolve
   alternate locations consistently and inventory components before filtering.
2. Validate ligand identities, charge/stereochemistry and receptor-frame poses.
   Audit contacting chains and peptide connectivity; repair numbering only when
   geometry supports the proposed connectivity.
3. PDBFixer performs structural completeness handling. PDB2PQR/PROPKA assigns
   protonation/hydrogens at pH 7.4. Missing/ambiguous segments are not automatically
   evidence of a trustworthy reconstructed binding site.
4. Graph/chemistry checks identify missing atoms, local bond problems, hydrogen
   clashes, nonstandard components and model requirements. Only registered,
   verified repair routes execute automatically, followed by re-audit. Retained
   ions/phosphoresidues/caps use explicit supported representations; this is not
   support for arbitrary metals or covalent adducts.
5. Apply the label-free target-level water/contact policy, maintaining consistent
   receptor treatment across the predefined ligand cohort. This is **not**
   unconditional deletion of every water or every hydrogen. The default water
   hydrogen orientation comes from PDB2PQR; source-water orientation is a separate,
   explicitly requested experimental option.
6. Prepare ligands with OpenMM/OpenFF and restrained k10 minimization
   (10 kJ mol^-1 angstrom^-2). Receptor/water geometry is fixed in this step;
   ligand atoms are minimized with the registered restraints. No default
   multi-conformer search, entropy correction or affinity-driven pose selection.
7. Build the shared full-cohort 5 angstrom pocket and fragments; run GAMESS
   FMO2-DFTB3/PCM with PIEDA. Full tPIE sums the protein-ligand fragment interaction
   entries; it is not isolated A/B/AB total-energy subtraction. PBSA is not added.
8. Require complete coverage and final convergence proofs, with registered
   label-blind rescue when eligible. Record any cohort-wide alternate numerical
   protocol; do not cherry-pick a different scientific model for each label.
9. Rank by `s = -tPIE`. The optional soft-cap computes
   `median + width * tanh((s - median)/width)`, where
   width = 3 * scaled_MAD and scaled_MAD = 1.4826 * MAD.
   Its strictly monotone form preserves ordinary within-target rank order except
   possible floating-point saturation/ties. Pearson may change; improvement is
   not guaranteed. Zero spread is reported as undefined, not divided by zero.

For recognized historical target names the packaged references use frozen
retrospective centers/scales. For a new target, the local transformation is
descriptive and cohort-dependent unless a reference is frozen and reused using
`--score-reference`. Do not compare these proxies as calibrated cross-target
affinities. Input affinity labels never enter preparation or the rescue decision.

## Underlying methods to acknowledge

Before publication, supply verified citations and licensing checks for GAMESS,
FMO/PIEDA, DFTB3/3OB, PCM, PDBFixer/OpenMM, PDB2PQR/PROPKA, RDKit,
OpenFF/NAGL and the Amber/solvent force-field assets. The package does not claim
authorship or redistribute binaries for those external methods.
