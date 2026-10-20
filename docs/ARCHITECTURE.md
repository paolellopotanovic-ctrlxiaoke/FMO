# Architecture

There is one public application entry point: `fmotpiefinal`, implemented in
`pipeline/cli.py`. Helpers do not duplicate the public pipeline.

| Layer | Responsibility | Main files |
| --- | --- | --- |
| Input | Component inventory, PDB/mmCIF and pose validation | `core/structure_input.py`, `pipeline/target_intake.py` |
| Chemistry graph | Connectivity, missing atoms, unsupported chemistry, repair plan | `core/receptor_chemical_graph.py`, `core/chemical_graph_policy.py`, `core/graph_preparation.py` |
| Preparation | Chain selection, structural repair, protonation, water policy | `pipeline/chain_selection.py`, `pipeline/structure_repair.py`, `pipeline/structure_treatment.py`, `pipeline/water_policy.py` |
| MM | Restrained ligand preparation, receptor/component models | `core/minimization.py`, `core/special_chemistry.py`, `core/special_dispatch.py` |
| QM | Pocket, fragments, decks, execution, convergence, PIEDA parsing | `core/bound_geometry.py`, `core/fragmentation.py`, `core/fmo_input.py`, `core/runner.py`, `core/fmo_convergence.py`, `core/parser.py` |
| Control | Complete-cohort execution and registered failure recovery | `pipeline/target_workflow.py`, `pipeline/target_stages.py`, `pipeline/single_case_rescue.py`, `pipeline/cohort_rescue.py` |
| Output | Eligibility checks, score transformation, ranking | `scoring.py` |

The orchestration layer retains separate stages because a full calculation needs
restartable, inspectable boundaries. In-process stage execution removes repeated
interpreter startup without changing stage arguments or bypassing gates.
`FMOTPIE_STAGE_EXECUTION=subprocess` provides the conservative subprocess route.

Legacy experiment names were replaced with descriptive module names. The exact
mapping and source/export hashes are in `SOURCE_PROVENANCE.json`. Only two live
bound-geometry helpers were extracted from the old three-state-energy module;
the unused energy/optimization experiments are not distributed. Function bodies
in the retained computation path are otherwise unchanged by the export, apart
from the temporary-file prefix and the protocol utility's default output label.

Advanced interface-intake/interface-scoring commands remain explicit audits; the
909-structure evidence concerns protein-small-molecule ranking, not a validated
general protein/RNA binding predictor. Adding a new chemistry model requires
MM, fragmentation/QM parameter, graph and convergence tests before enabling it.
