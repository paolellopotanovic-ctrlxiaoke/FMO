import importlib.util
from pathlib import Path

import pytest


SPEC = importlib.util.spec_from_file_location("matched_openfe", Path(__file__).parents[1] / "tools/run_matched_openfe.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_network_ignores_experimental_column(tmp_path):
    path = tmp_path / "edges.csv"
    path.write_text("Ligand 1,Ligand 2,ddG (kcal/mol)\na,b,not-a-number\nb,c,secret\n")
    assert MODULE.read_edges(path, ["a", "b", "c"]) == [("a", "b"), ("b", "c")]


@pytest.mark.parametrize("edges", ["a,b\n", "a,b\nb,a\nb,c\n", "a,z\nb,c\n"])
def test_invalid_network_fails_closed(tmp_path, edges):
    path = tmp_path / "edges.csv"
    path.write_text("Ligand 1,Ligand 2\n" + edges)
    with pytest.raises(ValueError):
        MODULE.read_edges(path, ["a", "b", "c"])


def test_installed_openfe_api_and_smoke_settings(tmp_path):
    openfe = pytest.importorskip("openfe")
    from openff.units import unit
    from rdkit import Chem
    from rdkit.Chem import AllChem

    solvent = openfe.SolventComponent(ion_concentration=0.15 * unit.molar)
    protocol = MODULE.make_protocol("smoke", "CPU")
    assert protocol.settings.protocol_repeats == 1
    assert protocol.settings.output_settings.output_indices == "all"
    molecules = []
    for smiles in ["CC", "CCC"]:
        molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(molecule, randomSeed=123)
        molecules.append(openfe.SmallMoleculeComponent.from_rdkit(molecule))
    mapping = openfe.LigandAtomMapping(molecules[0], molecules[1], {0: 0, 1: 1})
    mapping.to_json(tmp_path / "mapping.json")
    assert openfe.LigandAtomMapping.from_json(tmp_path / "mapping.json") == mapping
    state_a = openfe.ChemicalSystem({"solvent": solvent, "ligand": molecules[0]})
    state_b = openfe.ChemicalSystem({"solvent": solvent, "ligand": molecules[1]})
    transformation = openfe.Transformation(stateA=state_a, stateB=state_b, mapping=mapping, protocol=protocol)
    transformation.to_json(tmp_path / "transformation.json")
    assert len(transformation.create().protocol_units) > 0
