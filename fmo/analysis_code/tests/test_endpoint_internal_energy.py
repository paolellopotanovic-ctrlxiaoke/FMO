import sys
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))
from analyze_endpoint_internal_energy import evaluate, mean_difference, parameterized_context, physical_atom_maps


def test_double_difference_uses_within_molecule_reference():
    bound = [{"mean_kcal_mol": {"total": 105.}}, {"mean_kcal_mol": {"total": 107.}}]
    solution = [{"mean_kcal_mol": {"total": 101.}}]
    assert mean_difference(bound, solution)["total"] == pytest.approx(5.)


def test_physical_atom_maps_use_original_order_for_state_b():
    from rdkit import Chem

    first = Chem.AddHs(Chem.MolFromSmiles("C"))
    second = Chem.AddHs(Chem.MolFromSmiles("N"))
    component = lambda mol: SimpleNamespace(to_rdkit=lambda: mol)
    mapping = SimpleNamespace(componentA=component(first), componentB=component(second),
                              componentA_to_componentB={1: 1, 2: 2, 3: 3})
    atoms = [{"class": .25, "element": "C"}, *[{"class": .5, "element": "H"}] * 3,
             {"class": .25, "element": "H"}, {"class": .75, "element": "N"}]
    indices, _ = physical_atom_maps(atoms, mapping)
    assert indices[0].tolist() == [0, 1, 2, 3, 4]
    assert indices[10].tolist() == [5, 1, 2, 3]
    atoms[-1]["element"] = "O"
    with pytest.raises(ValueError, match="element"):
        physical_atom_maps(atoms, mapping)


def test_cached_forcefield_split_does_not_double_count(tmp_path):
    import numpy as np
    from rdkit import Chem

    molecule = Chem.MolFromSmiles("[C].[C]")
    for atom, charge in zip(molecule.GetAtoms(), [.1, -.1], strict=True):
        atom.SetDoubleProp("PartialCharge", charge)
    text = ('<ForceField><AtomTypes><Type name="ct" class="ct" element="C" mass="12"/></AtomTypes>'
            '<Residues><Residue name="LIG"><Atom name="C1" type="ct" charge="0.1"/>'
            '<Atom name="C2" type="ct" charge="-0.1"/></Residue></Residues>'
            '<NonbondedForce coulomb14scale="0.833333" lj14scale="0.5"><UseAttributeFromResidue name="charge"/>'
            '<Atom type="ct" sigma="0.3" epsilon="0.2"/></NonbondedForce></ForceField>')
    path = tmp_path / "cache.json"
    path.write_text(json.dumps({"openff-2.2.0": {"1": {"ffxml": text}}}))
    context, integrator, groups, _ = parameterized_context(path, molecule)
    try:
        result = evaluate(context, groups, np.array([[0., 0., 0.], [4., 0., 0.]]))
        expected_lj = 4 * .2 * ((.3 / .4) ** 12 - (.3 / .4) ** 6) / 4.184
        assert result["Nonbonded_LennardJones"] == pytest.approx(expected_lj)
        assert result["total"] == pytest.approx(result["NonbondedForce"])
        assert result["NonbondedForce"] == pytest.approx(result["Nonbonded_LennardJones"] + result["Nonbonded_Coulomb"])
    finally:
        del context, integrator
