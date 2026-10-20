from pathlib import Path
from sophosqm_baseline.gamess import write_job_manifest, GamessConfig
from sophosqm_baseline.fmo_input import render_fmo2_pieda
from sophosqm_baseline.fragmentation import Fragment, _prepared_residue_charge
from sophosqm_baseline.structure import Atom, Ligand, Protein


def test_job_manifest(tmp_path: Path):
    path = write_job_manifest(tmp_path, GamessConfig("gamess", "job.inp"), "s.json", "f.json", {"x": 1})
    assert path.exists()


def test_fmo_deck_has_explicit_fragments():
    protein = Protein((Atom(0, "O", (0.0, 0.0, 0.0)), Atom(1, "H", (0.0, 0.7, 0.0))))
    ligand = Ligand("LIG", (Atom(0, "O", (2.0, 0.0, 0.0)), Atom(1, "H", (2.0, 0.7, 0.0))))
    fragments = [
        Fragment(1, "PRO", "protein_residue", (0, 1), charge=1),
        Fragment(2, "LIG", "ligand", (0, 1), charge=1),
    ]
    deck = render_fmo2_pieda(protein, ligand, fragments,
                             dftb_parameter_dir="/parameters/3OB-3-1")
    assert "UNITS=ANGS" in deck
    assert "$FMOPRP NPRINT=9 IPIEDA=1 MAXIT=200 CONV=1.0E-7 NGUESS=2" in deck
    assert " CNVDMP=50 MCONV(2)=530 MCONV(4)=530 $END" in deck
    assert "NFRAG=2 NBODY=2" in deck
    fmo_xyz = deck.split(" $FMOXYZ", 1)[1].split(" $END", 1)[0]
    assert sum(line.lstrip().startswith("O") and "  8 " in line for line in fmo_xyz.splitlines()) >= 2
    assert "INDAT(1)=0 1 -2 0 3 -4 0" in deck
    assert "$FMOBND" not in deck
    assert max(map(len, deck.splitlines())) <= 80


def test_fmo_deck_rejects_odd_electron_singlet_fragment():
    protein = Protein((Atom(0, "C", (0.0, 0.0, 0.0)),))
    ligand = Ligand("LIG", (Atom(0, "H", (2.0, 0.0, 0.0)),))
    fragments = [Fragment(1, "PRO", "protein_residue", (0,)), Fragment(2, "LIG", "ligand", (0,))]
    try:
        render_fmo2_pieda(protein, ligand, fragments,
                          dftb_parameter_dir="/parameters/3OB-3-1")
    except ValueError as error:
        assert "electrons" in str(error)
    else:
        raise AssertionError("odd-electron singlet fragments must fail before GAMESS")


def test_dftb3_deck_selects_gamess_dftb_backend():
    protein = Protein((Atom(0, "O", (0.0, 0.0, 0.0)), Atom(1, "H", (0.0, 0.7, 0.0))))
    ligand = Ligand("LIG", (Atom(0, "O", (2.0, 0.0, 0.0)), Atom(1, "H", (2.0, 0.7, 0.0))))
    fragments = [Fragment(1, "PRO", "protein_residue", (0, 1), charge=1),
                 Fragment(2, "LIG", "ligand", (0, 1), charge=1)]
    deck = render_fmo2_pieda(protein, ligand, fragments,
                             dftb_parameter_dir="/parameters/3OB-3-1")
    assert "$BASIS GBASIS=DFTB $END" in deck
    assert "$DFTB SCC=.TRUE. DFTB3=.TRUE. DAMPXH=.TRUE. DAMPEX=4.00" in deck
    assert "$DFT DC=.TRUE. IDCVER=4 $END" in deck
    assert "$PCM SOLVNT=WATER IEF=-10 ICOMP=0 ICAV=1 IDISP=1 IFMO=-1 $END" in deck

    explicit = render_fmo2_pieda(protein, ligand, fragments,
                                  dftb_parameter_dir="/parameters/3OB-3-1")
    assert "$DFTBSK" in explicit
    assert "/parameters/3OB-3-1/H-H.skf" in explicit
    assert "/parameters/3OB-3-1/Br-Br.skf" not in explicit
    assert max(map(len, explicit.splitlines())) <= 80


def test_dftb_deck_can_include_explicit_fragment_bonds():
    protein = Protein((Atom(0, "C", (0.0, 0.0, 0.0)), Atom(1, "H", (0.0, 1.0, 0.0))))
    ligand = Ligand("LIG", (Atom(0, "C", (2.0, 0.0, 0.0)), Atom(1, "H", (2.0, 1.0, 0.0))))
    fragments = [Fragment(1, "PRO", "protein_residue", (0, 1), charge=0),
                 Fragment(2, "LIG", "ligand", (0, 1), charge=0)]
    deck = render_fmo2_pieda(protein, ligand, fragments,
                             dftb_parameter_dir="/parameters", fragment_bonds=[(-1, 3)])
    assert "$FMOHYB" in deck
    assert "-1 3 HOP_C" in deck
    assert "ICHARG(1)=0, 0" in deck


def test_gamess_log_classification_prioritizes_scf_failure():
    from sophosqm_baseline.runner import classify_gamess_log
    text = "SCF IS UNCONVERGED, TOO MANY ITERATIONS\nEXECUTION OF GAMESS TERMINATED -ABNORMALLY-"
    assert classify_gamess_log(0, text) == "FAILED_SCF_CONVERGENCE"
    assert classify_gamess_log(1, "") == "FAILED_GAMESS_RUNTIME"
    assert classify_gamess_log(0, "EXECUTION OF GAMESS TERMINATED NORMALLY") == "SUCCEEDED"


def test_facio_like_partition_moves_carbonyl_and_adds_caps():
    from sophosqm_baseline.fragmentation import make_facio_like_fragments
    protein = Protein((
        Atom(0, "N", (0.0, 0.0, 0.0), "A", "ALA", "1", "ATOM", "N"),
        Atom(1, "C", (1.4, 0.0, 0.0), "A", "ALA", "1", "ATOM", "CA"),
        Atom(2, "C", (2.5, 0.0, 0.0), "A", "ALA", "1", "ATOM", "C"),
        Atom(3, "O", (3.5, 0.0, 0.0), "A", "ALA", "1", "ATOM", "O"),
        Atom(4, "N", (3.7, 0.0, 0.0), "A", "GLY", "2", "ATOM", "N"),
        Atom(5, "C", (5.1, 0.0, 0.0), "A", "GLY", "2", "ATOM", "CA"),
        Atom(6, "C", (6.2, 0.0, 0.0), "A", "GLY", "2", "ATOM", "C"),
        Atom(7, "O", (7.2, 0.0, 0.0), "A", "GLY", "2", "ATOM", "O"),
    ))
    ligand = Ligand("LIG", (Atom(0, "C", (10.0, 0.0, 0.0)),))
    prepared, fragments, metadata = make_facio_like_fragments(protein, ligand, include_waters=False)
    assert metadata["protocol"] == "SophosQM-like FACIO approximation"
    assert len(metadata["cap_records"]) == 2
    assert {atom.residue_name for atom in prepared.atoms if atom.residue_name in {"ACE", "NME"}} == {"ACE", "NME"}
    ace_atoms = [atom for atom in prepared.atoms if atom.residue_name == "ACE"]
    nme_atoms = [atom for atom in prepared.atoms if atom.residue_name == "NME"]
    assert len(ace_atoms) == 6 and len(nme_atoms) == 6
    assert any(fragment.kind == "protein_residue" and
               any(prepared.atoms[index].residue_name == "ACE" for index in fragment.atom_indices)
               for fragment in fragments)
    assert 2 not in fragments[0].atom_indices
    assert 2 in fragments[1].atom_indices
    assert metadata["hop_pairs_gamess_indices"]

    prepared_again, fragments_again, metadata_again = make_facio_like_fragments(
        prepared, ligand, include_waters=False)
    assert len([atom for atom in prepared_again.atoms if atom.residue_name in {"ACE", "NME"}]) == 12
    assert len(metadata_again["cap_records"]) == 2
    assert len(fragments_again) == len(fragments)


def test_nme_cap_replaces_c_terminal_oxt():
    from sophosqm_baseline.fragmentation import make_facio_like_fragments
    protein = Protein((
        Atom(0, "N", (0.0, 0.0, 0.0), "A", "ALA", "1", "ATOM", "N"),
        Atom(1, "C", (1.4, 0.0, 0.0), "A", "ALA", "1", "ATOM", "CA"),
        Atom(2, "C", (1.8, 1.3, 0.2), "A", "ALA", "1", "ATOM", "C"),
        Atom(3, "O", (1.2, 2.2, 0.6), "A", "ALA", "1", "ATOM", "O"),
        Atom(4, "O", (2.8, 1.4, -0.4), "A", "ALA", "1", "ATOM", "OXT"),
        Atom(5, "H", (-0.9, 0.2, 0.2), "A", "ALA", "1", "ATOM", "H"),
    ))
    ligand = Ligand("LIG", (Atom(20, "C", (10.0, 0.0, 0.0)),))
    _, fragments, metadata = make_facio_like_fragments(protein, ligand)
    fragment = fragments[0]
    assert 4 not in fragment.atom_indices
    assert metadata["nme_replaced_oxt_atom_indices"]["A:ALA:1"] == [4]
    nme_atoms = {
        tuple(record["atom_indices"]) for record in metadata["cap_records"]
        if record["type"] == "NME"
    }
    assert any(set(fragment.atom_indices) & set(atoms) for atoms in nme_atoms)


def test_openmm_histidine_protonation_is_restored_from_topology():
    atoms = [
        Atom(0, "N", (0.0, 0.0, 0.0), "A", "HIS", "907", "ATOM", name)
        for name in ("N", "HD1", "HE2")
    ]
    assert _prepared_residue_charge("HIS", atoms) == 1
    assert _prepared_residue_charge("HIS", atoms[:2]) == 0
