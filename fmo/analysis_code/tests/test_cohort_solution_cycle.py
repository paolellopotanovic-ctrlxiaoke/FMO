import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from run_cohort_solution_cycle import canonical_receptor, cycle_decks


def test_shared_receptor_normalizes_titles_but_preserves_physical_groups():
    first = "! ligand one\n $PCM IDISP=0 $END\n $DATA\nligand one\nC1\nC-1 6\n $END\n"
    second = first.replace("ligand one", "ligand two")
    assert canonical_receptor(first) == canonical_receptor(second)
    assert "$PCM IDISP=0 $END" in canonical_receptor(first)
    assert "C-1 6" in canonical_receptor(first)
    assert "$DATA\n\n" not in canonical_receptor(first)
    assert "$DATA\nMatched common receptor\nC1" in canonical_receptor(first)


def test_full_cycle_does_not_subtract_selective_totals_or_mismatch_solvent(tmp_path):
    source = """ $CONTRL RUNTYP=ENERGY SCFTYP=RHF $END
 $PCM SOLVNT=WATER IEF=-10 ICAV=1 IDISP=1 IFMO=-1 $END
 $FMOPRP IPIEDA=1 $END
 $DFT DC=.TRUE. IDCVER=4 $END
 $FMO NFRAG=2 NBODY=2 NACUT=0 MODMOL=1 MOLFRG(1)=2
 ICHARG(1)=0, 1
 FRGNAM(1)=REC, LIG
 INDAT(1)=0 1 -2 0 3 -4 0
 $END
 $FMOXYZ
 C 6 0 0 0
 H 1 1 0 0
 N 7 2 0 0
 H 1 3 0 0
 $END
 $DATA
 original
 C1
 C-1 6
 H-1 1
 N-1 7
 $END
"""
    decks = cycle_decks(source, tmp_path)
    assert all("MODMOL" not in deck and "IDISP=0" in deck for deck in decks.values())
    assert all("$DFT DC=.TRUE. IDCVER=4" in deck for deck in decks.values())
    assert "IPIEDA=1" in decks["complex"]
    assert "NFRAG=1 NBODY=1" in decks["ligand"]
    assert "ICHARG(1)=1" in decks["ligand"]
    assert "N 7 2 0 0" not in decks["receptor"]
