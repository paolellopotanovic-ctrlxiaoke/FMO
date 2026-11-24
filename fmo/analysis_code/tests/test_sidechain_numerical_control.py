from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from run_sidechain_numerical_control import strict_damping_control, strict_diis_control, tight_soscf_control


def test_strict_control_changes_only_registered_damping():
    source = " $SCF CONV=1.0E-8 $END\n $FMOPRP CONV=1.0E-7 CNVDMP=50 MCONV(2)=576 $END"
    candidate = strict_damping_control(source)
    assert candidate.replace("CNVDMP=0", "CNVDMP=50") == source
    assert strict_damping_control(candidate) == candidate


def test_ambiguous_or_missing_damping_is_rejected():
    for source in ("CONV=1.0E-7", "CNVDMP=50 CNVDMP=20"):
        with pytest.raises(ValueError):
            strict_damping_control(source)


def test_diis_retry_preserves_geometry_and_final_tolerances():
    source = "$SCF CONV=1.0E-8 DIIS=.FALSE. SOSCF=.TRUE. $END\n$FMOPRP CONV=1.0E-7 CNVDMP=50 MCONV(2)=576 MCONV(4)=576 $END\n$FMOXYZ\nC 6.0 1.2 3.4 5.6\n$END"
    candidate = strict_diis_control(source)
    assert candidate.split("$FMOXYZ")[1] == source.split("$FMOXYZ")[1]
    assert "CONV=1.0E-8" in candidate and "CONV=1.0E-7" in candidate
    assert "MCONV(2)=530" in candidate and "MCONV(4)=530" in candidate
    with pytest.raises(ValueError):
        strict_diis_control(candidate)


def test_tight_solver_only_tightens_inner_scf_tolerance():
    source = "$SCF CONV=1.0E-8 SOSCF=.TRUE. $END\n$FMOPRP CONV=1.0E-7 CNVDMP=50 $END\n$FMOXYZ\nC 6.0 1 2 3\n$END"
    candidate = tight_soscf_control(source)
    assert candidate.replace("CONV=1.0E-10", "CONV=1.0E-8").replace("CNVDMP=0", "CNVDMP=50") == source
    assert "CONV=1.0E-7" in candidate
