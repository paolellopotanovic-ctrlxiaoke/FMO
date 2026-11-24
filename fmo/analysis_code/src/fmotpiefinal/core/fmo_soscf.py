"""Explicit FMO density-SOSCF controls, without changing the physical model."""
import re


def enable_fmo_soscf(deck: str) -> str:
    before = "MCONV(2)=530 MCONV(4)=530"
    if deck.count(before) != 1 or len(re.findall(r"\bITYPMX\s*=\s*-1\b", deck)) != 1:
        raise ValueError("FMO SOSCF requires original convergers and ITYPMX=-1")
    if "DIIS=.FALSE. SOSCF=.TRUE." not in deck:
        raise ValueError("FMO SOSCF requires compatible SCF controls")
    return deck.replace(before, "MCONV(2)=576 MCONV(4)=576", 1)


def verify_fmo_soscf_log(path) -> None:
    standard, soscf = False, False
    with path.open(errors="replace") as handle:
        for line in handle:
            standard |= "USE GAMESS STANDARD SUBROUTINES FOR SCF" in line
            soscf |= "SOSCF IN EFFECT" in line
            if "DIIS CHARGE UPDATE" in line or "BROYDEN'S CHARGE MIXING" in line:
                raise ValueError("DFTB charge mixer overrides the requested density SOSCF")
    if not standard or not soscf:
        raise ValueError("actual density SOSCF is not verified in the GAMESS log")
