"""Bounded SOSCF startup controls that preserve final convergence accuracy."""


def soscf_onset_deck(deck: str, onset: float) -> str:
    if onset == 0.25:
        return deck
    before = " $SCF DIRSCF=.TRUE. NPUNCH=0 CONV=1.0E-8 DIIS=.FALSE. SOSCF=.TRUE. $END"
    if onset != 1.0 or deck.count(before) != 1 or "MCONV(2)=576 MCONV(4)=576" not in deck:
        raise ValueError("unregistered or incompatible SOSCF startup control")
    return deck.replace(before, " $SCF DIRSCF=.TRUE. NPUNCH=0 CONV=1.0E-8\n"
                        " DIIS=.FALSE. SOSCF=.TRUE. SOGTOL=1.0 $END", 1)
