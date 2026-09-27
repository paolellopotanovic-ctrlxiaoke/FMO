from pathlib import Path
from sophosqm_baseline.parser import aggregate_tpie, parse_pie, parse_fragment_names, find_selected_tie


def _pieda_row(i: int, j: int, total: float, es: float, ex: float, ct: float, di: float, sol: float = 0.0) -> str:
    return (
        f"{i:5d}{j:5d}{'C1':3s}{0:3d}{1.0:7.2f}{0.0:8.4f}"
        f"{0.0:10.3f}{0.0:9.3f}{total:10.3f}{es:10.3f}"
        f"{ex:9.3f}{ct:9.3f}{di:9.3f}{sol:9.3f}"
    ) + "\n"


def test_gamess_pieda_parser_and_tpie(tmp_path: Path):
    output = tmp_path / "gamess.out"
    output.write_text(
        "CONV\n" + "=" * 80 + "\n  1 LIG\n  2 PROA\n  3 PROB\n\n Close fragment pairs\n"
        + "-" * 105 + "\n"
        + _pieda_row(2, 1, -2.5, -2.0, -0.2, -0.1, -0.2)
        + _pieda_row(3, 1, 1.25, 1.0, 0.1, 0.05, 0.1)
        + _pieda_row(3, 2, -9.0, -8.0, -0.5, -0.2, -0.3)
        + "\n Total energy\n",
        encoding="utf-8",
    )
    records = parse_pie(output)
    assert aggregate_tpie(records, ligand_fragment=1) == -1.25
    assert find_selected_tie(output, ligand_fragment=1) is None
    assert records[0].electrostatic_kcal_mol == -2.0
    assert parse_fragment_names(output) == {1: "LIG", 2: "PROA", 3: "PROB"}


def test_parser_rejects_placeholder_records(tmp_path: Path):
    output = tmp_path / "fake.out"
    output.write_text("PIE(1, 3) = -2.50\n", encoding="utf-8")
    try:
        parse_pie(output)
    except ValueError as error:
        assert "PIEDA" in str(error)
    else:
        raise AssertionError("placeholder PIE text must not be treated as GAMESS output")


def test_parser_reads_standard_two_body_fmo_table(tmp_path: Path):
    output = tmp_path / "standard.out"
    output.write_text(
        "Two-body FMO properties.\n"
        "========================\n"
        "    I    J DL  Z    R   Q(I->J)  EIJ-EI-EJ dDIJ*VIJ    total     Ees      Eex    Ect+mix   Erc+di    Gsol\n"
        " ---------------------------------------------------------------------------------------------------------\n"
        "    2    1 C1  0   0.76  0.0282    -8.296   -0.455    -8.751    -9.826    5.191   -2.516   -1.599    0.000\n"
        "Total energy of the molecule: Ecorr (2)= -1.0\n",
        encoding="utf-8",
    )
    records = parse_pie(output)
    assert records[0].fragment_i == 2
    assert records[0].value_kcal_mol == -8.751
    assert records[0].electrostatic_kcal_mol == -9.826


def test_parser_reads_official_selected_tie(tmp_path: Path):
    output = tmp_path / "tie.out"
    output.write_text(
        " Total interaction (TIE) to selected fragments in kcal/mol\n"
        "    1(LIG     ): ES   -2.0000, SCF   -1.0000, final TIE   -1.2500\n",
        encoding="utf-8",
    )
    assert find_selected_tie(output, ligand_fragment=1) == -1.25
