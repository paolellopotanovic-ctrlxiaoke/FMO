"""Bound repeated SCF failures without accepting a partial cohort."""
import json
import re
from pathlib import Path


def first_failed_monomer(deck: str, log: Path) -> str | None:
    inventory = re.search(r'FRGNAM\(1\)\s*=\s*(.*?)\bINDAT\(1\)', deck, re.S)
    if inventory is None:
        return None
    names = [name.strip() for name in inventory[1].split(',') if name.strip()]
    current = None
    with log.open(errors='replace') as handle:
        for line in handle:
            header = re.search(r'\bRHF\s+(monomer|dimer)\s+(\d+)', line)
            if header:
                current = int(header[2]) if header[1] == 'monomer' else None
            if 'SCF IS UNCONVERGED' in line:
                return names[current - 1] if current is not None and 0 < current <= len(names) else None
    return None


def repeated_scf_failure(root: Path, failed_names: list[str]) -> dict | None:
    if len(failed_names) < 2:
        return None
    observations = {}
    for name in failed_names:
        folder = root / 'qm' / name
        result = folder / 'run/result.json'
        if not result.is_file():
            continue
        record = json.loads(result.read_text())
        if record.get('status') != 'FAILED_SCF_CONVERGENCE':
            continue
        budget = record.get('numerical_budget', {})
        if budget.get('reason') == 'OUTER_FMO_RESIDUAL_CYCLE_BUDGET':
            signature = f"OUTER_FMO_PERIOD_{budget['period']}"
            observations.setdefault(signature, []).append(name)
            continue
        deck, log = folder / 'input.inp', folder / 'run/gamess.log'
        if not deck.is_file() or not log.is_file():
            continue
        fragment = first_failed_monomer(deck.read_text(), log)
        if fragment is not None:
            observations.setdefault(fragment, []).append(name)
    repeated = sorted((fragment, names) for fragment, names in observations.items() if len(names) >= 2)
    if not repeated:
        return None
    fragment, names = repeated[0]
    return {'schema': 'repeated-monomer-scf-cohort-gate-v1', 'fragment': fragment,
            'failed_cases': names, 'ranking_eligible': False,
            'action': 'STOP_NEW_SUBMISSIONS_AND_REQUIRE_WHOLE_COHORT_RESCUE',
            'interpretation': 'Repeated first explicit failed monomer, not unique physical root cause'}
