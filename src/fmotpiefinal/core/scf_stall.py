"""Bound a stalled numerical attempt; never accept an unconverged result."""
import re
import math


class MonomerSCFBudget:
    def __init__(self, consecutive_cycles=3):
        if consecutive_cycles < 3:
            raise ValueError('at least three consecutive cycles required')
        self.limit = consecutive_cycles
        self.current = None
        self.failed = set()
        self.streaks = {}
        self.last_iteration = None
        self.evidence = None
        self.energy_residuals = []

    def consume(self, line):
        header = re.search(r'\bRHF\s+(monomer|dimer)\s+(\d+)', line)
        if header:
            self.current = int(header[2]) if header[1] == 'monomer' else None
        if 'SCF IS UNCONVERGED' in line and self.current is not None:
            self.failed.add(self.current)
        cycle = re.search(r'FMO iter\s*(\d+)\.(\d+),', line)
        if cycle is None:
            return self.evidence
        iteration, layer = int(cycle[1]), int(cycle[2])
        if layer != 1 or (self.last_iteration is not None and iteration != self.last_iteration + 1):
            self.streaks = {}
            self.energy_residuals = []
        self.streaks = {fragment: self.streaks.get(fragment, 0) + 1 for fragment in self.failed}
        self.failed = set()
        self.current = None
        self.last_iteration = iteration
        residual = re.search(r'maxE\s*([0-9.EeDd+-]+)', line)
        if layer == 1 and residual is not None:
            value = float(residual[1].replace('D', 'E').replace('d', 'e'))
            self.energy_residuals.append(value)
            self.energy_residuals = self.energy_residuals[-24:]
            if len(self.energy_residuals) == 24 and min(self.energy_residuals) > 1e-7:
                for period in (2, 3, 4, 6):
                    reference = self.energy_residuals[-period:]
                    if all(math.isclose(value, reference[index % period], rel_tol=0.05, abs_tol=1e-10)
                           for index, value in enumerate(self.energy_residuals)):
                        self.evidence = {'schema': 'bounded-outer-fmo-attempt-v1',
                            'reason': 'OUTER_FMO_RESIDUAL_CYCLE_BUDGET', 'period': period,
                            'observed_cycles': 24, 'last_fmo_iteration': iteration,
                            'minimum_energy_residual_hartree': min(self.energy_residuals),
                            'requires_numerical_retry': True, 'ranking_eligible': False}
                        break
        repeated = sorted(fragment for fragment, count in self.streaks.items() if count >= self.limit)
        if layer == 1 and repeated:
            self.evidence = {'schema': 'bounded-monomer-scf-attempt-v1',
                'reason': 'REPEATED_MONOMER_SCF_FAILURE_BUDGET', 'fragments': repeated,
                'consecutive_cycles': self.limit, 'last_fmo_iteration': iteration,
                'requires_numerical_retry': True, 'ranking_eligible': False}
        return self.evidence
