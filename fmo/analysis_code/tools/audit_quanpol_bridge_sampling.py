"""Audit all pilot evaluations, acceptance algebra, and saved coordinate mappings."""
import argparse
import json
from pathlib import Path

import numpy as np

from audit_quanpol_mm_parity import read_sections
from run_matched_openfe import digest, write_json
from sample_quanpol_bridge import log_acceptance
from validate_quanpol_embedded_forces import parse_embedded


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    report = {"status": "RUNNING", "states": {}, "full_free_energy_converged": False, "labels_used": False}
    for name in ("b0", "b10", "s0", "s10"):
        directory = args.root / name
        protocol = json.loads((directory / "protocol.json").read_text())
        if protocol["status"] != "PILOT_COMPLETE_NOT_CONVERGED_FREE_ENERGY":
            raise ValueError(f"{name} pilot incomplete")
        force_gate = json.loads((directory / "target_force_validation.json").read_text())
        if not force_gate["passed"]:
            raise ValueError("composed target force gate failed")
        rows = [json.loads(line) for line in (directory / "samples.jsonl").read_text().splitlines()]
        moves = protocol["moves_per_window"] + protocol["burnin_moves_per_window"]
        if len(rows) != moves * len(protocol["lambdas"]):
            raise ValueError("missing or extra samples")
        energies = [json.loads(path.read_text()) for path in sorted(directory.glob("bridge_*.energy.json"))]
        if len(energies) != len(rows) + 1:
            raise ValueError("missing QM evaluations")
        for ordinal, energy in enumerate(energies):
            if energy["call"] != ordinal:
                raise ValueError("QM call sequence mismatch")
            stem = directory / f"bridge_{name}_{ordinal:05d}"
            deck, log = stem.with_suffix(".inp"), stem.with_suffix(".stdout")
            if digest(deck) != energy["deck_sha256"] or digest(log) != energy["log_sha256"]:
                raise ValueError("QM input/output checksum mismatch")
            text = log.read_text()
            if "DENSITY CONVERGED" not in text:
                raise ValueError("explicit SCF convergence marker missing")
            parsed = parse_embedded(text)
            if parsed["n_physical"] != protocol["shell_particles"]:
                raise ValueError("physical atom inventory changed")
            if abs(parsed["energy_kcal_mol"] - energy["qm_kcal_mol"]) > 1e-8:
                raise ValueError("parsed QM energy mismatch")
            if abs(energy["qm_kcal_mol"] - energy["mm_shell_kcal_mol"] - energy["delta_kcal_mol"]) > 1e-8:
                raise ValueError("QM-minus-MM energy closure failed")
        previous = energies[0]["delta_kcal_mol"]
        last_accepted = 0
        snapshot_errors = []
        initial = np.load(directory / "initial.npz")
        indices = initial["shell_to_system"]
        for ordinal, row in enumerate(rows):
            if row["window"] != ordinal // moves or row["move"] != ordinal % moves or row["lambda"] != protocol["lambdas"][ordinal // moves]:
                raise ValueError("window/move/coupling sequence mismatch")
            proposal = energies[ordinal + 1]["delta_kcal_mol"]
            if abs(proposal - row["proposal_delta_kcal_mol"]) > 1e-8 or abs(previous - row["previous_delta_kcal_mol"]) > 1e-8:
                raise ValueError("proposal/current-state energy mapping mismatch")
            expected = log_acceptance(row["mm_integration_energy_error_kj_mol"], 0., proposal - previous, row["lambda"], protocol["temperature_kelvin"])
            if abs(expected - row["log_acceptance"]) > 1e-8 or expected >= 0 and not row["accepted"]:
                raise ValueError("Metropolis acceptance algebra mismatch")
            current = proposal if row["accepted"] else previous
            if abs(current - row["current_delta_kcal_mol"]) > 1e-8:
                raise ValueError("rejected state was not retained")
            if row["accepted"]:
                last_accepted = ordinal + 1
            frame = directory / f"frame_w{row['window']:02d}_m{row['move']:04d}.npz"
            if frame.exists():
                data = np.load(frame)
                actual = initial["shell_reference_angstrom"] + 10 * (data["positions_nm"][indices] - initial["full_reference_nm"][indices])
                deck = directory / f"bridge_{name}_{last_accepted:05d}.inp"
                coordinates = np.array([entry[2:5] for entry in read_sections(deck.read_text())["COORDINATES"]], dtype=float)
                error = float(np.max(np.abs(actual - coordinates)))
                snapshot_errors.append(error)
                if error > 1e-8 or abs(float(data["delta_kcal_mol"]) - current) > 1e-8:
                    raise ValueError("saved full-system frame/accepted QM shell mapping mismatch")
            previous = current
        report["states"][name] = {"passed": True, "qm_evaluations": len(energies), "normal_scf_converged_logs": len(energies),
                                  "samples_including_rejections": len(rows), "saved_full_system_frames_checked": len(snapshot_errors),
                                  "max_snapshot_shell_coordinate_error_angstrom": max(snapshot_errors),
                                  "acceptance_algebra_and_rejected_state_retention_passed": True,
                                  "composed_target_force_gate_passed": True}
    report["status"] = "PASSED_PILOT_INTEGRITY_NOT_FREE_ENERGY_CONVERGENCE"
    write_json(args.root / "integrity_audit.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
