#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

if [[ -z "${PY:-}" && -x "$ROOT/envs/fmo87/bin/python" ]]; then
  PY="$ROOT/envs/fmo87/bin/python"
fi
PY="${PY:-python}"
if [[ -x "$ROOT/envs/fmo87/bin/python" ]]; then
  export PATH="$ROOT/envs/fmo87/bin:$PATH"
fi
WORK="${WORK:-$ROOT/runs/fmo87}"
STRUCTURE_ROOT="${STRUCTURE_ROOT:-$ROOT/data/fmo87_structures}"
GAMESS_DIR="${GAMESS_DIR:-$ROOT/engines/gamess}"
GAMESS_LIB="${GAMESS_LIB:-$ROOT/engines/gamess/runtime/lib}"
DFTB_PARAMETER_DIR="${DFTB_PARAMETER_DIR:-$GAMESS_DIR/auxdata/DFTB/3OB-3-1}"
BOLTZ_CACHE="${BOLTZ_CACHE:-$ROOT/engines/boltz-cache}"
BOLTZ_PY="${BOLTZ_PY:-$PY}"
MSA_DIR="${MSA_DIR:-$ROOT/data/boltz_msa}"
SOPHOSQM_SI="${SOPHOSQM_SI:-$ROOT/data/references/sophosqm_si.pdf}"
TIMEOUT="${TIMEOUT:-1800}"
BOOTSTRAP_REPETITIONS="${BOOTSTRAP_REPETITIONS:-10000}"
FMO_JOBS="${FMO_JOBS:-4}"
NO_WATER_WORK="${NO_WATER_WORK:-$ROOT/runs/fmo87/no_water_preparation}"
TARGETS=(CDK2 JNK1 P38 TYK2)

usage() {
  cat <<'USAGE'
Usage: ./run_pipeline.sh STAGE

Stages:
  verify          Run package regression tests and CLI checks
  check-inputs    Validate all frozen case SDF/PDB inputs and source hashes
  cohort          Freeze the exact 87-compound cohort
  prepare         Prepare CDK2/JNK1/P38/TYK2 proteins
  fmo             Minimize complexes and run GAMESS FMO2-DFTB3/PIEDA
  fmo-parallel    Launch FMO in FMO_JOBS resumable shards, then collect them
  boltz-inputs    Build Boltz-2 YAML inputs and manifest
  download-boltz  Download official Boltz-2 model artifacts into engines/boltz-cache
  boltz           Run Boltz-2 affinity inference
  audit           Strict-audit FMO and Boltz-2 outputs
  merge           Merge audited FMO and Boltz-2 outputs
  evaluate        Compute paper-comparison metrics and bootstrap CIs
  calibrate       Run SophosQM-style calibration ablation
  prepare-no-water Remove waters before PDB2PQR/PROPKA
  fmo-no-water    Re-minimize and run FMO with no explicit water
  audit-no-water  Strict-audit the full no-explicit-water branch
  evaluate-no-water Merge/evaluate the no-water branch against baseline Boltz-2
  no-water        Run all preparation-level no-water stages
  all             Run stages 1-9 sequentially

All benchmark inputs, MSA files, GAMESS, and its native runtime are bundled.
Run scripts/bootstrap_engines.sh first. Boltz-2 package installation is bundled
by that bootstrap; its official multi-GB checkpoints are downloaded by
download-boltz because individual weights exceed GitHub's file-size limit.
USAGE
}

need_env() {
  local name=$1
  if [[ -z "${!name:-}" ]]; then
    echo "error: environment variable $name is required for this stage" >&2
    exit 2
  fi
}

need_file() {
  if [[ ! -f "$1" ]]; then
    echo "error: required file does not exist: $1" >&2
    exit 2
  fi
}

stage_cohort() {
  need_env STRUCTURE_ROOT
  "$PY" "$ROOT/scripts/build_fmo87_cohort.py" \
    --structure-root "$STRUCTURE_ROOT" \
    --output "$ROOT/data/fmo87_cohort.csv"
}

stage_check_inputs() {
  need_env STRUCTURE_ROOT
  "$PY" "$ROOT/scripts/validate_inputs.py" \
    --cohort "$ROOT/data/fmo87_cohort.csv" \
    --structure-root "$STRUCTURE_ROOT" \
    --output "$WORK/input_validation.json"
}

stage_prepare() {
  mkdir -p "$WORK/prepared_proteins"
  "$PY" "$ROOT/scripts/prepare_protein.py" \
    --input "$STRUCTURE_ROOT/cdk2_protein.pdb" \
    --output "$WORK/prepared_proteins/CDK2_prepared.pdb"
  "$PY" "$ROOT/scripts/prepare_protein.py" \
    --input "$STRUCTURE_ROOT/jnk1_manual_flips_protein.pdb" \
    --output "$WORK/prepared_proteins/JNK1_prepared.pdb"
  "$PY" "$ROOT/scripts/prepare_protein.py" \
    --input "$STRUCTURE_ROOT/p38_protein.pdb" \
    --output "$WORK/prepared_proteins/P38_prepared.pdb"
  "$PY" "$ROOT/scripts/prepare_protein.py" \
    --input "$STRUCTURE_ROOT/tyk2_protein.pdb" \
    --output "$WORK/prepared_proteins/TYK2_prepared.pdb"
}

stage_fmo() {
  need_env GAMESS_DIR
  need_env GAMESS_LIB
  need_env DFTB_PARAMETER_DIR
  for target in "${TARGETS[@]}"; do
    "$PY" "$ROOT/scripts/run_fmo87.py" \
      --cohort "$ROOT/data/fmo87_cohort.csv" \
      --structure-root "$STRUCTURE_ROOT" \
      --prepared-proteins "$WORK/prepared_proteins" \
      --output "$WORK/fmo_rebuild/targets/$target" \
      --only-target "$target" \
      --gamess-dir "$GAMESS_DIR" \
      --library-dir "$GAMESS_LIB" \
      --dftb-parameter-dir "$DFTB_PARAMETER_DIR" \
      --timeout "$TIMEOUT" --stage all
  done
}

stage_fmo_parallel() {
  need_env GAMESS_DIR
  need_env GAMESS_LIB
  need_env DFTB_PARAMETER_DIR
  if (( FMO_JOBS < 1 )); then
    echo "error: FMO_JOBS must be >= 1" >&2
    exit 2
  fi
  mkdir -p "$WORK/fmo_parallel_logs"
  pids=()
  for index in $(seq 0 $((FMO_JOBS - 1))); do
    "$PY" "$ROOT/scripts/run_fmo87.py" \
      --cohort "$ROOT/data/fmo87_cohort.csv" \
      --structure-root "$STRUCTURE_ROOT" \
      --prepared-proteins "$WORK/prepared_proteins" \
      --output "$WORK/fmo_shards/$index" \
      --shard-index "$index" --shard-count "$FMO_JOBS" \
      --gamess-dir "$GAMESS_DIR" \
      --library-dir "$GAMESS_LIB" \
      --dftb-parameter-dir "$DFTB_PARAMETER_DIR" \
      --timeout "$TIMEOUT" --stage all \
      >"$WORK/fmo_parallel_logs/$index.log" 2>&1 &
    pids+=("$!")
  done
  failed=0
  for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then
      echo "error: FMO shard $index failed; see $WORK/fmo_parallel_logs/$index.log" >&2
      failed=1
    fi
  done
  if (( failed )); then exit 1; fi
  "$PY" "$ROOT/scripts/collect_fmo_shards.py" \
    --shard-root "$WORK/fmo_shards" \
    --output-root "$WORK/fmo_rebuild/targets" \
    --require-complete
}

stage_boltz_inputs() {
  need_env MSA_DIR
  "$PY" "$ROOT/scripts/build_boltz_fmo87_inputs.py" \
    --cohort "$ROOT/data/fmo87_cohort.csv" \
    --structure-root "$STRUCTURE_ROOT" \
    --output "$WORK/boltz_inputs" \
    --msa-dir "$MSA_DIR"
}

stage_download_boltz() {
  need_file "$ROOT/scripts/download_boltz_models.py"
  "$BOLTZ_PY" "$ROOT/scripts/download_boltz_models.py" \
    --cache "$BOLTZ_CACHE"
}

stage_boltz() {
  need_env BOLTZ_CACHE
  for target in "${TARGETS[@]}"; do
    "$BOLTZ_PY" "$ROOT/scripts/run_boltz_fmo87.py" \
    --manifest "$WORK/boltz_inputs/manifest.json" \
      --out-dir "$WORK/boltz_predictions/$target" \
      --cache "$BOLTZ_CACHE" \
      --python "$BOLTZ_PY" \
      --only-target "$target" \
      --accelerator gpu --devices 1 --seed 0 \
      --recycling-steps 3 \
      --sampling-steps-affinity 200 \
      --diffusion-samples-affinity 5 \
      --timeout "$TIMEOUT"
  done
}

stage_audit() {
  "$PY" "$ROOT/scripts/audit_fmo87_production.py" \
    --root "$WORK/fmo_rebuild/targets" \
    --prepared-proteins "$WORK/prepared_proteins" \
    --output-dir "$WORK/final" --require-complete
  "$PY" "$ROOT/scripts/audit_boltz87_production.py" \
    --root "$WORK/boltz_predictions" \
    --cache "$BOLTZ_CACHE" \
    --output-dir "$WORK/final" --require-complete
}

stage_prepare_no_water() {
  need_env STRUCTURE_ROOT
  mkdir -p "$NO_WATER_WORK/prepared_proteins"
  for spec in \
    "CDK2 cdk2_protein.pdb" \
    "JNK1 jnk1_manual_flips_protein.pdb" \
    "P38 p38_protein.pdb" \
    "TYK2 tyk2_protein.pdb"; do
    read -r target input <<<"$spec"
    "$PY" "$ROOT/scripts/prepare_protein.py" \
      --input "$STRUCTURE_ROOT/$input" \
      --output "$NO_WATER_WORK/prepared_proteins/${target}_prepared.pdb" \
      --drop-water
  done
}

stage_fmo_no_water() {
  need_env GAMESS_DIR
  need_env GAMESS_LIB
  need_env DFTB_PARAMETER_DIR
  for target in "${TARGETS[@]}"; do
    "$PY" "$ROOT/scripts/run_fmo87.py" \
      --cohort "$ROOT/data/fmo87_cohort.csv" \
      --structure-root "$STRUCTURE_ROOT" \
      --prepared-proteins "$NO_WATER_WORK/prepared_proteins" \
      --output "$NO_WATER_WORK/fmo_rebuild/targets/$target" \
      --only-target "$target" \
      --explicit-water-policy removed_before_preparation \
      --gamess-dir "$GAMESS_DIR" \
      --library-dir "$GAMESS_LIB" \
      --dftb-parameter-dir "$DFTB_PARAMETER_DIR" \
      --timeout "$TIMEOUT" --stage all
  done
}

stage_audit_no_water() {
  "$PY" "$ROOT/scripts/audit_fmo87_production.py" \
    --root "$NO_WATER_WORK/fmo_rebuild/targets" \
    --prepared-proteins "$NO_WATER_WORK/prepared_proteins" \
    --output-dir "$NO_WATER_WORK/audit" \
    --require-complete \
    --explicit-water-policy removed_before_preparation
}

stage_evaluate_no_water() {
  fmo_targets=()
  for target in "${TARGETS[@]}"; do
    fmo_targets+=("$NO_WATER_WORK/fmo_rebuild/targets/$target")
  done
  "$PY" "$ROOT/scripts/merge_fmo87_outputs.py" \
    --fmo-targets "${fmo_targets[@]}" \
    --boltz-targets \
      "$WORK/boltz_predictions/CDK2" \
      "$WORK/boltz_predictions/JNK1" \
      "$WORK/boltz_predictions/P38" \
      "$WORK/boltz_predictions/TYK2" \
    --output-dir "$NO_WATER_WORK/evaluation_inputs"
  "$PY" "$ROOT/scripts/evaluate_fmo87.py" \
    --fmo "$NO_WATER_WORK/evaluation_inputs/fmo_features.csv" \
    --boltz "$NO_WATER_WORK/evaluation_inputs/boltz_predictions.csv" \
    --boltz-input-manifest "$WORK/boltz_inputs/manifest.json" \
    --output-dir "$NO_WATER_WORK/evaluation" \
    --bootstrap-repetitions "$BOOTSTRAP_REPETITIONS"
}

stage_merge() {
  fmo_targets=()
  boltz_targets=()
  for target in "${TARGETS[@]}"; do
    fmo_targets+=("$WORK/fmo_rebuild/targets/$target")
    boltz_targets+=("$WORK/boltz_predictions/$target")
  done
  "$PY" "$ROOT/scripts/merge_fmo87_outputs.py" \
    --fmo-targets "${fmo_targets[@]}" \
    --boltz-targets "${boltz_targets[@]}" \
    --output-dir "$WORK/final"
}

stage_evaluate() {
  "$PY" "$ROOT/scripts/evaluate_fmo87.py" \
    --fmo "$WORK/final/fmo_features.csv" \
    --boltz "$WORK/final/boltz_predictions.csv" \
    --boltz-input-manifest "$WORK/boltz_inputs/manifest.json" \
    --output-dir "$WORK/final" \
    --bootstrap-repetitions "$BOOTSTRAP_REPETITIONS"
}

stage_calibrate() {
  need_file "$SOPHOSQM_SI"
  "$PY" "$ROOT/scripts/ablate_sophosqm_calibration.py" \
    --compounds "$WORK/final/compounds.csv" \
    --sophosqm-si "$SOPHOSQM_SI" \
    --output-dir "$WORK/final"
}

stage_verify() {
  "$PY" -m compileall -q "$ROOT/src" "$ROOT/scripts" "$ROOT/tests"
  (cd "$ROOT" && "$PY" -m pytest -q)
  for script in "$ROOT"/scripts/*.py; do
    "$PY" "$script" --help >/dev/null
  done
  need_file "$ROOT/data/fmo87_cohort.csv"
  need_file "$ROOT/reference_results/final/metrics.json"
  need_file "$ROOT/reference_results/final/production_audit.json"
  need_file "$ROOT/reference_results/final/boltz_production_audit.json"
  "$PY" - "$ROOT" <<'PY'
import csv
import hashlib
import json
import sys
from pathlib import Path

import yaml

root = Path(sys.argv[1])
config = yaml.safe_load((root / "configs/fmo87_production.yaml").read_text())
assert config["benchmark"]["targets"] == {"CDK2": 16, "JNK1": 21, "P38": 34, "TYK2": 16}
assert config["benchmark"]["cohort_sha256"] == "19edfb8e2ccf3ab1029d848e4062d70cb35532d2b4fee63674f160b2a292997b"
assert config["protein_preparation"] == {
    "pdb2pqr": "3.6.1", "forcefield": "AMBER", "titration": "PROPKA 3.5.1",
    "ph": 7.4, "hydrogen_bond_optimization": True, "water_policy": "retain_crystal_waters",
}
assert config["geometry"]["ligand_restraint_k_kj_mol_angstrom2"] == 100
assert config["geometry"]["protein_coordinates"] == "free"
assert config["fmo"]["pocket_cutoff_angstrom"] == 5.0
assert config["fmo"]["method"] == "FMO2-DFTB3/PIEDA"
assert config["fmo"]["fragmentation"] == "shifted_facio_like_ca_c_prime"
assert config["fmo"]["terminal_caps"] == "full_neutral_ACE_NME"
assert config["fmo"]["internal_boundaries"] == "FMOBND_HOP_C"
assert config["boltz2"]["seed"] == 0
assert config["boltz2"]["recycling_steps"] == 3
assert config["boltz2"]["sampling_steps_affinity"] == 200
assert config["boltz2"]["diffusion_samples_affinity"] == 5
cohort_path = root / "data/fmo87_cohort.csv"
with cohort_path.open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle))
expected = {"CDK2": 16, "JNK1": 21, "P38": 34, "TYK2": 16}
actual = {target: sum(row["target"] == target for row in rows) for target in expected}
assert len(rows) == 87 and actual == expected
assert all(int(row["formal_charge"]) == 0 for row in rows)
assert hashlib.sha256(cohort_path.read_bytes()).hexdigest() == "19edfb8e2ccf3ab1029d848e4062d70cb35532d2b4fee63674f160b2a292997b"
for name in ("production_audit.json", "boltz_production_audit.json"):
    audit = json.loads((root / "reference_results/final" / name).read_text())
    assert audit["strict_complete"] is True
    assert audit["n_compounds"] == 87 and audit["n_passed"] == 87
metrics = json.loads((root / "reference_results/final/metrics.json").read_text())["metrics"]
assert metrics["paper_reference"]["FMO"] == {"pearson_r": 0.55, "kendall_tau": 0.38}
assert metrics["paper_reference"]["Boltz-2"]["pearson_r"] == 0.66
assert metrics["paper_reference"]["Boltz-2"]["kendall_tau"] == 0.48
PY
  echo "package verification passed"
}

case "${1:-help}" in
  cohort) stage_cohort ;;
  check-inputs) stage_check_inputs ;;
  prepare) stage_prepare ;;
  fmo) stage_fmo ;;
  fmo-parallel) stage_fmo_parallel ;;
  boltz-inputs) stage_boltz_inputs ;;
  download-boltz) stage_download_boltz ;;
  boltz) stage_boltz ;;
  audit) stage_audit ;;
  merge) stage_merge ;;
  evaluate) stage_evaluate ;;
  calibrate) stage_calibrate ;;
  verify) stage_verify ;;
  prepare-no-water) stage_prepare_no_water ;;
  fmo-no-water) stage_fmo_no_water ;;
  audit-no-water) stage_audit_no_water ;;
  evaluate-no-water) stage_evaluate_no_water ;;
  no-water)
    stage_check_inputs
    stage_prepare_no_water
    stage_fmo_no_water
    stage_audit_no_water
    stage_evaluate_no_water
    ;;
  all)
    stage_check_inputs
    stage_cohort
    stage_prepare
    stage_fmo
    stage_boltz_inputs
    stage_download_boltz
    stage_boltz
    stage_audit
    stage_merge
    stage_evaluate
    stage_calibrate
    ;;
  help|-h|--help) usage ;;
  *) usage; exit 2 ;;
esac
