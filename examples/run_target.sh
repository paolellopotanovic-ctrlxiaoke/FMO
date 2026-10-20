#!/usr/bin/env bash
set -euo pipefail
: "${PROTEIN:?Set PROTEIN to the receptor PDB/mmCIF}"
: "${LIGANDS:?Set LIGANDS to the label-free ligand CSV}"
: "${POSES:?Set POSES to receptor-frame pose SDF}"
: "${GAMESS_DIR:?Set GAMESS_DIR}"
: "${LIBRARY_DIR:?Set LIBRARY_DIR}"
: "${DFTB_PARAMETER_DIR:?Set DFTB_PARAMETER_DIR}"
fmotpiefinal blind-run \
  --target "${TARGET:-EXAMPLE}" --protein "$PROTEIN" \
  --ligands "$LIGANDS" --poses "$POSES" \
  --gamess-dir "$GAMESS_DIR" --library-dir "$LIBRARY_DIR" \
  --dftb-parameter-dir "$DFTB_PARAMETER_DIR" \
  --workers "${WORKERS:-1}" --ncores "${NCORES:-1}" \
  --output "${OUTPUT:-runs/example}"
