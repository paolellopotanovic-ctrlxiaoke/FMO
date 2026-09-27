#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_DIR="${FMO87_ENV_DIR:-$ROOT/envs/fmo87}"
PYTHON="${FMO87_PYTHON:-}"

if [[ -z "$PYTHON" ]]; then
  for candidate in python3.12 python3.11 python3; do
    if command -v "$candidate" >/dev/null 2>&1; then
      PYTHON="$candidate"
      break
    fi
  done
fi
[[ -n "$PYTHON" ]] || { echo "error: no Python interpreter found" >&2; exit 2; }

if [[ ! -x "$ENV_DIR/bin/python" ]]; then
  "$PYTHON" -m venv "$ENV_DIR"
fi
"$ENV_DIR/bin/python" -m pip install --upgrade pip wheel
"$ENV_DIR/bin/python" -m pip install -r "$ROOT/requirements-lock.txt"
"$ENV_DIR/bin/python" -m pip install -e "$ROOT/engines/boltz[cuda]"
"$ENV_DIR/bin/python" - "$ROOT" <<'PY'
from importlib.metadata import version
from pathlib import Path
import sys

import boltz
import openmm
import rdkit
import scipy
import yaml

assert version("boltz") == "2.2.1"
assert (Path(sys.argv[1]) / "engines" / "boltz" / "SOURCE_COMMIT").read_text().strip() == "b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc"
print("engine imports passed", version("boltz"), openmm.__version__, rdkit.__version__, scipy.__version__)
PY

echo "FMO87 Python environment is ready: $ENV_DIR/bin/python"
echo "Bundled GAMESS: $ROOT/engines/gamess/rungms"
echo "Bundled Boltz source commit: b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc"
