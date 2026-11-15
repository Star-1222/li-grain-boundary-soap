#!/usr/bin/env bash
# setup_soap_env.sh -- create and verify the Python environment for the SOAP
# analysis package.
#
# [NEW 2026-08-21] Part of the 09_SOAP_Expansion_20260821 deliverable.
#
# Usage (from the 09_SOAP_Expansion_20260821 folder):
#     bash environment/setup_soap_env.sh [VENV_PATH]
#
# VENV_PATH defaults to ../.soap_venv (i.e. next to this deliverable folder,
# inside the project root, so it is not nested in the deliverable itself).
#
# What it does:
#   1. creates a virtualenv (skipped if it already exists);
#   2. installs requirements_soap.txt;
#   3. runs verify_soap_env.py inside the new environment;
#   4. optionally runs the full pytest suite when RUN_TESTS=1.
#
# Environment variables:
#   PYTHON      interpreter to build the venv from   (default: python3 or python)
#   VENV_PATH   same as the positional argument      (default: ../.soap_venv)
#   RUN_TESTS   1 -> also run  python -m pytest tests/  (default: 0)

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"

VENV_PATH="${1:-${VENV_PATH:-$ROOT/../.soap_venv}}"
PYTHON="${PYTHON:-}"
if [ -z "$PYTHON" ]; then
    if command -v python3 >/dev/null 2>&1; then PYTHON=python3; else PYTHON=python; fi
fi

echo "[setup] root      : $ROOT"
echo "[setup] venv path : $VENV_PATH"
echo "[setup] python    : $($PYTHON --version 2>&1)"

# --- 1. create venv ---------------------------------------------------------
if [ -f "$VENV_PATH/pyvenv.cfg" ]; then
    echo "[setup] venv already exists, reusing it"
else
    "$PYTHON" -m venv "$VENV_PATH"
    echo "[setup] created venv"
fi

# --- 2. activate + install --------------------------------------------------
if [ -f "$VENV_PATH/Scripts/activate" ]; then
    # shellcheck disable=SC1091
    source "$VENV_PATH/Scripts/activate"     # Windows (Git Bash / MSYS)
else
    # shellcheck disable=SC1091
    source "$VENV_PATH/bin/activate"         # Linux / macOS
fi

python -m pip install --upgrade pip
python -m pip install -r "$HERE/requirements_soap.txt"

# --- 3. verify --------------------------------------------------------------
python "$HERE/verify_soap_env.py"

# --- 4. optional test suite -------------------------------------------------
if [ "${RUN_TESTS:-0}" = "1" ]; then
    echo "[setup] running full test suite"
    (cd "$ROOT" && python -m pytest tests/ -v)
fi

echo "[setup] done. Activate the environment with:"
if [ -f "$VENV_PATH/Scripts/activate" ]; then
    echo "    source $VENV_PATH/Scripts/activate"
else
    echo "    source $VENV_PATH/bin/activate"
fi
