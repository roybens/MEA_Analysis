#!/bin/bash
set -euo pipefail

REPO_DIR="/MEA_Analysis"
BRANCH="main"


echo "=== Updating MEA_Analysis before run ==="
cd "$REPO_DIR"
git fetch origin "$BRANCH" --depth=1 || true
git reset --hard origin/"$BRANCH" || true


echo "=== Running pipeline ==="
exec python3 "$REPO_DIR/IPNAnalysis/run_pipeline_driver.py" "$@"