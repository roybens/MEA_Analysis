#!/bin/bash
set -euo pipefail

REPO_DIR="/MEA_Analysis"
BRANCH="main"


echo "=== Updating MEA_Analysis before run ==="
cd "$REPO_DIR"
git fetch origin "$BRANCH" --depth=1 || true
git reset --hard origin/"$BRANCH" || true

echo "=== RUNNING PRE-PIPELINE HOOKS ==="
#excute docker patches before running the pipeline
if [ -f "$REPO_DIR/dockers/spikesorter/docker_patches.sh" ]; then
    echo "Executing docker patches..."
    chmod +x "$REPO_DIR/dockers/spikesorter/docker_patches.sh"
    bash "$REPO_DIR/dockers/spikesorter/docker_patches.sh"
fi


echo "=== Running pipeline ==="
exec python3 "$REPO_DIR/IPNAnalysis/run_pipeline_driver.py" "$@"