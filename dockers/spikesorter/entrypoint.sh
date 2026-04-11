#!/bin/bash
# entrypoint.sh — lightweight container init
# Runs any docker_patches.sh present in the repo, then exec's the CMD.
# The CMD defaults to the pipeline driver but can be overridden:
#
#   docker run mea-analysis bash                              # interactive shell
#   docker run mea-analysis python3 IPNAnalysis/script.py    # custom script
#   docker run mea-analysis python3 IPNAnalysis/run_pipeline_driver.py /data --sorter kilosort4
#
set -euo pipefail

REPO_DIR="/MEA_Analysis"

echo "=== RUNNING PRE-PIPELINE HOOKS ==="
if [ -f "$REPO_DIR/dockers/spikesorter/docker_patches.sh" ]; then
    echo "Executing docker patches..."
    chmod +x "$REPO_DIR/dockers/spikesorter/docker_patches.sh"
    bash "$REPO_DIR/dockers/spikesorter/docker_patches.sh"
fi

exec "$@"