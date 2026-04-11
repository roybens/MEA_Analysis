#!/usr/bin/env bash
# build.sh — convenience wrapper to build the MEA Analysis Docker image.
#
# Usage:
#   ./build.sh [TAG]
#
# Examples:
#   ./build.sh                  # builds as mea-analysis:latest
#   ./build.sh mea-analysis:v1  # builds with a custom tag
#
# The build context must be the repository root (this script enforces that).
# Pass DATA_DIR to override the default /data mount when running afterward:
#   DATA_DIR=/mnt/lab-data docker compose run --rm mea-analysis

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TAG="${1:-mea-analysis:latest}"

echo "=== Building MEA Analysis image: ${TAG} ==="
echo "    Build context : ${REPO_ROOT}"
echo "    Dockerfile    : dockers/spikesorter/Dockerfile"
echo ""

docker build \
  --file "${REPO_ROOT}/dockers/spikesorter/Dockerfile" \
  --tag  "${TAG}" \
  "${REPO_ROOT}"

echo ""
echo "=== Build complete: ${TAG} ==="
echo ""
echo "Quick start:"
echo "  Run pipeline:    docker run --gpus all -v /your/data:/data ${TAG}"
echo "  Interactive shell: docker run --gpus all -v /your/data:/data -it ${TAG} bash"
echo ""
echo "Or use docker compose:"
echo "  DATA_DIR=/your/data docker compose run --rm mea-analysis"
