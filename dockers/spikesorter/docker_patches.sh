#`=/bin/bash
#missing modules in missing_modules.txt will be installed before running the pipeline
if [ -f "$REPO_DIR/dockers/spikesorter/missing_modules.txt" ]; then
    echo "Installing missing python modules..."
    pip install -r "$REPO_DIR/dockers/spikesorter/missing_modules.txt"
fi