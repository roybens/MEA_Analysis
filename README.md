# MEA_Analysis

End-to-end pipeline for neuronal spike sorting and network burst analysis on **Maxwell Biosystems MEA** recordings.

Built on [SpikeInterface](https://github.com/SpikeInterface/spikeinterface) with Kilosort4 as the default sorter.

---

## Repository Layout

```
MEA_Analysis/
├── IPNAnalysis/          ← Main production pipeline (start here)
├── NetworkAnalysis/      ← MATLAB-based network analysis tools
├── NeuronClassification/ ← UMAP-based neuron classification
├── Archive/              ← Legacy scripts (not maintained)
├── MEAProcessingLibrary/ ← Installable Python processing library
├── MaxwellBiosystemsDeviceInterface/  ← Device control scripts
├── StimulationAnalysis/  ← Stimulation experiment analysis
└── dockers/              ← Docker build configs for spike sorters
```

## Getting Started

The active pipeline lives in [`IPNAnalysis/`](IPNAnalysis/README.md).  
See the **[IPNAnalysis README](IPNAnalysis/README.md)** for:

- Architecture overview (driver → routine → config_loader)
- Pipeline stages (preprocessing, sorting, analyzer, reports)
- Configuration system (`mea_config.json`)
- Typical workflows and full CLI reference

### Quick Start

```bash
# Generate a config template
python IPNAnalysis/config_loader.py mea_config.json

# Dry run on a directory
python IPNAnalysis/run_pipeline_driver.py /data/experiment --config mea_config.json --dry

# Full batch run
python IPNAnalysis/run_pipeline_driver.py /data/experiment --config mea_config.json

# Single well
python IPNAnalysis/mea_analysis_routine.py /data/exp/run_001/Network/data.raw.h5 \
  --well well000 --rec rec0001 --config mea_config.json
```

## Dependencies

- Python ≥ 3.9
- [SpikeInterface](https://github.com/SpikeInterface/spikeinterface)
- Kilosort4 (GPU recommended, ≥ 8 GB VRAM)
- See `requirements.txt` for full list
  - use pip install requirements.txt for install dependencies

## Docker

The repo provides a Docker image (CUDA 12, SpikeInterface, Kilosort4, Maxwell HDF5 plugin) for containerised spike sorting and full-pipeline runs.

### Build the image

From the repo root (requires [Docker](https://docs.docker.com/get-docker/) and an NVIDIA GPU with drivers + [nvidia-docker](https://github.com/NVIDIA/nvidia-docker)):

```bash
docker build -t mea-spikesorter -f dockers/spikesorter/Dockerfile .
```

**Note:** The `dockers/spikesorter/Dockerfile` clones `roybens/MEA_Analysis` from GitHub during the build. This means the image is built from the remote repository (typically the latest `main`), **not** from your local working tree or branch. Local, unpushed changes in this checkout will not be included in the image. To build an image from your own fork or branch, adjust the clone URL/branch in the Dockerfile before running the `docker build` command above.

### Run the pipeline in Docker

Mount your data and config, and pass pipeline arguments after the image name. The container entrypoint runs `run_pipeline_driver.py` and will attempt to update the repo to the latest `main` before each run; if the update fails (for example, when offline or if the remote is unavailable), it continues using the existing checkout.

```bash
docker run --gpus all -it --rm \
  -v /path/to/your/experiment:/data/experiment \
  -v /path/to/your/output:/output \
  -v /path/to/mea_config.json:/config/mea_config.json:ro \
  mea-spikesorter /data/experiment --config /config/mea_config.json --output-dir /output
```

Adjust `-v` paths and the positional `/data/experiment` and `--output-dir` to match your layout. Use `--dry` to preview without processing.

### Use a Docker image for sorting only

When running the pipeline on the host (not inside the container), you can run the spike-sorting step inside Docker by passing the image name:

```bash
python IPNAnalysis/run_pipeline_driver.py /data/experiment --config mea_config.json --docker mea-spikesorter
```

Or set `"docker_image": "mea-spikesorter"` in the `sorting` section of `mea_config.json`.

---

### Docker basics (Ubuntu)

Use these steps to try a Docker environment with an Ubuntu image.

**Pull the Ubuntu image:**

```bash
docker pull ubuntu
```

(On Linux you may need `sudo docker pull ubuntu`. For a specific version: `docker pull ubuntu:20.04`.)

**Run an Ubuntu container and get a shell:**

The `-it` options run the container in interactive mode with a TTY so you can use the container’s shell:

```bash
docker run -it ubuntu
```

Or explicitly request a shell:

```bash
docker run -it --rm ubuntu /bin/bash
```

`--rm` removes the container when you exit.

---

### Installing Python, pip, and sudo in an Ubuntu container

Inside the Ubuntu container:

```bash
apt update && apt upgrade -y
apt install -y sudo
apt install -y python3-pip
pip install spikeinterface[full]
```

---

### Creating a new image from a container

If you install packages in a running container and want to save that state as a new image:

1. Exit the container: `exit`
2. List containers: `docker ps -a`
3. Copy the **container ID** of the one you just exited
4. Create an image from it:  
   `docker commit [container-id] [image-name]`  
   (e.g. `docker commit abc123 ubuntu-rohan`)
5. Run the new image:  
   `docker run -it ubuntu-rohan /bin/bash`
6. Inspect installed packages: `pip freeze`

---

### Creating a Dockerfile

Example Dockerfile that builds an image with SpikeInterface:

```dockerfile
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y sudo python3-pip

RUN python3 -m pip install spikeinterface

CMD ["bash"]
```

From the folder where the Dockerfile is saved:

```bash
docker build -t ubuntu-spikeinterface .
docker run -it --rm ubuntu-spikeinterface
```

---

### Running a Docker image on the lab server

Run the spike-sorter image with GPU access, mounted lab disks, and memory limits so Kilosort runs reliably. Opens an interactive bash shell; use it to run the pipeline manually or for debugging.

```bash
docker run -it --gpus all --rm --entrypoint bash \
  --shm-size=1g --ulimit memlock=-1 --ulimit stack=67108864 \
  -v /mnt/disk15tb/:/mnt/disk15tb -v /mnt/disk20tb:/mnt/disk20tb \
  -p 8884:8884 \
  mandarmp/benshalom:v1
```

- `--gpus all`: use all GPUs; `--entrypoint bash`: get a shell instead of running the pipeline entrypoint.
- `--shm-size=1g`, `memlock=-1`, `stack=67108864`: increase shared memory and stack limits for the sorter.
- `-v`: mount lab data disks; `-p 8884:8884`: expose port (e.g. for Jupyter). Replace image name if you use a different one.

---

## License

Pending

