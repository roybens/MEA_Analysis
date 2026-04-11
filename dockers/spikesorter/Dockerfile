FROM nvidia/cuda:12.8.0-cudnn-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTORCH_ALLOC_CONF="expandable_segments:True"
ENV TORCH_CUDA_ARCH_LIST="8.6;8.9;9.0"
ENV CUDA_DEVICE_MAX_CONNECTIONS=32

# Maxwell HDF5 plugin path
ENV HDF5_PLUGIN_PATH=/opt/hdf5/plugins

# system deps
# python3-is-python provides the `python` command as an alias for python3
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-dev python3-pip python3-is-python \
    git wget build-essential \
    hdf5-tools libhdf5-dev \
    libsm6 libxext6 libxrender1 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# install Maxwell HDF5 compression plugin
RUN mkdir -p /opt/hdf5/plugins && \
    wget -O /opt/hdf5/plugins/libcompression.so \
    "https://share.mxwbio.com/d/7f2d1e98a1724a1b8b35/files/?p=%2FLinux%2Flibcompression.so&dl=1" && \
    chmod 755 /opt/hdf5/plugins/libcompression.so

# install Python deps
# Install PyTorch with CUDA support first; the default PyPI index only ships
# CPU-only wheels, so kilosort's transitive torch dep would silently run on CPU.
# torch 2.9.1 ships cu126 wheels (CUDA 12.6 is ABI-compatible with 12.8 drivers).
COPY dockers/spikesorter/requirements.txt /tmp/requirements.txt
RUN pip3 install --upgrade pip setuptools wheel && \
    pip3 install torch==2.9.1 --extra-index-url https://download.pytorch.org/whl/cu126 && \
    pip3 install -r /tmp/requirements.txt && \
    pip3 cache purge

# Bundle the repo code into the image at build time.
# For development, override with: docker run -v $(pwd):/MEA_Analysis ...
COPY . /MEA_Analysis

WORKDIR /MEA_Analysis

# Declare the canonical data mount point
VOLUME ["/data"]

COPY dockers/spikesorter/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
# Default: run the pipeline driver against /data.
# Override by passing a different command, e.g.:
#   docker run ... mea-analysis bash
#   docker run ... mea-analysis python3 IPNAnalysis/some_script.py
CMD ["python3", "IPNAnalysis/run_pipeline_driver.py", "/data"]