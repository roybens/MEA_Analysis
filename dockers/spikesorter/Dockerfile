FROM nvidia/cuda:12.8.0-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTORCH_ALLOC_CONF="expandable_segments:True"
ENV TORCH_CUDA_ARCH_LIST="8.6;8.9;9.0"
ENV CUDA_DEVICE_MAX_CONNECTIONS=32

# system deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-dev python3-pip \
    git wget build-essential \
    hdf5-tools libhdf5-dev \
    libsm6 libxext6 libxrender1 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /

# install Python deps
COPY dockers/spikesorter/requirements.txt /requirements.txt
RUN pip3 install --upgrade pip setuptools wheel && \
    pip3 install -r /requirements.txt && \
    pip3 cache purge

# baseline clone at build time
RUN git clone https://github.com/roybens/MEA_Analysis.git /MEA_Analysis

# runtime updater
COPY dockers/spikesorter/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
