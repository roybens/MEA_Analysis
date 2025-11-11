# --- Base image: CUDA 12.8 runtime for Ubuntu 22.04 (matches your driver stack) ---
FROM nvidia/cuda:12.8.0-runtime-ubuntu22.04

# --- System dependencies ---
RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-venv python3-dev \
    git build-essential hdf5-tools libhdf5-dev \
    libopenblas-dev liblapack-dev ffmpeg wget curl vim \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

# --- Copy your Python environment spec ---
COPY requirements.txt .

# --- Install dependencies (matches your si_ks4 venv) ---
RUN pip install --upgrade pip setuptools wheel \
 && pip install -r requirements.txt

# --- Environment configuration ---
ENV PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    CUDA_VISIBLE_DEVICES=0 \
    PATH="/usr/local/cuda/bin:${PATH}"

# --- Copy your analysis scripts ---
COPY . /workspace

# --- Default entry ---
CMD ["bash"]
