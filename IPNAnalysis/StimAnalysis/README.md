# Stimulation Assay Analysis Pipeline

This folder contains Python tools and Jupyter notebooks for analyzing single-neuron stimulation experiments from multielectrode array (MEA) recordings. The analysis extracts spike times, filters artifacts, visualizes neuronal activity across different stimulation phases, and compares pre- and post-stimulation spike behavior.

## Folder Structure

```
StimAnalysis/
├── StimAnalysis.py                # Main class for running full spike analysis
├── stim_helper_functions.py      # Helper functions for filtering, plotting, waveform extraction
├── exp1_stim_analysis.ipynb      # Example notebook for trial 1
├── exp2_stim_analysis.ipynb      # Example notebook for trial 2
├── exp1106.ipynb                 # Additional experiments
├── exp1113.ipynb
├── royExp.ipynb
└── README.md                     # This file
```

## Requirements

Core dependencies include:
- spikeinterface
- numpy, scipy, matplotlib
- pandas, h5py
- scikit-learn

## Getting Started

### 1. Initialize the Analysis

```python
from StimAnalysis import StimulationAnalysis

analysis = StimulationAnalysis(
    file_path="data/exp1106.raw.h5",   # Path to your .h5 recording file
    stim_frequency=4,                  # Stimulation frequency in Hz
    recording_electrode=73,           # Electrode used for recording
    stim_electrode=141,               # Electrode used for stimulation
    artifact_electrode=96,            # Optional: electrode where artifact is visible
    spike_threshold=9,                # Spike detection threshold
    peak_sign="neg"                   # Options: "neg", "pos", or "both"
)
```

### 2. Run the Full Analysis Pipeline

```python
analysis.run_full_analysis(
    stim_start=60,                    # Time (in seconds) when stimulation begins
    stim_length=120                   # Duration (in seconds) of the stimulation period
)
```

This will:
- Detect spikes and remove artifact overlaps
- Plot stimulation and recording traces
- Plot spike counts across pre-, during-, and post-stimulation phases

## Features

- Artifact-aware spike detection
- Spike waveform extraction and t-SNE visualization
- Pre- vs. post-stimulation spike comparison
- ISI and fano factor analysis
- Modular, customizable methods for deeper analysis

## Custom Analysis

You can call additional methods such as:

```python
analysis.plot_stim_traces(...)
analysis.extract_spike_waveforms()
analysis.plot_spike_counts(...)
analysis.isi()
analysis.calculate_fano_factor(...)
```

## Author

Nathaniel Maffly  
UC Davis – MIND Institute
