# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

**MEA_Analysis** is an end-to-end neuronal spike sorting and network burst analysis pipeline for Maxwell Biosystems microelectrode array (MEA) recordings. It processes raw HDF5 data through preprocessing, spike sorting (Kilosort4), waveform analysis, quality curation, and network burst detection.

**Status**: Active production pipeline. Primary development on `main` and `dev_branch`. Current branch: `parameter-detection`.

**Key Authors**: Roy Benshaalom Lab members (Adam Weiner, Mandar Patil, Yuxin Ren, others).

## Architecture Overview

### Two-Tier Processing Model

The pipeline uses a **driver-worker pattern**:

1. **run_pipeline_driver.py** (Orchestrator)
   - Entry point for batch processing
   - Scans directories or single HDF5 files
   - Discovers all recordings and wells in each file
   - Maps recordings to wells using `recording_map` dictionary
   - Launches subprocess workers via `mea_analysis_routine.py` for each well
   - Supports `--dry` mode to preview without processing
   - Handles reference file filtering by assay type

2. **mea_analysis_routine.py** (Core Pipeline)
   - Contains `MEAPipeline` class: the main processing engine
   - Processes a single well (identified by `--well` and `--rec` flags)
   - Can run independently or be called by driver
   - **4-Stage Pipeline**:
     - **Preprocessing**: Highpass filter (300 Hz), local common median reference, float32 conversion, binary cache
     - **Spike Sorting**: Kilosort4 (default), SpikeInterface integration, Docker support
     - **Analyzer**: Template computation, quality metrics (firing rate, presence ratio, ISI violations, amplitude)
     - **Reports**: Waveform visualizations, probe locations, burst analysis, automatic curation
   - Checkpoint-based resumption: saves JSON state per well, resumes on crash, skips completed stages
   - Can re-analyze bursts on existing spike times via `--reanalyze-bursts`

3. **config_loader.py** (Shared Configuration)
   - Three-level priority: CLI flag → config file → hardcoded defaults
   - `build_extra_args()`: constructs subprocess argument strings for driver
   - Run directly to generate config template: `python config_loader.py mea_config.json`

### Data Flow

```
Input File (data.raw.h5)
  ├─ Multiple Recordings (rec0001, rec0002, ...)
  └─ Multiple Wells per Recording (well000-well005)
       │
       ├─ [Preprocessing] → binary/ (cached preprocessed recording)
       ├─ [Sorting] → sorter_output/ (Kilosort4 results)
       ├─ [Analyzer] → analyzer_output/ (waveforms, templates, QC metrics)
       └─ [Reports] → raster plots, burst statistics, curated units
```

### Key Classes and Modules

**MEAPipeline** (`mea_analysis_routine.py`)
- `__init__()`: Metadata parsing, checkpoint loading, logger setup
- `_parse_metadata()`: Extracts project/date/chip/run from file paths
- `run_preprocessing()`: Loads recording, filters, caches to binary
- `run_sorting()`: Runs Kilosort4, removes excess spikes
- `run_analyzer()`: Computes templates, waveforms, quality metrics
- `generate_reports()`: Plots, curates units, exports to Phy
- `_run_burst_analysis()`: Network burst detection on spike times
- Checkpoint methods: `_load_checkpoint()`, `_save_checkpoint()`, `should_skip()`

**Burst Detection** (`parameter_free_burst_detector.py`)
- `compute_network_bursts()`: Parameter-free adaptive burst detection
  - Biological calibration from ISI distributions
  - Per-unit ISI bursts + population firing rate signal
  - Gaussian smoothing + adaptive thresholding
  - Synchrony metrics + network burst merging
  - Returns burst intervals, participation metrics, spike counts

**Plotting & Analysis** (`meaplotter.py`, `helper_functions.py`)
- `MEAPlotter`: Publication-grade figure generation with configurable styling
- Helper functions: peak detection, raster plotting, network burst visualization
- Supports multiple plot modes (`separate`, `merged`) and raster sorting options

### Configuration System

**Priority Chain**: CLI flag > config file > hardcoded default

**Config Sections** (in `mea_config.json`):
```json
{
  "io": {
    "output_dir": null,
    "checkpoint_dir": null,
    "export_to_phy": false,
    "clean_up": false
  },
  "sorting": {
    "sorter": "kilosort4",
    "docker_image": null
  },
  "filtering": {
    "reference_file": null,
    "assay_types": ["network today", "network today/best"]
  },
  "plotting": {
    "plot_mode": "separate",
    "raster_sort": "none",
    "plot_debug": false
  },
  "curation": {
    "no_curation": false,
    "quality_thresholds": {
      "presence_ratio": 0.75,
      "rp_contamination": 0.15,
      "firing_rate": 0.05,
      "amplitude_median": -20,
      "amplitude_cv_median": 0.5
    }
  }
}
```

**CLI-Only Flags** (never in config):
- `--well`, `--rec`: Per-file identifiers
- `--debug`, `--dry`, `--force-restart`, `--reanalyze-bursts`, `--skip-spikesorting`: Runtime decisions

## Supporting Modules

### MEAProcessingLibrary
- Under construction (`setup.py` present but minimal)
- Intended as reusable installable package for spike processing utilities
- Currently mostly integration layer around SpikeInterface

### StimulationAnalysis
- Single-neuron stimulation experiment analysis
- Artifact-aware spike detection, spike waveform extraction
- Pre- vs post-stimulation spike comparison
- Main class: `StimulationAnalysis`

### Connections
- Network connectivity analysis using Lasso-based models (pyUOI)
- `classifier.py`: Unit classification
- `UoI_Lasso*.py`: Various Union of Intersections implementations
- `plot_connections.py`: Connection visualization

### NetworkAnalysis
- MATLAB-based legacy tools for network analysis
- **InHouseWebBasedGrapher**: Dash-based interactive visualization
- **PlottingFunctions**: Activity and burst property plotting
- **ActivityQuality**: Entropy and LDA-based analysis

### MaxwellBiosystemsDeviceInterface
- Hardware control for Maxwell Biosystems MEA systems
- Real-time stimulation and oscilloscope visualization
- C++ streamer (`mxw_streamer.cpp`) for live data acquisition
- Python-only alternative: `recorder.py` + `emulator.py`
- See `HOWTO.md` for detailed setup and calibration workflow

### Organoid & WildtypeSegregation
- Specialized analysis notebooks and MATLAB scripts
- Organoid-specific MEA analysis tools

### GUI
- MATLAB-based GUI for analysis control
- Reference manual and parameter exploration tools

## Common Workflows

### 1. Generate Config Template
```bash
cd IPNAnalysis
python config_loader.py /path/to/mea_config.json
# Edit the file with your settings
```

### 2. Dry Run (Preview Without Processing)
```bash
python IPNAnalysis/run_pipeline_driver.py /data/experiment \
  --config mea_config.json \
  --dry
```

### 3. Full Batch Processing
```bash
python IPNAnalysis/run_pipeline_driver.py /data/experiment \
  --config mea_config.json \
  --output-dir /results
```

### 4. Single Well Processing
```bash
python IPNAnalysis/mea_analysis_routine.py /data/exp/run_001/Network/data.raw.h5 \
  --well well000 \
  --rec rec0001 \
  --config mea_config.json
```

### 5. Re-analyze Bursts on Existing Spike Times
```bash
python IPNAnalysis/mea_analysis_routine.py /path/to/data.raw.h5 \
  --well well000 \
  --config mea_config.json \
  --reanalyze-bursts \
  --plot-mode merged \
  --raster-sort firing_rate
```

### 6. Override Config at Runtime
```bash
python IPNAnalysis/run_pipeline_driver.py /data/experiment \
  --config mea_config.json \
  --sorter mountainsort5 \
  --plot-debug \
  --no-curation
```

### 7. Skip Spike Sorting (Detection Only)
```bash
python IPNAnalysis/run_pipeline_driver.py /data/experiment \
  --config mea_config.json \
  --skip-spikesorting
```

### 8. Resume After Crash
```bash
# Just re-run the same command — checkpoints handle resumption automatically
python IPNAnalysis/run_pipeline_driver.py /data/experiment \
  --config mea_config.json

# Force full restart ignoring checkpoints
python IPNAnalysis/run_pipeline_driver.py /data/experiment \
  --config mea_config.json \
  --force-restart
```

## HPC Batch Submission

### Single Well Job
```bash
sbatch sbatch.sh
# Runs with 1 GPU, 32 CPUs, 1.5h timeout
# Edit INPUT_PATH, OUTPUT_ROOT in the script
```

### Parallel 4-Well Job (Maxwell MaxTwo plates)
```bash
sbatch sbatch_parallel.sh /path/to/run_directory
# Launches 4 wells in parallel on 4 GPUs
# Auto-detects .raw.h5 file
# Each well gets 1 GPU + 8 CPUs
```

### Re-analyze Bursts in Batch
```bash
sbatch sbatch_parallel_reanalyzebursts.sh /path/to/output_dir
```

**Environment Variables Set in SBATCH Scripts**:
- `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`: Prevents OOM with large batches
- `HDF5_PLUGIN_PATH=/pscratch/sd/m/mpatil1/hdf5_plugin`: Maxwell HDF5 compression plugin
- `OMP_NUM_THREADS`: Parallel processing threads per task

## Docker Support

**Containerized Spike Sorting** (for reproducibility)

```bash
# Build image
cd dockers/spikesorter
docker build -t mandarmp/benshalomlab_spikesorter:latest .

# Use in SBATCH
sbatch --image=mandarmp/benshalomlab_spikesorter:latest sbatch.sh

# Or pass via config
python run_pipeline_driver.py /data/experiment \
  --config mea_config.json \
  --docker mandarmp/benshalomlab_spikesorter:latest
```

**Dockerfile**:
- Base: `nvidia/cuda:12.8.0-runtime-ubuntu22.04`
- Installs SpikeInterface, Kilosort4, PyTorch, HDF5 plugin
- Entrypoint: `entrypoint.sh` (pulls latest MEA_Analysis code at runtime)

## Output Directory Structure

```
<output_dir>/<project>/<date>/<chip>/<run_id>/<well_id>/
├── binary/                          # Preprocessed recording cache
├── sorter_output/                   # Kilosort4 results
├── analyzer_output/                 # Waveforms, templates, QC metrics
├── checkpoints/                     # Resume state JSON per well
├── raster_burst_plot.svg            # Full recording raster + network burst
├── raster_burst_plot_30s.svg        # 30s zoom
├── raster_burst_plot_60s.svg        # 60s zoom
├── network_results.json             # Burst statistics, spike counts, participation
├── spike_times.npy                  # Spike times per unit (dict)
├── qm_unfiltered.xlsx               # Quality metrics (all units)
├── qm_curated.xlsx                  # Quality metrics (after curation)
├── tm_unfiltered.xlsx               # Template metrics
├── rejection_log.xlsx               # Rejected units and rejection reasons
├── locations_unfiltered.pdf         # All unit probe locations
├── waveforms_grid.pdf               # Waveform overview grid
├── run_id_well_id_pipeline.log      # Per-well processing log
└── phy_output/ (if --export-to-phy) # Phy-compatible spike viewer format
```

**Key Output Files**:
- `network_results.json`: Burst intervals, firing rates, network metrics (needed for downstream analysis)
- `spike_times.npy`: Dict of unit IDs → spike time arrays (seconds)
- `*_curated.xlsx`: Filtered units post-quality thresholding (for publication)

## Dependencies

### Core
- **spikeinterface** (0.103.0): Spike sorting and waveform analysis framework
- **kilosort** (4.1.1): Default spike sorter (GPU-accelerated)
- **neo** (0.14.3): Data I/O and standardization
- **h5py** (3.13.0): HDF5 file reading (Maxwell format)

### Scientific Computing
- numpy, scipy, scikit-learn, pandas
- matplotlib, seaborn: Plotting
- PyTorch (2.9.1): GPU acceleration for Kilosort4

### Hardware & Visualization
- PyQt6, PySide6: GUI frameworks
- pyqtgraph: Fast plotting for real-time display
- umap-learn: Dimensionality reduction for neuron classification
- wandb: Experiment tracking (optional)

### Development Notes
- Python ≥ 3.9 required
- CUDA 12.8 recommended (GPU required for Kilosort4)
- Minimum 8 GB VRAM for Kilosort4

## Testing & Validation

**No formal test suite exists yet.** Validation is done via:
1. Dry runs: `--dry` flag to preview without processing
2. Per-well debugging: `--debug` flag for verbose logging
3. Checkpoint validation: Ensure checkpoint files are valid JSON
4. Manual curation: MATLAB GUI or Phy viewer for spike inspection

**Debugging a Single Well**:
```bash
python IPNAnalysis/mea_analysis_routine.py /path/to/data.raw.h5 \
  --well well000 \
  --rec rec0001 \
  --config mea_config.json \
  --debug \
  --force-restart  # ignore checkpoint
```

## Development Workflow

### Git Branches
- **main**: Stable production branch
- **dev_branch**: Active development
- **parameter-detection**: Current development (WIP)
- **adamwea/issue###**: Issue-specific branches for fixes
- **copilot/###**: AI-assisted feature branches

### Common Git Tasks
```bash
# View recent changes
git log --oneline --max-count=20

# Diff between branches
git diff main parameter-detection -- IPNAnalysis/

# Create feature branch
git checkout -b feature/your-feature dev_branch

# Commit with author attribution
git commit -m "Description

Co-Authored-By: Name <email@example.com>"
```

### Remote
- **Origin**: https://github.com/roybens/MEA_Analysis.git

## Key Insights for Development

1. **Metadata Parsing**: The pipeline infers project structure from file paths:
   ```
   <project>/<date>/<chip>/<run_id>/Network/data.raw.h5
   ```
   Fallback to `.metadata` file if present (ConfigParser format).

2. **Checkpoint System**: Essential for long-running jobs. State is JSON with enum-based stage tracking. Always test checkpoint resumption when modifying pipeline stages.

3. **Memory Management**: Kilosort4 requires careful VRAM management:
   - `expandable_segments:True` prevents fragmentation
   - Auto-detects GPU (logs warning if CPU-only)
   - Adjusts batch size based on available VRAM (≥14GB → high; <14GB → low)

4. **Recording Types**: Pipeline handles multiple formats:
   - Maxwell HDF5 (`.h5`): via `si.read_maxwell()`
   - NWB (`.nwb`): via `si.read_nwb()`
   - Binary/folders: via `si.load_extractor()`

5. **Burst Detection**: Parameter-free, biological-calibration-based:
   - ISI distribution sets adaptive bin size
   - Population firing rate signal + unit participation
   - Adaptive thresholding prevents false positives on quiet recordings
   - Network merge gaps prevent over-fragmentation

6. **Quality Thresholds**: Configurable per project. Defaults assume standard network cultures. Organoid/high-noise recordings may need tuning.

7. **Parallel Processing**: Driver launches subprocesses (one per well). Each worker uses multi-threading for preprocessing (OMP_NUM_THREADS=16). SBATCH scripts coordinate GPU allocation across wells.

## Code Structure Highlights

### IPNAnalysis/ (Main Pipeline)
- `mea_analysis_routine.py` (906 lines): MEAPipeline class + main()
- `run_pipeline_driver.py` (346 lines): Driver orchestrator + file discovery
- `config_loader.py` (171 lines): Config resolution + template generation
- `parameter_free_burst_detector.py` (390 lines): Burst detection algorithm
- `helper_functions.py` (667 lines): Peak detection, plotting, file utilities
- `meaplotter.py` (1242 lines): Publication-grade figure generation
- `gaussianNetworkBursts.py` (257 lines): Legacy Gaussian-based burst model
- `scalebury.py` (87 lines): Scale bar overlay for plots

### MEAProcessingLibrary/
- Intended as installable package for reusable utilities
- Currently minimal; mostly a wrapper layer

### Other Modules
- See sections above for StimulationAnalysis, Connections, NetworkAnalysis, etc.

## Troubleshooting

### GPU Out of Memory
- Reduce `batch_size` in Kilosort4 parameters
- Increase `cluster_downsampling` (reduce cluster precision)
- Reduce `max_cluster_subset` cap
- Set `expandable_segments:True` in PYTORCH_CUDA_ALLOC_CONF

### HDF5 Plugin Errors
- Ensure `HDF5_PLUGIN_PATH` is set to correct libcompression.so path
- For Docker: path is `/opt/hdf5/plugins`
- For local: path is `/pscratch/sd/m/mpatil1/hdf5_plugin`

### Kilosort4 Crashes on Large Artifacts
- Ensure preprocessing uses highpass (300 Hz) + CMR (local, median)
- Float32 conversion prevents signal crushing
- Check for unsigned→signed conversion

### Checkpoint Resume Issues
- Delete checkpoint JSON if corrupted: `rm <checkpoint_file>.json`
- Re-run with `--force-restart` to ignore checkpoints
- Check logs for stage that failed

### Missing Wells / Empty Recording
- Pipeline skips empty wells gracefully
- Check HDF5 structure: `h5py.File(path, 'r').keys()` should show `recordings/` or `wells/`
- Verify `--rec` and `--well` parameters match HDF5 groups

## Useful Commands

```bash
# Inspect HDF5 file structure
python -c "import h5py; f=h5py.File('/path/data.raw.h5'); print(list(f.keys())); f.close()"

# Preview what would be processed
python IPNAnalysis/run_pipeline_driver.py /data/experiment --config mea_config.json --dry

# Check a checkpoint state
python -c "import json; print(json.dumps(json.load(open('/path/checkpoint.json')), indent=2))"

# Tail logs during processing
tail -f /results/path/*/run_id_well_id_pipeline.log

# Count spikes per unit post-curation
python -c "import numpy as np; spikes=np.load('/path/spike_times.npy', allow_pickle=True).item(); [print(u, len(spikes[u])) for u in sorted(spikes.keys())]"
```

