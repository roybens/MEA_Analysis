# MEA Data Processing Pipeline

This is a production-grade **MEA (Microelectrode Array) processing pipeline** for neuronal spike sorting, analysis, and burst detection. It processes neuronal data recorded from Maxwell Biosystems MEA systems with comprehensive preprocessing, spike sorting using Kilosort, waveform extraction, and neuronal network burst analysis.

## Architecture Overview

### Two-Tier Design

**1. run_pipeline_driver.py — Orchestrator**
- Scans directories or processes single HDF5 files
- Discovers all recordings and wells in each data file
- Maps recordings to their associated wells using `recording_map` dictionary
- Launches subprocesses for each recording-well combination via `mea_analysis_routine.py`
- Handles reference filtering, dry-runs, and logging
- Supports batch processing with checkpoint-based resumption
- Accepts a `--config` JSON file for project-level defaults; CLI flags always override config

**2. mea_analysis_routine.py — Core Pipeline (`MEAPipeline` class)**
- Per-well worker executing the full processing pipeline
- Stages: Preprocessing → Sorting → Analysis → Reports
- Checkpoint-based resumption (resume on crash, skip completed stages)
- Metadata parsing and intelligent output directory structuring
- Can be invoked independently of the driver for single-well processing

**3. config_loader.py — Shared Configuration**
- Shared by both scripts
- Resolves values in priority order: CLI flag → config file → hardcoded default
- `build_extra_args()` constructs subprocess argument strings for the driver
- Run directly (`python config_loader.py`) to generate a config template

## HDF5 File Structure

The pipeline expects HDF5 files with the following hierarchical structure:

```
file.h5
├── recordings/
│   ├── rec0001/
│   │   ├── well000
│   │   ├── well001
│   │   └── ...
│   ├── rec0002/
│   │   ├── well000
│   │   └── ...
│   └── ...
└── wells/ (legacy support)
```

The pipeline uses a `recording_map` to efficiently map all recordings to their corresponding wells without keeping the entire file open.

## Processing Pipeline Stages

| Stage | Method | Purpose |
|-------|--------|---------|
| **Preprocessing** | `run_preprocessing()` | Highpass filter (300 Hz), local common median reference, float32 conversion, binary cache |
| **Spike Sorting** | `run_sorting()` | Kilosort4 (default), SpikeInterface integration, Docker support for reproducibility |
| **Analyzer** | `run_analyzer()` | Template computation, quality metrics (firing rate, presence ratio, ISI violations, amplitude) |
| **Reports** | `generate_reports()` | Waveform visualizations, probe locations, burst analysis, automatic curation |

## Key Features

✅ **Multi-Recording Support** — Handle multiple recordings per HDF5 file  
✅ **Checkpointing** — Resume from failures; state saved as JSON  
✅ **Logging** — Per-well logs + driver-level logs with timestamps  
✅ **Quality Curation** — Automatic unit filtering via configurable thresholds  
✅ **Containerization** — Docker image support for Kilosort reproducibility  
✅ **Flexible I/O** — HDF5 (primary), placeholders for NWB and raw binary formats  
✅ **Burst Detection** — Parameter-free adaptive network burst detection  
✅ **Visualization** — Waveforms, probe maps, rasters, network burst overlays  
✅ **Config System** — JSON config file with CLI override support  

## Supporting Modules

- **config_loader.py** — Config file loading, CLI/config/default resolution, subprocess arg builder
- **helper_functions.py** — Peak detection, file discovery, raster and network plotting utilities
- **parameter_free_burst_detector.py** — Adaptive network burst detection with:
  - Per-unit ISI bursts (biological calibration)
  - Population firing rate signal computation
  - Adaptive thresholding + squelch logic
  - Network synchrony metrics
- **spikeMatrix.py** — Spike raster representation and matrix operations
- **gaussianNetworkBursts.py** — Gaussian-based burst modeling

## Directory Structure & Metadata

The pipeline infers project structure from file paths:
```
<project>/<date>/<chip>/<run_id>/Network/data.raw.h5
```

Output organized as:
```
<output_dir>/<project>/<date>/<chip>/<run_id>/well000/
  ├── binary/                    (preprocessed recording cache)
  ├── sorter_output/             (kilosort4 outputs)
  ├── analyzer_output/           (waveforms, templates, quality metrics)
  ├── raster_burst_plot.svg      (full recording raster + network burst)
  ├── raster_burst_plot_30s.svg  (30s zoom)
  ├── raster_burst_plot_60s.svg  (60s zoom)
  ├── network_results.json       (burst statistics)
  ├── spike_times.npy            (spike times per unit)
  ├── metrics_curated.xlsx       (quality metrics post-curation)
  ├── rejection_log.xlsx         (rejected units and reasons)
  ├── waveforms_grid.pdf         (waveform overview)
  ├── locations_unfiltered.pdf   (all unit probe locations)
  └── checkpoints/               (resume state per well)
```

## Configuration System

The pipeline uses a three-level priority chain:

```
CLI flag  →  config file (mea_config.json)  →  hardcoded defaults
```

**Generate a config template:**
```bash
python config_loader.py mea_config.json
```

Edit `mea_config.json` for your project, then pass it to either script via `--config`. CLI flags always override config file values. `--well` and `--rec` are always CLI-only (per-file identifiers, not project settings).

**Config sections:**

| Section | Keys |
|---------|------|
| `io` | `output_dir`, `checkpoint_dir`, `export_to_phy`, `clean_up` |
| `sorting` | `sorter`, `docker_image` |
| `filtering` | `reference_file`, `assay_types` (driver only) |
| `plotting` | `plot_mode`, `raster_sort`, `plot_debug` |
| `curation` | `no_curation`, `quality_thresholds` |

## Typical Workflow

### 1. Set up config for a new project
```bash
python config_loader.py mea_config.json
# edit mea_config.json — set output_dir, sorter, thresholds
```

### 2. Dry run to verify what would be processed
```bash
python run_pipeline_driver.py /data/experiment --config mea_config.json --dry
```

### 3. Full batch run
```bash
python run_pipeline_driver.py /data/experiment --config mea_config.json
```

### 4. Override config at runtime
```bash
# config says kilosort4, override for this run
python run_pipeline_driver.py /data/experiment --config mea_config.json --sorter mountainsort5
```

### 5. Single well processing
```bash
python mea_analysis_routine.py /data/exp/run_001/Network/data.raw.h5 \
  --well well000 --rec rec0001 \
  --config mea_config.json \
  --debug
```

### 6. Skip sorting, detection only
```bash
python run_pipeline_driver.py /data/experiment --config mea_config.json --skip-spikesorting
```

### 7. Re-analyze bursts on existing spike times
```bash
python mea_analysis_routine.py /data/file.h5 --well well000 \
  --config mea_config.json \
  --reanalyze-bursts \
  --plot-mode merged \
  --raster-sort firing_rate \
  --plot-debug
```

### 8. Resume after crash (automatic via checkpoint)
```bash
# just re-run the same command — checkpoints handle resumption automatically
python run_pipeline_driver.py /data/experiment --config mea_config.json

# force full restart ignoring checkpoints
python run_pipeline_driver.py /data/experiment --config mea_config.json --force-restart
```

## Command-Line Reference

### run_pipeline_driver.py

| Group | Argument | Description |
|-------|----------|-------------|
| positional | `path` | File or directory path containing MEA data |
| input/output | `--config` | Path to config JSON (CLI flags always override) |
| input/output | `--output-dir` | Output directory for all results |
| input/output | `--checkpoint-dir` | Checkpoint directory (default: output-dir/checkpoints) |
| input/output | `--export-to-phy` | Export results to Phy format |
| input/output | `--clean-up` | Remove intermediate files after processing |
| filtering | `--reference` | Excel file to filter runs by assay type |
| filtering | `--type` | Assay types to include |
| sorting | `--sorter` | Spike sorter to use (default: kilosort4) |
| sorting | `--docker` | Docker image for containerized sorting |
| sorting | `--skip-spikesorting` | Spike detection only, skip full sorting |
| plotting | `--plot-mode` | `separate` or `merged` (default: separate) |
| plotting | `--raster-sort` | `none`, `firing_rate`, `location_y`, `unit_id` |
| plotting | `--plot-debug` | Overlay burst/superburst intervals on plot |
| curation | `--no-curation` | Skip automatic unit curation |
| curation | `--params` | JSON string or file with quality thresholds |
| run control | `--force-restart` | Ignore checkpoints, restart from scratch |
| run control | `--reanalyze-bursts` | Re-run burst analysis on existing spike times |
| run control | `--dry` | Preview what would run without processing |
| run control | `--debug` | Enable verbose logging |

### mea_analysis_routine.py

Same groups and flags as the driver, minus `--dry` and the filtering group, plus:

| Group | Argument | Description |
|-------|----------|-------------|
| input/output | `--well` | Well ID to process (e.g. well000) — **required** |
| input/output | `--rec` | Recording name in HDF5 file (default: rec0000) |

## Notes

- `--well` and `--rec` are always CLI-only — they identify a specific file/recording and are never set in config
- `--debug`, `--dry`, `--force-restart`, `--reanalyze-bursts`, `--skip-spikesorting` are CLI-only run control flags — they represent one-off decisions and are never set in config
- Everything else can be set in `mea_config.json` and overridden per-run from CLI