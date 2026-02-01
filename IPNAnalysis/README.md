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

**2. mea_analysis_routine.py — Core Pipeline (`MEAPipeline` class)**
- Per-well worker executing the full processing pipeline
- Stages: Preprocessing → Sorting → Analysis → Reports
- Checkpoint-based resumption (resume on crash, skip completed stages)
- Metadata parsing and intelligent output directory structuring

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
| **Preprocessing** | `run_preprocessing()` | Bandpass filter (300–6000 Hz), common average reference, optional Zarr compression |
| **Spike Sorting** | `run_sorting()` | Kilosort4 (default), SpikeInterface integration, Docker support for reproducibility |
| **Analyzer** | `run_analyzer()` | Template computation, quality metrics (firing rate, SNR, ISI violations) |
| **Reports** | `generate_reports()` | Waveform visualizations, probe locations, burst analysis, automatic curation |

## Key Features

✅ **Multi-Recording Support** — Handle multiple recordings per HDF5 file  
✅ **Checkpointing** — Resume from failures; state saved as JSON  
✅ **Logging** — Per-well logs + driver-level logs with timestamps  
✅ **Quality Curation** — Automatic unit filtering via configurable thresholds (ISI violations, SNR, firing rate)  
✅ **Containerization** — Docker image support for Kilosort reproducibility  
✅ **Flexible I/O** — HDF5 (primary), placeholders for NWB and raw binary formats  
✅ **Burst Detection** — Parameter-free adaptive network burst detection  
✅ **Visualization** — Waveforms, probe maps, rasters, metrics  

## Supporting Modules

- **helper_functions.py** — Peak detection, file discovery, directory management utilities
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
AnalyzedData/<project>/<date>/<chip>/<run_id>/well000/
  ├── kilosort4__rec0000/        (sorter outputs)
  ├── waveforms__rec0000/        (extracted waveforms)
  ├── reports/                   (plots, metrics)
  └── checkpoints/               (resume state)
```

## Typical Workflow

### Process Entire Directory Tree with Docker
```bash
python IPNAnalysis/run_pipeline_driver.py /data/experiment \
  --docker si-98-ks4-maxwell \
  --sorter kilosort4 \
  --output-dir /results
```

### Resume on Failure (Automatic via Checkpoint)
```bash
python IPNAnalysis/run_pipeline_driver.py /data/experiment \
  --docker si-98-ks4-maxwell \
  --force-restart
```

### Single File Processing
```bash
python IPNAnalysis/mea_analysis_routine.py /data/exp/run_001/Network/data.raw.h5 \
  --rec rec0001 \
  --well well000 \
  --debug
```

### Dry Run (Preview Processing Steps)
```bash
python IPNAnalysis/run_pipeline_driver.py /data/experiment --dry
```

### Re-analyze Bursts Only
```bash
python IPNAnalysis/run_pipeline_driver.py /data/experiment \
  --skip-spikesorting \
  --reanalyze-bursts
```

## Command-Line Arguments Reference

| Argument | Type | Description |
|----------|------|-------------|
| `path` | str | File or directory path containing MEA data |
| `--reference` | str | Excel file to filter runs by assay type |
| `--type` | list | Assay types to include (default: ['network today', 'network today/best']) |
| `--params` | str | JSON file or string with quality thresholds |
| `--docker` | str | Docker image name for containerized sorting |
| `--skip-spikesorting` | flag | Skip spike sorting stage |
| `--force-restart` | flag | Restart even if checkpoint complete |
| `--sorter` | str | Sorter algorithm to use (default: kilosort4) |
| `--debug` | flag | Debug mode with verbose logging |
| `--dry` | flag | Dry run only (preview without processing) |
| `--clean-up` | flag | Clear sorter outputs etc. |
| `--export-to-phy` | flag | Export results to Phy format |
| `--checkpoint-dir` | str | Checkpoint directory (default: output-dir/checkpoints) |
| `--no-curation` | flag | Skip automatic curation |
| `--output-dir` | str | Output directory for results |
| `--reanalyze-bursts` | flag | Re-analyze bursts even if present |

## Class Descriptions and Key Functions

![UML Diagram](image.png)

### Logger
- **Attributes:**
  - `logger`
  - `logging`
- **Methods:**
  - `setup()`: Configures the logger.
  - `log_info(msg)`: Logs informational messages.
  - `log_error(msg)`: Logs error messages.

### DataHandler
- **Attributes:**
  - `file_path`
  - `stream_id`
- **Methods:**
  - `get_data_maxwell()`: Reads Maxwell MEA recordings.
  - `save_to_zarr()`: Saves and compresses the file into the Zarr format.

### Preprocessor
- **Methods:**
  - `preprocess()`: Performs bandpass filtering and common average referencing on the recorded data to enhance signal quality.
  - `get_channel_recording_stats()`: Retrieves channel recording statistics such as sampling frequency, number of channels, and total recording duration.

### KilosortHandler
- **Methods:**
  - `run_kilosort()`: Runs the Kilosort algorithm to detect and sort spikes.
  - `get_kilosort_result()`: Retrieves the result from a specified Kilosort output folder.

### WaveformHandler
- **Methods:**
  - `extract_waveforms()`: Extracts waveform snippets around detected spikes.
  - `get_waveforms_result()`: Loads waveform results, optionally with the associated recording and sorting.

### QualityMetrics
- **Methods:**
  - `get_quality_metrics()`: Computes various quality metrics such as firing rates, SNR, and inter-spike interval violations.
  - `remove_violated_units()`: Filters out units that do not meet specific quality criteria.

### TemplateHandler
- **Methods:**
  - `get_template_metrics()`: Computes metrics related to the waveform templates.
  - `get_temp_similarity_matrix()`: Computes a similarity matrix for the templates.
  - `merge_similar_templates()`: Merges similar templates based on a similarity score.
  - `remove_similar_templates()`: Identifies and removes redundant templates based on a similarity score.

### ChannelAnalyzer
- **Methods:**
  - `get_unique_templates_channels()`: Analyzes units and their corresponding extremum channels, removes units that correspond to the same channel, and keeps the highest amplitude one.
  - `get_channel_locations_mapping()`: Retrieves channel location mappings.

### Visualization
- **Methods:**
  - `plot_waveforms()`: Visualizes the spatial distribution of units on the MEA probe and saves waveform plots.

### Utilities
- **Methods:**
  - `find_connected_components()`: Finds connected components in a graph, used for identifying transitive relationships between templates.

### Helper
- **Methods:**
  - `isexists_folder_not_empty()`: Checks if a folder exists and is not empty.
  - `empty_directory()`: Empties the contents of a directory.

### Main
- **Methods:**
  - `main()`: Handles command-line arguments and determines whether to process a single file or a directory of files.
  - `process_block(file_path, time_in_s=300, stream_id='well000', recnumber=0, sorting_folder=f"{BASE_FILE_PATH}/Sorting_Intermediate_files", clear_temp_files=True)`: Processes a block of data, including preprocessing, spike sorting, waveform extraction, and saving results.
  - `routine_sequential(file_path, number_of_configurations, time_in_s)`: Processes multiple recording configurations sequentially.
  - `routine_parallel(file_path, number_of_configurations, time_in_s)`: Processes multiple recording configurations in parallel using multiprocessing.

## Data Flow and Interactions

1. **Logger Setup:**
   - The `Logger` class is instantiated and configured to log messages to `application.log`.

2. **Data Handling:**
   - The `DataHandler` class reads Maxwell MEA recordings and can save them in the Zarr format.

3. **Preprocessing:**
   - The `Preprocessor` class preprocesses the data by performing bandpass filtering and common average referencing. It also retrieves channel recording statistics.

4. **Spike Sorting:**
   - The `KilosortHandler` class runs the Kilosort algorithm to detect and sort spikes. The results are retrieved from the specified output folder.

5. **Waveform Extraction and Analysis:**
   - The `WaveformHandler` class extracts waveform snippets around detected spikes and loads waveform results.
   - The `QualityMetrics` class computes quality metrics and filters out units that do not meet specific criteria.

6. **Template Analysis:**
   - The `TemplateHandler` class computes metrics related to waveform templates, identifies and merges similar templates, and removes redundant templates.

7. **Channel Analysis:**
   - The `ChannelAnalyzer` class analyzes units and their corresponding extremum channels, ensuring only unique templates are retained.

8. **Visualization:**
   - The `Visualization` class visualizes the spatial distribution of units on the MEA probe and saves waveform plots.

9. **Utilities and Helper Functions:**
   - The `Utilities` class finds connected components for identifying transitive relationships between templates.
   - The `Helper` class provides utility functions to check and empty directories.

10. **Main Execution:**
    - The `Main` class handles the command-line arguments and processes data either sequentially or in parallel. It coordinates the entire data processing pipeline, from reading the data to saving the final results.

## Conclusion

This script provides a comprehensive pipeline for MEA data analysis, from preprocessing and spike sorting to detailed waveform analysis and visualization. It utilizes advanced sorting algorithms and quality metric computations to ensure high-quality data analysis, suitable for neuroscientific research.
