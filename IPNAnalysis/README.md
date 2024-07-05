# MEA Data Processing Pipeline

This script is designed to process neuronal data recorded from MEA (Microelectrode Arrays) systems. The main functionality includes data preprocessing, spike sorting using Kilosort, waveform extraction, and various post-processing steps to analyze the neuronal signals.

## Class Descriptions and Key Functions

![UML Diagram](/mnt/disk15tb/rohan/MEA_Analysis/IPNAnalysis/image.png)

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
