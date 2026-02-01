# High-Density MEA Analysis Pipeline Pseudocode

## 1. Global Setup & Entry Point (`main`)

**Goal:** Parse command line arguments, setup logging, and initiate the processing block.

```text
PROGRAM Main:
    IMPORT libraries (SpikeInterface, PyTorch, Plotting, etc.)
    DEFINE Global Configuration (Job settings, Paths)

    FUNCTION Start():
        // 1. Parse Command Line Arguments
        INPUT: FilePath, WellID, Parameters, DockerImage, Flags (Resume, Cleanup, etc.)
        
        // 2. Setup Environment
        DETERMINE Output_Directory and Checkpoint_Directory
        PARSE Metadata from input file path (Project, Date, ChipID, RunID)
        INITIALIZE Logger (Log file specific to run & well)
        
        LOAD Curation Thresholds from JSON (if provided)

        // 3. Execute Pipeline
        TRY:
            CALL process_block(FilePath, WellID, Thresholds, Directories...)
        CATCH Error:
            LOG Error details and Traceback
            EXIT Program with Failure
```


## 2. Process Block logic
```text
FUNCTION process_block(File, Well, Thresholds, Config):

    // --- Initialization ---
    DETECT Input Format (Maxwell .h5 or Binary)
    CREATE Output Directory Structure
    INITIALIZE Checkpoint System
    
    IF Checkpoint.ShouldSkip():
        PERFORM Resource Cleanup
        RETURN Success

    // --- STAGE 0: Preprocessing ---
    IF Checkpoint.Stage < SORTING_COMPLETE:
        LOG "Starting Preprocessing"
        
        LOAD Raw Recording (Maxwell H5)
        
        // SOTA Preprocessing Steps
        APPLY Filter: Highpass (300Hz)
        APPLY Reference: Local Common Median Reference (CMR) or Global fallback
        CONVERT Data Type: Float32 -> Int16 (Scaling)
        
        SAVE Preprocessed Chunk to Binary/Zarr format (for fast I/O)
        RELOAD Preprocessed Recording
    
    // --- STAGE 1: Spike Sorting ---
    IF Checkpoint.Stage < SORTING_COMPLETE:
        LOG "Starting Spike Sorting"
        
        DEFINE Sorter Parameters (Kilosort4 specific: nblocks=0 for 2D, thresholds)
        
        IF UseDocker:
            RUN Sorter inside Container
        ELSE:
            RUN Sorter Locally
            
        // Post-Sort Cleaning
        REMOVE Empty Units
        REMOVE Excess Spikes
        REMOVE Duplicates (0.1ms window)
        
        SAVE Checkpoint (Stage = SORTING_COMPLETE)

    // --- STAGE 2: Analyzer & Extensions ---
    IF Checkpoint.Stage < ANALYZER_COMPLETE:
        LOG "Computing Sorting Analyzer"
        
        // 1. Sparsity
        ESTIMATE Sparsity (Radius method, ~50um)
        
        // 2. Create Analyzer Object
        CREATE SortingAnalyzer (links Recording + Sorting)
        
        // 3. Compute Extensions (Features needed for curation)
        COMPUTE Missing Extensions:
            - Waveforms
            - Templates
            - Noise Levels
            - Quality Metrics (SNR, ISI violation, etc.)
            - Unit Locations (Monopolar Triangulation)
            
        SAVE Checkpoint (Stage = ANALYZER_COMPLETE)

    // --- STAGE 3: Reports, Curation & Export ---
    IF Checkpoint.Stage < REPORTS_COMPLETE:
        LOG "Generating Reports"
        
        // 1. Export Raw Metrics
        EXPORT "quality_metrics_unfiltered.xlsx"
        
        // 2. Curation (Filtering)
        IF NoCuration Flag is SET:
            KEEP All Units
        ELSE:
            APPLY Automatic Curation (Thresholds: SNR, FiringRate, PresenceRatio)
            LOG Rejected Units to Excel
        
        // 3. Visualization
        GENERATE Plots:
            - Probe Map (All Units)
            - Probe Map (Curated Units)
            - Waveform PDFs (Grid view of waveforms)
        
        // 4. Burst Analysis
        EXTRACT Spike Trains for Good Units
        DETECT Bursts (ISI Threshold method)
        COMPUTE Network Metrics (Synchrony, Burst Rate)
        GENERATE Raster Plots
        SAVE "network_data.json"
        
        // 5. Phy Export (Optional)
        IF ExportToPhy Flag is SET:
            CONVERT Data to Phy format for manual review
            
        SAVE Checkpoint (Stage = REPORTS_COMPLETE)

    // --- Cleanup ---
    IF Cleanup Flag is SET:
        DELETE Intermediate Binary Files
        DELETE Intermediate Sorter Outputs
        FREE Memory (RAM/GPU)

    RETURN Success
```
## 3. Helper fucntions
```text
FUNCTION preprocess(recording):
    // Used in Stage 0
    INPUT: Raw Recording
    Step 1: Unsigned to Signed conversion
    Step 2: Highpass Filter (remove LFP drift)
    Step 3: Common Median Reference (Local radius 250um)
    Step 4: Scale to Int16
    OUTPUT: Processed Recording

FUNCTION automatic_curation(metrics, thresholds):
    // Used in Stage 3
    FOR EACH Unit in metrics:
        CHECK Rules:
            - Is Presence Ratio > Threshold?
            - Is ISI Contamination < Threshold?
            - Is Firing Rate > Threshold?
            - Is Amplitude < Threshold?
        IF Fails any rule:
            Mark as Rejected
            Log Reason
    RETURN Cleaned_Metrics, Rejection_Log
```
## 4. Checkpointing

```text
CLASS PipelineCheckpoint:
    PROPERTIES: RunID, WellID, CurrentStage, CheckpointFile

    FUNCTION Init():
        IF CheckpointFile exists:
            LOAD State from JSON
        ELSE:
            CREATE new State (Stage = NOT_STARTED)

    FUNCTION ShouldSkip():
        IF User requests Resume AND Stage == REPORTS_COMPLETE:
            RETURN True (Skip processing)
        IF User requests ForceRestart:
            DELETE CheckpointFile
            RETURN False
        RETURN False

    FUNCTION Save(NewStage):
        UPDATE State with NewStage and Timestamp
        WRITE State to JSON file
```