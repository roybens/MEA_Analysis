#general imports
import os
import platform
import re
import logging

#spikeinterface imports
import spikeinterface
import spikeinterface.sorters as ss
import spikeinterface.full as si

#local imports
from file_selector import main as file_selector_main
#from AxonReconPipeline.src.archive.extract_raw_assay_data import main as extract_raw_assay_data_main
from generate_sorting_objects import main as generate_sorting_objects
from axon_trace_classes import axon_trace_objects
from reconstruct_axons import reconstruct_axon_trace_objects
import mea_processing_library as MPL
import helper_functions as helper

#Logger Setup
#Create a logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Create handlers
stream_handler = logging.StreamHandler()  # logs to console
#file_handler = logging.FileHandler('file.log')  # logs to a file

# Set level of handlers
stream_handler.setLevel(logging.DEBUG)
#file_handler.setLevel(logging.ERROR)

# Add handlers to the logger
logger.addHandler(stream_handler)
#logger.addHandler(file_handler)

# Create formatters and add it to handlers
#formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - %(module)s.%(funcName)s')
stream_handler.setFormatter(formatter)

def clear_terminal():
    if platform.system() == "Windows":
        os.system('cls')
    else:
        os.system('clear')

def exclude_non_continuous_data(informed_h5_dirs):
    #Test for continuity (spiking data only turned off) in the .h5 files
    #Exclude non-continuous files
    continuous_h5_dirs = {}
    i=0
    for dir in informed_h5_dirs:
        h5_file_path = dir['h5_file_path']
        if MPL.test_continuity(h5_file_path, verbose = True):
            #copy informed_h5_dir to continuous_h5_dirs
            continuous_h5_dirs[i] = dir
            i=i+1
            print(f"Data for {h5_file_path} is continuous.")
    return continuous_h5_dirs

def main(debug_mode=False, debug_folders=None):
    # Clear the terminal
    clear_terminal()

    ##1. Select the folders to process, check for continuity, and extract the .h5 file paths
    logger.info("Selecting folders to process, checking for continuity, and extracting the .h5 file paths...")
    # Select the folders to process
    selected_folders = file_selector_main(pre_selected_folders= debug_folders, debug_mode=debug_mode)
    #Extract all .h5 files from the selected folders
    h5_dirs = MPL.extract_raw_h5_filepaths(selected_folders)
    #Extract chipIDs and record dates from the .h5 file paths
    informed_h5_dirs = MPL.extract_recording_details(h5_dirs)
    #Test for continuity (spiking data only turned off) in the .h5 files
    #Exclude non-continuous files
    continuous_h5_dirs = exclude_non_continuous_data(informed_h5_dirs)

    ##2. Process the .h5 files through MEA processing pipeline, save pre-processed data: 
    #recordings, spikesortins, waveforms, and templates
    logger.info("Executing MEA processing pipeline. This will process the .h5 files and save pre-processed data including recordings, spike sortings, waveforms, and templates.")
    
    #Extract merged recordings and sortings from each .h5 file
    allowed_scan_types = ["ActivityScan", "AxonTracking", "Network"]
    merged_recordings_by_scan = []
    merged_sortings_by_scan = []
    for dir in continuous_h5_dirs.values():
        h5_file_path = dir['h5_file_path']
        scan_type = dir['scanType']
        if scan_type in allowed_scan_types:
            logger.info(f"Processing {h5_file_path}...")
            merged_recordings_by_stream, merged_sorting_by_stream = generate_sorting_objects(h5_file_path, MaxWell_ID=2)
            merged_recordings_by_scan.append(merged_recordings_by_stream)
            merged_sortings_by_scan.append(merged_sorting_by_stream)
        else:
            logger.error("Error: scan type not recognized.")

    ##3. Extract waveforms
    #Extract recording details from the .h5 file path    
    for dir in continuous_h5_dirs.values():
        recording_details = MPL.extract_recording_details(dir['h5_file_path'])
        date = recording_details[0]['date']
        chip_id = recording_details[0]['chipID']
        scan_type = recording_details[0]['scanType']     
        # - Save spike sorting objects.
        logger.info(f"Extracting waveforms")
        # Check if recordings folder exists in ./data/temp_data, if not create it.
        waveforms_dir = './AxonReconPipeline/data/temp_data/waveforms'
        if not os.path.exists(waveforms_dir):
            os.makedirs(waveforms_dir)
        waveforms_by_stream = []
        for scan_num in range(len(merged_sortings_by_scan)):
            for stream_num in range(len(merged_sortings_by_scan[scan_num])):
                #debug
                if stream_num > 0: break
                #debug
                merged_recording = merged_recordings_by_scan[scan_num][stream_num]
                merged_sorting = merged_sortings_by_scan[scan_num][stream_num]
                stream_name = f"Well#{stream_num+1}"
                relative_path = waveforms_dir + f"/{date}/{chip_id}/{scan_type}/{stream_name}"
                # try: 
                #     waveforms = si.load_waveforms(relative_path)
                # except:
                tot_units = len(merged_sorting.get_unit_ids())
                units_per_extraction = tot_units/10
                we = MPL.generate_waveform_extractor_unit_by_unit(
                    merged_recording,
                    merged_sorting,folder = relative_path, 
                    n_jobs = 8,
                    units_per_extraction = units_per_extraction, 
                    sparse = False, 
                    fresh_extractor = False,
                    load_if_exists = True)
                waveforms_by_stream.append(we)
    
    ##4. Extract templates
    #axon_trace_precursors = axon_trace_objects(recording_dir[well_name], dirs)
    #reconstruct_axon_trace_objects(axon_trace_precursors)          
    logger.info("done.")
    
if __name__ == "__main__":
    debug_mode = True
    #testing Full Activity Scan and Network Scan...
    selected_folders_test1 = [        
        "/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/ActivityScan/000023",
        "/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/Network/000024",
        #testing continuity test with non-continuous data
        #"/mnt/disk20tb/PrimaryNeuronData/Maxone/SPTAN1_1/230725/16657/ActivityScan",
        ]
    #...vs. Axon Tracking Scan Alone
    selected_folders_test2 = [  
        "/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/AxonTracking/000026",
        ]  
    if debug_mode:
        main(debug_mode, debug_folders=selected_folders_test1)
        main(debug_mode, debug_folders=selected_folders_test2)
    else:
        main()