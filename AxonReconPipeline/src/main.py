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
from generate_recording_and_sorting_objects import main as generate_recording_and_sorting_objects
from axon_trace_classes import axon_trace_objects
from reconstruct_axons import reconstruct_axon_trace_objects
import mea_processing_library as MPL
import helper_functions as helper
import spike_sorting as hp #hornauerp, https://github.com/hornauerp/axon_tracking/blob/main/axon_tracking/spike_sorting.py

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

def generate_merged_recs_and_sorts_by_scan(continuous_h5_dirs, 
                                           allowed_scan_types=["ActivityScan", "AxonTracking", "Network"], 
                                           recordings_dir = './AxonReconPipeline/data/temp_data/recordings', 
                                           spikesorting_dir = './AxonReconPipeline/data/temp_data/sortings'):
    merged_recordings_by_scan = []
    merged_sortings_by_scan = []
    h5_file_paths = [dir['h5_file_path'] for dir in continuous_h5_dirs.values()]
    for dir in continuous_h5_dirs.values():
        h5_file_path = dir['h5_file_path']
        scan_type = dir['scanType']
        if scan_type in allowed_scan_types:
            logger.info(f"Processing {h5_file_path}...")
            merged_recordings_by_stream, merged_sorting_by_stream = generate_recording_and_sorting_objects(h5_file_paths, 
                                                                                                           h5_file_path, 
                                                                                                           recordings_dir = recordings_dir, 
                                                                                                           spikesorting_dir = spikesorting_dir)
            merged_recordings_by_scan.append(merged_recordings_by_stream)
            merged_sortings_by_scan.append(merged_sorting_by_stream)
        else:
            logger.error("Error: scan type not recognized.")
    return merged_recordings_by_scan, merged_sortings_by_scan

def extract_sortings_per_stream(
                continuous_h5_dirs, 
                allowed_scan_types=[
                                    #"ActivityScan", 
                                    "AxonTracking", 
                                    #"Network"
                                    ], 
                recordings_dir = './AxonReconPipeline/data/temp_data/recordings',): 
    
    """ Extract sorting objects per stream from sliced and cocatenated recordings (most/all fixed electrodes during axontracking scan)
    """
    
    #merged_recordings_by_scan = []
    #merged_sortings_by_scan = []
    h5_file_paths = [dir['h5_file_path'] for dir in continuous_h5_dirs.values() if dir['scanType'] in allowed_scan_types]
    recording_details = MPL.extract_recording_details(h5_file_paths[0])
    date = recording_details[0]['date']
    chip_id = recording_details[0]['chipID']
    scanType = recording_details[0]['scanType']
    run_id = recording_details[0]['runID']
    
    #Sorter params
    sorter_params = si.get_default_sorter_params(si.Kilosort2_5Sorter)
    sorter_params['n_jobs'] = -1
    sorter_params['detect_threshold'] = 7
    sorter_params['minFR'] = 0.01
    sorter_params['minfr_goodchannels'] = 0.01
    sorter_params['keep_good_only'] = False
    sorter_params['do_correction'] = False
    #
    verbose = True
    spikesorting_dir = './AxonReconPipeline/data/temp_data/sortings'
    spikesorting_root = spikesorting_dir+f'/{date}/{chip_id}/{scanType}/{run_id}'
    merged_sorting_list_by_stream = None
    
    h5_file_paths = [dir['h5_file_path'] for dir in continuous_h5_dirs.values() if dir['scanType'] == "AxonTracking"]
    sorting_list_by_stream = hp.sort_recording_list(h5_file_paths,
                                                            save_root=spikesorting_root,
                                                            #save_path_changes= 5,
                                                            sorter = 'kilosort2_5',
                                                            sorter_params = sorter_params,
                                                            verbose = verbose,
                                                            scan_merge = False) 
    return sorting_list_by_stream

def extract_waveforms_by_scan(continuous_h5_dirs, merged_recordings_by_scan, merged_sortings_by_scan, waveforms_dir='./AxonReconPipeline/data/temp_data/waveforms'):
    waveforms_by_scan = []
    for scan_num, (dir, merged_sortings) in enumerate(zip(continuous_h5_dirs.values(), merged_sortings_by_scan)):
        recording_details = MPL.extract_recording_details(dir['h5_file_path'])
        date = recording_details[0]['date']
        chip_id = recording_details[0]['chipID']
        scan_type = recording_details[0]['scanType'] 
        run_id = recording_details[0]['runID']    

        if not os.path.exists(waveforms_dir):
            os.makedirs(waveforms_dir)

        #debug
        if scan_num > 0: break
        #debug

        for stream_num in range(len(merged_sortings)):
            #debug
            if stream_num > 0: break
            #debug
            waveforms_by_stream = []
            merged_recording = merged_recordings_by_scan[scan_num][stream_num]
            merged_sorting = merged_sortings_by_scan[scan_num][stream_num]
            stream_name = f"Well#{stream_num+1}"
            relative_path = waveforms_dir + f"/{date}/{chip_id}/{scan_type}/{run_id}/{stream_name}"

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
        waveforms_by_scan.append(waveforms_by_stream)
    return waveforms_by_scan

def post_process_waveforms(continuous_h5_dirs, we_by_scan, good_we_dir='./AxonReconPipeline/data/temp_data/good_waveforms'):
    good_we_by_scan = []
    for scan_num, (dir, we) in enumerate(zip(continuous_h5_dirs.values(), we_by_scan)):
        recording_details = MPL.extract_recording_details(dir['h5_file_path'])
        date = recording_details[0]['date']
        chip_id = recording_details[0]['chipID']
        scan_type = recording_details[0]['scanType'] 
        run_id = recording_details[0]['runID']    

        if not os.path.exists(good_we_dir):
            os.makedirs(good_we_dir)

        for stream_num, we_stream in enumerate(we):
            #debug
            if stream_num > 0: break
            #debug
            good_we_by_stream = []
            we = we_by_scan[scan_num][stream_num]
            stream_name = f"Well#{stream_num+1}"
            relative_path = good_we_dir + f"/{date}/{chip_id}/{scan_type}/{run_id}/{stream_name}"
            logger.info(f"Post-processing waveforms. Saving conserved waveforms and metrics to {relative_path}")
            good_we = MPL.postprocess_waveforms(we,
                                                relative_path,
                                                get_quality =True,
                                                remove_violations = True,
                                                remove_similar = True)
            good_we_by_stream.append(good_we)
        good_we_by_scan.append(good_we_by_stream)
    return good_we_by_scan

def compute_templates_by_scan(continuous_h5_dirs, we_by_scan, template_dir='./AxonReconPipeline/data/temp_data/templates'):
    templates_by_scan = []
    for scan_num, dir in enumerate(continuous_h5_dirs.values()):
        recording_details = MPL.extract_recording_details(dir['h5_file_path'])
        date = recording_details[0]['date']
        chip_id = recording_details[0]['chipID']
        scan_type = recording_details[0]['scanType'] 
        run_id = recording_details[0]['runID']    

        if not os.path.exists(template_dir):
            os.makedirs(template_dir)

        #debug
        if scan_num > 0: break
        #debug

        for stream_num in range(len(we_by_scan[scan_num])):
            #debug
            if stream_num > 0: break
            #debug
            templates_by_stream = []
            we = we_by_scan[scan_num][stream_num]
            stream_name = f"Well#{stream_num+1}"
            relative_path = template_dir + f"/{date}/{chip_id}/{scan_type}/{run_id}/{stream_name}"
            templates = we.get_all_templates()
            templates_by_stream.append(templates)
        templates_by_scan.append(templates_by_stream)
    return templates_by_scan

def main(debug_mode=False, debug_folders=None):
    # Clear the terminal
    clear_terminal()

    ##1. Select the folders to process, check for continuity, and extract the .h5 file paths
    logger.info("Selecting folders to process, checking for continuity, and extracting the .h5 file path:")
    # Select the folders to process
    selected_folders = file_selector_main(pre_selected_folders= debug_folders, debug_mode=debug_mode)
    logger.info(f"Selected folders: {selected_folders}")
    #Extract all .h5 files from the selected folders
    h5_dirs = MPL.extract_raw_h5_filepaths(selected_folders)
    #Extract chipIDs and record dates from the .h5 file paths
    informed_h5_dirs = MPL.extract_recording_details(h5_dirs)
    #Test for continuity (spiking data only turned off) in the .h5 files
    #Exclude non-continuous files
    continuous_h5_dirs = exclude_non_continuous_data(informed_h5_dirs)
    logger.info(f"Continuous data found in {len(continuous_h5_dirs)} files.")

    ##2. Extract sorting objects per stream from sliced and cocatenated recordings (most/all fixed electrodes during axontracking scan)
    logger.info("Loading and merging the recordings as needed from the .h5 files. Generating spike sorting objects and merging as needed.")
    allowed_scan_types = [
        #"ActivityScan", Dense activity scan only needed to define locations of electrodes in network and axon tracking scans
        "AxonTracking", 
        #"Network"
        ]
    logger.info(f"Allowed scan types: {allowed_scan_types}")
    recordings_dir = './AxonReconPipeline/data/temp_data/recordings'
    sorting_list_by_stream = extract_sortings_per_stream(continuous_h5_dirs, allowed_scan_types, recordings_dir)



    ##2. Perform spikesorting. Load and merge the recordings as needed from the .h5 files. Generate spike sorting objects and merge as needed: 
    spikesorting_dir = './AxonReconPipeline/data/temp_data/sortings'
    merged_recordings_by_scan, merged_sortings_by_scan = generate_merged_recs_and_sorts_by_scan(continuous_h5_dirs, 
                                                                                                allowed_scan_types, 
                                                                                                recordings_dir = recordings_dir,
                                                                                                spikesorting_dir = spikesorting_dir)

    ##3. Extract waveforms
    #Extract recording details from the .h5 file path
    waveforms_dir='./AxonReconPipeline/data/temp_data/waveforms'    
    we_by_scan = extract_waveforms_by_scan(continuous_h5_dirs, 
                                            merged_recordings_by_scan, 
                                            merged_sortings_by_scan, 
                                            waveforms_dir=waveforms_dir)
                
    ##4. Post-process waveforms
    good_waveforms_dir='./AxonReconPipeline/data/temp_data/good_waveforms'
    good_we_by_scan = post_process_waveforms(continuous_h5_dirs, 
                                            we_by_scan, 
                                            good_we_dir=good_waveforms_dir)

    
    
    ##5. Extract templates
    templates_dir='./AxonReconPipeline/data/temp_data/templates'
    templates_by_scan = compute_templates_by_scan(continuous_h5_dirs, 
                                                  we_by_scan, 
                                                  templates_dir=templates_dir)
    
    
    #axon_trace_precursors = axon_trace_objects(recording_dir[well_name], dirs)
    #reconstruct_axon_trace_objects(axon_trace_precursors)          
    logger.info("done.")
    
if __name__ == "__main__":
    debug_mode = True
    #testing Full Activity Scan and Network Scan...
    selected_folders_test1 = [        
        #"/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/ActivityScan/000023",
        "/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/Network/000024",
        "/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/AxonTracking/000026",
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