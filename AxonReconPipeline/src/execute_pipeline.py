#general imports
import os
import platform
import re
import logging
from pprint import pprint

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
import template_extraction as te_hp

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

def extract_templates_hp(continuous_h5_dirs, 
                         sorting_lists_by_scan,
                         allowed_scan_types = ["AxonTracking"], 
                         templates_dir='./AxonReconPipeline/data/temp_data/templates',
                         stream_select = None,
                         just_load_waveforms = False):
    te_params = dict()
    te_params['align_cutout'] = True #Align waveforms by max waveform peak
    te_params['upsample'] = 2 #Factor by which to upsample waveforms
    te_params['rm_outliers'] = True #Check if outliers should be removed
    te_params['n_jobs'] = 16 #Number of cores to use for waveform extraction
    te_params['n_neighbors'] = 10 #Number of neighbors for outlier detection
    te_params['peak_cutout'] = 2 #Looking for peak +- this value around the expected peak (removing minor offsets)
    te_params['overwrite_wf'] = False #Flag if waveform extraction should be repeated (e.g. different cutouts)
    te_params['overwrite_tmp'] = True #Flag if templates should be recalculated if already existing

    qc_params = dict()
    qc_params['min_n_spikes'] = 500 #Minimum number of spikes to be detected for a unit for template extraction to take place
    qc_params['exclude_mua'] = True #Exclude units that were labelled multi unit activity by kilosort

    for dir in continuous_h5_dirs.values():
        h5_file_path = dir['h5_file_path']
        sorting_list_by_stream = sorting_lists_by_scan[0] #for debugging
        scantype = dir['scanType']
        if scantype in allowed_scan_types:
        #logger.info(f"Processing {h5_file_path}...")
            te_hp.extract_templates_from_sorting_dict(sorting_list_by_stream, h5_file_path, qc_params, te_params, stream_select = stream_select, just_load_waveforms = just_load_waveforms)

def get_merged_templates_by_unit_by_scan(continuous_h5_dirs, sorting_lists_by_scan, allowed_scan_types = ["AxonTracking"], 
                                         templates_dir='./AxonReconPipeline/data/temp_data/templates', stream_select = None, just_load_waveforms = False):
    templates_dir='./AxonReconPipeline/data/temp_data/templates'    
    templates_by_scan = extract_templates_hp(continuous_h5_dirs, 
                                      sorting_lists_by_scan, 
                                      allowed_scan_types = allowed_scan_types,
                                      templates_dir=templates_dir,
                                      stream_select = stream_select, just_load_waveforms= just_load_waveforms)
    return templates_by_scan

def generate_multirecordings_by_stream_by_scan(continuous_h5_dirs, 
                       allowed_scan_types=["AxonTracking"],
                       recording_dir = './AxonReconPipeline/data/temp_data/recordings/', 
                       stream_select = None):
    
    multi_recs_by_scan = []
    common_el_by_scan = []
    rec_list = [dir['h5_file_path'] for dir in continuous_h5_dirs.values()]
    for rec_path in rec_list:
        h5 = h5py.File(rec_path)
        stream_ids = list(h5['wells'].keys())
        if stream_select is not None:
            stream_ids = stream_ids[stream_select]
            if isinstance(stream_ids, str):
                stream_ids = [stream_ids]
        multirecs_by_stream = []
        common_el_by_stream = []
        recording_details = MPL.extract_recording_details(rec_path)
        scanType = recording_details[0]['scanType']
        if scanType in allowed_scan_types:
            for stream_id in stream_ids:
                multirecording, common_el = hp.concatenate_recording_slices(rec_path, 
                                                                    stream_id,
                                                                    recording_dir = recording_dir,
                                                                    #stream_break = stream_break
                                                                    )
                multirecs_by_stream.append(multirecording)
                common_el_by_stream.append(common_el)
                # if stream_break is not None:
                #     formatted = 'well{:03}'.format(stream_break)
                #     if stream_id == formatted:
                #         break
            multi_recs_by_scan.append(multirecs_by_stream)
            common_el_by_scan.append(common_el_by_stream)
    return multi_recs_by_scan, common_el_by_scan

import matplotlib.pyplot as plt

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

def select_folders_and_get_continuous_h5_dirs(debug_mode=False, pre_selected_folders=None):
    logger.info("Selecting folders to process, checking for continuity, and extracting the .h5 file path:")
    # Select the folders to process
    selected_folders = file_selector_main(pre_selected_folders= pre_selected_folders, debug_mode=debug_mode)
    logger.info(f"Selected folders: {selected_folders}")
    #Extract all .h5 files from the selected folders
    h5_dirs = MPL.extract_raw_h5_filepaths(selected_folders)
    #Extract chipIDs and record dates from the .h5 file paths
    informed_h5_dirs = MPL.extract_recording_details(h5_dirs)
    #Test for continuity (spiking data only turned off) in the .h5 files
    #Exclude non-continuous files
    continuous_h5_dirs = exclude_non_continuous_data(informed_h5_dirs)
    logger.info(f"Continuous data found in {len(continuous_h5_dirs)} files.")
    return continuous_h5_dirs

def main():
    debug_mode = True
    #testing Full Activity Scan and Network Scan...
    pre_selected_folders = [        
        #"/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/ActivityScan/000023",
        #"/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/Network/000024",
        #"/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/AxonTracking/000026",
        #testing continuity test with non-continuous data
        #"/mnt/disk20tb/PrimaryNeuronData/Maxone/SPTAN1_1/230725/16657/ActivityScan",

        #new comparison 09Feb24
        #"/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/FolicAcid_AxonTracking_05022024/FolicAcid_AxonTrack_Comparison/240205/M05506/AxonTracking/000003"
        #"/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/FolicAcid_AxonTracking_05022024/FolicAcid_AxonTrack_Comparison/240205/M07038/AxonTracking"
        ]
    #Organoid effort 31Mar24
    pre_selected_folders = [
        #reconstructions not working here...
        #'/mnt/ben-shalom_nas/irc/media/harddrive8tb/IMR90_T3_03262024_AP/IMR90_T3_03262024_AP/240328/16719/AxonTracking/000018'
        #same chip, fuller recording
        #'/mnt/ben-shalom_nas/irc/media/harddrive8tb/IMR90_T3_03262024_AP/IMR90_T3_03262024_AP/240328/16719/AxonTracking/000017',
        
        #Other organoid chips
        #'/mnt/ben-shalom_nas/irc/media/harddrive8tb/IMR90_T3_03262024_AP/IMR90_T3_03262024_AP/240326/18832/AxonTracking/000009',
        
        #Primary neuron Folic Acid
        #new comparison 09Feb24
        #"/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/FolicAcid_AxonTracking_05022024/FolicAcid_AxonTrack_Comparison/240205/M05506/AxonTracking/000003",
        #"/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/FolicAcid_AxonTracking_05022024/FolicAcid_AxonTrack_Comparison/240205/M07038/AxonTracking/000006"
        ]

    KCNT1_effort = [
        '/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/KCNT1_T4_C1_04122024/KCNT1_T4_C1_04122024/240503/M08034/AxonTracking/000082'
    ]
    continuous_h5_dirs = select_folders_and_get_continuous_h5_dirs(pre_selected_folders = KCNT1_effort, debug_mode = debug_mode)
    pprint(continuous_h5_dirs)

    #Set allowed scan types, axontracking alone is currently the best option
    allowed_scan_types = [
        #"ActivityScan", Dense activity scan only needed to define locations of electrodes in network and axon tracking scans
        "AxonTracking", 
        #"Network"
        ]
    recordings_dir = './AxonReconPipeline/data/temp_data/recordings'
    #only get sortings for well001
    stream_select = 0
    ##2. Generate multirecording objects for each stream for spikesorting and waveform extraction steps
    multi_recs_by_scan, common_el_by_scan = generate_multirecordings_by_stream_by_scan(continuous_h5_dirs,
                                                                                        allowed_scan_types,
                                                                                        recordings_dir,
                                                                                        stream_select = stream_select,
                                                                                        )