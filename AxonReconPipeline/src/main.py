#general imports
import os
import platform
import re
import logging
import os, time, h5py

#spikeinterface imports
import spikeinterface
import spikeinterface.sorters as ss
import spikeinterface.full as si
from spikeinterface.widgets import plot_probe_map

#local imports
from file_selector import main as file_selector_main
#from AxonReconPipeline.src.archive.extract_raw_assay_data import main as extract_raw_assay_data_main
from generate_recording_and_sorting_objects import main as generate_recording_and_sorting_objects
from axon_trace_classes import axon_trace_objects
from reconstruct_axons import reconstruct_axon_trace_objects
import mea_processing_library as MPL
import lib_helper_functions as helper
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
                h5_file_path, 
                allowed_scan_types=[
                                    #"ActivityScan", 
                                    "AxonTracking", 
                                    #"Network"
                                    ], 
                recordings_dir = './AxonReconPipeline/data/temp_data/recordings',
                sortings_dir = './AxonReconPipeline/data/temp_data/sortings',
                stream_select = None): 
    
    """ Extract sorting objects per stream from sliced and cocatenated recordings (most/all fixed electrodes during axontracking scan)
    """
    
    #merged_recordings_by_scan = []
    #merged_sortings_by_scan = []
    #h5_file_paths = [dir['h5_file_path'] for dir in continuous_h5_dirs.values() if dir['scanType'] in allowed_scan_types]
    recording_details = MPL.extract_recording_details(h5_file_path)
    date = recording_details[0]['date']
    chip_id = recording_details[0]['chipID']
    scanType = recording_details[0]['scanType']
    run_id = recording_details[0]['runID']
    
    #Sorter params - dont really need to define this here, remove later
    sorter_params = si.get_default_sorter_params(si.Kilosort2_5Sorter)
    sorter_params['n_jobs'] = -1
    sorter_params['detect_threshold'] = 7
    sorter_params['minFR'] = 0.01
    sorter_params['minfr_goodchannels'] = 0.01
    sorter_params['keep_good_only'] = False
    sorter_params['do_correction'] = False
    #
    verbose = True
    #spikesorting_dir = './AxonReconPipeline/data/temp_data/sortings'
    #spikesorting_dir = os.path.abspath(spikesorting_dir)
    #spikesorting_dir = str(spikesorting_dir)
    spikesorting_root = sortings_dir+f'/{date}/{chip_id}/{scanType}/{run_id}'
    merged_sorting_list_by_stream = None
    
    #h5_file_paths = [dir['h5_file_path'] for dir in continuous_h5_dirs.values() if dir['scanType'] == "AxonTracking"]
    sorting_list_by_stream = hp.sort_recording_list(h5_file_path,
                                                    save_root=spikesorting_root,
                                                    #save_path_changes= 5,
                                                    sorter = 'kilosort2_5',
                                                    sorter_params = sorter_params,
                                                    verbose = verbose,
                                                    #scan_merge = False,
                                                    stream_select = stream_select) 
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

def get_list_of_sortings_by_scan(
        continuous_h5_dirs, 
        recordings_dir = './AxonReconPipeline/data/temp_data/recordings', 
        sortings_dir = './AxonReconPipeline/data/temp_data/sortings',   
        allowed_scan_types=["AxonTracking"], 
        stream_select = None):
    logger.info("Loading and merging the recordings as needed from the .h5 files. Generating spike sorting objects and merging as needed.")    
    logger.info(f"Allowed scan types: {allowed_scan_types}")
    recordings_dir = './AxonReconPipeline/data/temp_data/recordings'
    sorting_lists_by_scan = []
    for dir in continuous_h5_dirs.values():
        h5_file_path = dir['h5_file_path']
        scantype = dir['scanType']
        if scantype in allowed_scan_types:
            #logger.info(f"Processing {h5_file_path}...")
            sorting_list_by_stream = extract_sortings_per_stream(h5_file_path, 
                                                                 allowed_scan_types=allowed_scan_types, 
                                                                 recordings_dir=recordings_dir,
                                                                 sortings_dir = sortings_dir,  
                                                                 stream_select=stream_select)
            sorting_lists_by_scan.append(sorting_list_by_stream)
    return sorting_lists_by_scan

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

#from probe_map import ProbeMapWidget

#plot locations of common electrodes
def plot_electrodes(recs_ext_by_scan, recording_dir = None, auto_plot = True, name = None):
    #fig, axs = plt.subplots(len(multi_recs_by_scan), 1, figsize=(10, 10))
    
    widgets_by_scan = []
    for scan_i, multirec_by_stream in enumerate(recs_ext_by_scan):   
        widgets_by_stream = []
        for stream_i, rec_ext in enumerate(multirec_by_stream):
            widget = plot_probe_map(rec_ext)

            #extract recording details
            attributes = ['file_path', 'folder_path', 'recording_list', 'recording_list_parent_folder']
            rec_ext_path = None
            for attr in attributes:
                try:
                    if attr == 'recording_list':
                        try: rec_ext_path = rec_ext.recording_list[0]._parent_recording._kwargs['file_path']
                        except: rec_ext_path = rec_ext.recording_list[0]._parent_recording._kwargs['recording']._kwargs['file_path']
                    elif attr == 'recording_list_parent_folder':
                        rec_ext_path = rec_ext._kwargs['recording_list'][0]._kwargs['parent_recording']._kwargs['folder_path']
                        #get parent folder instead of file path
                        rec_ext_path = os.path.dirname(rec_ext_path)
                    else:
                        rec_ext_path = rec_ext._kwargs[attr]
                    break
                except KeyError:
                    continue

            if rec_ext_path is None:
                raise ValueError("Could not extract the folder path from the recording object.")

            if rec_ext_path is None:
                raise ValueError("Could not extract the folder path from the recording object.")
            recording_details = MPL.extract_recording_details(rec_ext_path)
            date = recording_details[0]['date']
            chip_id = recording_details[0]['chipID']
            scan_type = recording_details[0]['scanType']
            run_id = recording_details[0]['runID']
            try: stream_id = multirec_by_stream[0].stream_name
            except: 
                stream_id = multirec_by_stream[0].recording_list[0]._kwargs['parent_recording']._kwargs['folder_path']
                stream_id = os.path.dirname(stream_id)
                stream_id = os.path.basename(stream_id)
            
            #fig title
            widget.ax.set_title(f"{date}-{chip_id}-{stream_id}")          
            
            #save to recording_dir if not None`
            if recording_dir is not None:              
                recording_dir = recording_dir + f"/{date}/{chip_id}/{scan_type}/{run_id}/{stream_id}"
                if not os.path.exists(recording_dir):
                    os.makedirs(recording_dir)
                #plt.savefig(recording_dir+f"/probe_map_{stream_id}.png")
                if name is not None:
                    widget.figure.savefig(recording_dir+f"/{name}_probe_map_{stream_id}.png", dpi=600)
                else: 
                    widget.figure.savefig(recording_dir+f"/probe_map_{stream_id}.png", dpi=600)
            
            #append to widgets_by_stream
            widgets_by_stream.append(widget)
        #append to widgets_by_scan
        widgets_by_scan.append(widgets_by_stream)
            
    if auto_plot:
        #plt.show()
        plt.tight_layout()
        plt.show()

    return widgets_by_scan

def main():
    # Clear the terminal
    clear_terminal()

    ##1. Select the folders to process, check for continuity, and extract the .h5 file paths
    continuous_h5_dirs = select_folders_and_get_continuous_h5_dirs()

    ##2. Generate multirecording objects for each stream for spikesorting and waveform extraction steps
    multi_recs_by_scan, common_el_by_scan = generate_multirecordings_by_stream_by_scan(continuous_h5_dirs)

    
    ##3. Extract sorting objects per stream from sliced and cocatenated recordings (most/all fixed electrodes during axontracking scan)
    sorting_lists_by_scan = get_list_of_sortings_by_scan(continuous_h5_dirs)

    ##4. Extract Waveforms
                                                                                        

    ##5. Extract Templates
    #Extract recording details from the .h5 file path
    templates_by_scan = get_merged_templates_by_unit_by_scan(continuous_h5_dirs, sorting_lists_by_scan)
    
            
    logger.info("done.")
    
if __name__ == "__main__":
    main()