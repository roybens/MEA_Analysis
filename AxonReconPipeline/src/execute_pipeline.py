# General imports
import os
import re
import logging
from pprint import pprint
import h5py
import numpy as np
import matplotlib.pylab as plt
from scipy import io
import networkx as nx
from pathlib import Path
import json
import sys
import shutil

# Spikeinterface imports
import spikeinterface.sorters as ss
import spikeinterface.full as si

# Local imports
from file_selector import main as file_selector_main
from axon_trace_classes import axon_trace_objects
from reconstruct_axons import reconstruct_axon_trace_objects
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from MEAProcessingLibrary import mea_processing_library as MPL
import helper_functions as helper
import spike_sorting as hp #hornauerp, https://github.com/hornauerp/axon_tracking/blob/main/axon_tracking/spike_sorting.py
import template_extraction as te_hp

# Axon trace algo imports
import MEAutility as mu
from axon_velocity import *
from axon_velocity.models import load_cell
from axon_velocity.evaluation import *
from probeinterface import plotting as plotting_probe

# Logger Setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler()  # logs to console
stream_handler.setLevel(logging.DEBUG)
logger.addHandler(stream_handler)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - %(module)s.%(funcName)s')
stream_handler.setFormatter(formatter)

def extract_waveform_from_segments(h5_file_path, stream_id, segment_sorting, sel_unit_id, save_root, peak_cutout=2, align_cutout=True, upsample=2, rm_outliers=True, n_jobs=16, n_neighbors=10, overwrite_wf=False, overwrite_tmp = False):
    #full_path = segment_sorting._kwargs['recording_list'][0]._parent_recording._kwargs['recording']._kwargs['file_path']
    try:
        full_path = segment_sorting._kwargs['recording_or_recording_list'][0]._kwargs['parent_recording']._kwargs['recording']._kwargs['file_path']
    except:
        full_path = segment_sorting._kwargs['recording_or_recording_list'][0]._parent_recording._kwargs['folder_path']
        #segment_sorting._kwargs['recording_list'][0]._kwargs['parent_recording']._kwargs['folder_path']
        #get parent folder instead of file path
        full_path = os.path.dirname(full_path)
    
    if save_root is None:
        #save_root = os.path.dirname(full_path) 
        recording_details = MPL.extract_recording_details(full_path)
        date = recording_details[0]['date']
        chip_id = recording_details[0]['chipID']
        scanType = recording_details[0]['scanType']
        run_id = recording_details[0]['runID']
        save_root = f'./AxonReconPipeline/data/temp_data/waveforms/{date}/{chip_id}/{scanType}/{run_id}/{stream_id}'
    waveforms = te_hp.extract_waveforms(h5_file_path, segment_sorting, stream_id, save_root, n_jobs, overwrite_wf)
    return waveforms
def extract_all_waveforms(h5_file_path, stream_id, segment_sorting, pos, te_params, save_root=None):
    sel_unit_ids = segment_sorting.get_unit_ids()

    try:
        full_path = segment_sorting._kwargs['recording_or_recording_list'][0]._kwargs['parent_recording']._kwargs['recording']._kwargs['file_path']
    except:
        full_path = segment_sorting._kwargs['recording_or_recording_list'][0]._parent_recording._kwargs['folder_path']
        #segment_sorting._kwargs['recording_list'][0]._kwargs['parent_recording']._kwargs['folder_path']
        #get parent folder instead of file path
        full_path = os.path.dirname(full_path)
    
    if save_root is None:
        #save_root = os.path.dirname(full_path) 
        recording_details = MPL.extract_recording_details(full_path)
        date = recording_details[0]['date']
        chip_id = recording_details[0]['chipID']
        scanType = recording_details[0]['scanType']
        run_id = recording_details[0]['runID']
        save_root = f'./AxonReconPipeline/data/temp_data/waveforms/{date}/{chip_id}/{scanType}/{run_id}/{stream_id}'
    waveform_save_path = os.path.join(save_root, 'waveforms')
    if not os.path.exists(waveform_save_path):
        os.makedirs(waveform_save_path)
        
    for sel_unit_id in sel_unit_ids: 
        #template_save_file = os.path.join(template_save_path, str(sel_unit_id) + '.npy')
                ## debug
        # if sel_unit_id != 203: 
        # ##
        break
        
    if te_params['overwrite_tmp']:
        try:
            extract_waveform_from_segments(h5_file_path, stream_id, segment_sorting, sel_unit_id, save_root, **te_params)
            #grid = convert_to_grid(template_matrix, pos)
            #np.save(template_save_file, merged_template)
        except Exception as e:
            print(f'{e}')
def extract_waveforms_from_sorting_dict(sorting_list, h5_file_path, qc_params={}, te_params={}, stream_break = None):
    #rec_list = list(sorting_dict.keys())
    
    rec_list = h5_file_path
    if isinstance(h5_file_path, str):
        rec_list = [h5_file_path]    

    for rec_path in rec_list:
        #sorting_list = sorting_dict[rec_path]
        #sorting_list = sorting_dict
        #sorting_list.sort()

        for sorting in sorting_list:
        #for sorting_path in sorting_list:
            #sorting = si.KiloSortSortingExtractor(sorting_path)
            sorting_path = sorting._kwargs['folder_path']
            # #Temporary debug
            # #if 'Network' in stream_id, replace with 'AxonTracking'
            # if 'Network' in sorting_path:
            #     sorting_path = sorting_path.replace('Network', 'AxonTracking')
            stream_id = [p for p in sorting_path.split('/') if p.startswith('well')][0] #Find out which well this belongs to
            
            # print(stream_id)
            # #debug:
            # stream_number = int(stream_id[-1])
            # if stream_number > 0:
            #     continue
            # #debug end
            
            #rec_names, common_el, pos = ss.find_common_electrodes(rec_path, stream_id)
            multirecording, pos = hp.concatenate_recording_slices(rec_path, stream_id)

            #duration = int(h5['assay']['inputs']['record_time'][0].decode('UTF-8')) * n_recs #In case we want to use firing rate as criterion
            cleaned_sorting = te_hp.select_units(sorting, **qc_params)
            cleaned_sorting = si.remove_excess_spikes(cleaned_sorting, multirecording) #Relevant if last spike time == recording_length
            cleaned_sorting.register_recording(multirecording)
            segment_sorting = si.SplitSegmentSorting(cleaned_sorting, multirecording)
            extract_all_waveforms(h5_file_path, stream_id, segment_sorting, pos, te_params, save_root=None)

            if stream_break is not None:
                formatted = 'well{:03}'.format(stream_break)
                if stream_id == stream_break:
                    break
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
def get_list_of_sortings_by_scan(
        continuous_h5_dirs, 
        recordings_dir = './AxonReconPipeline/data/temp_data/recordings', 
        sortings_dir = './AxonReconPipeline/data/temp_data/sortings',   
        allowed_scan_types=["AxonTracking"], 
        stream_select = None):
    logger.info("Loading and merging the recordings as needed from the .h5 files. Generating spike sorting objects and merging as needed.")    
    logger.info(f"Allowed scan types: {allowed_scan_types}")
    #recordings_dir = './AxonReconPipeline/data/temp_data/recordings'
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
def extract_templates_hp(continuous_h5_dirs, 
                         sorting_lists_by_scan,
                         allowed_scan_types = ["AxonTracking"], 
                         templates_dir='./AxonReconPipeline/data/temp_data/templates',
                         stream_select = None,
                         just_load_waveforms = False):
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
def main(h5_parent_dirs, allowed_scan_types=['AxonTracking'], stream_select = None):

    '''Select folders and get continuous h5 filepaths'''
    continuous_h5_file_info = select_folders_and_get_continuous_h5_dirs(pre_selected_folders = h5_parent_dirs, debug_mode = debug_mode)

    '''Generate multirecordings objects by stream for spikesorting and waveform extraction steps'''
    recordings_dir = './AxonReconPipeline/data/temp_data/recordings'
    multi_recs_by_scan, common_el_by_scan = generate_multirecordings_by_stream_by_scan(
        continuous_h5_file_info,
        allowed_scan_types,
        recordings_dir,
        stream_select = stream_select,
        )
    
    '''Extract sorting objects per stream from sliced and cocatenated recordings (most/all fixed electrodes during axontracking scan)'''
    sortings_dir = './AxonReconPipeline/data/temp_data/sortings'
    sorting_lists_by_scan = get_list_of_sortings_by_scan(
        continuous_h5_file_info, 
        recordings_dir = recordings_dir,
        sortings_dir = sortings_dir, 
        allowed_scan_types = allowed_scan_types,
        stream_select = stream_select,
        )
    
    '''Extract Waveforms from sorting objects'''
    for i, sorting_list_by_stream in enumerate(sorting_lists_by_scan):
        extract_waveforms_from_sorting_dict(
            sorting_list_by_stream, 
            continuous_h5_file_info[i]['h5_file_path'], 
            qc_params=qc_params, 
            te_params=te_params, 
            stream_break = None)

    '''Extract Templates'''
    templates_dir='./AxonReconPipeline/data/temp_data/templates'  
    templates_by_scan = get_merged_templates_by_unit_by_scan(continuous_h5_file_info, 
                                                            sorting_lists_by_scan,
                                                            allowed_scan_types=allowed_scan_types,
                                                            templates_dir = templates_dir,
                                                            stream_select = stream_select, 
                                                            just_load_waveforms = load_waveforms)

    #Define Algorithm Params
    params = get_default_graph_velocity_params()
    options = ['axon_velocity', 'organoid']
    option = options[0]
    if option == 'axon_velocity':
        # change params (fig6 Buccino axon_velocity paper example)
        params['detect_threshold'] = 0.01
        params['kurt_threshold'] = 0.1
        params['peak_std_threshold'] = 0.8
        params['upsample'] = 5
        params['neighbor_radius'] = 100
        params['r2_threshold'] = 0.8

    elif option == 'organoid':
        ## Modulating params for organoid data

        # General
        params['upsample'] = 5 #1  # Upsample waveforms by this factor, upsample = 1 means no upsampling
        params['min_selected_points'] = 20 #30  # Minimum number of points to be selected for graph propagation, typically 30
        params['verbose'] = True  # Verbose output

        # Channel selection
        params['detect_threshold'] = 0.01  # Detection threshold (with respect to channel featuring maximal signal)
        params['detection_type'] = 'relative'  # Whether to use an 'absolute' or 'relative' detection threshold
        params['kurt_threshold'] = 0.05 #0.3  # Kurtosis threshold below which a channel is discarded
        params['peak_std_threshold'] = 0.5 #1  # Peak time standard deviation threshold in ms below which a channel is discarded
        params['init_delay'] = 0 #0.1  # Initial delay in seconds (with respect to maximum channel) below which a channel is discarded
        params['peak_std_distance'] = 50 #30  # Distance in µm to select channel neighborhood to compute peak time standard deviation
        params['remove_isolated'] = False #True  # If True, isolated channels are removed from selection

        # Graph
        params['init_amp_peak_ratio'] = 0.2  # Scalar value that weighs the contribution of the amplitude and the peak latency for h_init (ω_init in equation (2))
        params['max_distance_for_edge'] = 200 #100  # Maximum distance in µm between channels to create a graph edge
        params['max_distance_to_init'] = 500 #200  # Maximum distance in µm between a channel and the init_channel to create a graph edge
        params['n_neighbors'] = 5 #3  # Maximum number of edges that one channel can connect to below which an axonal branch is discarded
        params['distance_exp'] = 2  # Exponent for distance computation (e in equation (3))
        params['edge_dist_amp_ratio'] = 0.5 #0.3  # Relative weight between distance and amplitude to select neighbor nodes for graph edges

        # Axonal reconstruction
        params['min_path_length'] = 50 #100  # Minimum axon path length in µm to include an axonal branch
        params['min_path_points'] = 4 #5  # Minimum number of channels in an axon path to include an axonal branch
        params['neighbor_radius'] = 150 #100  # Radius in µm to exclude neighboring channels around an identified path
        params['min_points_after_branching'] = 3  # Minimum number of points after a branching to avoid pruning

        # Path cleaning/velocity estimation
        params['mad_threshold'] = 8  # Threshold in median absolute deviations on the fit error to consider points as outliers in the velocity estimation
        params['split_paths'] = True  # If True, the final path splitting step is enabled
        params['max_peak_latency_for_splitting'] = 0.5  # If a jump in the peak latencies of a path exceeds this value, the path can be split in sub-paths
        params['r2_threshold'] = 0.7 #0.9  # R^2 threshold for velocity linear fit below which an axon branch is discarded
        params['r2_threshold_for_outliers'] = 0.9 #0.98  # R^2 threshold below which outliers are detected and removed
        params['min_outlier_tracking_error'] = 75 #50  # Tracking error in µm above which a point can be considered an outlier and removed
    
    def get_extremum(template, locations):
        """
        Get the extremum of the template.
        """
        # Get the extremum index along the first dim of the template, aligning with channel locations
        extremum_idx = np.unravel_index(np.argmax(np.abs(template)), template.shape)[0]
        extremum = locations[extremum_idx]
        return extremum


    stream_select = 0
    templates_dir='./AxonReconPipeline/data/temp_data/templates'
    for i, continuous_h5_dir in continuous_h5_file_info.items():
        recording_details = MPL.extract_recording_details(continuous_h5_dir['h5_file_path'])
        date = recording_details[0]['date']
        chip_id = recording_details[0]['chipID']
        scanType = recording_details[0]['scanType']
        run_id = recording_details[0]['runID']
        h5 = h5py.File(continuous_h5_dir['h5_file_path'], 'r')
        stream_ids = list(h5['wells'].keys())
        for stream_id in stream_ids:
            if stream_select is not None:
                if stream_id != f'well{stream_select:03}':
                    continue
                else:
                    print(f'Processing stream {stream_id}')
                    break

        temp_dir = os.path.dirname(templates_dir)
        sorting_dir = Path(temp_dir) /'sortings'/ date / chip_id / scanType / run_id / stream_id / 'sorter_output'
        sorting = ss.Kilosort2_5Sorter._get_result_from_folder(sorting_dir)
        recon_dir = Path(temp_dir) /'reconstructions'/ date / chip_id / scanType / run_id / stream_id
        successful_recons = {}  
        unit_ids = sorting.get_unit_ids()
        successful_recons[str(recon_dir)] = {}
        successful_recons[str(recon_dir)]["successful_units"] = {}
        units = []
        branch_ids = []
        velocities = []
        path_lengths = []
        r2s = []
        extremums =[]
        for unit_id in unit_ids:

            print(f'\nProcessing unit {unit_id}')
            # debug
            #if unit_id != 51: continue
            #
                    
            merged_template_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}.npy'
            merged_channel_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}_channels.npy'
            merged_template_filled_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}_filled.npy'
            merged_channel_filled_dir = Path(templates_dir) / date / chip_id / scanType / run_id / stream_id / 'templates' / f'{unit_id}_channels_filled.npy'
            try:              

                merged_template = np.load(merged_template_dir)
                merged_channel_loc = np.load(merged_channel_dir)
                merged_template_filled = np.load(merged_template_filled_dir)
                merged_channel_filled_loc = np.load(merged_channel_filled_dir)

                #merged_template = np.load('/mnt/disk15tb/adam/git_workspace/MEA_Analysis/AxonReconPipeline/data/temp_data/templates/240118/M06844/AxonTracking/000026/well001/partial_templates/seg3_unit17.npy')
                #merged_channel_loc = np.load('/mnt/disk15tb/adam/git_workspace/MEA_Analysis/AxonReconPipeline/data/temp_data/templates/240118/M06844/AxonTracking/000026/well001/partial_templates/seg3_unit17_channels.npy')
            except:
                print(f'No merged template found for {unit_id}')
                continue

            ## Channels
            template = merged_template[-1]
            transformed_template = template.T
            #mirror matrix about x-axis to match maxwell orientation
            #transformed_template = np.flip(transformed_template, axis=0)
            
            template_filled = merged_template_filled[-1]
            transformed_template_filled = template_filled.T
            #mirror matrix about x-axis to match maxwell orientation\
            #transformed_template_filled = np.flip(transformed_template_filled, axis=0)
            
            ## Locations
            trans_loc = merged_channel_loc[-1]
            #multiply all y values by -1 to match maxwell orientation
            trans_loc = np.array([[loc[0], loc[1]*-1] for loc in trans_loc])

            trans_loc_filled = merged_channel_filled_loc[-1]
            #multiply all y values by -1 to match maxwell orientation
            trans_loc_filled = np.array([[loc[0], loc[1]*-1] for loc in trans_loc_filled])
            #mirror matrix about x-axis to match maxwell orientation
            #trans_loc_filled = np.flip(trans_loc_filled, axis = 0)

            ##
            gtrs = {}
            fs = 10000
            plot_dir = Path(recon_dir) / f'unit_{unit_id}' / 'axon_tracking_plots'
            if os.path.exists(plot_dir):
                #print(f'Plot directory already exists for unit {unit_id}, skipping...')
                #continue
                shutil.rmtree(plot_dir)
            key = 0

            # Generate amplitude_map figure
            try:
                fig_amp = plt.figure(figsize=(10, 5))
                ax_amp = fig_amp.add_subplot(111)
                #set_trace()
                ax_amp = plot_amplitude_map(transformed_template_filled, 
                                            trans_loc_filled, 
                                            log=False, 
                                            ax=ax_amp, 
                                            cmap="PRGn", 
                                            colorbar=True, 
                                            colorbar_orientation="horizontal")
                ax_amp.set_title(f"Amplitude", fontsize=20)

                #Save figure
                fig_name = "amplitude_map.png"
                fig_path = plot_dir / fig_name
                if not os.path.exists(plot_dir):
                    os.makedirs(plot_dir)
                plt.savefig(fig_path, dpi=300, bbox_inches='tight')
                plt.close()

                #accept_units_amp_map.append(key)
            except Exception as e:
                #reject_units_amp_map.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate amplitude_map, error: {e}")

            # Generate peak latency map figure
            try:
                fig_peaks = plt.figure(figsize=(10, 5))
                ax_peaks = fig_peaks.add_subplot(111)
                ax_peaks = plot_peak_latency_map(transformed_template_filled, 
                                                trans_loc_filled, fs=fs, log=False, 
                                                ax=ax_peaks, colorbar=True, 
                                                colorbar_orientation="horizontal")
                ax_peaks.set_title(f"Peak latency", fontsize=20)

                #Save figure
                fig_name = "peak_latency_map.png"
                fig_path = plot_dir / fig_name
                plt.savefig(fig_path, dpi=300, bbox_inches='tight')
                plt.close()

                #accept_units_peak_latency_map.append(key)
            except Exception as e:
                #reject_units_peak_latency_map.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate peak_latency_map, error: {e}")

            # Select channels and generate selected channels figures
            def plot_selected_channels(fig_title, selected_channels, locations, filename):
                fig = plt.figure(figsize=(10, 5))
                ax = fig.add_subplot(111)
                ax.set_title(fig_title, fontsize=20)
                ax.plot(locations[:, 0], locations[:, 1], marker=".", color="grey", ls="", alpha=0.2)
                ax.plot(locations[selected_channels, 0], 
                        locations[selected_channels, 1], marker=".", color="k", ls="", alpha=0.5)
                ax.axis("off")
                fig_name = filename
                fig_path = plot_dir / fig_name
                plt.savefig(fig_path, dpi=300, bbox_inches='tight')
                plt.close()
            try:
                # Build axon tracking object
                verbose_bool = True
                gtr0 = None #initialize gtr0
                gtr0 = GraphAxonTracking(transformed_template, trans_loc, fs, 
                                        #verbose=verbose_bool, 
                                        **params)           

                # Select channels
                gtr0.select_channels()                

                # Use the function for each block
                plot_selected_channels(f"Selected after detection threshold: {gtr0._detect_threshold} ms", 
                                    np.array(list(gtr0._selected_channels_detect)), gtr0.locations, "detection_threshold.png")
                plot_selected_channels(f"Selected after detection threshold: {gtr0._kurt_threshold} ms", 
                                    np.array(list(gtr0._selected_channels_kurt)), gtr0.locations, "kurt_threshold.png")
                plot_selected_channels(f"Selected after peak std threshold: {gtr0._peak_std_threhsold} ms", 
                                        np.array(list(gtr0._selected_channels_peakstd)), gtr0.locations, "peak_std_threshold.png")
                plot_selected_channels(f"Selected after init delay threshold: {gtr0._init_delay} ms", 
                                        np.array(list(gtr0._selected_channels_init)), gtr0.locations, "init_delay_threshold.png")
                plot_selected_channels("Selected after all thresholds", gtr0.selected_channels, gtr0.locations, "all_thresholds.png")

                #if all thresholds are passed, compute graph propagation velocity
                try:
                    gtr0 = compute_graph_propagation_velocity(transformed_template, trans_loc, fs, 
                        #verbose=False, 
                        **params)
                except Exception as e:
                    print(f"unit {unit_id} failed to compute_graph_propagation_velocity, error: {e}")
                #accept_units_selected_channels.append(key)
            except Exception as e:
                #reject_units_selected_channels.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate selected_channels, error: {e}")              
            
            # Generate axon analytics
            try:

                for i, br in enumerate(gtr0.branches):
                    path = br["channels"]
                    velocity = br["velocity"]
                    r2 = br["r2"]
                    length = gtr0.compute_path_length(path)
                    units.append(unit_id)
                    branch_ids.append(i)
                    velocities.append(velocity)
                    path_lengths.append(length)
                    r2s.append(r2)
                    extremums.append(get_extremum(transformed_template, trans_loc))

                #accept_units_axon_reconstruction_heuristics.append(key)
            except Exception as e:
                #reject_units_axon_reconstruction_heuristics.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate axon_analytics, error: {e}")

            # Generate axon_reconstruction_heuristics figure
            try:
                gtr0.build_graph()
                gtr0.find_paths()
                gtr0._verbose = 1

                fig_graph = plt.figure(figsize=(10, 7))
                fig_graph = gtr0.plot_graph(node_search_labels=False, fig=fig_graph, cmap_nodes="viridis", cmap_edges="YlGn")

                fig_name = "axon_reconstruction_heuristics.png"
                fig_path = plot_dir / fig_name
                plt.savefig(fig_path, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"unit {unit_id} successfully generated axon_reconstruction_heuristics")

                #accept_units_axon_reconstruction_heuristics.append(key)
            except Exception as e:
                #reject_units_axon_reconstruction_heuristics.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate axon_reconstruction_heuristics, error: {e}")

            # # Generate axon_reconstruction figures
            # try:
            #     fig_graph = plt.figure(figsize=(12, 7))
            #     #ax_raw = fig_graph.add_subplot(111)                
            #     #axpaths_raw = ax_raw
            #     axpaths_raw = gtr0.plot_branches(cmap="tab20", fig = fig_graph,
            #                                      #plot_bp=True, plot_neighbors=True, 
            #                                         #plot_full_template=True, 
            #                                      #ax=axpaths_raw
            #                                      )
            #     #axpaths_raw.legend(fontsize=6)
            #     #axpaths_raw.set_title("All Branches")

            #     fig_name = "axon_reconstruction_all.png"
            #     fig_path = plot_dir / fig_name
            #     plt.savefig(fig_path, dpi=600, bbox_inches='tight')
            #     plt.close()
                
            #     successful_recons[str(recon_dir)]["successful_units"][str(unit_id)] = {}
            #     #successful_recons[str(recon_dir)]["successful_units"][str(unit_id)]["stream_id"] = stream_id
            #     #successful_recons[str(recon_dir)]["successful_units"][str(unit_id)]["fig_path"] = fig_path
            #     print(f"unit {unit_id} successfully generated axon_reconstruction")
            #     #successful_recons[continuous_h5_dir]["template"] = template
            #     #successful_recons[continuous_h5_dir]["merged_channel_loc"] = merged_channel_loc
            #     # gtrs[key] = compute_graph_propagation_velocity(template.T, merged_channel_loc[i], fs, 
            #     #                                         verbose=False, **params)

            #     #accept_units_axon_reconstruction.append(key)
            #     #flag = True
            #     #break
            # except Exception as e:
            #     #reject_units_axon_reconstruction.append(key)
            #     plt.close()
            #     print(f"unit {unit_id} failed to generate all axon_reconstruction, error: {e}")
            
            # Generate axon_reconstruction figures
            try:
                fig_graph = plt.figure(figsize=(12, 7))
                ax_raw = fig_graph.add_subplot(111)                
                axpaths_raw = ax_raw
                axpaths_raw = gtr0.plot_clean_branches(cmap="tab20", plot_bp=True, 
                                                    #plot_neighbors=True, 
                                                    plot_full_template=True, ax=axpaths_raw)
                axpaths_raw.legend(fontsize=6)
                axpaths_raw.set_title("Clean Branches")

                fig_name = "axon_reconstruction_clean.png"
                fig_path = plot_dir / fig_name
                plt.savefig(fig_path, dpi=600, bbox_inches='tight')
                plt.close()
                
                successful_recons[str(recon_dir)]["successful_units"][str(unit_id)] = {}
                #successful_recons[str(recon_dir)]["successful_units"][str(unit_id)]["stream_id"] = stream_id
                #successful_recons[str(recon_dir)]["successful_units"][str(unit_id)]["fig_path"] = fig_path
                print(f"unit {unit_id} successfully generated clean axon_reconstruction")
                #successful_recons[continuous_h5_dir]["template"] = template
                #successful_recons[continuous_h5_dir]["merged_channel_loc"] = merged_channel_loc
                # gtrs[key] = compute_graph_propagation_velocity(template.T, merged_channel_loc[i], fs, 
                #                                         verbose=False, **params)

                #accept_units_axon_reconstruction.append(key)
                #flag = True
                #break
            except Exception as e:
                #reject_units_axon_reconstruction.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate clean axon_reconstruction, error: {e}")

            # Generate axon_reconstruction figures
            try:
                fig_graph = plt.figure(figsize=(12, 7))
                ax_raw = fig_graph.add_subplot(111)                
                axpaths_raw = ax_raw
                axpaths_raw = gtr0.plot_raw_branches(cmap="tab20", plot_bp=True, plot_neighbors=True, 
                                                    plot_full_template=True, ax=axpaths_raw)
                axpaths_raw.legend(fontsize=6)
                axpaths_raw.set_title("Raw Branches")

                fig_name = "axon_reconstruction_raw.png"
                fig_path = plot_dir / fig_name
                plt.savefig(fig_path, dpi=600, bbox_inches='tight')
                plt.close()
                
                successful_recons[str(recon_dir)]["successful_units"][str(unit_id)] = {}
                #successful_recons[str(recon_dir)]["successful_units"][str(unit_id)]["stream_id"] = stream_id
                #successful_recons[str(recon_dir)]["successful_units"][str(unit_id)]["fig_path"] = fig_path
                print(f"unit {unit_id} successfully generated axon_reconstruction")
                #successful_recons[continuous_h5_dir]["template"] = template
                #successful_recons[continuous_h5_dir]["merged_channel_loc"] = merged_channel_loc
                # gtrs[key] = compute_graph_propagation_velocity(template.T, merged_channel_loc[i], fs, 
                #                                         verbose=False, **params)

                #accept_units_axon_reconstruction.append(key)
                #flag = True
                #break
            except Exception as e:
                #reject_units_axon_reconstruction.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate raw axon_reconstruction, error: {e}")
                
            # Generate axon_reconstruction_velocities figure
            try:
                fvel, ax_vel = plt.subplots(figsize=(9, 15))
                ax_vel.axis("off")
                fvel = gtr0.plot_velocities(cmap="tab20", plot_outliers=True, fig=fvel, markersize=12, markersize_out=18, fs=20)
                plt.xticks(fontsize=16)
                plt.yticks(fontsize=16)

                fig_name = "axon_reconstruction_velocities.png"
                fig_path = plot_dir / fig_name
                plt.savefig(fig_path, dpi=600, bbox_inches='tight')
                plt.close()
                print(f"unit {unit_id} successfully generated axon_reconstruction_velocities")

                #accept_units_axon_reconstruction_velocities.append(key)
            except Exception as e:
                #reject_units_axon_reconstruction_velocities.append(key)
                plt.close()
                print(f"unit {unit_id} failed to generate axon_reconstruction_velocities, error: {e}")

        import pandas as pd
        df_mea1k = pd.DataFrame({"unit_ids": units, "unit location": extremums, "branch_id": branch_ids, "velocity": velocities,
                            "length": path_lengths, "r2": r2s})
        df_mea1k.to_csv(Path(recon_dir) / "axon_analytics.csv", index=False)
        
        # def convert_keys_to_int(d):
        #     if isinstance(d, dict):
        #         return {int(k) if isinstance(k, np.int64) else k: convert_keys_to_int(v) for k, v in d.items()}
        #     elif isinstance(d, list):
        #         return [convert_keys_to_int(v) for v in d]
        #     else:
        #         return d
        with open(f'{recon_dir}/successful_recons.json', 'w') as f:
            json.dump(successful_recons, f)

if __name__ == "__main__":
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
    
    #Set allowed scan types, axontracking alone is currently the best option
    allowed_scan_types = [
        #"ActivityScan", Dense activity scan only needed to define locations of electrodes in network and axon tracking scans
        "AxonTracking", 
        #"Network"
        ]

    h5_parent_dirs = KCNT1_effort
    
    '''options'''
    # homo
    stream_select = 0
    # hetero
    stream_select = 3
    load_waveforms = True #Load waveforms from disk instead of extracting them, if possible
    main(h5_parent_dirs, allowed_scan_types = allowed_scan_types, stream_select = stream_select)