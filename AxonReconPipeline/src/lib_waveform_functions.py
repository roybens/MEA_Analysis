''' This file contains waveform-related functions that are used in the main pipeline functions. '''

'''imports'''
import sys
import os
import h5py
import spikeinterface.full as si
import spikeinterface.sorters as ss
import shutil
import numpy as np
from tqdm import tqdm
from pathlib import Path
import time

''' Local imports '''
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from MEAProcessingLibrary import mea_processing_library as MPL
import lib_sorting_functions as sorter

'''logging setup'''
import logging
logger = logging.getLogger(__name__) #Create a logger
logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler()  # Create handlers, logs to console
stream_handler.setLevel(logging.DEBUG) # Set level of handlers
logger.addHandler(stream_handler) # Add handlers to the logger
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - %(module)s.%(funcName)s') # Create formatters and add it to handlers
stream_handler.setFormatter(formatter)

'''default settings'''
default_recording_dir = './AxonReconPipeline/data/temp_data/recordings'
default_sorting_dir = './AxonReconPipeline/data/temp_data/sortings'
default_n_jobs = 4

def extract_waveforms_of_units_by_segment(h5_file_path, stream_id, segment_sorting, save_root, peak_cutout=2, align_cutout=True, upsample=2, rm_outliers=True, n_jobs=16, n_neighbors=10, overwrite_wf=False, overwrite_tmp = False):
    def get_assay_information(rec_path):
        h5 = h5py.File(rec_path)
        pre, post, well_id = -1, -1, 0
        while pre <= 0 or post <= 0: #some failed axon trackings give negative trigger_post values, so we try different wells
            well_name = list(h5['wells'].keys())[well_id]
            rec_name = list(h5['wells'][well_name].keys())[well_id]
            pre = h5['wells'][well_name][rec_name]['groups']['routed']['trigger_pre'][0]
            post = h5['wells'][well_name][rec_name]['groups']['routed']['trigger_post'][0]
            well_id += 1
            
        return [pre, post]
    def extract_waveforms(h5_file_path, segment_sorting, stream_id, save_root, n_jobs, overwrite_wf):
        try:
            save_root = segment_sorting._kwargs['recording_or_recording_list'][0]._kwargs['parent_recording']._kwargs['recording']._kwargs['file_path']
            save_root = os.path.dirname(save_root)
        except:
            save_root = segment_sorting._kwargs['recording_or_recording_list'][0]._parent_recording._kwargs['folder_path']
            #segment_sorting._kwargs['recording_list'][0]._kwargs['parent_recording']._kwargs['folder_path']
            #get parent folder instead of file path
            save_root = os.path.dirname(save_root)
            #replace 'recordings' with waveforms
            save_root = save_root.replace('recordings', 'waveforms')
        full_path = h5_file_path
        cutout = [x / (segment_sorting.get_sampling_frequency()/1000) for x in get_assay_information(full_path)] #convert cutout to ms
        h5 = h5py.File(full_path)
        rec_names = list(h5['wells'][stream_id].keys())
        #logger.info(f'Extracting waveforms. Stream: {stream_id}, {len(rec_names)} segments')
        
        for sel_idx, rec_name in enumerate(rec_names):
            wf_path = os.path.join(save_root, 'waveforms', 'seg' + str(sel_idx))        
            rec = si.MaxwellRecordingExtractor(full_path,stream_id=stream_id,rec_name=rec_name)
            chunk_size = np.min([10000, rec.get_num_samples()]) - 100 #Fallback for ultra short recordings (too little activity)
            rec_centered = si.center(rec, chunk_size=chunk_size)            
            seg_sort = si.SelectSegmentSorting(segment_sorting, sel_idx)
            seg_sort = si.remove_excess_spikes(seg_sort, rec_centered)
            seg_sort.register_recording(rec_centered) 
            #Try loading waveforms, if they don't exist or overwrite_wf is True, extract them
            if not overwrite_wf:
                try: 
                    seg_we = si.WaveformExtractor.load(wf_path, with_recording=True, sorting = seg_sort)                
                    # seg_we = si.load_waveforms(wf_path, with_recording = True, sorting = seg_sort)
                    if len(seg_sort.get_unit_ids()) == len(seg_we.unit_ids):
                        logger.info(f'Waveforms loaded from {wf_path} already exist')
                        logger.info(f'Loaded waveforms from units in Stream: {stream_id}, Seg: {sel_idx}')
                        #logger.info(f'Units: {seg_sort.get_unit_ids()}')
                        continue
                        #break
                    else:
                        raise ValueError('Unit IDs do not match')    
                except:
                    if os.path.exists(wf_path):
                        shutil.rmtree(wf_path)
                    pass

            if not os.path.exists(wf_path) or overwrite_wf:                           
                logger.info(f'Extracting waveforms to {wf_path}')
                os.makedirs(wf_path, exist_ok=True)
                seg_we = si.WaveformExtractor.create(rec_centered, seg_sort,
                                                    wf_path, 
                                                    allow_unfiltered=True,
                                                    remove_if_exists=True)
                seg_we.set_params(ms_before=cutout[0], ms_after=cutout[1], return_scaled = True)
                #logger.info(f'Extracting waveforms. Stream: {stream_id}, Unit: {sel_unit_id}, Seg: {sel_idx}, {n_jobs} jobs')
                logger.info(f'Extracting waveforms from units in Stream: {stream_id}, Seg: {sel_idx}')
                #logger.info(f'Units: {seg_sort.get_unit_ids()}')
                seg_we.run_extract_waveforms(n_jobs=n_jobs)
    
    # #full_path = segment_sorting._kwargs['recording_list'][0]._parent_recording._kwargs['recording']._kwargs['file_path']
    # try: full_path = segment_sorting._kwargs['recording_or_recording_list'][0]._kwargs['parent_recording']._kwargs['recording']._kwargs['file_path']
    # except:
    #     full_path = segment_sorting._kwargs['recording_or_recording_list'][0]._parent_recording._kwargs['folder_path']
    #     #segment_sorting._kwargs['recording_list'][0]._kwargs['parent_recording']._kwargs['folder_path']
    #     #get parent folder instead of file path
    #     full_path = os.path.dirname(full_path)

    assert save_root is not None, 'save_root must be provided'
    full_path = save_root.replace('waveforms', 'recordings')
    
    if save_root is None:
        #save_root = os.path.dirname(full_path) 
        recording_details = MPL.extract_recording_details(full_path)
        date = recording_details[0]['date']
        chip_id = recording_details[0]['chipID']
        scanType = recording_details[0]['scanType']
        run_id = recording_details[0]['runID']
        save_root = f'./AxonReconPipeline/data/temp_data/waveforms/{date}/{chip_id}/{scanType}/{run_id}/{stream_id}'
    waveforms = extract_waveforms(h5_file_path, segment_sorting, stream_id, save_root, n_jobs, overwrite_wf)
    return waveforms
def extract_waveforms_from_stream(h5_file_path, stream_id, segment_sorting, pos, te_params, save_root=None):
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
        
    #for sel_unit_id in sel_unit_ids:         
    #if te_params['overwrite_tmp']:
    try:
        extract_waveforms_of_units_by_segment(h5_file_path, stream_id, segment_sorting, save_root, **te_params)
        #grid = convert_to_grid(template_matrix, pos)
        #np.save(template_save_file, merged_template)
    except Exception as e:
        print(f'{e}')
def select_units(sorting, min_n_spikes=500, exclude_mua=True):
    if exclude_mua:
        ks_label = sorting.get_property('KSLabel')
        mua_idx = ks_label == 'mua'
    else:
        mua_idx = np.full((sorting.get_num_units(),), False, dtype='bool')

    
    n_spikes = [len(sorting.get_unit_spike_train(x)) for x in sorting.get_unit_ids()]

    
    bad_n_spikes_idx = np.array(n_spikes) < min_n_spikes
    bad_idx = mua_idx | bad_n_spikes_idx
    bad_id = [i for i, x in enumerate(bad_idx) if x]
    
    cleaned_sorting = sorting.remove_units(bad_id)
    
    return cleaned_sorting
def extract_all_waveforms(continuous_h5_info, sorting_lists_by_scan = None, qc_params={}, te_params={}, stream_break = None):
    #rec_list = list(sorting_dict.keys())
    #h5_file_path = h5_file_path['h5_file_path']

    allowed_scan_types = ['AxonTracking']
    if stream_break is not None: stream_select = stream_break
    stream_select = stream_break
    if sorting_lists_by_scan is None: sorting_lists_by_scan = sorter.spikesort_recordings(continuous_h5_info, sortings_dir = default_sorting_dir, allowed_scan_types = allowed_scan_types, stream_select = stream_select,)

    for h5 in continuous_h5_info.values():
        h5_file_path = h5['h5_file_path']

        for i, sorting_list_by_stream in enumerate(sorting_lists_by_scan):
        
            rec_list = h5_file_path
            if isinstance(h5_file_path, str):
                rec_list = [h5_file_path]    

            for rec_path in rec_list:
                #sorting_list = sorting_dict[rec_path]
                #sorting_list = sorting_dict
                #sorting_list.sort()

                for sorting in sorting_list_by_stream:
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
                    multirecording, pos = sorter.concatenate_recording_slices(rec_path, stream_id)

                    #duration = int(h5['assay']['inputs']['record_time'][0].decode('UTF-8')) * n_recs #In case we want to use firing rate as criterion
                    cleaned_sorting = select_units(sorting, **qc_params)
                    cleaned_sorting = si.remove_excess_spikes(cleaned_sorting, multirecording) #Relevant if last spike time == recording_length
                    cleaned_sorting.register_recording(multirecording)
                    segment_sorting = si.SplitSegmentSorting(cleaned_sorting, multirecording)
                    extract_waveforms_from_stream(h5_file_path, stream_id, segment_sorting, pos, te_params, save_root=None)

                    if stream_break is not None:
                        formatted = 'well{:03}'.format(stream_break)
                        if stream_id == stream_break:
                            break