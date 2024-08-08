''' This file contains template-related functions that are used in the main pipeline functions. '''

'''imports'''
import sys
import os
import h5py
import spikeinterface.full as si
import numpy as np
from tqdm import tqdm

''' Local imports '''
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from MEAProcessingLibrary import mea_processing_library as MPL
import AxonReconPipeline.src.lib_sorting_functions as sorter
import AxonReconPipeline.src.lib_waveform_functions as waveformer

''' Template merging functions '''
from AxonReconPipeline.src.func_merge_templates import merge_templates
from AxonReconPipeline.src.func_merge_templates import merge_templates_concurrent 

default_n_jobs = 4

'''functions'''
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_template_segment(sel_unit_id, rec_name, waveforms, save_root, sel_idx, logger):
    # Get the waveform extractor object
    try:
        seg_we = waveforms[rec_name]['waveforms']
    except KeyError: 
        if logger is not None: logger.error(f'{rec_name} does not exist')
        else: print(f'{rec_name} does not exist')
        return None, None
    
    # Get the template for the selected unit
    try: template = seg_we.get_template(sel_unit_id)
    except Exception as e: 
        if logger is not None: logger.error(f'Could not get templates for {sel_unit_id} in {rec_name}: {e}')
        else: print(f'Could not get templates for {sel_unit_id} in {rec_name}: {e}')
        return None, None
    
    # Check if the template is empty
    if np.isnan(template).all():
        if logger is not None: logger.warning(f'Unit {sel_unit_id} in segment {sel_idx} is empty. Skipping.')
        else: print(f'Unit {sel_unit_id} in segment {sel_idx} is empty. Skipping.')
        return None, None
    
    # Get the channel locations
    locs = seg_we.get_channel_locations()
    dir_path = os.path.join(save_root, 'template_segments')
    file_name = f'seg{sel_idx}_unit{sel_unit_id}.npy'
    channel_loc_file_name = f'seg{sel_idx}_unit{sel_unit_id}_channels.npy'
    channel_loc_save_file = os.path.join(dir_path, channel_loc_file_name)
    template_save_file = os.path.join(dir_path, file_name)
    
    # Warnings for NaN values
    if not np.isnan(template).all() and np.isnan(template).any():
        if logger is not None: logger.warning(f'Unit {sel_unit_id} in segment {sel_idx} has NaN values')
        else: print(f'Unit {sel_unit_id} in segment {sel_idx} has NaN values')
            
    return (rec_name, {'template': template, 'path': template_save_file}), (rec_name, {'channel_locations': locs, 'path': channel_loc_save_file})

def extract_template_segments(sel_unit_id, h5_path, stream_id, waveforms, save_root=None, logger=None, max_workers=12):
    if logger is not None: logger.info(f'Extracting template segments for unit {sel_unit_id}')
    else: print(f'Extracting template segments for unit {sel_unit_id}')
    full_path = h5_path
    h5 = h5py.File(full_path)
    rec_names = list(h5['wells'][stream_id].keys())
    
    template_segments = {}
    channel_locations = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(process_template_segment, sel_unit_id, rec_name, waveforms, save_root, sel_idx, logger)
            for sel_idx, rec_name in enumerate(rec_names)
        ]
        for future in as_completed(futures):
            try:
                template_result, loc_result = future.result()
                if template_result and loc_result:
                    rec_name_t, template_data = template_result
                    rec_name_l, loc_data = loc_result
                    template_segments[rec_name_t] = template_data
                    channel_locations[rec_name_l] = loc_data
                    if logger is not None: logger.debug(f'Processed segment {rec_name_t} for unit {sel_unit_id}')
                    else: print(f'Processed segment {rec_name_t} for unit {sel_unit_id}')
            except Exception as e:
                if logger is not None: logger.error(f'Error processing unit {sel_unit_id} in segment {rec_name_l}: {e}')
                else: print(f'Error processing unit {sel_unit_id} in segment {rec_name_l}: {e}')
                continue
    if logger is not None: logger.info(f'Finished extracting {len(template_segments)} template segments for unit {sel_unit_id}')
    else: print(f'Finished extracting {len(template_segments)} template segments for unit {sel_unit_id}')
    return template_segments, channel_locations  

def merge_template_segments(unit_segments, channel_locations, logger=None):
    if logger is not None: logger.info(f'Merging partial templates')
    else: print(f'Merging partial templates')        
    template_list = [tmp['template'] for rec_name, tmp in unit_segments.items()]
    channel_locations_list = [ch_loc['channel_locations'] for rec_name, ch_loc in channel_locations.items()]
    merged_template, merged_channel_loc, merged_template_filled, merged_channel_locs_filled = merge_templates(template_list, channel_locations_list, logger=logger)
    
    # aw20July2024 - This doesnt currently work. TODO: Fix this later. not important.
    # merged_template, merged_channel_loc, merged_template_filled, merged_channel_locs_filled = merge_templates_concurrent(template_list, channel_locations_list, logger=logger)
    
    merged_template = merged_template[0]
    merged_channel_loc = merged_channel_loc[0]
    merged_template_filled = merged_template_filled[0]
    merged_channel_locs_filled = merged_channel_locs_filled[0]
    return merged_template, merged_channel_loc, merged_template_filled, merged_channel_locs_filled

def get_time_derivative(merged_template, merged_template_filled, unit='seconds', sampling_rate=10000, axis=0):
    """
    Computes the time derivative of the input arrays along the specified axis.
    
    Parameters:
    - merged_template: numpy.ndarray
        The first input array where the y-axis describes channels and x-axis describes voltage(time) signals.
    - merged_template_filled: numpy.ndarray
        The second input array where the y-axis describes channels and x-axis describes voltage(time) signals.
    - unit: str, optional
        The time unit for the derivative, either 'seconds' or 'ms'. Default is 'seconds'.
    - sampling_rate: int, optional
        The sampling rate in samples per second. Default is 10000.
    - axis: int, optional
        The axis along which to compute the derivative. Default is 0 (assuming time signals along y-axis).
    
    Returns:
    - d_merged_template: numpy.ndarray
        The time derivative of merged_template.
    - d_merged_template_filled: numpy.ndarray
        The time derivative of merged_template_filled.
    """
    
    if unit == 'seconds':
        delta_t = 1.0 / sampling_rate
    elif unit == 'ms':
        delta_t = 1.0 / (sampling_rate / 1000)
    else:
        raise ValueError("Invalid unit. Choose either 'seconds' or 'ms'.")

    d_merged_template = np.diff(merged_template, axis=axis) / delta_t
    d_merged_template_filled = np.diff(merged_template_filled, axis=axis) / delta_t
    return d_merged_template, d_merged_template_filled

def extract_merged_templates(h5_path, stream_id, segment_sorting, waveforms, te_params, save_root=None, unit_limit=None, logger=None, template_bypass=False, unit_select=None, **temp_kwargs):
    def get_existing_unit_ids():
        template_files = os.listdir(template_save_path)
        sel_unit_ids = []
        for f in template_files:
            try: sel_unit_ids.append(int(f.split('.')[0]))
            except ValueError: continue  # Skip files that don't start with an integer          
        return sel_unit_ids
    
    def add_to_dict(unit_segments, channel_locations, merged_template, merged_channel_loc, merged_template_filled, 
                    merged_channel_locs_filled, template_save_file, channel_loc_save_file, template_save_file_fill, 
                    channel_loc_save_file_fill, d_merged_template, d_merged_template_filled):
        unit_templates[sel_unit_id] = {
            'template_segments': unit_segments, 
            'channel_locations': channel_locations,
            'merged_template_path': template_save_file,
            'merged_template': merged_template,
            'dvdt_merged_template': d_merged_template,
            'merged_channel_loc_path': channel_loc_save_file,
            'merged_channel_loc': merged_channel_loc,
            'merged_template_filled_path': template_save_file_fill,
            'dvdt_merged_template_filled': d_merged_template_filled,
            'merged_template_filled': merged_template_filled,
            'merged_channel_locs_filled_path': channel_loc_save_file_fill,
            'merged_channel_locs_filled': merged_channel_locs_filled, 
        }
    
    if logger is not None: logger.info(f'Extracting merged templates for {h5_path}')
    else: print(f'Extracting merged templates for {h5_path}')        

    h5_details = MPL.extract_recording_details(h5_path)
    h5_details = h5_details[0]
    date = h5_details['date']
    chip_id = h5_details['chipID']
    scanType = h5_details['scanType']
    run_id = h5_details['runID']
    if save_root is not None: template_save_path = save_root+f'/{date}/{chip_id}/{scanType}/{run_id}/{stream_id}'
    if template_bypass: 
        sel_unit_ids = get_existing_unit_ids()
        assert len(sel_unit_ids) > 0, 'No existing templates found. Generating New Templates.'    
    else: sel_unit_ids = segment_sorting.get_unit_ids()
    if not os.path.exists(template_save_path): os.makedirs(template_save_path)
        
    unit_templates = {}
    unit_count = 0
    for sel_unit_id in tqdm(sel_unit_ids):
        if unit_select is not None and sel_unit_id != unit_select: continue # manage unit selection option
        
        template_save_file = os.path.join(template_save_path, str(sel_unit_id) + '.npy')
        dvdt_template_save_file = os.path.join(template_save_path, str(sel_unit_id) + '_dvdt.npy')
        channel_loc_save_file = os.path.join(template_save_path, str(sel_unit_id) + '_channels.npy')
        template_save_file_fill = os.path.join(template_save_path, str(sel_unit_id) + '_filled.npy')
        dvdt_template_save_file_fill = os.path.join(template_save_path, str(sel_unit_id) + '_dvdt_filled.npy')
        channel_loc_save_file_fill = os.path.join(template_save_path, str(sel_unit_id) + '_channels_filled.npy')

        try:
            assert te_params.get('load_merged_templates', False) == True, 'load_merged_templates is set to False. Generating New Templates.'
            assert te_params.get('overwrite_tmp', False) == False, 'overwrite_tmp is set to True. Generating New Templates.'
            
            if logger is not None: logger.info(f'Loading merged template for unit {sel_unit_id}')
            else: print(f'Loading merged template for unit {sel_unit_id}')                
            merged_template = np.load(template_save_file)
            dvdt_template_save_file = np.load(dvdt_template_save_file)
            merged_channel_loc = np.load(channel_loc_save_file)
            merged_template_filled = np.load(template_save_file_fill)
            dvdt_template_save_file_fill = np.load(dvdt_template_save_file_fill)
            merged_channel_locs_filled = np.load(channel_loc_save_file_fill)
            add_to_dict(
                "loading segments is not currently supported", 
                "loading segment channels is not currently supported", 
                merged_template, merged_channel_loc, merged_template_filled, merged_channel_locs_filled,
                template_save_file, channel_loc_save_file, template_save_file_fill, channel_loc_save_file_fill,
                dvdt_template_save_file, dvdt_template_save_file_fill
            )
            unit_count += 1
            if unit_limit is not None and unit_count >= unit_limit: break
            continue 
        except Exception as e: 
            if logger is not None: logger.warning(f'loading template for unit {sel_unit_id}:\n{e}. Generating New Templates.')
            else: print(f'Error loading template for unit {sel_unit_id}:\n{e}. Generating New Templates.'); pass 

        try:
            if logger is not None: logger.info(f'Extracting template segments for unit {sel_unit_id}')
            else: print(f'Extracting template segments for unit {sel_unit_id}')
                
            n_jobs = te_params.get('n_jobs', 8)
            unit_segments, channel_locations = extract_template_segments(sel_unit_id, h5_path, stream_id, waveforms, save_root=save_root, logger=logger, max_workers=n_jobs)            
            if logger is not None: logger.info(f'Merging partial templates for unit {sel_unit_id}')
            else: print(f'Merging partial templates for unit {sel_unit_id}')
                
            merged_template, merged_channel_loc, merged_template_filled, merged_channel_locs_filled = merge_template_segments(unit_segments, channel_locations, logger=logger)
            dvdt_merged_template, dvdt_merged_template_filled = get_time_derivative(merged_template, merged_template_filled)
            
            add_to_dict(
                unit_segments, channel_locations, merged_template, merged_channel_loc, merged_template_filled, 
                merged_channel_locs_filled, template_save_file, channel_loc_save_file, template_save_file_fill, 
                channel_loc_save_file_fill, dvdt_merged_template, dvdt_merged_template_filled
            )
            
            if te_params.get('save_merged_templates', False):
                np.save(template_save_file, merged_template)
                np.save(dvdt_template_save_file, dvdt_merged_template)
                np.save(channel_loc_save_file, merged_channel_loc)
                np.save(template_save_file_fill, merged_template_filled)
                np.save(dvdt_template_save_file_fill, dvdt_merged_template_filled)
                np.save(channel_loc_save_file_fill, merged_channel_locs_filled)
            if logger is not None: logger.info(f'Merged template saved to {save_root}/templates')
            else: print(f'Merged template saved to {save_root}/templates')
            unit_count += 1
            if unit_limit is not None and unit_count >= unit_limit: break                
        except Exception as e:
            if logger is not None: logger.error(f'Unit {sel_unit_id} encountered the following error: {e}')
            else: print(f'Unit {sel_unit_id} encountered the following error:\n {e}')
    return unit_templates

def extract_templates(multirec, sorting, waveforms, h5_path, stream_id, save_root=None, te_params={}, qc_params={}, unit_limit=None, logger=None, template_bypass=False, **temp_kwargs):
    if logger is not None: logger.info(f'Extracting templates for {h5_path}')
    else: print(f'Extracting templates for {h5_path}')
    if template_bypass:
        try: 
            unit_templates = extract_merged_templates(h5_path, stream_id, None, None, te_params, save_root=save_root, unit_limit=unit_limit, logger=logger, template_bypass=True, **temp_kwargs)
            return unit_templates
        except Exception as e: 
            if logger is not None: logger.error(f'Error loading templates via bypass:\n{e}')
            else: print(f'Error loading templates via bypass:{e}')
    cleaned_sorting = waveformer.select_units(sorting, **qc_params)
    cleaned_sorting = si.remove_excess_spikes(cleaned_sorting, multirec) 
    cleaned_sorting.register_recording(multirec)
    segment_sorting = si.SplitSegmentSorting(cleaned_sorting, multirec)
    unit_templates = extract_merged_templates(h5_path, stream_id, segment_sorting, waveforms, te_params, save_root=save_root, unit_limit=unit_limit, logger=logger, **temp_kwargs)
    return unit_templates