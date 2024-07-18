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
import lib_waveform_functions as waveformer

'''default settings'''
default_recording_dir = './AxonReconPipeline/data/temp_data/recordings'
default_sorting_dir = './AxonReconPipeline/data/temp_data/sortings'
default_template_dir = './AxonReconPipeline/data/temp_data/templates'  

default_n_jobs = 4

'''functions'''
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_rec_name(sel_unit_id, rec_name, waveforms, save_root, sel_idx, logger):
    
    # Get the waveform extractor object
    try: seg_we = waveforms[rec_name]['waveforms']
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
    # if '14' in rec_name:
    #     print()
    locs = seg_we.get_channel_locations()
    dir_path = os.path.join(save_root, 'template_segments')
    file_name = f'seg{sel_idx}_unit{sel_unit_id}.npy'
    channel_loc_file_name = f'seg{sel_idx}_unit{sel_unit_id}_channels.npy'
    channel_loc_save_file = os.path.join(dir_path, channel_loc_file_name)
    template_save_file = os.path.join(dir_path, file_name)
    
    #Warnings for NaN values
    if not np.isnan(template).all() and np.isnan(template).any():
        if logger is not None: logger.warning(f'Unit {sel_unit_id} in segment {sel_idx} has NaN values')
        else: print(f'Unit {sel_unit_id} in segment {sel_idx} has NaN values')
            
    return (rec_name, {'template': template, 'path': template_save_file}), (rec_name, {'channel_locations': locs, 'path': channel_loc_save_file})

def extract_template_segments(sel_unit_id, h5_path, stream_id, waveforms, save_root=None, logger=None, max_workers=12):
    full_path = h5_path
    h5 = h5py.File(full_path)
    rec_names = list(h5['wells'][stream_id].keys())
    
    template_segments = {}
    channel_locations = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(process_rec_name, sel_unit_id, rec_name, waveforms, save_root, sel_idx, logger)
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
                #debug
                #break

    return template_segments, channel_locations

def fill_template(merged_templates, merged_channel_loc, merged_count_at_channel_by_sample, logger=None):
    def get_pitch():
        distances = []
        for i in range(len(merged_channel_loc[-1])):
            for j in range(i + 1, len(merged_channel_loc[-1])):
                distances.append(np.linalg.norm(np.array(merged_channel_loc[-1][i]) - np.array(merged_channel_loc[-1][j])))
                min_distance = np.min(distances)
                if min_distance == 17.5:
                    return min_distance
        assert min_distance % 17.5 == 0, "Minimum distance is not a multiple of 17.5 um."
        return min_distance
    
    logger.info('Filling in missing channels...')
    x_min = np.min([loc[0] for loc in merged_channel_loc[-1]])
    x_max = np.max([loc[0] for loc in merged_channel_loc[-1]])
    y_min = np.min([loc[1] for loc in merged_channel_loc[-1]])
    y_max = np.max([loc[1] for loc in merged_channel_loc[-1]])

    min_distance = get_pitch()

    x_grid = np.arange(x_min, x_max + min_distance, min_distance)
    y_grid = np.arange(y_min, y_max + min_distance, min_distance)
    xx, yy = np.meshgrid(x_grid, y_grid)
    grid = np.array([xx.flatten(), yy.flatten()]).T

    existing_channels_set = set(map(tuple, merged_channel_loc[-1]))
    expected_channels_set = set(map(tuple, grid))
    missing_channel_set = [ch for idx, ch in enumerate(expected_channels_set) if tuple(ch) not in existing_channels_set]
    assert len(missing_channel_set) + len(merged_channel_loc[-1]) == len(grid), "Missing channels do not match."
    assert len(missing_channel_set) + len(merged_channel_loc[-1]) <= 26400, "Too many channels for Maxwell MEAs."

    logger.info(f'Processing {len(missing_channel_set)} unrecorded channels')
    new_size = merged_templates[-1][0].shape[0] + len(missing_channel_set)
    new_merged_templates = [np.array([np.zeros(new_size, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
    for l, t in enumerate(merged_templates[-1]):
        new_array = np.zeros(new_size, dtype=t.dtype)
        new_array[:t.shape[0]] = t
        new_merged_templates[-1][l] = new_array

    num_channels = new_merged_templates[-1][0].shape[0]
    num_samples = len(merged_templates[-1])
    new_merged_count_at_channel_by_sample = [np.array([np.zeros(num_channels, dtype='float32')] * num_samples, dtype='float32')]
    for sample_idx, merge_counts_in_sample in enumerate(merged_count_at_channel_by_sample[-1]):
        for channel_idx, count_num in enumerate(merged_count_at_channel_by_sample[-1][sample_idx]):
            new_merged_count_at_channel_by_sample[-1][sample_idx][channel_idx] = count_num

    same_size = merged_channel_loc[-1][0].shape[0]
    new_length = len(merged_channel_loc[-1]) + len(missing_channel_set)
    new_merged_channels = [np.array([np.zeros(same_size, dtype='float64')] * new_length, dtype='float64')]
    for l, c in enumerate(merged_channel_loc[-1]):
        new_tuple = np.zeros(same_size, dtype=c.dtype)
        new_tuple = c
        new_merged_channels[-1][l] = new_tuple

    for idx, j in enumerate(missing_channel_set):
        position = merged_templates[-1][0].shape[0] + idx
        channel_loc_to_append = tuple(missing_channel_set[idx])
        tuple_to_append = np.array(channel_loc_to_append)
        new_merged_channels[-1][position] = tuple_to_append

    for footprint in new_merged_templates[-1]:
        assert np.array(footprint[-len(missing_channel_set)]).all() == 0, "Filler channels are not zeroed out."
    for footprint in new_merged_count_at_channel_by_sample[-1]:
        assert np.array(footprint[-len(missing_channel_set)]).all() == 0, "Filler channels are not zeroed out."
    for ch in missing_channel_set:
        assert ch in new_merged_channels[-1], "Missing channels are not in the channel locations."

    merged_templates = new_merged_templates
    merged_channel_loc = new_merged_channels
    merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample

    return merged_templates, merged_channel_loc, merged_count_at_channel_by_sample
   
def merge_templates(templates, channel_locations, logger=None):
    merged_templates = []
    merged_channel_loc = []
    verbose = True

    if logger is not None: logger.debug('Starting merge...')
    else: print('Starting merge...')
        
    test_channel_measure = 0
    for i in range(len(templates)):
        channel_loc = channel_locations[i]

        if len(merged_templates) == 0:
            incoming_template = templates[i]
            if np.isnan(incoming_template).all(): continue
            assert np.isnan(incoming_template).any() == False, "First template contains NaN values."

            merged_templates.append(templates[i])
            merged_channel_loc.append(channel_loc)
            test_channel_measure += len(channel_loc)
            num_channels = merged_templates[-1][0].shape[0]
            merged_count_at_channel_by_sample = [np.array([np.zeros(num_channels, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
            new_merged_count_at_channel_by_sample = [np.array([np.zeros(num_channels, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
            for sample_idx, merge_counts_in_sample in enumerate(merged_count_at_channel_by_sample[-1]):
                for channel_idx, count_num in enumerate(merged_count_at_channel_by_sample[-1][sample_idx]):
                    new_merged_count_at_channel_by_sample[-1][sample_idx][channel_idx] = count_num + 1
            merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample 
        else:
            if logger is not None:
                logger.debug(f'Merging matching channels...')
                logger.debug(f'Total merged channels: {len(merged_channel_loc[-1])}')            
                logger.debug(f'Total channels: {test_channel_measure}')  
            else:
                print(f'Merging matching channels...')
                print(f'Total merged channels: {len(merged_channel_loc[-1])}')            
                print(f'Total channels: {test_channel_measure}')  
                
            new_template = []
            new_channel_loc = []
            incoming_template = templates[i]
            if logger is not None:
                logger.debug(f'Processing {len(channel_loc)} channels')
            else:
                print(f'Processing {len(channel_loc)} channels')

            existing_channels_set = set(map(tuple, merged_channel_loc[-1]))
            incoming_channels_set = set(map(tuple, channel_loc))
            overlapping_channels = existing_channels_set.intersection(incoming_channels_set)
            
            idx_incoming_channels_in_overlap = [idx for idx, ch in enumerate(channel_loc) if tuple(ch) in overlapping_channels]
            idx_existing_channels_in_overlap = [idx for idx, ch in enumerate(merged_channel_loc[-1]) if tuple(ch) in overlapping_channels]
            idx_new_incoming_channels = [idx for idx, ch in enumerate(channel_loc) if tuple(ch) not in overlapping_channels]
            assert len(idx_incoming_channels_in_overlap) + len(idx_new_incoming_channels) == len(channel_loc), "Indices do not match."
            test_channel_measure += len(idx_new_incoming_channels)

            channels_already_merged = merged_channel_loc[-1][idx_existing_channels_in_overlap]
            merging_channels = channel_loc[idx_incoming_channels_in_overlap]
            merging_channels_ordered = [ch for idx, ch in enumerate(channels_already_merged) if tuple(ch) in merging_channels]
            merging_channels_ordered = np.array(merging_channels_ordered, dtype='float64')
            assert len(channels_already_merged) == len(merging_channels_ordered), "The number of channels does not match."
            for idx in range(len(channels_already_merged)):
                assert np.array_equal(channels_already_merged[idx], merging_channels_ordered[idx]), "Matching channels are not the same."
            shared_ordered_channels = merging_channels_ordered                   
           
            if logger is not None: logger.debug(f'Processing {len(channels_already_merged)} matching channels')
            else: print(f'Processing {len(channels_already_merged)} matching channels')
                
            if idx_existing_channels_in_overlap or idx_incoming_channels_in_overlap:               
                if logger is not None: logger.debug(f'Processing matching channels across footprints, template {i+1}/{len(templates)}...')
                else: print(f'Processing matching channels across footprints, template {i+1}/{len(templates)}...')     
                    
                merged_sample_at_channel = merged_templates[-1][:, idx_existing_channels_in_overlap]
                sample_to_be_merged = incoming_template[:, idx_incoming_channels_in_overlap]

                merge_count = merged_count_at_channel_by_sample[-1][:, idx_existing_channels_in_overlap]             
                total_samples_at_channel = merge_count + 1
                total_samples_at_channel[np.isnan(sample_to_be_merged)] = total_samples_at_channel[np.isnan(sample_to_be_merged)] - 1
                sample_to_be_merged[np.isnan(sample_to_be_merged)] = 0                               
                weighted_average = np.zeros_like(total_samples_at_channel)
                non_zero_mask = total_samples_at_channel != 0
                weighted_average[non_zero_mask] = (merged_sample_at_channel[non_zero_mask] * merge_count[non_zero_mask] + sample_to_be_merged[non_zero_mask]) / total_samples_at_channel[non_zero_mask]
                assert not np.isnan(weighted_average).any(), "Weighted average contains NaN values."
                
                merged_templates[-1][:, idx_existing_channels_in_overlap] = weighted_average
                merged_count_at_channel_by_sample[-1][:, idx_existing_channels_in_overlap] = total_samples_at_channel
                        
            if idx_new_incoming_channels:         
                if logger is not None: logger.debug(f'Processing {len(idx_new_incoming_channels)} non-matching channels')
                else: print(f'Processing {len(idx_new_incoming_channels)} non-matching channels')
                    
                new_size = merged_templates[-1][0].shape[0] + len(idx_new_incoming_channels)
                new_merged_templates = [np.array([np.zeros(new_size, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
                for l, t in enumerate(merged_templates[-1]):
                    new_array = np.zeros(new_size, dtype=t.dtype)
                    new_array[:t.shape[0]] = t
                    new_merged_templates[-1][l] = new_array

                num_channels = new_merged_templates[-1][0].shape[0]
                new_merged_count_at_channel_by_sample = [np.array([np.zeros(num_channels, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
                for sample_idx, merge_counts_in_sample in enumerate(merged_count_at_channel_by_sample[-1]):
                    for channel_idx, count_num in enumerate(merged_count_at_channel_by_sample[-1][sample_idx]):
                        new_merged_count_at_channel_by_sample[-1][sample_idx][channel_idx] = count_num

                same_size = merged_channel_loc[-1][0].shape[0]
                new_length = len(merged_channel_loc[-1]) + len(idx_new_incoming_channels)
                new_merged_channels = [np.array([np.zeros(same_size, dtype='float64')] * new_length, dtype='float64')]
                for l, c in enumerate(merged_channel_loc[-1]):
                    new_tuple = np.zeros(same_size, dtype=c.dtype)
                    new_tuple = c
                    new_merged_channels[-1][l] = new_tuple

                for k in tqdm(range(len(templates[i])), desc=f'Processing non-matching channels in each footprint, template {i+1}/{len(templates)}'):
                    for idx, j in enumerate(idx_new_incoming_channels):
                        position = merged_templates[-1][0].shape[0] + idx
                        sample_to_append = templates[i][k][j]
                        if np.isnan(sample_to_append):
                            sample_to_append = 0
                            new_merged_count_at_channel_by_sample[-1][k][position] = 0
                        else:
                            new_merged_count_at_channel_by_sample[-1][k][position] += 1
                        new_merged_templates[-1][k][position] = sample_to_append               
                        channel_loc_to_append = tuple(channel_loc[j])
                        tuple_to_append = np.array(channel_loc_to_append)
                        new_merged_channels[-1][position] = tuple_to_append          

                merged_templates = new_merged_templates
                merged_channel_loc = new_merged_channels
                merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample                   
      
            else:
                raise Exception("Something funky is going on")
    if logger is not None: logger.debug('done merging')
    else: print('done merging')
        
    if logger is not None: logger.debug(f'Total merged channels: {len(merged_channel_loc[-1])}')
    else: print(f'Total merged channels: {len(merged_channel_loc[-1])}')
        
    for i in range(len(merged_templates[-1])): assert len(merged_templates[-1][i]) <= 26400
    assert len(merged_channel_loc[-1]) <= 26400

    merged_templates_filled, merged_channel_loc_filled, merged_count_at_channel_by_sample_filled = fill_template(merged_templates, merged_channel_loc, merged_count_at_channel_by_sample, logger=logger)

    if logger is not None: logger.debug(f'Total merged channels after fill: {len(merged_channel_loc_filled[-1])}')
    else: print(f'Total merged channels after fill: {len(merged_channel_loc_filled[-1])}')
        
    for i in range(len(merged_templates_filled[-1])): assert len(merged_templates_filled[-1][i]) <= 26400
    assert len(merged_channel_loc_filled[-1]) <= 26400
    return merged_templates, merged_channel_loc, merged_templates_filled, merged_channel_loc_filled  

def merge_template_segments(unit_segments, channel_locations, logger=None):
    if logger is not None: logger.info(f'Merging partial templates')
    else: print(f'Merging partial templates')        
    template_list = [tmp['template'] for rec_name, tmp in unit_segments.items()]
    channel_locations_list = [ch_loc['channel_locations'] for rec_name, ch_loc in channel_locations.items()]
    merged_template, merged_channel_loc, merged_template_filled, merged_channel_locs_filled = merge_templates(template_list, channel_locations_list, logger=logger)
    merged_template = merged_template[0]
    merged_channel_loc = merged_channel_loc[0]
    merged_template_filled = merged_template_filled[0]
    merged_channel_locs_filled = merged_channel_locs_filled[0]
    return merged_template, merged_channel_loc, merged_template_filled, merged_channel_locs_filled

def extract_merged_templates(h5_path, stream_id, segment_sorting, waveforms, te_params, save_root=None, unit_limit=None, logger=None, template_bypass=False):
    
    def get_existing_unit_ids():
        template_files = os.listdir(template_save_path)
        sel_unit_ids = []
        for f in template_files:
            try: sel_unit_ids.append(int(f.split('.')[0]))
            except ValueError: continue  # Skip files that don't start with an integer          
        return sel_unit_ids
    
    def add_to_dict(unit_segments, channel_locations, merged_template, merged_channel_loc, merged_template_filled, merged_channel_locs_filled, 
                    template_save_file, channel_loc_save_file, template_save_file_fill, channel_loc_save_file_fill):
        unit_templates[sel_unit_id] = {
            'template_segments': unit_segments, 
            'channel_locations': channel_locations,
            'merged_template_path': template_save_file,
            'merged_template': merged_template,
            'merged_channel_loc_path': channel_loc_save_file,
            'merged_channel_loc': merged_channel_loc,
            'merged_template_filled_path': template_save_file_fill,
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

        template_save_file = os.path.join(template_save_path, str(sel_unit_id) + '.npy')
        channel_loc_save_file = os.path.join(template_save_path, str(sel_unit_id) + '_channels.npy')
        template_save_file_fill = os.path.join(template_save_path, str(sel_unit_id) + '_filled.npy')
        channel_loc_save_file_fill = os.path.join(template_save_path, str(sel_unit_id) + '_channels_filled.npy')

        try:
            #if te_params.get('load_merged_templates', False):
            assert te_params.get('load_merged_templates', False) == True, 'load_merged_templates is set to False. Generating New Templates.'
            assert te_params.get('overwrite_tmp', False) == False, 'overwrite_tmp is set to True. Generating New Templates.'
            
            if logger is not None: logger.info(f'Loading merged template for unit {sel_unit_id}')
            else: print(f'Loading merged template for unit {sel_unit_id}')                
            merged_template = np.load(template_save_file)
            merged_channel_loc = np.load(channel_loc_save_file)
            merged_template_filled = np.load(template_save_file_fill)
            merged_channel_locs_filled = np.load(channel_loc_save_file_fill)
            add_to_dict("loading segments is not currently supported", 
                        "loading segment channels is not currently supported", 
                        merged_template, merged_channel_loc, merged_template_filled, merged_channel_locs_filled,
                        template_save_file, channel_loc_save_file, template_save_file_fill, channel_loc_save_file_fill)
            unit_count += 1
            if unit_limit is not None and unit_count >= unit_limit: break
            continue 
            #else: raise Exception('load_merged_templates is set to False. Generating New Templates.')
        except Exception as e: 
            if logger is not None: logger.info(f'Error loading template for unit {sel_unit_id}:\n{e}. Generating New Templates.')
            else: print(f'Error loading template for unit {sel_unit_id}:\n{e}. Generating New Templates.'); pass 

        try:
            if logger is not None: logger.info(f'Extracting template segments for unit {sel_unit_id}')
            else: print(f'Extracting template segments for unit {sel_unit_id}')
                
            n_jobs = te_params.get('n_jobs', 8)
            unit_segments, channel_locations = extract_template_segments(sel_unit_id, h5_path, stream_id, waveforms, save_root=save_root, logger=logger, max_workers=n_jobs)            
            if logger is not None: logger.info(f'Merging partial templates for unit {sel_unit_id}')
            else: print(f'Merging partial templates for unit {sel_unit_id}')
                
            merged_template, merged_channel_loc, merged_template_filled, merged_channel_locs_filled = merge_template_segments(unit_segments, channel_locations, logger=logger)
            add_to_dict(unit_segments, channel_locations, merged_template, merged_channel_loc, merged_template_filled, merged_channel_locs_filled,
                        template_save_file, channel_loc_save_file, template_save_file_fill, channel_loc_save_file_fill)
            if te_params.get('save_merged_templates', False):
                np.save(template_save_file, merged_template)
                np.save(channel_loc_save_file, merged_channel_loc)
                np.save(template_save_file_fill, merged_template_filled)
                np.save(channel_loc_save_file_fill, merged_channel_locs_filled)
            if logger is not None: logger.info(f'Merged template saved to {save_root}/templates')
            else: print(f'Merged template saved to {save_root}/templates')
            unit_count += 1
            if unit_limit is not None and unit_count >= unit_limit: break                
        except Exception as e:
            if logger is not None: logger.error(f'Unit {sel_unit_id} encountered the following error: {e}')
            else: print(f'Unit {sel_unit_id} encountered the following error:\n {e}')               
        

    return unit_templates

def extract_templates(multirec, sorting, waveforms, h5_path, stream_id, save_root=None, te_params={}, qc_params={}, unit_limit=None, logger=None, template_bypass=False):
    if logger is not None: logger.info(f'Extracting templates for {h5_path}')
    else:print(f'Extracting templates for {h5_path}')
    if template_bypass:
        try: 
            unit_templates = extract_merged_templates(h5_path, stream_id, None, None, te_params, save_root=save_root, unit_limit=unit_limit, logger=logger, template_bypass=True)
            return unit_templates
        except Exception as e: 
            if logger is not None: logger.error(f'Error loading templates via bypass:\n{e}')
            else: print(f'Error loading templates via bypass:{e}')
    cleaned_sorting = waveformer.select_units(sorting, **qc_params)
    cleaned_sorting = si.remove_excess_spikes(cleaned_sorting, multirec) 
    cleaned_sorting.register_recording(multirec)
    segment_sorting = si.SplitSegmentSorting(cleaned_sorting, multirec)
    unit_templates = extract_merged_templates(h5_path, stream_id, segment_sorting, waveforms, te_params, save_root=save_root, unit_limit=unit_limit, logger=logger)
    return unit_templates