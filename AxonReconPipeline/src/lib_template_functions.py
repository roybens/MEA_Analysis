''' This file contains template-related functions that are used in the main pipeline functions. '''

'''imports'''
import sys
import os
import h5py
import spikeinterface.full as si
#import spikeinterface.sorters as ss
#import shutil
import numpy as np
from tqdm import tqdm
#from pathlib import Path
#import time

#from axon_tracking import spike_sorting as ss
#import spike_sorting as ss

''' Local imports '''
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from MEAProcessingLibrary import mea_processing_library as MPL
import lib_sorting_functions as sorter
import lib_waveform_functions as waveformer
#import reconstruct_axons as ra

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
default_template_dir='./AxonReconPipeline/data/temp_data/templates'  

default_n_jobs = 4

'''functions'''
def extract_template_segments_for_unit(h5_file_path, stream_id, segment_sorting, sel_unit_id, save_root, peak_cutout=2, align_cutout=True, upsample=2, rm_outliers=True, n_jobs=16, n_neighbors=10, overwrite_wf=False, overwrite_tmp = True, just_load_waveforms = False):
    #full_path = segment_sorting._kwargs['recording_list'][0]._parent_recording._kwargs['recording']._kwargs['file_path']
    #full_path = segment_sorting._kwargs['recording_or_recording_list'][0]._kwargs['parent_recording']._kwargs['recording']._kwargs['file_path']
    full_path = h5_file_path
    # cutout = get_assay_information(full_path)
    # if align_cutout:
    #     wf_length = np.int16((sum(cutout) - 2*peak_cutout) * upsample) #length of waveforms after adjusting for potential peak alignments
    #     template_matrix = np.full([wf_length, 26400], np.nan)
    # else:
    #     wf_length = np.int16(sum(cutout) * upsample)
    #     template_matrix = np.full([wf_length, 26400], np.nan)
        

        
    
    # if not just_load_waveforms:
    #     extract_waveforms(h5_file_path, segment_sorting, stream_id, save_root, n_jobs, overwrite_wf)

    h5 = h5py.File(full_path)
    rec_names = list(h5['wells'][stream_id].keys())
    #Find a way to average common electrodes
    # Initialize count_matrix with the same shape as template_matrix, filled with zeros
    #count_matrix = np.zeros(template_matrix.shape)
    template_list = []
    channel_locations_list = []
    #replace templates for waveforms in save_root
    if 'templates' in save_root:
        wf_load_dir = save_root.replace('templates', 'waveforms')
        
    # aw = True
    # if aw:
    for sel_idx, rec_name in enumerate(rec_names): 
        try: rec = si.MaxwellRecordingExtractor(full_path,stream_id=stream_id,rec_name=rec_name)
        except Exception as e: 
            print(f'Error loading recording {rec_name} in segment {sel_idx}:\n{e}')
            continue
        #els = rec.get_property("contact_vector")["electrode"]
        seg_sort = si.SelectSegmentSorting(segment_sorting, sel_idx)
        seg_we = si.load_waveforms(os.path.join(wf_load_dir, 'waveforms', 'seg' + str(sel_idx)), sorting = seg_sort)
        template = seg_we.get_template(sel_unit_id)
        if np.isnan(template).all():
            logger.info(f'Unit {sel_unit_id} in segment {sel_idx} is empty')
            continue
        channel_locations = seg_we.get_channel_locations()
        # Define the directory and filename
        dir_path = os.path.join(save_root, 'partial_templates')
        file_name = f'seg{sel_idx}_unit{sel_unit_id}.npy'
        channel_loc_file_name = f'seg{sel_idx}_unit{sel_unit_id}_channels.npy'
        # Create the directory if it doesn't exist
        os.makedirs(dir_path, exist_ok=True)
        # Combine the directory and filename to get the full path
        channel_loc_save_file = os.path.join(dir_path, channel_loc_file_name)
        template_save_file = os.path.join(dir_path, file_name)
        
        if not np.isnan(template).all() and np.isnan(template).any():
            logger.info(f'Unit {sel_unit_id} in segment {sel_idx} has NaN values')
            #template = np.nan_to_num(template)    
        np.save(template_save_file, template)
        np.save(channel_loc_save_file, channel_locations)
        logger.info(f'Partial template saved to {template_save_file}')
        channel_locations_list.append(channel_locations)
        template_list.append(template)

    return template_list, channel_locations_list
def fill_template(merged_templates, merged_channel_loc, merged_count_at_channel_by_sample):
    def get_pitch():
        # Between any two adjacent channels, there is some minimmum distance, which we can use to determine if a channel is missing
        # Get distances between channels and find the minimum distance
        distances = []
        for i in range(len(merged_channel_loc[-1])):
            for j in range(i+1, len(merged_channel_loc[-1])):
                distances.append(np.linalg.norm(np.array(merged_channel_loc[-1][i]) - np.array(merged_channel_loc[-1][j])))
                min_distance = np.min(distances)
                if min_distance == 17.5: #specific to maxwell MEAs, 17.5 um
                    return min_distance
        #assert min distance is a mult of 17.5um
        assert min_distance % 17.5 == 0, "Minimum distance is not a multiple of 17.5 um."
        return min_distance
    
    ## Fill in missing channels
    print('Filling in missing channels...')
    # Find min and max x and y values in channel locations
    x_min = np.min([loc[0] for loc in merged_channel_loc[-1]])
    x_max = np.max([loc[0] for loc in merged_channel_loc[-1]])
    y_min = np.min([loc[1] for loc in merged_channel_loc[-1]])
    y_max = np.max([loc[1] for loc in merged_channel_loc[-1]])

    # Between any two adjacent channels, there is some minimmum distance, which we can use to determine if a channel is missing
    # Get distances between channels and find the minimum distance
    min_distance = get_pitch() #basically, expcted channel resolution

    # Create a grid of channels
    x_grid = np.arange(x_min, x_max+min_distance, min_distance)
    y_grid = np.arange(y_min, y_max+min_distance, min_distance)
    xx, yy = np.meshgrid(x_grid, y_grid)
    grid = np.array([xx.flatten(), yy.flatten()]).T

    #get channel locations present in grid, missing from channel locations
    existing_channels_set = set(map(tuple, merged_channel_loc[-1]))
    expected_channels_set = set(map(tuple, grid))
    missing_channel_set = [ch for idx, ch in enumerate(expected_channels_set) if tuple(ch) not in existing_channels_set]
    assert len(missing_channel_set) + len(merged_channel_loc[-1]) == len(grid), "Missing channels do not match."
    assert len(missing_channel_set) + len(merged_channel_loc[-1]) <= 26400, "Too many channels for Maxwell MEAs."

    # Fill in missing channels               
    
    #print total merged channels so far                
    #print(f'Total merged channels: {len(merged_channel_loc[-1])}')
    print(f'Processing {len(missing_channel_set)} unrecorded channels')
                        ## Copy existing templates over
    # Calculate size of new sub-arrays
    #new_size = merged_templates[-1][0].shape[0] + 1 #get the size of any sample in merged templates, add 1
    new_size = merged_templates[-1][0].shape[0] + len(missing_channel_set)
    # Calculate new size, Create pre-sized merged_template array
    new_merged_templates = [np.array([np.zeros(new_size, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
    #Allocated updated data to new_merged_templates
    for l, t in enumerate(merged_templates[-1]):
        # Allocate new array
        new_array = np.zeros(new_size, dtype=t.dtype)  
        # Copy existing data
        new_array[:t.shape[0]] = t
        # Save to list
        new_merged_templates[-1][l] = new_array
    #Update new value at sample at channel
    # Save back to merged templates  
    #merged_templates = new_merged_templates

    #update array to count how many times a channel is merged
    num_samples = new_merged_templates[-1][0].shape[0]
    new_merged_count_at_channel_by_sample = [np.array([np.zeros(num_samples, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
    for sample_idx, merge_counts_in_sample in enumerate(merged_count_at_channel_by_sample[-1]):
        for channel_idx, count_num in enumerate (merged_count_at_channel_by_sample[-1][sample_idx]):
            new_merged_count_at_channel_by_sample[-1][sample_idx][channel_idx] = count_num          
    #merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample

    ## Copy existing Channels over
    # Calculate size of new sub-arrays and create space for new channel location
    same_size = merged_channel_loc[-1][0].shape[0] #get the size of any channel location tuple
    new_length = len(merged_channel_loc[-1]) + len(missing_channel_set)
    new_merged_channels = [np.array([np.zeros(same_size, dtype='float64')] * new_length, dtype='float64')]
    for l, c in enumerate(merged_channel_loc[-1]):
        # Allocate new array
        new_tuple = np.zeros(same_size, dtype=c.dtype)  
        # Copy existing data
        new_tuple = c
        # Save to list
        new_merged_channels[-1][l] = new_tuple
    #Update new value at sample at channel
    #merged_channel_loc = new_merged_channel   

    ###Update Values
    for idx, j in enumerate(missing_channel_set):
        position = merged_templates[-1][0].shape[0] + idx
        #new_merged_count_at_channel_by_sample[-1][k][position] = 0
        channel_loc_to_append = tuple(missing_channel_set[idx])   
        tuple_to_append = np.array(channel_loc_to_append)  
        new_merged_channels[-1][position] = tuple_to_append          

    #validate
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
def merge_templates(templates, channel_locations):
    ###Merge
    merged_templates = []
    merged_channel_loc = []
    verbose = True

    print('Starting merge...')
    test_channel_measure = 0
    for i in range(len(templates)):
        #print(f'Processing template {i+1}/{len(templates)}')
        #template key
        #rec_num = rec_nums[key]
        #are_keys_sequential(channel_locations)
        channel_loc = channel_locations[i]

        if len(merged_templates) == 0:
            # First template, so just copy over
            incoming_template = templates[i]
            #continue until first template contains no NaN values
            #missing channels will be handled by fill_template
            if np.isnan(incoming_template).all(): continue
            assert np.isnan(incoming_template).any() == False, "First template contains NaN values."

            merged_templates.append(templates[i])
            merged_channel_loc.append(channel_loc)
            test_channel_measure += len(channel_loc)            
            #prepare array to count how many times samples at a channel are merged
            num_samples = merged_templates[-1][0].shape[0]
            # Calculate new size, Create pre-sized merged_template array
            merged_count_at_channel_by_sample = [np.array([np.zeros(num_samples, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
            new_merged_count_at_channel_by_sample = [np.array([np.zeros(num_samples, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
            #Add 1 indicating that first template represents 1 legit value for ready to be merged 
            for sample_idx, merge_counts_in_sample in enumerate(merged_count_at_channel_by_sample[-1]):
                for channel_idx, count_num in enumerate (merged_count_at_channel_by_sample[-1][sample_idx]):
                    new_merged_count_at_channel_by_sample[-1][sample_idx][channel_idx] = count_num+1
            merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample 
        else:
            if verbose: print(f'Merging matching channels...')
            print(f'Total merged channels: {len(merged_channel_loc[-1])}')            
            if verbose: print(f'Total channels: {test_channel_measure}')  
            # Not first template, need to merge
            new_template = []
            new_channel_loc = []
            incoming_template = templates[i]
            if verbose: print(f'Processing {len(channel_loc)} channels')
            
            

            # Find matching and non-matching indices between channel_loc and merged_channel_loc
            # This assumes that channel locations are unique and can be identified with a set operation
            existing_channels_set = set(map(tuple, merged_channel_loc[-1]))
            incoming_channels_set = set(map(tuple, channel_loc))
            overlapping_channels = existing_channels_set.intersection(incoming_channels_set)
            
            idx_incoming_channels_in_overlap = [idx for idx, ch in enumerate(channel_loc) if tuple(ch) in overlapping_channels]
            idx_existing_channels_in_overlap = [idx for idx, ch in enumerate(merged_channel_loc[-1]) if tuple(ch) in overlapping_channels]
            idx_new_incoming_channels = [idx for idx, ch in enumerate(channel_loc) if tuple(ch) not in overlapping_channels]
            assert len(idx_incoming_channels_in_overlap) + len(idx_new_incoming_channels) == len(channel_loc), "Indices do not match."
            test_channel_measure += len(idx_new_incoming_channels)

            #validate
            channels_already_merged = merged_channel_loc[-1][idx_existing_channels_in_overlap]
            merging_channels = channel_loc[idx_incoming_channels_in_overlap]
            #these should be the same, but merging channels may be out of order. Bring merging channels into the same order as channels_already_merged
            merging_channels_ordered = [ch for idx, ch in enumerate(channels_already_merged) if tuple(ch) in merging_channels]
            merging_channels_ordered = np.array(merging_channels_ordered, dtype='float64')
            assert len(channels_already_merged) == len(merging_channels_ordered), "The number of channels does not match."
            # Now, we should compare the actual channel locations to ensure they are identical
            # Since we already have matched indices, we should use them to compare the channels directly
            for idx in range(len(channels_already_merged)):
                assert np.array_equal(channels_already_merged[idx], merging_channels_ordered[idx]), "Matching channels are not the same."
            shared_ordered_channels = merging_channels_ordered                   
           
            if verbose: print(f'Processing {len(channels_already_merged)} matching channels')
            if idx_existing_channels_in_overlap or idx_incoming_channels_in_overlap:               
                print(f'Processing matching channels across footprints, template {i+1}/{len(templates)}...')            
                
                #k = matching_indices
                
                merged_sample_at_channel = merged_templates[-1][:, idx_existing_channels_in_overlap]
                #sample_to_be_merged = incoming_template[:, k]
                sample_to_be_merged = incoming_template[:, idx_incoming_channels_in_overlap]

                # Compute weighted average
                merge_count = merged_count_at_channel_by_sample[-1][:, idx_existing_channels_in_overlap]             
                total_samples_at_channel = merge_count+1
                #uncount nan values
                total_samples_at_channel[np.isnan(sample_to_be_merged)] = total_samples_at_channel[np.isnan(sample_to_be_merged)] - 1
                sample_to_be_merged[np.isnan(sample_to_be_merged)] = 0  #turn any nans into zeros                               
                # Calculate the weighted average
                weighted_average = np.zeros_like(total_samples_at_channel)
                non_zero_mask = total_samples_at_channel != 0
                weighted_average[non_zero_mask] = (merged_sample_at_channel[non_zero_mask] * merge_count[non_zero_mask] + sample_to_be_merged[non_zero_mask]) / total_samples_at_channel[non_zero_mask]
                assert not np.isnan(weighted_average).any(), "Weighted average contains NaN values."
                
                # Update merged templates
                merged_templates[-1][:, idx_existing_channels_in_overlap] = weighted_average
                merged_count_at_channel_by_sample[-1][:, idx_existing_channels_in_overlap] = total_samples_at_channel
                        
            if idx_new_incoming_channels:
                #if verbose: print(f'Non-matching indices: {non_matching_indices}')
                #update array to count how many times a channel is merged
                # num_samples = merged_templates[-1][0].shape[0]
                # # Calculate new size
                # new_merged_count_at_channel_by_sample = [np.array([np.zeros(num_samples, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
                # for sample_idx, merge_counts_in_sample in enumerate(merged_count_at_channel_by_sample[-1]):
                #     for channel_idx, count_num in enumerate (merged_count_at_channel_by_sample[-1][sample_idx]):
                #         new_merged_count_at_channel_by_sample[-1][sample_idx][channel_idx] = count_num          
                # #merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample                
                
                #print total merged channels so far                
                #print(f'Total merged channels: {len(merged_channel_loc[-1])}')
                print(f'Processing {len(idx_new_incoming_channels)} non-matching channels')
                                    ## Copy existing templates over
                # Calculate size of new sub-arrays
                #new_size = merged_templates[-1][0].shape[0] + 1 #get the size of any sample in merged templates, add 1
                new_size = merged_templates[-1][0].shape[0] + len(idx_new_incoming_channels)
                # Calculate new size, Create pre-sized merged_template array
                new_merged_templates = [np.array([np.zeros(new_size, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
                #Allocated updated data to new_merged_templates
                for l, t in enumerate(merged_templates[-1]):
                    # Allocate new array
                    new_array = np.zeros(new_size, dtype=t.dtype)  
                    # Copy existing data
                    new_array[:t.shape[0]] = t
                    # Save to list
                    new_merged_templates[-1][l] = new_array
                #Update new value at sample at channel
                # Save back to merged templates  
                #merged_templates = new_merged_templates

                #update array to count how many times a channel is merged
                num_samples = new_merged_templates[-1][0].shape[0]
                new_merged_count_at_channel_by_sample = [np.array([np.zeros(num_samples, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
                for sample_idx, merge_counts_in_sample in enumerate(merged_count_at_channel_by_sample[-1]):
                    for channel_idx, count_num in enumerate (merged_count_at_channel_by_sample[-1][sample_idx]):
                        new_merged_count_at_channel_by_sample[-1][sample_idx][channel_idx] = count_num          
                #merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample

                ## Copy existing Channels over
                # Calculate size of new sub-arrays and create space for new channel location
                same_size = merged_channel_loc[-1][0].shape[0] #get the size of any channel location tuple
                new_length = len(merged_channel_loc[-1]) + len(idx_new_incoming_channels)
                new_merged_channels = [np.array([np.zeros(same_size, dtype='float64')] * new_length, dtype='float64')]
                for l, c in enumerate(merged_channel_loc[-1]):
                    # Allocate new array
                    new_tuple = np.zeros(same_size, dtype=c.dtype)  
                    # Copy existing data
                    new_tuple = c
                    # Save to list
                    new_merged_channels[-1][l] = new_tuple
                #Update new value at sample at channel
                #merged_channel_loc = new_merged_channels  
                for k in tqdm(range(len(templates[i])), desc=f'Processing non-matching channels in each footprint, template {i+1}/{len(templates)}'):
                #for k in range(len(templates[i])):
                    #if verbose: print(f'Processing sample {k+1}/{len(templates[i])}')
                    #footprint k
                    #channel_to_index = {tuple(existing_channel): idx for idx, existing_channel in enumerate(merged_channel_loc[-1])}        
                
                    ###Update Values
                    for idx, j in enumerate(idx_new_incoming_channels):
                        position = merged_templates[-1][0].shape[0] + idx
                        ### No match, append
                        sample_to_append = templates[i][k][j]
                        if np.isnan(sample_to_append):
                            sample_to_append = 0
                            new_merged_count_at_channel_by_sample[-1][k][position] = 0
                        else:
                            new_merged_count_at_channel_by_sample[-1][k][position] += 1
                        # Now we can extend merged_templates[0]
                        new_merged_templates[-1][k][position] = sample_to_append               
                        channel_loc_to_append = tuple(channel_loc[j])   
                        tuple_to_append = np.array(channel_loc_to_append)  
                        new_merged_channels[-1][position] = tuple_to_append          

                merged_templates = new_merged_templates
                merged_channel_loc = new_merged_channels
                merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample                   
      
            else:
                raise Exception("Something funky is going on")
    print('done merging')
    # with open(os.path.join(plot_dir, 'merged_template_save.pkl'), 'wb') as f:
    #     pickle.dump(merged_templates, f)
    # with open(os.path.join(plot_dir, 'merged_channel_loc_save.pkl'), 'wb') as f:
    #     pickle.dump(merged_channel_loc, f)

    print(f'Total merged channels: {len(merged_channel_loc[-1])}')
    for i in range(len(merged_templates[-1])):
        assert len(merged_templates[-1][i]) <= 26400
    assert len(merged_channel_loc[-1]) <= 26400

    merged_templates_filled, merged_channel_loc_filled, merged_count_at_channel_by_sample_filled = fill_template(merged_templates, merged_channel_loc, merged_count_at_channel_by_sample)

    print(f'Total merged channels after fill: {len(merged_channel_loc_filled[-1])}')
    for i in range(len(merged_templates_filled[-1])):
        assert len(merged_templates_filled[-1][i]) <= 26400
    assert len(merged_channel_loc_filled[-1]) <= 26400
    
    return merged_templates, merged_channel_loc, merged_templates_filled, merged_channel_loc_filled  
def merge_template_segments(template_list, channel_locations_list):
    logger.info(f'Merging partial templates')
    merged_template, merged_channel_loc, merged_template_filled, merged_channel_loc_filled = merge_templates(template_list, channel_locations_list)
    return merged_template, merged_channel_loc, merged_template_filled, merged_channel_loc_filled    
def extract_templates_segments(h5_file_path, stream_id, segment_sorting, pos, te_params, save_root=None, just_load_waveforms = False, quick_load_templates = False):
    sel_unit_ids = segment_sorting.get_unit_ids()
    if save_root is None:
        #save_root = os.path.dirname(full_path) 
        recording_details = MPL.extract_recording_details(h5_file_path)
        date = recording_details[0]['date']
        chip_id = recording_details[0]['chipID']
        scanType = recording_details[0]['scanType']
        run_id = recording_details[0]['runID']
        save_root = f'./AxonReconPipeline/data/temp_data/templates/{date}/{chip_id}/{scanType}/{run_id}/{stream_id}'
    #template_save_path = os.path.join(save_root, 'waveforms')
    # if not os.path.exists(waveform_save_path):
    #     os.makedirs(waveform_save_path)
    template_save_path = os.path.join(save_root, 'templates')
    if not os.path.exists(template_save_path):
        os.makedirs(template_save_path)
        
    for sel_unit_id in tqdm(sel_unit_ids):

        template_save_file = os.path.join(template_save_path, str(sel_unit_id) + '.npy')
        channel_loc_save_file = os.path.join(template_save_path, str(sel_unit_id) + '_channels.npy')
        template_save_file_fill = os.path.join(template_save_path, str(sel_unit_id) + '_filled.npy')
        channel_loc_save_file_fill = os.path.join(template_save_path, str(sel_unit_id) + '_channels_filled.npy')

        '''Check if template files already exist for this unit and overwrite if necessary'''
        merge_template = True
        try: assert os.path.isfile(template_save_file) and os.path.isfile(channel_loc_save_file) and os.path.isfile(template_save_file_fill) and os.path.isfile(channel_loc_save_file_fill) == True, f'Template files already exist for unit {sel_unit_id}'
        except: merge_template = False #load_template = True

        '''Attempt to load merged template if it exists'''
        if merge_template and not te_params['overwrite_tmp']: #try_load_merged_template()
            try: 
                logger.info(f'Loading merged template for unit {sel_unit_id}')
                merged_template = np.load(template_save_file)
                merged_channel_loc = np.load(channel_loc_save_file)
                merged_template_filled = np.load(template_save_file_fill)
                merged_channel_locs_filled = np.load(channel_loc_save_file_fill)
                continue #if load is successful, skip to next unit
            except Exception as e: 
                logger.info(f'Error loading template for unit {sel_unit_id}:\n{e}')
                if quick_load_templates: raise e #if quick_load_templates is True, raise error and force normal template extraction
                pass #if load is unsuccessful, continue to extract and merge templates

        '''Extract and merge templates for this unit'''
        #if te_params['overwrite_tmp']:
        try:
            logger.info(f'Extracting template segments for unit {sel_unit_id}')
            template_list, channel_locations_list = extract_template_segments_for_unit(h5_file_path, stream_id, segment_sorting, sel_unit_id, save_root, just_load_waveforms = just_load_waveforms, **te_params)
            logger.info(f'Merging partial templates for unit {sel_unit_id}')
            merged_template, merged_channel_loc, merged_template_filled, merged_channel_locs_filled = merge_template_segments(template_list, channel_locations_list)
            #grid = convert_to_grid(template_matrix, pos)
            np.save(channel_loc_save_file, merged_channel_loc)
            np.save(template_save_file, merged_template)
            np.save(channel_loc_save_file_fill, merged_channel_locs_filled)
            np.save(template_save_file_fill, merged_template_filled)
            logger.info(f'Merged template saved to {save_root}/templates')
        except Exception as e:
            print(f'Unit {sel_unit_id} encountered the following error:\n {e}')
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
def extract_templates_from_sortings(sorting_list, h5_file_path, qc_params={}, te_params={}, stream_select = None, just_load_waveforms = False, quick_load_templates = False):
    #rec_list = list(sorting_dict.keys())
    if isinstance(h5_file_path, str): rec_list = [h5_file_path]    

    for rec_path in rec_list:
        for sorting in sorting_list:
            sorting_path = sorting._kwargs['folder_path']
            stream_id = [p for p in sorting_path.split('/') if p.startswith('well')][0] #Find out which well this belongs to
            if stream_select is not None:
                #stream_select = 3  # for example
                stream_select_str = f'well{stream_select:03}'
                #stream_ids = stream_ids[stream_select]
                if stream_select_str not in stream_id:
                    continue
            
            #rec_names, common_el, pos = ss.find_common_electrodes(rec_path, stream_id)
            multirecording, pos = sorter.concatenate_recording_slices(rec_path, stream_id)

            #duration = int(h5['assay']['inputs']['record_time'][0].decode('UTF-8')) * n_recs #In case we want to use firing rate as criterion
            cleaned_sorting = select_units(sorting, **qc_params)
            cleaned_sorting = si.remove_excess_spikes(cleaned_sorting, multirecording) #Relevant if last spike time == recording_length
            cleaned_sorting.register_recording(multirecording)
            segment_sorting = si.SplitSegmentSorting(cleaned_sorting, multirecording)
            extract_templates_segments(h5_file_path, stream_id, segment_sorting, pos, te_params, save_root=None, just_load_waveforms = just_load_waveforms, quick_load_templates = quick_load_templates)   
def extract_all_templates(continuous_h5_info, sorting_lists_by_scan = None, qc_params={}, te_params={}, allowed_scan_types = ["AxonTracking"], templates_dir='./AxonReconPipeline/data/temp_data/templates', stream_select = None, just_load_waveforms = False, quick_load_templates = False):
    
    if sorting_lists_by_scan is None: sorting_lists_by_scan = sorter.spikesort_recordings(continuous_h5_info, sortings_dir = default_sorting_dir, allowed_scan_types = allowed_scan_types, stream_select = stream_select,)
    if quick_load_templates is False: waveformer.extract_all_waveforms(continuous_h5_info, sorting_lists_by_scan = sorting_lists_by_scan, qc_params=qc_params, te_params=te_params, stream_break = stream_select)
    else: pass #skip trying waveform extraction/loading if quick_load_templates is True
    for dir in continuous_h5_info.values():
        h5_file_path = dir['h5_file_path']
        #sorting_list_by_stream = sorting_lists_by_scan[0] #for debugging
        scantype = dir['scanType']
        chipID = dir['chipID']
        runID = dir['runID']
        if scantype in allowed_scan_types:
            for sorting_list in sorting_lists_by_scan:
                sorting_path = sorting_list[0]._kwargs['folder_path']
                try: assert chipID and runID in sorting_path, f'ChipID and RunID not found in sorting path {sorting_path}'
                except: continue
                extract_templates_from_sortings(sorting_list, h5_file_path, qc_params, te_params, stream_select = stream_select, just_load_waveforms = just_load_waveforms, quick_load_templates = quick_load_templates)