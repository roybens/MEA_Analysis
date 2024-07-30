from AxonReconPipeline.src.lib_template_functions import *
import concurrent.futures
import numpy as np
from tqdm import tqdm

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
      
            # else:
            #     raise Exception("Something funky is going on")
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

def merge_templates_single(template_data): #aw20July2024 - This doesnt currently work. TODO: Fix this later. not important.
    templates, channel_locations, logger, start, end, part_merged_count_at_channel_by_sample = template_data
    merged_templates = []
    merged_channel_loc = []
    verbose = True

    if logger is not None: logger.debug('Starting merge...')
    else: print('Starting merge...')
        
    test_channel_measure = 0
    for i in range(start, end):
        channel_loc = channel_locations[i]

        if len(merged_templates) == 0:
            incoming_template = templates[i]
            if np.isnan(incoming_template).all(): continue
            assert np.isnan(incoming_template).any() == False, "First template contains NaN values."

            merged_templates.append(templates[i])
            merged_channel_loc.append(channel_loc)
            test_channel_measure += len(channel_loc)
            num_channels = merged_templates[-1][0].shape[0]
            if part_merged_count_at_channel_by_sample is None: merged_count_at_channel_by_sample = [np.array([np.zeros(num_channels, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
            else: merged_count_at_channel_by_sample = part_merged_count_at_channel_by_sample
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
                        
            elif idx_new_incoming_channels:         
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
        
    assert len(merged_channel_loc[-1]) == len(set(map(tuple, merged_channel_loc[-1]))), "Channels are not unique."
    for i in range(len(merged_templates[-1])): 
        assert len(merged_templates[-1][i]) <= 26400, f"Template {i} has more than 26400 channels."
    assert len(merged_channel_loc[-1]) <= 26400, "Too many channels for Maxwell MEAs."


    # merged_templates_filled, merged_channel_loc_filled, merged_count_at_channel_by_sample_filled = fill_template(merged_templates, merged_channel_loc, merged_count_at_channel_by_sample, logger=logger)

    # if logger is not None: logger.debug(f'Total merged channels after fill: {len(merged_channel_loc_filled[-1])}')
    # else: print(f'Total merged channels after fill: {len(merged_channel_loc_filled[-1])}')
        
    # for i in range(len(merged_templates_filled[-1])): assert len(merged_templates_filled[-1][i]) <= 26400
    # assert len(merged_channel_loc_filled[-1]) <= 26400
    return merged_templates, merged_channel_loc, merged_count_at_channel_by_sample #, merged_templates_filled, merged_channel_loc_filled

def merge_templates_concurrent(templates, channel_locations, logger=None): #aw20July2024 - This doesnt currently work. TODO: Fix this later. not important.
    def add_up_counts(merged_count_at_channel_by_sample, merged_channel_loc):
        final_merged_count_at_channel = [np.array([np.zeros(merged_channel_loc[-1].shape[0], dtype='float32')] * len(merged_count_at_channel_by_sample[-1]), dtype='float32')]
        for sample_idx, merge_counts_in_sample in enumerate(merged_count_at_channel_by_sample[-1]):
            for channel_idx, count_num in enumerate(merged_count_at_channel_by_sample[-1][sample_idx]):
                final_merged_count_at_channel[-1][sample_idx][channel_idx] = count_num
        return final_merged_count_at_channel
    
    num_threads = 4
    n = len(templates)
    step = n // num_threads
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        for i in range(num_threads):
            start = i * step
            end = (i + 1) * step if i != num_threads - 1 else n
            template_data = (templates, channel_locations, logger, start, end, None)
            futures.append(executor.submit(merge_templates_single, template_data))
        results = [future.result() for future in futures]
    
    # Combine results from all threads
    partial_merged_templates = results[0][0]
    partial_merged_channel_locs = results[0][1]
    partial_merged_count_at_channel_by_sample = results[0][2]
    # final_merged_templates_filled = results[0][2]
    # final_merged_channel_loc_filled = results[0][3]
    
    for result in results[0:]:
        partial_merged_templates += result[0]
        partial_merged_channel_locs += result[1]
        partial_merged_count_at_channel_by_sample += result[2]
        # final_merged_templates_filled += result[2]
        # final_merged_channel_loc_filled += result[3]

    part_merged_count_at_channel_by_sample = add_up_counts(partial_merged_count_at_channel_by_sample, partial_merged_channel_locs)
    template_data = (partial_merged_templates, partial_merged_channel_locs, logger, 0, len(partial_merged_templates), part_merged_count_at_channel_by_sample)
    final_merged_templates, final_merged_channel_loc, final_merged_count_at_channel_by_sample = merge_templates_single(template_data)
    final_merged_templates_filled, final_merged_channel_loc_filled, final_merged_count_at_channel_by_sample_filled = fill_template(final_merged_templates, final_merged_channel_loc, final_merged_count_at_channel_by_sample, logger=logger)

    if logger is not None: logger.debug(f'Total merged channels after fill: {len(final_merged_channel_loc_filled[-1])}')
    else: print(f'Total merged channels after fill: {len(final_merged_channel_loc_filled[-1])}')
        
    for i in range(len(final_merged_templates_filled[-1])): assert len(final_merged_templates_filled[-1][i]) <= 26400
    assert len(final_merged_channel_loc_filled[-1]) <= 26400   
    
    return final_merged_templates, final_merged_channel_loc, final_merged_templates_filled, final_merged_channel_loc_filled