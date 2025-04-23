### Initiate Axon trace algo
import numpy as np
import matplotlib.pylab as plt
import MEAutility as mu
from scipy.signal import resample_poly
from scipy.stats import kurtosis, linregress
from matplotlib import gridspec
from scipy import io
import numpy as np
import networkx as nx
from pathlib import Path
import sys
from pprint import pprint
from axon_velocity import *
from axon_velocity.models import load_cell
from axon_velocity.evaluation import *
import os
from pathlib import Path
import pickle
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from probeinterface import plotting as plotting_probe
import json

import concurrent.futures
from tqdm import tqdm
from multiprocessing import Manager

def debug():
    import h5py
    file_path = '/mnt/disk15tb/adam/git_workspace/random_variables_for_debug/templates.hdf5'
    templates_load = []
    with h5py.File(file_path, 'r') as h5f:
        for key in h5f.keys():
            templates_load.append(h5f[key][:])

    file_path = '/mnt/disk15tb/adam/git_workspace/random_variables_for_debug/channel_locations.hdf5'
    channel_load = []
    with h5py.File(file_path, 'r') as h5f:
        for key in h5f.keys():
            channel_load.append(h5f[key][:])

    templates = templates_load
    channel_locations = channel_load
    plot_dir = '/mnt/disk15tb/adam/git_workspace/random_variables_for_debug'

    return templates, channel_locations, plot_dir


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

def merge_templates(templates, channel_locations, plot_dir):
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

if __name__ == '__main__':
    templates, channel_locations, plot_dir = debug()
    merged_templates, merged_channel_loc = merge_templates(templates, channel_locations, plot_dir)
    #merged_templates, merged_channel_loc = fill_template(merged_templates, merged_channel_loc)
    # print(merged_templates)
    # print(merged_channel_loc)