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

def process_channel(args):
    i, k, j, channel = args[:4]
    channel_to_index,merging_template, merging_count_at_channel_by_sample = args[4:7]
    templates, merging_channel_loc_tuples = args[7:]

    # channel j
    #channel = channel_loc_tuples[j]
    found_idx = channel_to_index.get(channel)

    if found_idx is not None:
        # Found match, average values  
        merged_sample_at_channel = merging_template[k][found_idx]
        sample_to_be_merged = templates[i][k][j]

        # Choose to avg or replace based merge_count_at channel
        merge_count = merging_count_at_channel_by_sample[k][found_idx]

        if merge_count > 0:
            # Weighted Average - probably need to improve this method later
            total_samples_at_channel = merge_count + 1
            merged_sample_at_channel = (merged_sample_at_channel*merge_count + sample_to_be_merged) / (total_samples_at_channel)
        elif merge_count <= 0 and merged_sample_at_channel == 0:
            merged_sample_at_channel = sample_to_be_merged
        else:
            print("Something funky is going on")

        merging_template[k][found_idx] = merged_sample_at_channel
        merging_count_at_channel_by_sample[k][found_idx] += 1
    else:
        # No match, append
        sample_to_append = templates[i][k][j]
        #channel_loc_to_append = channel_loc[j] 

        # Extend merged_templates
        merging_template[k] = list(merging_template[k]) + [sample_to_append]

        # Extend merged_count_at_channel_by_sample, set to 1 representing 1 value at channel
        merging_count_at_channel_by_sample[k] = list(merging_count_at_channel_by_sample[k]) + [1]

        # Extend merged_channel_loc
        merging_channel_loc_tuples.append(channel)

    return merging_template, merging_channel_loc_tuples

def merge_templates(templates, channel_locations, plot_dir):
    ###Merge
    merged_template = []
    merged_channel_loc = []
    merging_template = None
    merging_count_at_channel_by_sample = None
    merging_channel_loc_tuples = None
    
    for i in range(len(templates)):
        #template key
        #rec_num = rec_nums[key]
        #are_keys_sequential(channel_locations)
        channel_loc = channel_locations[i]
        
        if len(merged_template) == 0:
            # First template, so just copy over
            merged_template.append(templates[i]) 
            merged_channel_loc.append(channel_loc)

            #prepare array to count how many times samples at a channel are merged
            num_samples = merged_template[-1][0].shape[0]
            # Calculate new size, Create pre-sized merged_template array
            merged_count_at_channel_by_sample = [np.array([np.zeros(num_samples, dtype='float32')] * len(merged_template[-1]), dtype='float32')]
            new_merged_count_at_channel_by_sample = [np.array([np.zeros(num_samples, dtype='float32')] * len(merged_template[-1]), dtype='float32')]
            #Add 1 indicating that first template represents 1 legit value for ready to be merged 
            new_merged_count_at_channel_by_sample[-1] += 1
            merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample 
        else:
            # Not first template, need to merge
            new_template = []
            new_channel_loc = []
            
            # #update array to count how many times a channel is merged
            # num_samples = merged_templates[-1][0].shape[0]
            # # Calculate new size
            # new_merged_count_at_channel_by_sample = [np.array([np.zeros(num_samples, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
            # for sample_idx, merge_counts_in_sample in enumerate(merged_count_at_channel_by_sample[-1]):
            #     for channel_idx, count_num in enumerate (merged_count_at_channel_by_sample[-1][sample_idx]):
            #         new_merged_count_at_channel_by_sample[-1][sample_idx][channel_idx] = count_num          
            # merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample  
            
            # Convert numpy arrays to tuples
            channel_loc_tuples = [tuple(channel) for channel in channel_loc]
            # make temporary list of lists template for merging
            if merging_template is None:
                merging_channel_loc_tuples = [tuple(channel) for channel in merged_channel_loc[-1]]
                merging_template = [arr.tolist() for arr in merged_template[-1]]
                merging_count_at_channel_by_sample = [arr.tolist() for arr in merged_count_at_channel_by_sample[-1]]            

            # Create a dictionary mapping channels to their indices
            channel_to_index = {channel: idx for idx, channel in enumerate(merging_channel_loc_tuples)}
            
            #merging_channel_loc = merged_channel_loc[-1]

            # for k in range(len(templates[i])):
            #     # sample k
            #     for j, channel in enumerate(channel_loc_tuples):
            #         # channel j
            #         found_idx = channel_to_index.get(channel)

            #         if found_idx is not None:
            #             # Found match, average values  
            #             merged_sample_at_channel = merging_template[k][found_idx]
            #             sample_to_be_merged = templates[i][k][j]

            #             # Choose to avg or replace based merge_count_at channel
            #             merge_count = merging_count_at_channel_by_sample[k][found_idx]

            #             if merge_count > 0:
            #                 # Weighted Average - probably need to improve this method later
            #                 total_samples_at_channel = merge_count + 1
            #                 merged_sample_at_channel = (merged_sample_at_channel*merge_count + sample_to_be_merged) / (total_samples_at_channel)
            #             elif merge_count <= 0 and merged_sample_at_channel == 0:
            #                 merged_sample_at_channel = sample_to_be_merged
            #             else:
            #                 print("Something funky is going on")

            #             merging_template[k][found_idx] = merged_sample_at_channel
            #             merging_count_at_channel_by_sample[k][found_idx] += 1
            #         else:
            #             # No match, append
            #             sample_to_append = templates[i][k][j]
            #             channel_loc_to_append = channel_loc[j] 

            #             # Extend merged_templates
            #             merging_template[k] = list(merging_template[k]) + [sample_to_append]

            #             # Extend merged_count_at_channel_by_sample, set to 1 representing 1 value at channel
            #             merging_count_at_channel_by_sample[k] = list(merging_count_at_channel_by_sample[k]) + [1]

            #             # Extend merged_channel_loc
            #             merging_channel_loc_tuples.append(channel)
                
            args = (channel_to_index, 
                    merging_template, merging_count_at_channel_by_sample, 
                    templates, merging_channel_loc_tuples)
            

            with Manager() as manager:
                merging_template = manager.list(merging_template)
                merging_channel_loc_tuples = manager.list(merging_channel_loc_tuples)

                with concurrent.futures.ProcessPoolExecutor() as executor:
                    for k in tqdm(range(len(templates[i])), desc='Processing templates'):
                        results = executor.map(process_channel, 
                                               [(i, k, j, channel) + args 
                                                #for k in range(len(templates[i])) 
                                                for j, channel in enumerate(channel_loc_tuples)])                        
    
    #Validation
    assert len(merged_template[-1]) == len(merging_template)
    for i in range(len(merging_template)):
        assert len(merging_template[i]) <= 26400
    assert len(merging_channel_loc_tuples) <= 26400          
    
    # Convert channel_loc to float64
    for i in range(merging_channel_loc_tuples):
        merging_channel_loc_tuples[i] = np.array(merging_channel_loc_tuples[i], dtype='float64')

    # Convert template to float32
    for i in range(len(merging_template)):
        merging_template[i] = np.array(merging_template[i], dtype='float32')

            
            # for k in range(len(templates[i])):
            #     #sample k
            #     for j in range(len(channel_loc)):
            #         #channel j
            #         channel = tuple(channel_loc[j])
            #         found_idx = None
            #         for idx, existing_channel in enumerate(merged_channel_loc[-1]):
            #             if tuple(existing_channel) == channel:
            #                 found_idx = idx
            #                 break

            #         if found_idx is not None:
            #             # Found match, average values  
            #             try: 
            #                 merged_sample_at_channel = merged_templates[-1][k][found_idx]
            #                 sample_to_be_merged = templates[i][k][j]
            #                 #Choose to avg or replace based merge_count_at channel
            #                 if merged_count_at_channel_by_sample[-1][k][found_idx] > 0:
            #                     #Weighted Average - probably need to improve this method later
            #                     merge_count = merged_count_at_channel_by_sample[-1][k][found_idx]
            #                     total_samples_at_channel = merge_count + 1
            #                     merged_sample_at_channel = (merged_sample_at_channel*merge_count + sample_to_be_merged) / (total_samples_at_channel)
            #                 elif merged_count_at_channel_by_sample[-1][k][found_idx] <= 0 and merged_sample_at_channel == 0:
            #                     merged_sample_at_channel = sample_to_be_merged
            #                 else:
            #                     raise Exception("Something funky is going on")
            #                 merged_templates[-1][k][found_idx] = merged_sample_at_channel
            #                 merged_count_at_channel_by_sample[-1][k][found_idx] += 1
            #                 #new_template.append((merged_sample_at_channel + sample_to_be_merged) / 2)
            #             except Exception as e: 
            #                 print(e)
            #         else:
            #             ### No match, append
            #             sample_to_append = templates[i][k][j]
            #             channel_loc_to_append = tuple(channel_loc[j]) 
                                        
            #             ## Copy existing templates over
            #             # Calculate size of new sub-arrays
            #             new_size = merged_templates[-1][0].shape[0] + 1 #get the size of any sample in merged templates, add 1
            #             # Calculate new size, Create pre-sized merged_template array
            #             new_merged_templates = [np.array([np.zeros(new_size, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
            #             #Allocated updated data to new_merged_templates
            #             for l, t in enumerate(merged_templates[-1]):
            #                 # Allocate new array
            #                 new_array = np.zeros(new_size, dtype=t.dtype)  
            #                 # Copy existing data
            #                 new_array[:t.shape[0]] = t
            #                 # Save to list
            #                 new_merged_templates[-1][l] = new_array

            #             #Update new value at sample at channel
            #             # Save back to merged templates  
            #             merged_templates = new_merged_templates
            #             # Now we can extend merged_templates[0]
            #             merged_templates[-1][k][-1] = sample_to_append

            #             #update array to count how many times a channel is merged
            #             num_samples = merged_templates[-1][0].shape[0]
            #             new_merged_count_at_channel_by_sample = [np.array([np.zeros(num_samples, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
            #             for sample_idx, merge_counts_in_sample in enumerate(merged_count_at_channel_by_sample[-1]):
            #                 for channel_idx, count_num in enumerate (merged_count_at_channel_by_sample[-1][sample_idx]):
            #                     new_merged_count_at_channel_by_sample[-1][sample_idx][channel_idx] = count_num          
            #             merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample          
            #             merged_count_at_channel_by_sample[-1][k][-1] += 1

            #             ## Copy existing Channels over
            #             # Calculate size of new sub-arrays and create space for new channel location
            #             same_size = merged_channel_loc[-1][0].shape[0] #get the size of any channel location tuple
            #             new_length = len(merged_channel_loc[-1]) + 1 #add space for new channel location
            #             new_merged_channels = [np.array([np.zeros(same_size, dtype='float64')] * new_length, dtype='float64')]
            #             for l, c in enumerate(merged_channel_loc[-1]):
            #                 # Allocate new array
            #                 new_tuple = np.zeros(same_size, dtype=c.dtype)  
            #                 # Copy existing data
            #                 new_tuple = c
            #                 # Save to list
            #                 new_merged_channels[-1][l] = new_tuple
            #             #Update new value at sample at channel
            #             tuple_to_append = np.array(channel_loc_to_append)   
            #             new_merged_channels[-1][-1] = tuple_to_append          
            #             merged_channel_loc = new_merged_channels

    # with open(os.path.join(plot_dir, 'merged_template_save.pkl'), 'wb') as f:
    #     pickle.dump(merged_templates, f)
    # with open(os.path.join(plot_dir, 'merged_channel_loc_save.pkl'), 'wb') as f:
    #     pickle.dump(merged_channel_loc, f)

    return merged_templates, merged_channel_loc