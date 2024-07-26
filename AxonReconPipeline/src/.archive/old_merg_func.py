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

def merge_templates(templates, channel_locations, plot_dir):
    ###Merge
    merged_templates = []
    merged_channel_loc = []
    verbose = False

    for i in range(len(templates)):
        print(f'Processing template {i+1}/{len(templates)}')
        #template key
        #rec_num = rec_nums[key]
        #are_keys_sequential(channel_locations)
        channel_loc = channel_locations[i]

        if len(merged_templates) == 0:
            # First template, so just copy over
            merged_templates.append(templates[i]) 
            merged_channel_loc.append(channel_loc)

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
            # Not first template, need to merge
            new_template = []
            new_channel_loc = []

            #update array to count how many times a channel is merged
            num_samples = merged_templates[-1][0].shape[0]
            # Calculate new size
            new_merged_count_at_channel_by_sample = [np.array([np.zeros(num_samples, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
            for sample_idx, merge_counts_in_sample in enumerate(merged_count_at_channel_by_sample[-1]):
                for channel_idx, count_num in enumerate (merged_count_at_channel_by_sample[-1][sample_idx]):
                    new_merged_count_at_channel_by_sample[-1][sample_idx][channel_idx] = count_num          
            merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample  

            for k in range(len(templates[i])):
                if verbose: print(f'Processing sample {k+1}/{len(templates[i])}')
                #sample k
                for j in range(len(channel_loc)):
                    if verbose: print(f'Processing channel {j+1}/{len(channel_loc)}')
                    #channel j
                    channel = tuple(channel_loc[j])
                    found_idx = None
                    for idx, existing_channel in enumerate(merged_channel_loc[-1]):
                        if tuple(existing_channel) == channel:
                            found_idx = idx
                            break

                    if found_idx is not None:
                        if verbose: print('found match, average')
                        # Found match, average values  
                        try: 
                            merged_sample_at_channel = merged_templates[-1][k][found_idx]
                            sample_to_be_merged = templates[i][k][j]
                            #Choose to avg or replace based merge_count_at channel
                            if merged_count_at_channel_by_sample[-1][k][found_idx] > 0:
                                #Weighted Average - probably need to improve this method later
                                merge_count = merged_count_at_channel_by_sample[-1][k][found_idx]
                                total_samples_at_channel = merge_count + 1
                                merged_sample_at_channel = (merged_sample_at_channel*merge_count + sample_to_be_merged) / (total_samples_at_channel)
                            elif merged_count_at_channel_by_sample[-1][k][found_idx] <= 0 and merged_sample_at_channel == 0:
                                merged_sample_at_channel = sample_to_be_merged
                            else:
                                raise Exception("Something funky is going on")
                            merged_templates[-1][k][found_idx] = merged_sample_at_channel
                            merged_count_at_channel_by_sample[-1][k][found_idx] += 1
                            #new_template.append((merged_sample_at_channel + sample_to_be_merged) / 2)
                        except Exception as e: 
                            print(e)
                    else:
                        if verbose: print('no match, append')
                        ### No match, append
                        sample_to_append = templates[i][k][j]
                        channel_loc_to_append = tuple(channel_loc[j]) 

                        ## Copy existing templates over
                        # Calculate size of new sub-arrays
                        new_size = merged_templates[-1][0].shape[0] + 1 #get the size of any sample in merged templates, add 1
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
                        merged_templates = new_merged_templates
                        # Now we can extend merged_templates[0]
                        merged_templates[-1][k][-1] = sample_to_append

                        #update array to count how many times a channel is merged
                        num_samples = merged_templates[-1][0].shape[0]
                        new_merged_count_at_channel_by_sample = [np.array([np.zeros(num_samples, dtype='float32')] * len(merged_templates[-1]), dtype='float32')]
                        for sample_idx, merge_counts_in_sample in enumerate(merged_count_at_channel_by_sample[-1]):
                            for channel_idx, count_num in enumerate (merged_count_at_channel_by_sample[-1][sample_idx]):
                                new_merged_count_at_channel_by_sample[-1][sample_idx][channel_idx] = count_num          
                        merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample          
                        merged_count_at_channel_by_sample[-1][k][-1] += 1

                        ## Copy existing Channels over
                        # Calculate size of new sub-arrays and create space for new channel location
                        same_size = merged_channel_loc[-1][0].shape[0] #get the size of any channel location tuple
                        new_length = len(merged_channel_loc[-1]) + 1 #add space for new channel location
                        new_merged_channels = [np.array([np.zeros(same_size, dtype='float64')] * new_length, dtype='float64')]
                        for l, c in enumerate(merged_channel_loc[-1]):
                            # Allocate new array
                            new_tuple = np.zeros(same_size, dtype=c.dtype)  
                            # Copy existing data
                            new_tuple = c
                            # Save to list
                            new_merged_channels[-1][l] = new_tuple
                        #Update new value at sample at channel
                        tuple_to_append = np.array(channel_loc_to_append)   
                        new_merged_channels[-1][-1] = tuple_to_append          
                        merged_channel_loc = new_merged_channels

    print('done merging')
    # with open(os.path.join(plot_dir, 'merged_template_save.pkl'), 'wb') as f:
    #     pickle.dump(merged_templates, f)
    # with open(os.path.join(plot_dir, 'merged_channel_loc_save.pkl'), 'wb') as f:
    #     pickle.dump(merged_channel_loc, f)

    return merged_templates, merged_channel_loc