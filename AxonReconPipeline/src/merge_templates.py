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


###

def merge_templates(templates, channel_locations, plot_dir):
    ###Merge
    merged_template = []
    merged_channel_loc = []
    merging_template = None
    merging_count_at_channel_by_sample = None
    merging_channel_loc_tuples = None
    verbose = False
    
    for i in tqdm(range(len(templates)), desc='Merging templates'):

        channel_loc = channel_locations[i]
        # Convert numpy arrays to tuples
        channel_loc_tuples = [tuple(channel) for channel in channel_loc]
        
        #initialize merging_template if it is None
        if merging_template is None:
            merged_template.append(templates[i]) 
            merged_channel_loc.append(channel_loc)
            num_samples = merged_template[-1][0].shape[0]
            merged_count_at_channel_by_sample = [np.array([np.zeros(num_samples, dtype='float32')] * len(merged_template[-1]), dtype='float32')]
            new_merged_count_at_channel_by_sample = [np.array([np.zeros(num_samples, dtype='float32')] * len(merged_template[-1]), dtype='float32')]
            merged_count_at_channel_by_sample = new_merged_count_at_channel_by_sample
            merging_channel_loc_tuples = [tuple(channel) for channel in merged_channel_loc[-1]]
            merging_template = [[val for val in arr] for arr in merged_template[-1]]
            merging_count_at_channel_by_sample = [arr.tolist() for arr in merged_count_at_channel_by_sample[-1]]
            #for any non-nan values, set count to 1 in merging_count_at_channel_by_sample
            for k in range(len(merging_template)):
                for j in range(len(merging_template[k])):
                    if not np.isnan(merging_template[k][j]):
                        merging_count_at_channel_by_sample[k][j] = 1
            assert not np.isnan(merging_count_at_channel_by_sample).any()
            continue

        #speed things up, check if all samples at channels are Nan
        # if np.isnan(templates[i]).all():
        #     continue

        # if np.isnan(templates[i]).any():
        #     print("Nans in template")             

        # Create a dictionary mapping channels to their indices
        channel_to_index = {channel: idx for idx, channel in enumerate(merging_channel_loc_tuples)}

        print(f"Merging template: {i}...")
        for k in range(len(templates[i])):                
            # sample k
            #speed things up, check if all samples at channels are Nan
            # if np.isnan(templates[i][k]).all():
            #     continue  
            for j, channel in enumerate(channel_loc_tuples):
                # channel j
                found_idx = channel_to_index.get(channel)

                if found_idx is not None:
                    if verbose: print("Found match, averaging")
                    # Found match, average values  
                    try:
                        merged_sample_at_channel = merging_template[k][found_idx]
                        sample_to_be_merged = templates[i][k][j]
                        #avoid merging Nan Values:
                        # if np.isnan(sample_to_be_merged):
                        #     continue
                    except:
                        print("Error")

                    # Choose to avg or replace based merge_count_at channel
                    try: merge_count = merging_count_at_channel_by_sample[k][found_idx]
                    except:
                        print("Error")

                    if not np.isnan(sample_to_be_merged):
                        if merge_count > 0:
                            # Weighted Average - probably need to improve this method later
                            total_samples_at_channel = merge_count + 1
                            merged_sample_at_channel = (merged_sample_at_channel*merge_count + sample_to_be_merged) / (total_samples_at_channel)
                        elif merge_count <= 0 and np.isnan(merged_sample_at_channel):
                            total_samples_at_channel = 1
                            merged_sample_at_channel = sample_to_be_merged
                        # elif not np.isnan(merged_sample_at_channel):
                        #     #merge_count = 1
                        #     total_samples_at_channel = 1
                        #     merged_sample_at_channel = (merged_sample_at_channel*merge_count + sample_to_be_merged) / (total_samples_at_channel)
                        else:
                            print("Something funky is going on")
                            print(merged_sample_at_channel)
                            print(sample_to_be_merged)
                            print(merge_count)
                            print(i, k, j)
                            #debugging``
                            break
                    else:
                        continue

                    merging_template[k][found_idx] = merged_sample_at_channel
                    merging_count_at_channel_by_sample[k][found_idx] += 1
                else:
                    if verbose: print("No match, appending")
                    # No match, append
                    sample_to_append = templates[i][k][j]
                    # #avoid appending Nan Values:
                    # if np.isnan(sample_to_append):
                    #     continue
                    channel_loc_to_append = channel_loc[j] 

                    # Extend merged_templates
                    merging_template[k] = list(merging_template[k]) + [sample_to_append]

                    # Extend merged_count_at_channel_by_sample, set to 1 representing 1 value at channel
                    merging_count_at_channel_by_sample[k] = list(merging_count_at_channel_by_sample[k]) + [1]

                    # Extend merged_channel_loc
                    merging_channel_loc_tuples.append(channel)

                    #Keep lists consistent length
                    # Find the maximum length of any list in merging_template
                    max_length = max(len(lst) for lst in merging_template)

                    # Extend all lists in merging_template to the maximum length
                    for lst in merging_template:
                        lst.extend([np.nan] * (max_length - len(lst)))

                    # Do the same for merging_count_at_channel_by_sample
                    for lst in merging_count_at_channel_by_sample:
                        lst.extend([0] * (max_length - len(lst)))

                    # Update dictionary mapping channels to their indices
                    channel_to_index = {chan: idx for idx, chan in enumerate(merging_channel_loc_tuples)}                     
    
    #Validation
    assert len(merged_template[-1]) == len(merging_template)
    for i in range(len(merging_template)):
        assert len(merging_template[i]) <= 26400
    assert len(merging_channel_loc_tuples) <= 26400
    # Assert that all elements in merging_count_at_channel_by_sample are not zero
    #assert not (merging_count_at_channel_by_sample == 0).any()
    for i in range(len(merging_count_at_channel_by_sample)):
        assert not (np.array(merging_count_at_channel_by_sample[i]) == 0).any()
    assert not np.isnan(merging_template).any()     
    

    
    # Convert channel_loc to float64 and make it a numpy array of numpy arrays
    new_merged_channel_loc = np.array([np.array(tuple, dtype='float64') for tuple in merging_channel_loc_tuples])
    merged_channel_loc[-1] = new_merged_channel_loc   

    # Convert template to float32 and make it a numpy array of numpy arrays
    #new_merged_template = np.array([np.array(template, dtype='float32') for template in merging_template])
    merged_template[-1] = np.array([np.array(template, dtype='float32') for template in merging_template])

    return merged_template, merged_channel_loc


    # remove_nans = False
    # if remove_nans:
    #     #Remove nanas from merged template (make function later)
    #     #At this point there will be nan values in columns (representing channels that were not present in any of the partial templates), remove those columns and correspoinding channel locations
    #     #Remove nan columns from merging_template and track column indices to remove from merging_channel_loc_tuplesimport numpy as np
    #     # Convert merging_template to a numpy array
    #     merging_template_nans = np.array(merging_template)

    #     # Find the indices of the columns where all values are NaN
    #     nan_cols = np.where(np.all(np.isnan(merging_template_nans), axis=0))[0]

    #     #verify each colum is all nans
    #     for i in nan_cols:
    #         assert np.all(np.isnan(merging_template_nans[:,i]))

    #     # Remove the NaN columns from merging_template
    #     merging_template_nans = np.delete(merging_template_nans, nan_cols, axis=1)

    #     # Remove the corresponding elements from merging_channel_loc_tuples
    #     merging_channel_loc_tuples_nan = [loc for i, loc in enumerate(merging_channel_loc_tuples) if i not in nan_cols]
    #     merging_channel_loc_tuples = merging_channel_loc_tuples_nan

    #     #Validation (make function later)
    #     assert np.all(np.isnan(merging_template_nans) == False)
    #     assert len(merged_template[-1]) == len(merging_template)
    #     for i in range(len(merging_template)):
    #         assert len(merging_template[i]) <= 26400
    #     assert len(merging_channel_loc_tuples) <= 26400

    #     #Convert back to list of lists
    #     merging_template = merging_template_nans.tolist()  