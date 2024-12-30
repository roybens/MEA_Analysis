
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


#Define Algorithm Params
params = get_default_graph_velocity_params()
# change params (fig6 Buccino axon_velocity paper example)
params['detect_threshold'] = 0.01
params['kurt_threshold'] = 0.1
params['peak_std_threshold'] = 0.8
params['upsample'] = 5
params['neighbor_radius'] = 100
params['r2_threshold'] = 0.8

# AW invention. I dont remember why I did this. 4Jan2023
params["theilsen_maxiter"] = 10000
pprint(params)

# Get the directory of the current script
script_dir = Path(os.path.dirname(os.path.realpath(__file__)))

# Get the parent directory (project root)
project_root = script_dir.parent

# Define the path to the reconstruction folder
reconstruction_dir = project_root / "data" / "temp_data" / "reconstruction"

# Create the reconstruction folder if it doesn't exist
reconstruction_dir.mkdir(parents=True, exist_ok=True)


def merge_templates(templates, channel_locations, plot_dir):
    ###Merge
    merged_templates = []
    merged_channel_loc = []
    
    for i in range(len(templates)):
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
                #sample k
                for j in range(len(channel_loc)):
                    #channel j
                    channel = tuple(channel_loc[j])
                    found_idx = None
                    for idx, existing_channel in enumerate(merged_channel_loc[-1]):
                        if tuple(existing_channel) == channel:
                            found_idx = idx
                            break

                    if found_idx is not None:
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

    with open(os.path.join(plot_dir, 'merged_template_save.pkl'), 'wb') as f:
        pickle.dump(merged_templates, f)
    with open(os.path.join(plot_dir, 'merged_channel_loc_save.pkl'), 'wb') as f:
        pickle.dump(merged_channel_loc, f)

    return merged_templates, merged_channel_loc

def generate_axon_trace_figure(template, channel_locations, fs, params, verbose_bool = False):

   
    #build 
    gtr = GraphAxonTracking(template[:,:].T, 
                            channel_locations, 
                            fs, 
                            verbose=verbose_bool, 
                            **params)
    gtr.select_channels()
    gtr.build_graph()
    gtr.find_paths()
    gtr._verbose = 0
    #gtr.clean_paths(remove_outliers=True)
    
    ##plot heuristics
    fig_graph = plt.figure(figsize=(12, 7))
    #ax1 = fig_graph.add_subplot(1,1,1)
    #fig_graph, figax = plt.subplots(1, 5)
    # fig_graph = gtr.plot_graph(node_search_labels=False, 
    #                             fig=fig_graph, 
    #                             cmap_nodes="viridis", 
    #                             cmap_edges="YlGn")   
    

    ## Create GridSpec
    #gs = plt.GridSpec(1, 6)
    #gs.set_width_ratios([2,0.25,2,0.25,2,2])
    
    ## Add new subplot
    ax_raw = fig_graph.add_subplot()  

    ##trace - raw
    #fpaths_raw, axpaths_raw = plt.subplots(figsize=(7, 10))
    axpaths_raw = ax_raw
    axpaths_raw = gtr.plot_raw_branches(
                                        cmap="tab20", 
                                        plot_bp=True, 
                                        plot_neighbors=True, 
                                        plot_full_template=True,
                                        ax=axpaths_raw)
    axpaths_raw.legend(fontsize=6)
    axpaths_raw.set_title("Raw Branches")

    ## Add new subplot
    #ax_clean = fig_graph.add_subplot(gs[:, -1])  

    ##trace - clean
    #fpaths_raw, axpaths_raw = plt.subplots(figsize=(7, 10))
    # axpaths_clean = ax_clean
    # axpaths_clean = gtr.plot_clean_branches(
    #                                     #cmap="tab20", 
    #                                     plot_bp=True, 
    #                                     #plot_neighbors=True, 
    #                                     plot_full_template=True,
    #                                     ax=axpaths_clean)
    # axpaths_clean.legend(fontsize=6)
    # axpaths_clean.set_title("Clean Branches")

    ## Format
    # for i in range(len(fig_graph.axes)): fig_graph.axes[i].set_subplotspec(gs[0, i])
    # # Distribute columns evenly
    # plt.subplots_adjust(wspace=.75)
    
    return fig_graph

def unique_fig_name(true_path):
    # Split the path using 'harddrive8tb' as a separator
    try: parts = true_path.split('harddrive8tb')
    except: 
        print('fix unique name function given new data source')
        return

    # Get the second part (everything after 'harddrive8tb')
    unique_name = parts[1].strip('/')

    # Replace slashes with underscores (or any other desired character)
    unique_name = unique_name.replace('/', '_')

    # Add a file extension (e.g., '.png') if needed
    #unique_name += '.png'

    return unique_name

def reconstruct_axon_trace_objects(axon_trace_precursors):
    """
    This function reconstructs axon trace objects from precursors.
    It iterates over the units in the precursors, merges their templates,
    and generates figures for each of them.
    """

    # Extract the units from the precursors
    units_to_be_merged = axon_trace_precursors.units_to_be_merged

    #Initialize list to store gtr objects
    gtrs = {}

    # Initialize an empty dictionary to store the unique channel locations
    MEA_channels = []
    channel_counter = 0
    
    # Initialize lists to keep track of accepted and rejected units for each figure
    accept_units_amp_map = []
    reject_units_amp_map = []
    accept_units_peak_latency_map = []
    reject_units_peak_latency_map = []
    accept_units_selected_channels = []
    reject_units_selected_channels = []
    accept_units_axon_reconstruction_heuristics = []
    reject_units_axon_reconstruction_heuristics = []
    accept_units_axon_reconstruction = []
    reject_units_axon_reconstruction = []
    accept_units_axon_reconstruction_velocities = []
    reject_units_axon_reconstruction_velocities = []

    # Iterate over the keys of the dictionary, each representing a unit
    #flag = False
    for key in units_to_be_merged.keys():
        #if flag:
        #    break

        # Get the full path
        full_path = units_to_be_merged[key]['source_path']

        # Split the path into components
        path_components = full_path[0].split(os.sep)

        # Find the index of the "sorting" component
        sorting_index = path_components.index('sorting')

        # Get the path components after "sorting"
        path_after_sorting = os.sep.join(path_components[sorting_index + 1:])

        # Extract the run_key for the current unit
        run_key = units_to_be_merged[key]['run_key']

        # Create a string from the run_key array
        run_key_str = ''.join([
            rk[0] + 'S' + re.search(r'\d+', rk).group() if 'ActivityScan' in rk 
            else rk[0] + re.search(r'\d+', rk).group() if 'Network' in rk 
            else rk[0]  + 'X' + re.search(r'\d+', rk).group() if 'AxonTracking' in rk 
            else rk 
            for rk in run_key
        ])

        # Create a unique folder name based on the run_key and unit_id
        folder_name = f"Unit{key}_{run_key_str}"

        # Define the path to the plot folder
        plot_dir = reconstruction_dir / path_after_sorting / folder_name

        # Create the plot folder if it doesn't exist
        plot_dir.mkdir(parents=True, exist_ok=True)

        # Merge the templates and channel locations for the current unit
        #debug
        new_templates = units_to_be_merged[key]['sparse_unit_templates']
        new_channel_locations = units_to_be_merged[key]['sparse_channel_locations']        
        old_channels = units_to_be_merged[key]['channel_locations']
        old_templates = units_to_be_merged[key]['templates']
        #debug

        ## Create a dictionary to hold all the data
        data = {
            'new_templates': new_templates,
            'new_channel_locations': new_channel_locations,
            'old_channels': old_channels,
            'old_templates': old_templates
        }

        # Define the path to the JSON file
        json_file_path = plot_dir / 'data.json'

        # Write the data to the JSON file
        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                return json.JSONEncoder.default(self, obj)

        # Write the data to the JSON file
        with open(json_file_path, 'w') as json_file:
            json.dump(data, json_file, cls=NumpyEncoder)

        print(f'Data saved to {json_file_path}')
                
        #Build MEA_channel_locations
        # Iterate over the old_channels dictionary
        # Iterate over the old_channels dictionary
        for channel_locations_list in old_channels:
            for channel_location in channel_locations_list:
                # Convert the numpy array to a tuple and append it to MEA_channels
                MEA_channels.append(tuple(channel_location))

        # Convert MEA_channels to a set to keep only unique elements
        MEA_channels = set(MEA_channels)

        # Convert the set back to a list
        MEA_channels = list(MEA_channels)    

        templates = new_templates
        channel_locations = new_channel_locations

        # Load or generate merged templates
        load_if_available = True
        template_path = plot_dir / "merged_template_save.pkl"
        channel_path = plot_dir / "merged_channel_loc_save.pkl"

        try:
            if load_if_available and template_path.exists() and channel_path.exists():
                with open(template_path, 'rb') as f:
                    merged_templates = pickle.load(f)
                with open(channel_path, 'rb') as f:
                    merged_channel_loc = pickle.load(f)
            else:
                raise FileNotFoundError
        except FileNotFoundError:
            # If the pickle files don't exist or loading is not enabled, run the function to generate the data
            merged_templates, merged_channel_loc = merge_templates(templates, channel_locations, plot_dir)

        # Iterate over the merged templates
        for i, template in enumerate(merged_templates):
           
            # Extract the sampling frequency for the current unit
            fs = units_to_be_merged[key]['sampling_frequency'][i]        

            transformed_template = np.array(template).T
            
            # Generate amplitude_map figure
            try:
                fig_amp = plt.figure(figsize=(10, 5))
                ax_amp = fig_amp.add_subplot(111)
                ax_amp = plot_amplitude_map(transformed_template, merged_channel_loc[i], log=True, ax=ax_amp, cmap="PRGn", colorbar=True, colorbar_orientation="horizontal")
                ax_amp.set_title(f"Amplitude", fontsize=20)

                #Save figure
                fig_name = "amplitude_map.png"
                fig_path = plot_dir / fig_name
                plt.savefig(fig_path, dpi=300, bbox_inches='tight')
                plt.close()

                accept_units_amp_map.append(key)
            except Exception as e:
                reject_units_amp_map.append(key)
                print(f"unit {key} failed to generate amplitude_map, error: {e}")

            # Generate peak latency map figure
            try:
                fig_peaks = plt.figure(figsize=(10, 5))
                ax_peaks = fig_peaks.add_subplot(111)
                ax_peaks = plot_peak_latency_map(transformed_template, merged_channel_loc[i], fs=fs, log=False, ax=ax_peaks, colorbar=True, colorbar_orientation="horizontal")
                ax_peaks.set_title(f"Peak latency", fontsize=20)

                #Save figure
                fig_name = "peak_latency_map.png"
                fig_path = plot_dir / fig_name
                plt.savefig(fig_path, dpi=300, bbox_inches='tight')
                plt.close()

                accept_units_peak_latency_map.append(key)
            except Exception as e:
                reject_units_peak_latency_map.append(key)
                print(f"unit {key} failed to generate peak_latency_map, error: {e}")

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
                gtr0 = GraphAxonTracking(transformed_template, merged_channel_loc[i], fs, verbose=verbose_bool, **params)
                

                # Select channels
                gtr0.select_channels()                

                # Use the function for each block
                plot_selected_channels(f"Selected after detection threshold: {gtr0._detect_threshold} ms", 
                                    np.array(list(gtr0._selected_channels_detect)), gtr0.locations, "detection_threshold.png")
                plot_selected_channels(f"Selected after peak std threshold: {gtr0._peak_std_threhsold} ms", 
                                        np.array(list(gtr0._selected_channels_peakstd)), gtr0.locations, "peak_std_threshold.png")
                plot_selected_channels(f"Selected after init delay threshold: {gtr0._init_delay} ms", 
                                        np.array(list(gtr0._selected_channels_init)), gtr0.locations, "init_delay_threshold.png")
                plot_selected_channels("Selected after all thresholds", gtr0.selected_channels, gtr0.locations, "all_thresholds.png")
                accept_units_selected_channels.append(key)
            except Exception as e:
                reject_units_selected_channels.append(key)
                print(f"unit {key} failed to generate selected_channels, error: {e}")              

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

                accept_units_axon_reconstruction_heuristics.append(key)
            except Exception as e:
                reject_units_axon_reconstruction_heuristics.append(key)
                print(f"unit {key} failed to generate axon_reconstruction_heuristics, error: {e}")

            # Generate axon_reconstruction figure
            try:
                fig_graph = plt.figure(figsize=(12, 7))
                ax_raw = fig_graph.add_subplot(111)                
                axpaths_raw = ax_raw
                axpaths_raw = gtr0.plot_raw_branches(cmap="tab20", plot_bp=True, plot_neighbors=True, 
                                                    plot_full_template=True, ax=axpaths_raw)
                axpaths_raw.legend(fontsize=6)
                axpaths_raw.set_title("Raw Branches")

                fig_name = "axon_reconstruction.png"
                fig_path = plot_dir / fig_name
                plt.savefig(fig_path, dpi=600, bbox_inches='tight')
                plt.close()
                gtrs[key] = compute_graph_propagation_velocity(template.T, merged_channel_loc[i], fs, 
                                                        verbose=False, **params)

                accept_units_axon_reconstruction.append(key)
                #flag = True
                #break
            except Exception as e:
                reject_units_axon_reconstruction.append(key)
                print(f"unit {key} failed to generate axon_reconstruction, error: {e}")
             
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

                accept_units_axon_reconstruction_velocities.append(key)
            except Exception as e:
                reject_units_axon_reconstruction_velocities.append(key)
                print(f"unit {key} failed to generate axon_reconstruction_velocities, error: {e}")
            
    # Print the number of accepted and rejected units for each figure
    print(f"Amplitude map: {len(accept_units_amp_map)} accepted, {len(reject_units_amp_map)} rejected")
    print(f"Peak latency map: {len(accept_units_peak_latency_map)} accepted, {len(reject_units_peak_latency_map)} rejected")
    print(f"Selected channels: {len(accept_units_selected_channels)} accepted, {len(reject_units_selected_channels)} rejected")
    print(f"Axon reconstruction heuristics: {len(accept_units_axon_reconstruction_heuristics)} accepted, {len(reject_units_axon_reconstruction_heuristics)} rejected")
    print(f"Axon reconstruction: {len(accept_units_axon_reconstruction)} accepted, {len(reject_units_axon_reconstruction)} rejected")
    print(f"Axon reconstruction velocities: {len(accept_units_axon_reconstruction_velocities)} accepted, {len(reject_units_axon_reconstruction_velocities)} rejected")

    #Plot all reconstructed axons

    MEA = plotting.get_probe(MEA_channels)

    fig_mea1k, ax = plt.subplots(figsize=(10, 7))
    _ = plotting_probe.plot_probe(MEA, ax=ax, contacts_kargs={"alpha": 0.1}, probe_shape_kwargs={"alpha": 0.1})
    ax.axis("off")

    i = 0
    i_sel = 0
    cmap = "tab20"
    cm = plt.get_cmap(cmap)
    for i, gtr in gtrs.items():        
        
        color = cm(i / len(gtrs))
        lw = 1
        alpha = 1
        zorder = 1
        try:
            #if len(gtr.branches) > 0:
            ax.plot(gtr.locations[gtr.init_channel, 0], gtr.locations[gtr.init_channel, 1], 
                    marker="o", markersize=5, color=color, alpha=alpha, zorder=zorder)

            # for visualization purposes, plot raw branches
            for b_i, path in enumerate(gtr._paths_raw):
                if b_i == 0:
                    ax.plot(gtr.locations[path, 0], gtr.locations[path, 1], marker="", color=color,
                            lw=lw, alpha=alpha, zorder=zorder, label=i)
                else:
                    ax.plot(gtr.locations[path, 0], gtr.locations[path, 1], marker="", color=color,
                            lw=lw, alpha=alpha, zorder=zorder)
        except:
            pass


    #ax.plot([0, 500], [1900, 1900], color="k", marker="|")
    #ax.text(20, 1950, "500$\mu$m", color="k", fontsize=18)
    ax.set_title("")

    # So that this is comparable to scope, Invert the y-axis to flip the plot vertically
    ax.invert_yaxis()

    fig_name = "all_axon_reconstruction.png"
    plot_dir = reconstruction_dir / path_after_sorting
    fig_path = plot_dir / fig_name
    plt.savefig(fig_path, dpi=600, bbox_inches='tight')
