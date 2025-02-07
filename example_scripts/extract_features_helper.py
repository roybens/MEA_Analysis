# use this script to process experimental data and exract features to be used developing network models.
# primarily we aim to get experimental data into a form that is similar to the simulated data 
# - such that the same analysis can be applied to both.
import os
import glob
import MEA_Analysis.MEAProcessingLibrary.mea_processing_library as mea
# NOTE: to import MEA_Analysis this way, need to pip install like so:
# pip install -e **/MEA_Analysis # ** = whatever path your MEA_Analysis copy is in.

import numpy as np
import matplotlib.pyplot as plt
# =============================================================================
def extract_network_features(
    raw_data_paths, 
    sorted_data_dir = None, 
    output_dir = None, 
    stream_select=None, 
    plot=True, 
    conv_params=None,
    mega_params=None,
    ):
    
    #init
    assert output_dir is not None, f"Error: output_dir is None"
    output_dir = os.path.abspath(output_dir)
    #assert os.path.exists(output_dir), f"Error: output_dir does not exist. {output_dir}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
                                                                                   
    #assert sorted_data_dir is not None, f"Error: sorted_data_dir is None"  
    raw_data_paths = [os.path.abspath(raw_data_path) for raw_data_path in raw_data_paths]
    sorted_data_dir = os.path.abspath(sorted_data_dir)
    
    #assert conv_params is defined, this is required for network analysis
    assert conv_params is not None, f"Error: conv_params is None - must be provided for network analysis"
    assert mega_params is not None, f"Error: mega_params is None - must be provided for network analysis"
    # mega_params = conv_params.mega_params
    # conv_params = conv_params.conv_params
    
    
    #iterate through raw data paths and get all h5 files
    h5_paths = []
    for raw_data_path in raw_data_paths:
        if os.path.isdir(raw_data_path):
            h5_paths.extend(glob.glob(os.path.join(raw_data_path, '**/*.h5'), recursive=True))
        else:
            h5_paths.append(raw_data_path)
                
    # iterate through sorter_output_dir and get all well output directories
    sorted_output_folders = []
    for root, dirs, files in os.walk(sorted_data_dir):
        # for dir in dirs:
        #     sorted_data_paths.append(os.path.join(root, dir))
        if root.endswith('sorter_output'):
            sorted_output_folders.append(os.path.join(root))
     
    # get recording details from each h5 file, check for match where all details match in sorted_output_folders
    # this way, we pair up recordings with their corresponding sorting output
    # get paired data objects - network analysis requires both recording and sorting objects
    def get_data_obj_pairs():
        well_data_list = []
        path_pairs = []
        for h5_path in h5_paths:
            # get recording details
            recording_details = mea.extract_recording_details(raw_data_path)[0]     # NOTE: this works for a list of dirs or a single dir - but treats single dir 
                                                                                        # as a list of a single dir
            # h5_file_path = recording_details['h5_file_path']
            # runID = recording_details['runID']
            # scanType = recording_details['scanType']
            # chipID = recording_details['chipID']
            # date = recording_details['date']
            # projectName = recording_details['projectName'] 
            
            #remove h5_file_path from recording_details - this wont match, this is the old path
            h5_path = recording_details.pop('h5_file_path')
            
            for sorted_output_folder in sorted_output_folders:         
                
                # check if all details match, long form for debugging
                # elements = [recording_details[key] for key in recording_details.keys()]
                # elements_in_sorted_output_folder = [element in sorted_output_folder for element in elements]
                # found = all(elements_in_sorted_output_folder)
                
                # shortform
                found = all([f'/{recording_details[key]}/' in sorted_output_folder for key in recording_details.keys()])
                
                if found:
                    path_pairs.append((h5_path, sorted_output_folder))
                    # get recording and sorting objects
                    #print(f"Found matching recording and sorting objects for {h5_path}")
                    
                    def load_recording_and_sorting_objects():
                        # this function expects to load one sorting_obj and one recording_obj for the related well.
                        # there may be any number of rec_segments in the recording_obj, but only one sorting_obj
                        
                        # get stream_select from sorted_output_folder path
                        # look for the word 'well' in the string. It will beb followed by three digits.
                        # get the int value of those digits. Should be 0-5
                        stream_select = int(sorted_output_folder.split('well')[1][:3])
                        wellid = f'well{str(0).zfill(2)}{stream_select}'
                        
                        # try to load sorting object
                        try: sort_obj = mea.load_kilosort2_results(sorted_output_folder)
                        except Exception as e:
                            print(f"Error: Could not load sorting object for {sorted_output_folder}")
                            print(e)
                            sort_obj = e # put error in sort_obj for debugging
                        
                        #try to load recording object
                        try: 
                            _, well_recs, _, _ = mea.load_recordings(h5_path, stream_select=stream_select)
                            rec_segments = well_recs[wellid]
                        except Exception as e:
                            print(f"Error: Could not load recording object for {h5_path}")
                            print(e)
                            #well_recs = e
                            rec_segments = e
                        
                        # return paired objects
                        return (rec_segments, sort_obj), recording_details                        
                        #print(f"Loaded recording and sorting objects for {h5_path}")
                    well_data, recording_details = load_recording_and_sorting_objects()
                    well_data_list.append(well_data)
                    
        return well_data_list, recording_details, path_pairs
    well_data_list, recording_details, path_pairs = get_data_obj_pairs()
    
    # iterate through data_obj_list and get network metrics for each pair
    for well_data in well_data_list:
        try:
            recording_segments, sort_obj = well_data    # NOTE: I think in general, for this whole process, we shouldnt need more than one
                                                        # recording segment per well - we just need to get attributes from the recording object
                                                        # TODO: if so, reduce the amount of data we're passing around to just the attributes we need.
            
            #if sort_obj is an error, skip
            if isinstance(sort_obj, Exception):
                #print(f"Error: Could not load sorting object for {sort_obj}")
                print(f"Error: Could not load sorting object.")
                print(f'Error details: {sort_obj}')
                print(f"Skipping well...")
                continue
            stream_id = recording_segments[0].stream_id
            stream_num = int(stream_id.split('well')[1][:3])
            kwargs = recording_details.copy()
            extract_network_metrics(
                recording_segments[0], # TODO: testing, only sending one seg.
                sort_obj, 
                stream_num,
                conv_params,
                mega_params,                
                # mega_params,
                # conv_params, 
                output_dir,             
                plot=plot, details=recording_details, **kwargs)    
        except Exception as e:
            print(f"Error: Could not get network metrics for:")
            from pprint import pprint
            #pprint('recording_segments:', recording_segments)
            #pprint(recording_segments)
            pprint(recording_details)
            pprint(sort_obj)
            print(e)
            continue         
    print('done')
        
def extract_network_metrics(
        recording_object, 
        sorting_object, 
        stream_num, 
        conv_params,
        mega_params,
        #convolution_params_path, 
        output_path,
        save_path=None,
        bursting_plot_path=None, bursting_fig_path=None, 
        plot=False, **kwargs):
    
    # 
    assert sorting_object is not None, f"Error: sorting_object is None"
    
    # import network analysis module
    from RBS_network_models.network_analysis import get_experimental_network_activity_metrics
    
    # get network metrics
    well_id = f'well{str(0).zfill(2)}{stream_num}'
    well_recording_segment = recording_object 
    try: 
        network_metrics = get_experimental_network_activity_metrics(sorting_object, well_recording_segment, conv_params, mega_params)
    except Exception as e:
        print(e)
        print(f"Error: Could not get network metrics for {well_id}")
    
    # get mega network metrics
    #mega_burst_conv_params = conv_params
    #mega_burst_conv_params.conv_params['binSize'] *= 5 # HACK: this is a hack to make the mega burst conv params bin size 10x larger
    #mega_burst_conv_params.conv_params['gaussianSigma'] *= 15 # HACK: this is a hack to make the mega burst conv params sigma 10x larger
    # mega_burst_conv_params = mega_params
    # mega_metrics = get_experimental_network_activity_metrics(sorting_object, well_recording_segment, mega_burst_conv_params)
    
    # add mega metrics to network metrics
    #network_metrics['mega_bursting_data'] = mega_metrics['bursting_data']
    
    # get recording details
    print("Saving network metrics as numpy...")
    recording_details = kwargs['details']
    projectName = recording_details['projectName']
    date = recording_details['date']
    chipID = recording_details['chipID']
    scanType = recording_details['scanType']
    runID = recording_details['runID']
    
    # save network metrics
    if save_path is None: 
        save_path = os.path.join(output_path, f"{projectName}_{date}_{chipID}_{scanType}_{runID}_network_metrics_well00{stream_num}.npy")
    # save_path = os.path.join(output_path, f"network_metrics_well00{stream_num}.npy")
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    np.save(save_path, network_metrics)
    
    # #save mega network metrics
    # mega_save_path = save_path.replace('.npy', '_mega.npy')
    # np.save(mega_save_path, mega_metrics)
    
    #plot network activity
    #print("Plotting network activity...")
    if plot:
        try: 
            print("Generating network summary plot...")
            # aw 2025-01-20 17:25:21 - I guess I'll just plot both for now. I like how mine looks, but additional context is nice for Roy.
            save_path_3p = save_path.replace('.npy', '_3p.pdf') 
            plot_network_metrics(
                network_metrics, 
                bursting_plot_path, 
                bursting_fig_path,
                save_path=save_path_3p,
                #mode = '2p',
                mode = '3p',
                )
            
            save_path_2p = save_path.replace('.npy', '_2p.pdf')
            plot_network_metrics(
                network_metrics, 
                bursting_plot_path, 
                bursting_fig_path,
                save_path=save_path_2p,
                mode = '2p',
                )  
        except Exception as e:
            print(e)
            print(f"Error: Could not plot network activity for {well_id}")  

    # return network metrics and save path
    print(f"Saved network metrics to {save_path}")
    return network_metrics, save_path
    
def plot_network_metrics(
    network_metrics,
    bursting_plot_path,
    bursting_fig_path,
    save_path=None,
    mode='2p',
    ):
    
    # TODO: blend with plot comparison plot? I think.
    #from RBS_network_models.network_analysis import plot_network_summary
    from .network_analysis import plot_network_summary
    plot_network_summary(network_metrics, bursting_plot_path, bursting_fig_path, save_path=save_path, mode=mode)   
    
# =============================================================================
# reference
def save_network_metrics(recording_object, sorting_object, stream_num, 
                         convolution_params_path, output_path,
                         bursting_plot_path=None, bursting_fig_path=None, 
                         plot=False, **kwargs):
    add_repo_root_to_sys_path()
    #import _external.RBS_network_simulation_optimization_tools.modules.analysis.analyze_network_activity as ana
    import RBS_network_models.developing.utils.analyze_network_activity as ana 
    #load recordings and sorting data                    
    #get recording segment - note: network scans should only have one recording segment
    well_id = f'well{str(0).zfill(2)}{stream_num}'
    #recording_object = kwargs['recording_object']
    well_recording_segment = recording_object[well_id]['recording_segments'][0]  #note: network scans should only have one recording segment
    
    #get network activity metrics
    #import conv_params from convolution_params_path
    conv_params = import_module_from_path(convolution_params_path)
    conv_params = conv_params.conv_params
    #sorting_object = kwargs['sorting_object']
    network_metrics = ana.get_experimental_network_activity_metrics(sorting_object, well_recording_segment, conv_params)
    
    # save network metrics
    print("Saving network metrics as numpy...")
    recording_details = kwargs['details']
    projectName = recording_details['projectName']
    date = recording_details['date']
    chipID = recording_details['chipID']
    scanType = recording_details['scanType']
    runID = recording_details['runID']
    save_path = os.path.join(output_path, f"{projectName}_{date}_{chipID}_{scanType}_{runID}_network_metrics_well00{stream_num}.npy")
    #save_path = os.path.join(output_path, f"network_metrics_well00{stream_num}.npy")
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    np.save(save_path, network_metrics)
    
    #plot network activity
    print("Plotting network activity...")
    if plot:
        fig, ax = plt.subplots(2, 1, figsize=(16, 9))
        if bursting_plot_path is None and bursting_fig_path is None:
            #modify save path
            raster_plot_path = save_path.replace('.npy', '_raster_plot.pdf')
            raster_fig_path = save_path.replace('.npy', '_raster_fig.pkl')
            raster_plot_path_png = save_path.replace('.npy', '_raster_plot.png')
            bursting_plot_path = save_path.replace('.npy', '_bursting_plot.pdf')
            bursting_fig_path = save_path.replace('.npy', '_bursting_fig.pkl')
            bursting_plot_path_png = save_path.replace('.npy', '_bursting_plot.png')
        # assert bursting_plot_path is not None, f"Error: bursting_plot_path is None"
        # assert bursting_fig_path is not None, f"Error: bursting_fig_path is None"
        # bursting_plot_path_png = save_path.replace('.npy', '_bursting_plot.png')

        #generate_network_bursting_plot(network_metrics, bursting_plot_path, bursting_fig_path)
        print("Generating raster plot...")
        spiking_data_by_unit = network_metrics['spiking_data']['spiking_data_by_unit']
        ax[0] = plot_raster_plot_experimental(ax[0], spiking_data_by_unit)       
        print("Generating network bursting plot...")
        bursting_ax = network_metrics['bursting_data']['bursting_summary_data']['ax']
        ax[1] = plot_network_bursting_experimental(ax[1], bursting_ax)
        
        #save plots
        plt.tight_layout()
        fig.savefig(bursting_plot_path) #save as pdf
        print(f"Network summary plot saved to {bursting_plot_path}")
        fig.savefig(bursting_plot_path_png, dpi=600) #save as png
        print(f"Network summary plot saved to {bursting_plot_path_png}")
    
    print(f"Saved network metrics to {save_path}")

def load_recording_and_sorting_object_tuples(recording_paths, sorting_output_parent_path, stream_select=None):
    paired_objects = []
    for recording_path in recording_paths:
    
        #prep sorting objects
        sorting_output_paths = check_for_sorting_objects_associated_to_recordings(
            recording_paths, 
            sorting_output_parent_path
            )
        
        recording_details = mea_lib.extract_recording_details(recording_path)[0]
        
        for sorting_output_path in sorting_output_paths:
            stream_sorting_object = load_spike_sorted_data(sorting_output_path)

            #load recordings and sorting data
            MaxID, recording_object, expected_well_count, rec_counts = mea_lib.load_recordings(recording_path, stream_select=stream_select)
            
            # Append the paired objects as a tuple
            paired_objects.append((recording_object, stream_sorting_object, recording_details))
            
    return paired_objects

def load_recording_and_sorting_object_tuples_dep(recording_paths, sorting_output_parent_path, stream_select=None):
    paired_objects = []
    for recording_path in recording_paths:
    
        #prep sorting objects
        sorting_output_paths = check_for_sorting_objects_associated_to_recordings(
            recording_paths, 
            sorting_output_parent_path
            )
        
        recording_details = mea_lib.extract_recording_details(recording_path)[0]
        
        for sorting_output_path in sorting_output_paths:
            stream_sorting_object = load_spike_sorted_data(sorting_output_path)

            #load recordings and sorting data
            MaxID, recording_object, expected_well_count, rec_counts = mea_lib.load_recordings(recording_path, stream_select=stream_select)
            
            # Append the paired objects as a tuple
            paired_objects.append((recording_object, stream_sorting_object, recording_details))
            
    return paired_objects

