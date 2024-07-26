'''Main script to execute the AxonReconPipeline'''

'''Imports'''
import shutil
import axon_velocity as av

# Local Python Function Libraries
import lib_helper_functions as helper
import lib_sorting_functions as sorter
import lib_waveform_functions as waveformer
import lib_template_functions as templater
import lib_axon_velocity_functions as axoner

# Logger Setup
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler()  # logs to console
stream_handler.setLevel(logging.DEBUG)
logger.addHandler(stream_handler)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - %(module)s.%(funcName)s')
stream_handler.setFormatter(formatter)

def attempt_load_merged_templates(continuous_h5_file_info, templates_dir, allowed_scan_types, stream_select = None):
    logger.info('Attempting to load merged templates from disk without validating template segments, waveforms, or spike sorting...')
    logger.info('WARNING: Only use this option if you are sure that the merged templates are correctly extracted and saved on disk.')
    try:
        templater.extract_all_templates(
            continuous_h5_file_info, 
            qc_params=qc_params, te_params=te_params, 
            allowed_scan_types=allowed_scan_types, 
            templates_dir = templates_dir, 
            stream_select = stream_select, 
            just_load_waveforms = load_waveforms,
            quick_load_templates = True
            )
        return True
    except Exception as e:
        logger.error(f'Failed to load merged templates from disk: {e}')
        logger.error('Continuing with normal template extraction...')
        return False   

def main(h5_parent_dirs, allowed_scan_types=['AxonTracking'], stream_select = None):

    '''Delete temp data if it exists'''
    temp_data_path = './AxonReconPipeline/data/temp_data'
    #if os.path.exists(temp_data_path): shutil.rmtree(temp_data_path, ignore_errors=True)
    
    '''Select folders and get continuous h5 filepaths'''
    continuous_h5_file_info = helper.select_folders_and_get_continuous_h5_dirs(pre_selected_folders = h5_parent_dirs, debug_mode = debug_mode, stream_select = stream_select) #TODO: add stream select here.
    
    '''Attempt to load temporary data from disk'''
    # quick_templates_loaded = False
    # if load_merged_templates: quick_templates_loaded = attempt_load_merged_templates(continuous_h5_file_info, templates_dir, allowed_scan_types, stream_select = stream_select)
    
    '''Execute Main Pipeline Steps, as needed'''
    #if quick_templates_loaded is False:
    '''Generate multirecordings objects for spikesorting and waveform extraction steps'''
    #sorter.build_multirecording_objects(continuous_h5_file_info, allowed_scan_types, recordings_dir, stream_select = stream_select, n_jobs = 8)
    
    '''Extract sorting objects per stream from sliced and cocatenated recordings (most/all fixed electrodes during axontracking scan)'''
    sorter.spikesort_recordings(continuous_h5_file_info, sortings_dir = sortings_dir, allowed_scan_types = allowed_scan_types, stream_select = stream_select,)
    
    '''Extract Waveforms from sorting objects'''
    waveformer.extract_all_waveforms(continuous_h5_file_info, qc_params=qc_params, te_params=te_params, stream_break = stream_select)
    if run_lean: shutil.rmtree(sortings_dir, ignore_errors=True)
    
    '''Get Templates from Waveform Objects'''
    templater.extract_all_templates(continuous_h5_file_info, qc_params=qc_params, te_params=te_params, allowed_scan_types=allowed_scan_types, templates_dir = templates_dir, stream_select = stream_select, just_load_waveforms = load_waveforms)   
    if run_lean: shutil.rmtree(waveforms_dir, ignore_errors=True)
    
    '''Analysis and Reconstruction Steps'''
    axoner.analyze_and_reconstruct(continuous_h5_file_info, templates_dir, params, stream_select = stream_select, n_jobs = 32)
    #if run_lean: shutil.rmtree(templates_dir, ignore_errors=True)

    '''Copy reconstruction data to destination directories'''
    #recon_dir = Path(temp_dir) / 'reconstructions'
    #helper.copy_reconstruction_data(h5_parent_dirs, destination_dirs, recon_dir)

if __name__ == "__main__":

    '''directories'''
    recordings_dir = './AxonReconPipeline/data/temp_data/recordings'
    sortings_dir = './AxonReconPipeline/data/temp_data/sortings'
    waveforms_dir = './AxonReconPipeline/data/temp_data/waveforms'
    templates_dir='./AxonReconPipeline/data/temp_data/templates'
    recon_dir = './AxonReconPipeline/data/temp_data/reconstructions'
    KCNT1_effort = [        
        #'/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/KCNT1_T4_C1_04122024/KCNT1_T4_C1_04122024/240503/M08034/AxonTracking/000082', #homo and het 
        #'/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/KCNT1_T3_C1_03122024/KCNT1_T3_C1_03122024/240409/M07037/AxonTracking/000095', #homo, het, wt       
        #'/mnt/disk20tb/KCNT1_T3_NeuronalScans/KCNT1_T3/KCNT1_T3/240328/M07037/AxonTracking/000055',  #data I'm using for refactoring, since NAS is out.
        
        '/mnt/disk20tb/PrimaryNeuronData/Maxtwo/KCNT1_T3_C1_03122024/KCNT1_T3_C1_03122024/240322/M07037/AxonTracking', #homo, het, wt  DIV 10
        #'/mnt/disk20tb/PrimaryNeuronData/Maxtwo/KCNT1_T3_C1_03122024/KCNT1_T3_C1_03122024/240328/M07037/AxonTracking', #homo, het, wt  DIV 16
        #'/mnt/disk20tb/PrimaryNeuronData/Maxtwo/KCNT1_T3_C1_03122024/KCNT1_T3_C1_03122024/240405/M07037/AxonTracking', #homo, het, wt  DIV 24
        #'/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/KCNT1_T3_C1_03122024/KCNT1_T3_C1_03122024/240409/M07037/AxonTracking', #homo, het, wt  DIV 28
        #'/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/KCNT1_T3_C1_03122024/KCNT1_T3_C1_03122024/240412/M07037/AxonTracking', #homo, het, wt  DIV 31

    ]
    h5_parent_dirs = KCNT1_effort
    
    #Check if all paths exist
    import os
    for i in range(len(h5_parent_dirs)): assert os.path.exists(h5_parent_dirs[i]), f"Path {h5_parent_dirs[i]} does not exist."

    #check cwd 
    print(f"Current working directory: {os.getcwd()}")

    #for each dir in h5_parent_dirs, get the name of the line, which should directly preceed the date folder (YYMMDD)
    destination_parent_dir = '/mnt/ben-shalom_nas/axonal_reconstruction_analyses/'
    destination_dirs = helper.build_destination_dirs(h5_parent_dirs, destination_parent_dir)

    '''runtime/loading options''' 
    debug_mode = True #avoid launching gui for file selection, use pre-selected folders 
    load_merged_templates = False #Load merged templates from disk instead of extracting them, if possible
    load_template_segments = True #Load templates from disk instead of extracting them, if possible
    load_waveforms = True #Load waveforms from disk instead of extracting them, if possible
    load_sortings = True #Load sortings from disk instead of extracting them, if possible
    run_lean = False #Delete sortings after waveform extraction to save disk space

    '''allowed scan types'''
    allowed_scan_types = [
        #"ActivityScan", Dense activity scan only needed to define locations of electrodes in network and axon tracking scans
        "AxonTracking",   #Set allowed scan types, axontracking alone is currently the best option
        #"Network"
        ]

    '''template extraction params'''
    te_params = dict()
    te_params['align_cutout'] = True #Align waveforms by max waveform peak
    te_params['upsample'] = 2 #Factor by which to upsample waveforms
    te_params['rm_outliers'] = True #Check if outliers should be removed
    te_params['n_jobs'] = 16 #Number of cores to use for waveform extraction
    te_params['n_neighbors'] = 10 #Number of neighbors for outlier detection
    te_params['peak_cutout'] = 2 #Looking for peak +- this value around the expected peak (removing minor offsets)
    te_params['overwrite_wf'] = False #Flag if waveform extraction should be repeated (e.g. different cutouts)
    te_params['overwrite_tmp'] = False #Flag if templates should be recalculated if already existing

    '''quality control params'''
    qc_params = dict()
    qc_params['min_n_spikes'] = 500 #Minimum number of spikes to be detected for a unit for template extraction to take place
    qc_params['exclude_mua'] = True #Exclude units that were labelled multi unit activity by kilosort

    '''reconstruction and analysis params'''
    #Define Algorithm Params
    params = av.get_default_graph_velocity_params()
    options = ['axon_velocity', 'organoid']
    option = options[0]
    if option == 'axon_velocity':
        # change params (fig6 Buccino axon_velocity paper example)
        params['detect_threshold'] = 0.01
        params['kurt_threshold'] = 0.1
        params['peak_std_threshold'] = 0.8
        params['upsample'] = 5
        params['neighbor_radius'] = 100
        params['r2_threshold'] = 0.8
    elif option == 'organoid':
        ## Modulating params for organoid data

        ## NOTE: These params didnt work very well to get organoid data. We concluded there are fundamental differences between organoid and MEA data that need to be addressed.

        # General
        params['upsample'] = 5 #1  # Upsample waveforms by this factor, upsample = 1 means no upsampling
        params['min_selected_points'] = 20 #30  # Minimum number of points to be selected for graph propagation, typically 30
        params['verbose'] = True  # Verbose output

        # Channel selection
        params['detect_threshold'] = 0.01  # Detection threshold (with respect to channel featuring maximal signal)
        params['detection_type'] = 'relative'  # Whether to use an 'absolute' or 'relative' detection threshold
        params['kurt_threshold'] = 0.05 #0.3  # Kurtosis threshold below which a channel is discarded
        params['peak_std_threshold'] = 0.5 #1  # Peak time standard deviation threshold in ms below which a channel is discarded
        params['init_delay'] = 0 #0.1  # Initial delay in seconds (with respect to maximum channel) below which a channel is discarded
        params['peak_std_distance'] = 50 #30  # Distance in µm to select channel neighborhood to compute peak time standard deviation
        params['remove_isolated'] = False #True  # If True, isolated channels are removed from selection

        # Graph
        params['init_amp_peak_ratio'] = 0.2  # Scalar value that weighs the contribution of the amplitude and the peak latency for h_init (ω_init in equation (2))
        params['max_distance_for_edge'] = 200 #100  # Maximum distance in µm between channels to create a graph edge
        params['max_distance_to_init'] = 500 #200  # Maximum distance in µm between a channel and the init_channel to create a graph edge
        params['n_neighbors'] = 5 #3  # Maximum number of edges that one channel can connect to below which an axonal branch is discarded
        params['distance_exp'] = 2  # Exponent for distance computation (e in equation (3))
        params['edge_dist_amp_ratio'] = 0.5 #0.3  # Relative weight between distance and amplitude to select neighbor nodes for graph edges

        # Axonal reconstruction
        params['min_path_length'] = 50 #100  # Minimum axon path length in µm to include an axonal branch
        params['min_path_points'] = 4 #5  # Minimum number of channels in an axon path to include an axonal branch
        params['neighbor_radius'] = 150 #100  # Radius in µm to exclude neighboring channels around an identified path
        params['min_points_after_branching'] = 3  # Minimum number of points after a branching to avoid pruning

        # Path cleaning/velocity estimation
        params['mad_threshold'] = 8  # Threshold in median absolute deviations on the fit error to consider points as outliers in the velocity estimation
        params['split_paths'] = True  # If True, the final path splitting step is enabled
        params['max_peak_latency_for_splitting'] = 0.5  # If a jump in the peak latencies of a path exceeds this value, the path can be split in sub-paths
        params['r2_threshold'] = 0.7 #0.9  # R^2 threshold for velocity linear fit below which an axon branch is discarded
        params['r2_threshold_for_outliers'] = 0.9 #0.98  # R^2 threshold below which outliers are detected and removed
        params['min_outlier_tracking_error'] = 75 #50  # Tracking error in µm above which a point can be considered an outlier and removed

    '''data selection options'''
    #stream_select = 0 #KCNT1, homo
    #stream_select = 3 #KCNT1, het
    #stream_select = 2 #KCNT1, wt
    stream_select = None

    '''RUN'''
    main(h5_parent_dirs, allowed_scan_types = allowed_scan_types, stream_select = stream_select)