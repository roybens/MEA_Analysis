# run_pipeline.py
import lib_helper_functions as helper
from axon_reconstructor import AxonReconstructor

# Set up kwargs with default reconstruction parameters from the paper
kwargs = {
    # runtime options
    'sorting_params': {
        'allowed_scan_types': ['AxonTracking'],
        'load_existing_sortings': True,
        'keep_good_only': False,
    },
    'te_params': {
        'load_merged_templates': True,
        'save_merged_templates': True,
    },
    'av_params': {
        'upsample': 1,  # Upsampling factor for template
        'min_selected_points': 30,  # Minimum number of selected points to run axon tracking
        'verbose': True,  # If True, the output is verbose

        # Channel selection
        'detect_threshold': 0.02,  # Detection threshold for amplitude
        'detection_type': 'relative',  # Detection threshold type
        'kurt_threshold': 0.3,  # Kurtosis threshold for noise filtering
        'peak_std_threshold': 1.0,  # Peak time standard deviation threshold
        'init_delay': 0.1,  # Initial delay threshold
        'peak_std_distance': 30.0,  # Distance in µm for peak time std calculation
        'remove_isolated': True,  # Remove isolated channels

        # Graph
        'init_amp_peak_ratio': 0.2,  # Weight for amplitude and peak latency in initial sorting
        'max_distance_for_edge': 100.0,  # Max distance in µm between channels for graph edges
        'max_distance_to_init': 200.0,  # Max distance to initial channel for graph edges
        'n_neighbors': 3,  # Max number of edges per channel in the graph
        'distance_exp': 2,  # Exponent for distance calculation
        'edge_dist_amp_ratio': 0.3,  # Weight between distance and amplitude for neighbor selection

        # Axonal reconstruction
        'min_path_length': 100.0,  # Minimum axon path length in µm
        'min_path_points': 5,  # Minimum number of channels in a path
        'neighbor_radius': 100.0,  # Radius in µm to exclude neighboring channels
        'min_points_after_branching': 3,  # Min number of points after a branching

        # Path cleaning/velocity estimation
        'mad_threshold': 8.0,  # MAD threshold for outlier detection
        'split_paths': True,  # Enable path splitting
        'max_peak_latency_for_splitting': 0.5,  # Max peak latency jump for splitting paths
        'r2_threshold': 0.9,  # R2 threshold for velocity fit
        'r2_threshold_for_outliers': 0.98,  # R2 threshold for outlier detection
        'min_outlier_tracking_error': 50.0,  # Min tracking error in µm for outlier detection
    },
    'save_reconstructor_object': True,
    'reconstructor_save_options': {
        'recordings': True, 
        'multirecordings': True, 
        'sortings': True,
        'waveforms': True,
        'templates': True
    },
    'reconstructor_load_options': {
        'load_reconstructor': True,
        'load_multirecs': True,
        'load_sortings': True,
        'load_wfs': True,
        'load_templates': True,
        'load_templates_bypass': False,
        'restore_environment': False,
    },
    'verbose': True,
    'debug_mode': True,
    'n_jobs': 1,
    'max_workers': 32,
    'logger_level': 'DEBUG',
    'run_lean': True,
}

h5_parent_dirs = ['/mnt/disk20tb/PrimaryNeuronData/Maxtwo/KCNT1_T3_C1_03122024/KCNT1_T3_C1_03122024/']

'''Run the pipeline normally'''
#reconstructor = AxonReconstructor(h5_parent_dirs, **kwargs)
# reconstructor.run_pipeline(**kwargs)

'''hacky way to run the pipeline on maxtwo, one well at a time to minimize size of temp data'''
h5_files = helper.get_list_of_h5_files(h5_parent_dirs, **kwargs)
for h5_file in h5_files:
    max_two_wells = 6 # 6 wells total
    max_two_wells_analyzed = 0
    while max_two_wells_analyzed < max_two_wells:
        kwargs['stream_select'] = max_two_wells_analyzed # should go 0 through 5
        reconstructor = AxonReconstructor([h5_file], **kwargs)
        reconstructor.run_pipeline(**kwargs)
        max_two_wells_analyzed += 1
