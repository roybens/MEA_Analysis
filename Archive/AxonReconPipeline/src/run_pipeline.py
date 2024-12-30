import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from MEAProcessingLibrary import mea_processing_library as MPL
import AxonReconPipeline.src.lib_helper_functions as helper
from AxonReconPipeline.src.axon_reconstructor import AxonReconstructor

# Set up kwargs with default reconstruction parameters from the paper
kwargs = {
    # runtime options
    'sorting_params': {
        'allowed_scan_types': ['AxonTracking'],
        'load_existing_sortings': True,
        'keep_good_only': False,
    },
    'te_params': {
        'load_merged_templates': False,
        'save_merged_templates': True,
        'time_derivative': True,
    },
    'av_params': {
        'upsample': 2,  # Upsampling factor to better capture finer details (unitless). Higher values increase temporal resolution.
        'min_selected_points': 20,  # Minimum number of selected points to include more potential signals (unitless). Determines the minimum dataset size for analysis.
        'verbose': False,  # Verbosity flag for debugging (boolean). If True, additional logs are printed.

        # Channel selection
        'detect_threshold': 0.01,  # Threshold for detecting signals, sensitive to smaller amplitudes (relative or absolute units, depending on 'detection_type'). Works with 'detection_type'.
        'detection_type': 'relative',  # Type of detection threshold ('relative' or 'absolute'). Determines if 'detect_threshold' is a relative or absolute value.
        'kurt_threshold': 0.2,  # Kurtosis threshold for noise filtering (unitless, kurtosis measure). Lower values allow more channels through noise filtering.
        'peak_std_threshold': 0.5,  # Peak time standard deviation threshold (milliseconds). Filters out channels with high variance in peak times.
        'init_delay': 0.05,  # Initial delay threshold to include faster signals (milliseconds). Minimum delay for considering a channel.
        'peak_std_distance': 20.0,  # Distance for calculating peak time standard deviation (micrometers). Defines the neighborhood for peak time calculation.
        'remove_isolated': True,  # Flag to remove isolated channels (boolean). If True, only connected channels are kept for analysis.

        # Graph
        'init_amp_peak_ratio': 0.2,  # Ratio between initial amplitude and peak time (unitless). Used to sort channels for path search.
        'max_distance_for_edge': 150.0,  # Maximum distance for creating graph edges (micrometers). Defines the maximum allowed distance between connected channels.
        'max_distance_to_init': 4000.0,  # Maximum distance to initial channel for creating edges (micrometers). Determines the initial connectivity range.
        'n_neighbors': 9,  # Maximum number of neighbors (edges) per channel (unitless). Enhances connectivity by increasing the number of edges.
        'distance_exp': 1.5,  # Exponent for distance calculation (unitless). Adjusts the influence of distance in edge creation.
        'edge_dist_amp_ratio': 0.2,  # Ratio between distance and amplitude for neighbor selection (unitless). Balances the importance of proximity and amplitude in selecting edges.

        # Axonal reconstruction
        'min_path_length': 80.0,  # Minimum path length to include shorter paths (micrometers). Ensures that only significant paths are considered.
        'min_path_points': 3,  # Minimum number of channels in a path (unitless). Defines the minimum size of a valid path.
        'neighbor_radius': 80.0,  # Radius to include neighboring channels (micrometers). Defines the search radius for neighbors.
        'min_points_after_branching': 2,  # Minimum points after branching to keep paths (unitless). Ensures that branches have enough data points.

        # Path cleaning/velocity estimation
        'mad_threshold': 10.0,  # Median Absolute Deviation threshold for path cleaning (unitless). Higher values allow more variability in the path.
        'split_paths': True,  # Flag to enable path splitting (boolean). If True, paths can be split for better velocity fit.
        'max_peak_latency_for_splitting': 0.7,  # Maximum peak latency jump for splitting paths (milliseconds). Allows capturing more variations by splitting paths at significant jumps.
        'r2_threshold': 0.75,  # R-squared threshold for velocity fit (unitless). Lower values include more paths with less perfect fits.
        'r2_threshold_for_outliers': 0.95,  # R-squared threshold for outlier detection (unitless). Defines the threshold below which the algorithm looks for outliers.
        'min_outlier_tracking_error': 40.0,  # Minimum tracking error to consider a channel as an outlier (micrometers). Sets the error tolerance for tracking outliers.
    },
    'analysis_options': {
        'generate_animation': True,
        'generate_summary': False,
    },
    'save_reconstructor_object': True,
    'reconstructor_save_options': {
        'recordings': True, 
        'multirecordings': True, 
        'sortings': True,
        'waveforms': False,
        'templates': True,
    },
    'reconstructor_load_options': {
        'load_reconstructor': True,
        
        #Only relevant if load_reconstructor is True:
        'load_multirecs': True,
        'load_sortings': True,
        'load_wfs': False,
        'load_templates': True,
        'load_templates_bypass': False,
        'restore_environment': False,
    },
    'verbose': True,
    'debug_mode': True,
    'n_jobs': 1,
    'max_workers': 8,
    'logger_level': 'DEBUG',
    'run_lean': True,
    
    #custom pathing
    'recon_dir':'./AxonReconPipeline/data/reconstructions_KCNT1',
    'reconstructor_dir' : './AxonReconPipeline/data/reconstructors_KCNT1',
    # 'log_file': 'axon_reconstruction.log',
    # 'error_log_file': 'axon_reconstruction_error.log',
}

h5_parent_dirs = [
    '/mnt/disk20tb/PrimaryNeuronData/Maxtwo/KCNT1_T4_C1_04122024/KCNT1_T4_C1_04122024/240416/M08034/AxonTracking/000015/data.raw.h5',
    '/mnt/disk20tb/PrimaryNeuronData/Maxtwo/KCNT1_T4_C1_04122024/KCNT1_T4_C1_04122024/240430/M08034/AxonTracking/000067/data.raw.h5',
    '/mnt/disk20tb/PrimaryNeuronData/Maxtwo/KCNT1_T4_C1_04122024/KCNT1_T4_C1_04122024/240426/M08034/AxonTracking/000055/data.raw.h5',
    '/mnt/disk20tb/PrimaryNeuronData/Maxtwo/KCNT1_T4_C1_04122024/KCNT1_T4_C1_04122024/240510/M08034/AxonTracking/000106/data.raw.h5',
    '/mnt/disk20tb/PrimaryNeuronData/Maxtwo/KCNT1_T4_C1_04122024/KCNT1_T4_C1_04122024/240514/M08034/AxonTracking/000120/data.raw.h5',
]

def run_pipeline(h5_parent_dirs, mode = 'normal', **kwargs):
    if mode == 'normal':
        '''Run the pipeline normally'''
        h5_files = helper.get_list_of_h5_files(h5_parent_dirs, **kwargs)
        reconstructor = AxonReconstructor(h5_files, **kwargs)
        reconstructor.run_pipeline(**kwargs)
    elif mode == 'lean':
        '''hacky way to run the pipeline on maxtwo, one well at a time to minimize size of temp data'''
        h5_files = helper.get_list_of_h5_files(h5_parent_dirs, **kwargs)
        for h5_file in h5_files:
            max_two_wells = 6 # 6 wells total
            max_two_wells_analyzed = 0
            while max_two_wells_analyzed < max_two_wells:
                
                #Get reconstructor ID
                h5_details = MPL.extract_recording_details([h5_file])
                date = h5_details[0]['date']
                chipID = h5_details[0]['chipID']
                runID = h5_details[0]['runID']        
                if 'stream_select' not in kwargs: kwargs['stream_select'] = max_two_wells_analyzed # should go 0 through 5
                else: # keep 'stream_select' in kwargs as it is
                    if max_two_wells_analyzed == kwargs['stream_select']: pass
                    else: max_two_wells_analyzed += 1; continue
                reconstructorID = f'{date}_{chipID}_{runID}_well00{kwargs["stream_select"]}'

                #Set log files to be unique for each reconstructor
                kwargs['log_file'] = f'{kwargs["reconstructor_dir"]}/{reconstructorID}_axon_reconstruction.log'
                kwargs['error_log_file'] = f'{kwargs["reconstructor_dir"]}/{reconstructorID}_axon_reconstruction_error.log'

                #Run the pipeline
                reconstructor = AxonReconstructor([h5_file], **kwargs)
                reconstructor.run_pipeline(**kwargs)
                max_two_wells_analyzed += 1