'''Setup Python environment '''
import setup_environment
import sys
import os
setup_environment.set_pythonpath()

''' Import the necessary modules '''
import argparse
from modules import mea_processing_library as MPL
from modules.axon_reconstructor import AxonReconstructor
from pprint import pprint

''' Helper functions '''

import h5py
import logging
from h5py.h5f import FileID

def log_hdf5_errors(file_path):
    try:
        with h5py.File(file_path, 'r') as f:
            fid = FileID(f.id.id)
            print(f"{file_path} opened successfully.")
            #check_fletcher32(file_path)
            print("HDF5 File ID:", fid)
    except Exception as e:
        print(f"Error opening {file_path}: {e}")
        #check_fletcher32(file_path)

def check_fletcher32(file_path):
    try:
        with h5py.File(file_path, 'r') as f:
            for name, dataset in f.items():
                if dataset.fletcher32:
                    print(f"Dataset {name} has checksum validation.")
    except Exception as e:
        print(f"Error checking Fletcher32 for {file_path}: {e}")

def parse_arguments():
    print("Parsing arguments...")
    print(sys.argv)
    parser = argparse.ArgumentParser(description="Run the axon reconstruction pipeline for a single well.")
    parser.add_argument("--plate_file", required=True, help="Path to the HDF5 file for the plate.")
    parser.add_argument("--output_dir", required=True, help="Output directory for reconstruction results.")
    parser.add_argument("--stream_select", type=int, required=True, help="Stream select for the well (0-5).")
    return parser.parse_args()

''' Main function '''
def main():
    args = parse_arguments()
    plate_file = args.plate_file
    output_dir = args.output_dir
    stream_select = args.stream_select

    ''' Parameters '''
    # Reconstructor parameters
    kwargs = {
        # runtime options
        'sorting_params': {
            'allowed_scan_types': ['AxonTracking'],
            'load_existing_sortings': True,
            'keep_good_only': False,
            #'use_gpu': False,
        },
        'te_params': {
            'load_merged_templates': True,
            'save_merged_templates': True,
            'time_derivative': True,
        },
        'av_params': {
            # 'upsample': 2,  # Upsampling factor to better capture finer details (unitless). Higher values increase temporal resolution.
            # 'min_selected_points': 20,  # Minimum number of selected points to include more potential signals (unitless). Determines the minimum dataset size for analysis.
            # 'verbose': False,  # Verbosity flag for debugging (boolean). If True, additional logs are printed.

            # # Channel selection
            # 'detect_threshold': 0.01,  # Threshold for detecting signals, sensitive to smaller amplitudes (relative or absolute units, depending on 'detection_type'). Works with 'detection_type'.
            # 'detection_type': 'relative',  # Type of detection threshold ('relative' or 'absolute'). Determines if 'detect_threshold' is a relative or absolute value.
            # 'kurt_threshold': 0.2,  # Kurtosis threshold for noise filtering (unitless, kurtosis measure). Lower values allow more channels through noise filtering.
            # 'peak_std_threshold': 0.5,  # Peak time standard deviation threshold (milliseconds). Filters out channels with high variance in peak times.
            # 'init_delay': 0.05,  # Initial delay threshold to include faster signals (milliseconds). Minimum delay for considering a channel.
            # 'peak_std_distance': 20.0,  # Distance for calculating peak time standard deviation (micrometers). Defines the neighborhood for peak time calculation.
            # 'remove_isolated': True,  # Flag to remove isolated channels (boolean). If True, only connected channels are kept for analysis.

            # # Graph
            # 'init_amp_peak_ratio': 0.2,  # Ratio between initial amplitude and peak time (unitless). Used to sort channels for path search.
            # 'max_distance_for_edge': 150.0,  # Maximum distance for creating graph edges (micrometers). Defines the maximum allowed distance between connected channels.
            # 'max_distance_to_init': 300.0,  # Maximum distance to initial channel for creating edges (micrometers). Determines the initial connectivity range.
            # #'max_distance_to_init': 4000.0,  # Maximum distance to initial channel for creating edges (micrometers). Determines the initial connectivity range.
            # 'n_neighbors': 9,  # Maximum number of neighbors (edges) per channel (unitless). Enhances connectivity by increasing the number of edges.
            # 'distance_exp': 1.5,  # Exponent for distance calculation (unitless). Adjusts the influence of distance in edge creation.
            # 'edge_dist_amp_ratio': 0.2,  # Ratio between distance and amplitude for neighbor selection (unitless). Balances the importance of proximity and amplitude in selecting edges.

            # # Axonal reconstruction
            # 'min_path_length': 80.0,  # Minimum path length to include shorter paths (micrometers). Ensures that only significant paths are considered.
            # 'min_path_points': 3,  # Minimum number of channels in a path (unitless). Defines the minimum size of a valid path.
            # 'neighbor_radius': 80.0,  # Radius to include neighboring channels (micrometers). Defines the search radius for neighbors.
            # 'min_points_after_branching': 2,  # Minimum points after branching to keep paths (unitless). Ensures that branches have enough data points.

            # # Path cleaning/velocity estimation
            # 'mad_threshold': 10.0,  # Median Absolute Deviation threshold for path cleaning (unitless). Higher values allow more variability in the path.
            # 'split_paths': True,  # Flag to enable path splitting (boolean). If True, paths can be split for better velocity fit.
            # 'max_peak_latency_for_splitting': 0.7,  # Maximum peak latency jump for splitting paths (milliseconds). Allows capturing more variations by splitting paths at significant jumps.
            # 'r2_threshold': 0.75,  # R-squared threshold for velocity fit (unitless). Lower values include more paths with less perfect fits.
            # 'r2_threshold_for_outliers': 0.95,  # R-squared threshold for outlier detection (unitless). Defines the threshold below which the algorithm looks for outliers.
            # 'min_outlier_tracking_error': 40.0,  # Minimum tracking error to consider a channel as an outlier (micrometers). Sets the error tolerance for tracking outliers.
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
            'load_templates': False,
            'load_templates_bypass': False,  # This is a new parameter that allows the user to bypass pre-processing steps and load the templates directly. 
                                            # Useful if there isn't any need to reprocess the templates.
            'restore_environment': False,
        },
        'verbose': True,
        'debug_mode': True,
        #'n_jobs': 64,
        'n_jobs': 128,
        #'max_workers': 64,
        'logger_level': 'DEBUG',
        'run_lean': True,
    }

    kwargs['project_name'] = None # Use project name to create subdirectories, if true the paths below can be relative
    kwargs['output_dir'] = output_dir # Output directory for the reconstructions when running on NERSC
    kwargs['mode'] = 'lean'
    kwargs['stream_select'] = stream_select

    # Process a single well file
    print(f"Processing plate {plate_file} stream {stream_select}...")
    h5_files = [plate_file]
    
    h5_details = MPL.extract_recording_details([plate_file])
    date = h5_details[0]['date']
    chipID = h5_details[0]['chipID']
    runID = h5_details[0]['runID']
    scanType = h5_details[0]['scanType']
    
    projectName = h5_details[0]['projectName']
    reconstructorID = f'{date}_{chipID}_{runID}_well00{stream_select}'

    wellID = f'well00{stream_select}'
    well_path = f'{kwargs["output_dir"]}/{projectName}/{date}/{chipID}/{scanType}/{runID}/{wellID}'
    kwargs['log_file'] = f'{well_path}/{wellID}_axon_reconstruction.log'
    kwargs['error_log_file'] = f'{well_path}/{wellID}_axon_reconstruction_error.log'
    kwargs['recordings_dir'] = os.path.join(well_path, 'recordings')
    kwargs['sortings_dir'] = os.path.join(well_path, 'sortings')
    kwargs['waveforms_dir'] = os.path.join(well_path, 'waveforms')
    kwargs['templates_dir'] = os.path.join(well_path, 'templates')
    kwargs['recon_dir'] = os.path.join(well_path, 'reconstructions')
    kwargs['reconstructor_dir'] = well_path
    
    print(f'Starting reconstruction for {reconstructorID}...')
    print(f'Log file: {kwargs["log_file"]}')
    print(f'Error log file: {kwargs["error_log_file"]}')
    print(f'Recordings dir: {kwargs["recordings_dir"]}')
    print(f'Sortings dir: {kwargs["sortings_dir"]}')
    print(f'Waveforms dir: {kwargs["waveforms_dir"]}')
    print(f'Templates dir: {kwargs["templates_dir"]}')
    print(f'Reconstructions dir: {kwargs["recon_dir"]}')
    print(f'Reconstructor dir: {kwargs["reconstructor_dir"]}')
    
    # # Run the pipeline
    # print("Reconstructor Object:")
    # pprint(vars(reconstructor))

    # print("\nReconstructor Parameters:")
    # pprint(kwargs)
    # import h5py
    # with h5py.File(plate_file, 'r') as f:
    #     print(list(f.keys()))
    
    #kwargs switches, only do spike sorting for now
    kwargs['concatenate_switch'] = True
    kwargs['sort_switch'] = True
    kwargs['waveforms_switch'] = True
    kwargs['template_switch'] = True
    kwargs['recon_switch'] = True
    
    #log the HDF5 file
    log_hdf5_errors(plate_file)
    
    try:
        reconstructor = AxonReconstructor(h5_files, **kwargs)
        reconstructor.run_pipeline(**kwargs)
    except Exception as e:
        print(f"Error processing plate {plate_file} stream {stream_select}: {e}")

#Specify arguments for testing
# sys.argv = [
#     'run_pipeline_HPC.py',  # Script name
#     '--plate_file', 
#     #'/pscratch/sd/a/adammwea/RBS_synology_rsync/B6J_DensityTest_10012024_AR/B6J_DensityTest_10012024_AR/B6J_DensityTest_10012024_AR/241021/M08029/AxonTracking/000097/data.raw.h5',
#     #'/pscratch/sd/a/adammwea/RBS_synology_rsync/B6J_DensityTest_10012024_AR/B6J_DensityTest_10012024_AR/B6J_DensityTest_10012024_AR/241021/M06804/AxonTracking/000091/data.raw.h5',
#     '/pscratch/sd/a/adammwea/workspace/xInputs/xRBS_input_data/B6J_DensityTest_10012024_AR/241021/M06804/AxonTracking/000091/data.raw.h5',
#     '--output_dir', 
#     '/pscratch/sd/a/adammwea/yThroughput/zRBS_axon_reconstruction_output', 
#     '--stream_select', '0'
# ]

if __name__ == "__main__":
    main()
    
    #to debug in login node, init shifter container with the following command:
    #shifter --image adammwea/axonkilo_docker:v7 /bin/bash
    
    #to run spikesorting as needed in interactive node with gpu:
    '''
    salloc -A m2043_g -q interactive -C gpu -t 04:00:00 --nodes=4 --gpus=4 --image=adammwea/axonkilo_docker:v7
    '''
    #not sure why I cant specify tasks-per-node, in the previous line