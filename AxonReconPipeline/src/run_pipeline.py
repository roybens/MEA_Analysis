# run_pipeline.py
from axon_reconstructor import AxonReconstructor

# Set up kwargs
kwargs = {
    # runtime options
    'sorting_params': {
        'allowed_scan_types': ['AxonTracking'],
        'load_existing_sortings': True,
        'keep_good_only': True,
    },
    'te_params': {
        'load_merged_templates': False,
        'save_merged_templates': True,
    },
    'verbose': True,
    'debug_mode': True,
    'n_jobs': 4,
    'max_workers': 32,
    #'stream_select': 0,  # 0 for first stream, i.e., first well
    'run_lean': True,
    #'logger_level': 'INFO',  # Specify the desired logger level here
    #'logger_level': 'WARNING',  # Specify the desired logger level here
    #'logger_level': 'ERROR',  # Specify the desired logger level here
    'logger_level': 'DEBUG',  # Specify the desired logger level here
    #'logger_level': 'CRITICAL',  # Specify the desired logger level here
    #'unit_limit': 10, # 1 for one unit per well, useful for testing
}


#h5_parent_dirs = ['/mnt/disk20tb/PrimaryNeuronData/Maxtwo/KCNT1_T3_C1_03122024/KCNT1_T3_C1_03122024/240322/M07037/AxonTracking/000025']
h5_parent_dirs = ['/mnt/disk20tb/PrimaryNeuronData/Maxtwo/KCNT1_T3_C1_03122024/KCNT1_T3_C1_03122024/']
reconstructor = AxonReconstructor(h5_parent_dirs, **kwargs)

# Run the pipeline
reconstructor.run_pipeline()