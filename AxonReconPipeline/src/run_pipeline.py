# run_pipeline.py
import lib_helper_functions as helper
from axon_reconstructor import AxonReconstructor

# Set up kwargs
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
        #'n_jobs': 32
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
        'load_reconstructor': False, # If True, will attempt to load reconstructor object from disk - all following options are relevant only if this is True
        'load_multirecs': True, # If True, will attempt to load multirecordings from if present in reconstructor object
        'load_sortings': True, # If True, will attempt to load recordings from if present in reconstructor object
        'load_wfs': True, # If True, will attempt to load waveforms from if present in reconstructor object
        'load_templates': True, # If True, will attempt to load templates from if present in reconstructor object
        'load_templates_bypass': False, # If True, will attempt to load templates from if present in reconstructor object without checking presence of waveforms, sortings, etc. 
                                        # Useful if sorting and waveform temp data have been delteted but templates are still available
                                        # Use with caution
        'restore_environment': False, # Potentially useful if loading reconstructor object in a different environment, untested.
    },
    'verbose': True,
    'debug_mode': True,
    'n_jobs': 4,
    'max_workers': 32,
    #'stream_select': 0,  # 0 for first stream, i.e., first well
    'run_lean': False, # If True, temp files will be deleted after each recording is processed
    #'logger_level': 'INFO',  # Specify the desired logger level here
    #'logger_level': 'WARNING',  # Specify the desired logger level here
    #'logger_level': 'ERROR',  # Specify the desired logger level here
    'logger_level': 'DEBUG',  # Specify the desired logger level here
    #'logger_level': 'CRITICAL',  # Specify the desired logger level here
    #'unit_limit': 10, # 1 for one unit per well, useful for testing
}

#h5_parent_dirs = ['/mnt/disk20tb/PrimaryNeuronData/Maxtwo/KCNT1_T3_C1_03122024/KCNT1_T3_C1_03122024/240322/M07037/AxonTracking/000025']
h5_parent_dirs = ['/mnt/disk20tb/PrimaryNeuronData/Maxtwo/KCNT1_T3_C1_03122024/KCNT1_T3_C1_03122024/']

'''Run the pipeline normally'''
#reconstructor = AxonReconstructor(h5_parent_dirs, **kwargs)
# reconstructor.run_pipeline(**kwargs)

'''hacky way to run the pipeline on maxtwo, one well at a time to minimze size of temp data'''
h5_files = helper.get_list_of_h5_files(h5_parent_dirs, **kwargs)
for h5_file in h5_files:
    max_two_wells = 6 # 6 wells total
    max_two_wells_analyzed = 0
    while max_two_wells_analyzed < max_two_wells:
        kwargs['stream_select'] = max_two_wells_analyzed #should go 0 through 5
        kwargs['stream_select'] = 1
        reconstructor = AxonReconstructor([h5_file], **kwargs)
        reconstructor.run_pipeline(**kwargs)
        max_two_wells_analyzed += 1