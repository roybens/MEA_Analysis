import os
import shutil
import logging
import spikeinterface.sorters as ss
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import AxonReconPipeline.src.lib_sorting_functions as sorter
import AxonReconPipeline.src.lib_waveform_functions as waveformer
import AxonReconPipeline.src.lib_template_functions as templater
from AxonReconPipeline.src.func_analyze_and_reconstruct import analyze_and_reconstruct
import spikeinterface.full as si
import subprocess
import json
import numpy as np
import dill
import h5py
from MEAProcessingLibrary import mea_processing_library as MPL

class SingleLevelFilter(logging.Filter):
    def __init__(self, level):
        super().__init__()
        self.level = level

    def filter(self, record):
        return record.levelno == self.level

class AxonReconstructor:
    def __init__(self, h5_parent_dirs, **kwargs):
        self.log_file = kwargs.get('log_file', 'axon_reconstruction.log')
        self.error_log_file = kwargs.get('error_log_file', 'axon_reconstruction_error.log')
        self.logger_level = kwargs.get('logger_level', 'INFO')
        self.logger = self.setup_logger()
        self.logger.info("Initializing AxonReconstructor")

        self.h5_parent_dirs = h5_parent_dirs
        self.allowed_scan_types = kwargs.get('allowed_scan_types', ['AxonTracking'])
        self.stream_select = kwargs.get('stream_select', None)
        self.unit_select = kwargs.get('unit_select', None)
        self.unit_limit = kwargs.get('unit_limit', None)
        self.save_reconstructor_object = kwargs.get('save_reconstructor_object', False)

        self.debug_mode = kwargs.get('debug_mode', True)
        self.load_existing_sortings = kwargs.get('load_existing_sortings', True)
        self.load_merged_templates = kwargs.get('load_merged_templates', False)
        self.load_template_segments = kwargs.get('load_template_segments', True)
        self.load_waveforms = kwargs.get('load_waveforms', True)
        self.run_lean = kwargs.get('run_lean', False)
        self.verbose = kwargs.get('verbose', True)

        self.recordings_dir = kwargs.get('recordings_dir', './AxonReconPipeline/data/temp_data/recordings')
        self.sortings_dir = kwargs.get('sortings_dir', './AxonReconPipeline/data/temp_data/sortings')
        self.waveforms_dir = kwargs.get('waveforms_dir', './AxonReconPipeline/data/temp_data/waveforms')
        self.templates_dir = kwargs.get('templates_dir', './AxonReconPipeline/data/temp_data/templates')
        self.recon_dir = kwargs.get('recon_dir', './AxonReconPipeline/data/reconstructions')
        self.reconstructor_dir = kwargs.get('reconstructor_dir', './AxonReconPipeline/data/reconstructors')

        self.n_jobs = kwargs.get('n_jobs', 1)
        self.max_workers = kwargs.get('max_workers', 16)
        assert self.max_workers % self.n_jobs == 0, "max_workers must be a multiple of n_jobs to avoid oversubscription."

        self.sorting_params = self.get_sorting_params(kwargs)
        self.te_params = self.get_te_params(kwargs)
        self.qc_params = self.get_qc_params(kwargs)
        self.av_params = self.get_av_params(kwargs)
        self.analysis_options = self.get_analysis_options(kwargs)
        self.reconstructor_save_options, self.reconstructor_load_options = self.get_reconstructor_options(kwargs)

        self.continuous_h5_file_info = None
        self.sorting_lists_by_scan = None
        self.waveforms = None
        self.templates = None
        self.reconstruction_results = None
        self.requirements = None
        self.git_version = None

        self.logger.debug("AxonReconstructor initialized with parameters: %s", self.__dict__)

    def setup_logger(self, prefix=None):
        logger_level = self.logger_level
        log_file = self.log_file
        error_log_file = self.error_log_file
        logger = logging.getLogger('axon_reconstructor')

        # Remove all existing handlers
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)

        # Clear the log files by opening them in write mode
        reconstructor_path = os.path.dirname(log_file)
        if prefix is not None:
            if not os.path.exists(reconstructor_path): os.makedirs(reconstructor_path)
            if os.path.exists(log_file): open(log_file, 'w').close()
            if os.path.exists(error_log_file): open(error_log_file, 'w').close()

        if prefix is None: formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(funcName)s - %(message)s')
        else: formatter = logging.Formatter(f'%(asctime)s - {prefix} - %(name)s - %(levelname)s - %(funcName)s - %(message)s')
        
        #logger.warning(os.environ['HDF5_FILE_PATH'])
        logger.setLevel(logger_level.upper())

        # Handlers for different log levels
        debug_handler = logging.FileHandler(error_log_file)
        debug_handler.setLevel(logging.DEBUG)
        debug_handler.addFilter(SingleLevelFilter(logging.DEBUG))
        debug_handler.setFormatter(formatter)
        logger.addHandler(debug_handler)

        info_handler = logging.FileHandler(log_file)
        info_handler.setLevel(logging.INFO)
        info_handler.addFilter(SingleLevelFilter(logging.INFO))
        info_handler.setFormatter(formatter)
        logger.addHandler(info_handler)

        warning_handler = logging.FileHandler(log_file)
        warning_handler.setLevel(logging.WARNING)
        warning_handler.addFilter(SingleLevelFilter(logging.WARNING))
        warning_handler.setFormatter(formatter)
        logger.addHandler(warning_handler)

        error_handler = logging.FileHandler(error_log_file)
        error_handler.setLevel(logging.ERROR)
        error_handler.addFilter(SingleLevelFilter(logging.ERROR))
        error_handler.setFormatter(formatter)
        logger.addHandler(error_handler)

        critical_handler = logging.FileHandler(error_log_file)
        critical_handler.setLevel(logging.CRITICAL)
        critical_handler.addFilter(SingleLevelFilter(logging.CRITICAL))
        critical_handler.setFormatter(formatter)
        logger.addHandler(critical_handler)

        # Stream handler for console output
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logger_level.upper())
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

        if prefix is None:
            logger.info("Logger initialized with level: %s", logger_level)
            logger.debug('This is a debug message')
            logger.info('This is an info message')
            logger.warning('This is a warning message')
            logger.error('This is an error message')
            logger.critical('This is a critical message')

        return logger

    def get_sorting_params(self, kwargs):
        default_params = ss.Kilosort2Sorter.default_params()
        default_params.update({
            'allowed_scan_types': self.allowed_scan_types,
            'stream_select': self.stream_select,
            'load_existing_sortings': self.load_existing_sortings,
            'sortings_dir': self.sortings_dir,
            'overwrite': False,
            'verbose': self.verbose,
            'n_jobs': self.n_jobs,
            'detect_threshold': 7,
            'minFR': 0.01,
            'minfr_goodchannels': 0.01,
            'keep_good_only': False,
            'do_correction': False,
            'chunk_duration': "1s",
            'progress_bar': True,
            'use_docker': True,
        })
        default_params.update(kwargs.get('sorting_params', {}))
        for key, value in default_params.items():
            for k, v in kwargs.get('sorting_params', {}).items():
                if k == key: default_params[key] = v
        return default_params

    def get_te_params(self, kwargs):
        default_params = {
            'align_cutout': True,
            'upsample': 2,
            'rm_outliers': True,
            'n_jobs': self.max_workers,
            'n_neighbors': 10,
            'peak_cutout': 2,
            'overwrite_wf': False,
            'overwrite_tmp': False,
            'load_merged_templates': True,
            'save_merged_templates': True,
            'time_derivative': False,
        }
        default_params.update(kwargs.get('te_params', {}))
        return default_params

    def get_qc_params(self, kwargs):
        default_params = {
            'min_n_spikes': 500,
            'exclude_mua': True
        }
        default_params.update(kwargs.get('qc_params', {}))
        for key, value in default_params.items():
            for k, v in kwargs.get('qc_params', {}).items():
                if k == key: default_params[key] = v
        return default_params

    def get_av_params(self, kwargs):
        av_params = {
            'upsample': 1,  # Upsampling factor for template
            'min_selected_points': 30,  # Minimum number of selected points to run axon tracking
            'verbose': False,  # If True, the output is verbose

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
        }
        av_params.update(kwargs.get('av_params', {}))
        return av_params
    
    def get_analysis_options(self, kwargs):
        analysis_options = {
            'generate_animation': False
        }
        analysis_options.update(kwargs.get('analysis_options', {}))
        return analysis_options

    def get_reconstructor_options(self, kwargs):
        reconstructor_save_options = {
            'recordings': False,
            'multirecordings': False,
            'sortings': False,
            'waveforms': False,
            'templates': True,
            'reconstructions': False
        }
        reconstructor_save_options.update(kwargs.get('reconstructor_save_options', {}))
        reconstructor_load_options = {
            'load_reconstructor': True,
            'load_multirecs': True,
            'load_sortings': True,
            'load_wfs': True,
            'load_templates': True,
            'load_templates_bypass': False,
            'restore_environment': False
        }
        reconstructor_load_options.update(kwargs.get('reconstructor_load_options', {}))
        return reconstructor_save_options, reconstructor_load_options

    def load_recordings(self):
        self.logger.info("Loading continuous h5 files")
        
        recordings = {}
        self.reconstructor_id = None        
        for dir in self.h5_parent_dirs:
            h5_file_paths = MPL.extract_raw_h5_filepaths(dir)
            for h5_path in h5_file_paths:
                h5_details = MPL.extract_recording_details(h5_path)
                date = h5_details[0]['date']
                chip_id = h5_details[0]['chipID']
                run_id = h5_details[0]['runID']
                scan_type = h5_details[0]['scanType']
                if scan_type in self.sorting_params['allowed_scan_types']:
                    try: device, recording_segments, stream_count, rec_counts = MPL.load_recordings(h5_path, stream_select=self.stream_select, logger=self.logger)
                    except Exception as e: 
                        self.logger.error(f"Error loading recording segments from {h5_path}: {e}")
                        self.logger.error(f"Skipping {h5_path}")
                        continue
                    if self.reconstructor_id is None: self.reconstructor_id = f"{date}_{chip_id}_{run_id}" #TODO: This is a hacky way to set the reconstructor_id, need to fix such that there is one reconstructor per recording
                    recordings[f"{date}_{chip_id}_{run_id}"] = {
                        'h5_path': h5_path,
                        'scanType': h5_details[0]['scanType'],
                        'date': date,
                        'chip_id': chip_id,
                        'run_id': run_id,
                        'device': device,
                        'streams': recording_segments,
                        'stream_count': stream_count,
                        'rec_counts': rec_counts
                }

        self.recordings = recordings
        self.logger.info("Completed loading .h5 files")

    def _validate_multirec(self, rec_key, stream_id, recording_segments):
        assert self.reconstructor_load_options['load_multirecs'], 'Load multirecs option is set to False. Generating new multirecording.'
        self.logger.info(f'Attempting to load multirec {rec_key} stream {stream_id}')
        assert 'multirecordings' in self.__dict__, "No multirecordings found in reconstructor. Generating new multirecording."
        assert self.multirecordings is not None, "No multirecordings found in reconstructor. Generating new multirecording."
        assert stream_id in self.multirecordings[rec_key]['streams'], f"No multirecordings found in reconstructor for {stream_id}. Generating new multirecording."
        assert self.multirecordings[rec_key]['streams'][stream_id]['multirecording'], f"No multirecording found in reconstructor for {stream_id}. Generating new multirecording."
        assert len(self.multirecordings[rec_key]['streams'][stream_id]['common_el'])>0, f"Empty common_el found in reconstructor for {stream_id}. Generating new multirecording."
        self.logger.info(f"Success! Using exisiting multirecording from reconstructor object for {stream_id}. Skipping concatenation.")
        return self.multirecordings[rec_key]['streams'][stream_id]
    
    def concatenate_recordings(self):        
        self.logger.info("Concatenating recordings")
        multirecordings = {}
        try: assert self.recordings, "No recordings found. Skipping concatenation."
        except Exception as e: self.logger.error(e); return
        for rec_key, recording in self.recordings.items():
            
            date = recording['date']
            chip_id = recording['chip_id']
            scanType = recording['scanType']
            run_id = recording['run_id']
            h5_path = recording['h5_path']
            h5 = h5py.File(h5_path)
            stream_ids = list(h5['wells'].keys())
            if self.stream_select is not None: stream_ids = [stream_ids[self.stream_select]] if isinstance(stream_ids[self.stream_select], str) else stream_ids[self.stream_select]
            streams = {}
            for stream_id in stream_ids:
                self.setup_logger(prefix=f'{rec_key}_{stream_id}')
                #test logger message
                self.logger.info(f"Concatenating recording segments for {date}_{chip_id}_{run_id} stream {stream_id}")
                recording_segments = recording['streams'][stream_id]['recording_segments']

                # Check if multirecording already exists in reconstructor object
                try: 
                    assert self.reconstructor_load_options['load_reconstructor'], 'Load reconstructor option from reconstructor is set to False. Concatenating new multirecording.'
                    multirecording_dict = self._validate_multirec(rec_key, stream_id, recording_segments)
                    streams[stream_id] = multirecording_dict
                # If not, generate new multirecording
                except AssertionError as e: 
                    self.logger.warning(e)
                    self.logger.info(f'Concatenating recording segments for {date}_{chip_id}_{run_id} stream {stream_id}')
                    multirecording, common_el, multirec_save_path = sorter.concatenate_recording_segments(
                        h5_path, recording_segments, stream_id, save_dir=self.recordings_dir, logger=self.logger)                     
                    streams[stream_id] = {
                        'multirecording': multirecording,
                        'common_el': common_el,
                        'multirec_save_path': multirec_save_path
                    } 

            multirecordings[f"{date}_{chip_id}_{run_id}"] = {
                'h5_path': h5_path, 
                'scanType': scanType,
                'date': date,
                'chip_id': chip_id,
                'run_id': run_id,
                'streams': streams
            }
        try: 
            self.update_nested_dict(self.multirecordings, multirecordings)
        except: self.multirecordings = multirecordings
        self.logger.info("Completed concatenation of recordings")
    
    def _validate_sort(self, rec_key, stream_id, streams):
        assert self.reconstructor_load_options['load_sortings'], 'Load existing sortings option is set to False. Generating new sorting.'
        self.logger.info(f'Attempting to load sorting {rec_key} stream {stream_id}')
        assert 'sortings' in self.__dict__, "No sortings found in reconstructor. Generating new sorting."
        assert self.sortings is not None, "No sortings found in reconstructor. Generating new sorting."
        assert stream_id in self.sortings[rec_key]['streams'], f"No sortings found in reconstructor for {stream_id}. Generating new sorting."
        assert self.sortings[rec_key]['streams'][stream_id]['sorting'], f"No sorting found in reconstructor for {stream_id}. Generating new sorting."
        assert os.path.exists(self.sortings[rec_key]['streams'][stream_id]['sorting_path']+'/sorter_output'), f"Sorting path not found for {stream_id}. Generating new sorting."
        self.logger.info(f"Success! Using exisiting sorting from reconstructor object for {stream_id}. Skipping spike sorting.")
        return self.sortings[rec_key]['streams'][stream_id]
    
    def _check_if_sorting_failed(self, rec_key, stream_id, streams):
        try: message = self.sortings[rec_key]['streams'][stream_id]['sorting']['message']
        except: return False # If no message found, older version of reconstructor object - try to load or generate new sorting as needed
        if 'Spike sorting' in message and 'failed' in message: return True
        else: return False
    
    def spikesort_recordings(self):
        self.logger.info("Starting spike sorting process")
        sortings = {}
        try: assert self.multirecordings, "No multirecordings found. Skipping spike sorting." 
        except Exception as e: self.logger.error(e); return 
        for rec_key, multirec in self.multirecordings.items():
            date = multirec['date']
            chip_id = multirec['chip_id']
            scanType = multirec['scanType']
            run_id = multirec['run_id']
            self.logger.info(f'Generating spike sorting objects for {date}_{chip_id}_{run_id}')
            spikesorting_root = os.path.join(self.sorting_params['sortings_dir'], f'{date}/{chip_id}/{scanType}/{run_id}')
            streams = {}
            for stream_id, stream in multirec['streams'].items():
                self.setup_logger(prefix=f'{rec_key}_{stream_id}')
                #test logger message
                self.logger.info(f"Spike sorting stream {stream_id}")
                if self.stream_select is not None and stream_id != 'well{:03}'.format(self.stream_select): continue

                # Check if sorting already exists in reconstructor object
                try:
                    assert self.reconstructor_load_options['load_reconstructor'], 'Load existing sortings option from reconstructor is set to False. Generating new sorting.'
                    sorting_dict = self._validate_sort(rec_key, stream_id, multirec['streams'])
                    if self._check_if_sorting_failed(rec_key, stream_id, multirec['streams']): 
                        self.logger.warning(f"Sorting previously failed for {rec_key} stream {stream_id}. Presumably it will fail again. Skipping well.")
                        continue
                    streams[stream_id] = sorting_dict
                # If not, generate new sorting
                except AssertionError as e:
                    self.logger.warning(e)
                    self.logger.info(f'Spike sorting stream {stream_id}')
                    mr = stream['multirecording']                    
                    sorting, stream_sort_path, message = sorter.sort_multirecording(mr, stream_id, save_root=spikesorting_root, sorting_params=self.sorting_params, logger=self.logger)
                    streams[stream_id] = {
                        'sorting_path': stream_sort_path,
                        'sorting': sorting,
                        'message': message
                    }

            sortings[f"{date}_{chip_id}_{run_id}"] = {
                'scanType': scanType,
                'date': date,
                'chip_id': chip_id,
                'run_id': run_id,
                'spiksorting_root': spikesorting_root,
                'streams': streams
            }
        try: 
            self.update_nested_dict(self.sortings, sortings)
        except: self.sortings = sortings
        self.logger.debug("Completed spike sorting")

    def _validate_wf(self, rec_key, stream_id, sorting):
        assert self.reconstructor_load_options['load_wfs'], 'Load existing waveforms option is set to False. Generating new waveforms.'
        self.logger.info(f'Attempting to load waveforms {rec_key} stream {stream_id}')
        assert 'waveforms' in self.__dict__, "No waveforms found in reconstructor. Generating new waveforms."
        assert self.waveforms is not None, "No waveforms found in reconstructor. Generating new waveforms."
        assert stream_id in self.waveforms[rec_key]['streams'], f"No waveforms found in reconstructor for {stream_id}. Generating new waveforms."
        assert self.waveforms[rec_key]['streams'][stream_id], f"No waveforms found in reconstructor for {stream_id}. Generating new waveforms."
        for rec_name, waveform_seg in self.waveforms[rec_key]['streams'][stream_id].items():
            assert os.path.exists(self.waveforms[rec_key]['streams'][stream_id][rec_name]['path']), f"Waveform path not found for {stream_id}, segment {rec_name}. Generating new waveforms as needed."
        self.logger.info(f"Success! Using exisiting waveforms from reconstructor object for {stream_id}. Skipping waveform extraction.")
        return self.waveforms[rec_key]['streams'][stream_id]
    
    def extract_waveforms(self):    
        waveforms = {}
        self.logger.info("Extracting waveforms")
        try: assert self.sortings, "No sortings found. Skipping waveform extraction."
        except Exception as e: self.logger.error(e); return
        for key, sorting in self.sortings.items():
            
            multirecs = self.multirecordings[key]
            streams = {}
            for stream_id, multirec in multirecs['streams'].items():
                self.setup_logger(prefix=f'{key}_{stream_id}')
                #test logger message
                self.logger.info(f"Extracting waveforms from stream {stream_id}")
                # Check if waveforms already exist in reconstructor object
                try:
                    assert self.reconstructor_load_options['load_reconstructor'], 'Load existing waveforms option is set to False. Generating new waveforms.'
                    stream_wfs = self._validate_wf(key, stream_id, sorting)
                    streams[stream_id] = stream_wfs
                except:
                    mr = multirec['multirecording']
                    sort = sorting['streams'][stream_id]['sorting']
                    if sort is None: 
                        self.logger.error(f"Sorting object not found for {key} stream {stream_id}. Skipping waveform extraction.")
                        continue
                    cleaned_sorting = waveformer.select_units(sort, **self.qc_params)
                    cleaned_sorting = si.remove_excess_spikes(cleaned_sorting, mr) # Relevant if last spike time == recording_length
                    cleaned_sorting.register_recording(mr)
                    segment_sorting = si.SplitSegmentSorting(cleaned_sorting, mr)
                    self.logger.info(f'Extracting waveforms from stream: {stream_id}')
                    h5_path = self.recordings[key]['h5_path']
                    wf_kwargs = {
                        'save_root': self.waveforms_dir,
                        'logger': self.logger,
                        'te_params': self.te_params,
                        'unit_limit': self.unit_limit,
                        'unit_select': self.unit_select,
                    }
                    stream_wfs = waveformer.extract_unit_waveforms(h5_path, stream_id, segment_sorting, **wf_kwargs) # TODO: make distinct wf_params
                    streams[stream_id] = stream_wfs

                if self.stream_select is not None:
                    formatted = 'well{:03}'.format(self.stream_select)
                    if stream_id == formatted:
                        break
            waveforms[key] = {
                'h5_path': multirecs['h5_path'],
                'scanType': multirecs['scanType'],
                'date': multirecs['date'],
                'chip_id': multirecs['chip_id'],
                'run_id': multirecs['run_id'],
                'waveforms_dir': self.waveforms_dir,
                'streams': streams
            }
        try: 
            self.update_nested_dict(self.waveforms, waveforms)
        except: self.waveforms = waveforms
        self.logger.info("Completed waveform extraction")

    def _validate_templates(self, rec_key, stream_id, datum):
        assert self.reconstructor_load_options['load_templates'], 'Load existing templates option is set to False. Generating new templates.'
        self.logger.info(f'Attempting to load templates {rec_key} stream {stream_id}')
        assert 'templates' in self.__dict__, "No templates found in reconstructor. Generating new templates."
        assert self.templates is not None, "No templates found in reconstructor. Generating new templates."
        assert stream_id in self.templates[rec_key]['streams'], f"No templates found in reconstructor for {stream_id}. Generating new templates."
        assert self.templates[rec_key]['streams'][stream_id], f"No templates found in reconstructor for {stream_id}. Generating new templates."
        unit_list = [unit for unit in self.templates[rec_key]['streams'][stream_id]['units'].keys()]
        assert len(unit_list) > 0, f"No units found in reconstructor for {stream_id}. Generating new templates."
        for unit, unit_template in self.templates[rec_key]['streams'][stream_id]['units'].items():
            assert isinstance(unit_template['merged_template'], np.ndarray), f"'merged_template' not found for stream_id {stream_id}, unit {unit}. Generating new templates as needed."
            assert isinstance(unit_template['merged_channel_loc'], np.ndarray), f"'merged_channel_loc' not found for stream_id {stream_id}, unit {unit}. Generating new templates as needed."
            assert isinstance(unit_template['dvdt_merged_template'], np.ndarray), f"'dvdt_merged_template' not found for stream_id {stream_id}, unit {unit}. Generating new templates as needed."
            assert isinstance(unit_template['merged_template_filled'], np.ndarray), f"'merged_template_filled' not found for stream_id {stream_id}, unit {unit}. Generating new templates as needed."
            assert isinstance(unit_template['dvdt_merged_template_filled'], np.ndarray), f"'dvdt_merged_template_filled' not found for stream_id {stream_id}, unit {unit}. Generating new templates as needed."
            assert isinstance(unit_template['merged_channel_locs_filled'], np.ndarray), f"'merged_channel_locs_filled' not found for stream_id {stream_id}, unit {unit}. Generating new templates as needed."
        self.logger.info(f"Success! Using exisiting templates from reconstructor object for {stream_id}. Skipping template extraction.")
        return self.templates[rec_key]['streams'][stream_id]
    
    def extract_templates(self, template_bypass=False):
        self.logger.info("Extracting templates")
        templates = {}
        # if template_bypass is False:
        try: 
            assert self.waveforms, "No waveforms found."
            data = self.waveforms
        except Exception as e: 
            self.logger.error(e)
            if template_bypass: 
                self.logger.warning("Template bypass is enabled. Will attempt to validate templates in loaded reconstructor without waveforms.")      
                data = self.recordings                 
                pass
            else: 
                logging.error("No waveforms found. Skipping template extraction.")
                return

        for key, datum in data.items():
            streams = {}
            for stream_id, wfs in datum['streams'].items():
                self.setup_logger(prefix=f'{key}_{stream_id}')
                #test logger message
                self.logger.info(f"Extracting templates from stream {stream_id}")
                if self.stream_select is not None and stream_id != 'well{:03}'.format(self.stream_select): continue

                # Check if templates already exist in reconstructor object
                try:
                    assert self.reconstructor_load_options['load_reconstructor'], 'Load existing templates option from reconstructor is set to False. Generating new templates.'
                    stream_templates = self._validate_templates(key, stream_id, datum)
                    streams[stream_id] = stream_templates
                except Exception as e:
                    if template_bypass:
                        self.logger.error(f"Error encountered during template validation: {e}")
                        self.logger.error(f"Error will be raised and pipeline will return to non-bypass mode.")
                        raise e
                    self.logger.warning(e)
                    self.logger.info(f'Extracting templates from stream: {stream_id}')
                    h5_path = self.recordings[key]['h5_path']
                    sorting = self.sortings[key]['streams'][stream_id]['sorting']
                    multirec = self.multirecordings[key]['streams'][stream_id]['multirecording']
                    temp_kwargs = {
                        'unit_select': self.unit_select,
                    }
                    unit_templates = templater.extract_templates(
                        multirec, sorting, wfs, h5_path, stream_id, save_root=self.templates_dir, 
                        te_params=self.te_params, qc_params=self.qc_params, unit_limit=self.unit_limit, 
                        logger=self.logger, template_bypass=template_bypass, **temp_kwargs)
                    streams[stream_id] = {'units': unit_templates}
            templates[key] = {
                'h5_path': datum['h5_path'],
                'scanType': datum['scanType'],
                'date': datum['date'],
                'chip_id': datum['chip_id'],
                'run_id': datum['run_id'],
                'templates_dir': self.templates_dir,
                'streams': streams
            }
        try: 
            self.update_nested_dict(self.templates, templates)
        except: self.templates = templates
        self.logger.debug("Completed template extraction")

    def analyze_and_reconstruct(self, load_existing_templates=False):
        self.logger.info("Reconstructing Axonal Morphology")
        try: assert self.templates, "No templates found. Skipping analysis and reconstruction."
        except AssertionError as e: self.logger.error(e); return
        recon_kwargs = {
            'unit_select': self.unit_select,
        }
        analyze_and_reconstruct(
            self.templates,
            recon_dir=self.recon_dir,
            params=self.av_params,
            analysis_options=self.analysis_options,
            stream_select=self.stream_select,
            n_jobs=self.max_workers,
            logger=self.logger,
            **recon_kwargs,
        )
        self.logger.debug("Completed analysis and reconstruction")

    def save_reconstructor(self):
        self.logger.info("Saving reconstructor object")
        for rec_key, recording in self.recordings.items():
            h5_path = recording['h5_path']
            h5 = h5py.File(h5_path)
            stream_ids = list(h5['wells'].keys())
            if self.stream_select is not None: stream_ids = [stream_ids[self.stream_select]] if isinstance(stream_ids[self.stream_select], str) else stream_ids[self.stream_select]
            for stream_id in stream_ids:
                # Ensure the directory exists
                if not os.path.exists(self.reconstructor_dir):
                    os.makedirs(self.reconstructor_dir)

                self.logger.info(f"Saving reconstructor object stream: {stream_id}")

                # Capture pip requirements
                self.requirements = subprocess.check_output(['pip', 'freeze']).decode().splitlines()

                # Capture git version
                self.git_version = subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode().strip()

                # Optionally remove excess data
                if not self.reconstructor_save_options['recordings']: del self.recordings
                if not self.reconstructor_save_options['multirecordings']: del self.multirecordings
                if not self.reconstructor_save_options['sortings']: del self.sortings
                if not self.reconstructor_save_options['waveforms']: del self.waveforms
                if not self.reconstructor_save_options['templates']: del self.templates

                # Save the entire object using dill        
                dill_file_path = os.path.join(self.reconstructor_dir, f'{rec_key}_{stream_id}.dill')
                with open(dill_file_path, 'wb') as f:
                    dill.dump(self, f)

                self.logger.info(f"Reconstructor object saved to {dill_file_path}")

    def update_nested_dict(self, target, source):
        """
        Recursively updates a nested dictionary `target` with values from `source`.
        """
        for key, value in source.items():
            if isinstance(value, dict):
                if target.get(key) is None:
                    target[key] = {}
                elif not isinstance(target[key], dict):
                    raise AttributeError(f"Expected dict at target[{key}], but got {type(target[key])}")
                self.update_nested_dict(target[key], value)
            else:
                target[key] = value
    
    def load_reconstructor(self):
        self.logger.info("Loading reconstructor object")
        for rec_key, recording in self.recordings.items():
            h5_path = recording['h5_path']
            h5 = h5py.File(h5_path)
            stream_ids = list(h5['wells'].keys())
            if self.stream_select is not None: stream_ids = [stream_ids[self.stream_select]] if isinstance(stream_ids[self.stream_select], str) else stream_ids[self.stream_select]
            for stream_id in stream_ids:
                dill_file_path = os.path.join(self.reconstructor_dir, f'{rec_key}_{stream_id}.dill')
                try:
                    assert os.path.exists(dill_file_path), f"Reconstructor object not found at {dill_file_path}. Generating new reconstructor object."
                    self.logger.info(f"Attempting to load reconstructor object from {dill_file_path}")
                    with open(dill_file_path, 'rb') as f:
                        loaded_obj = dill.load(f)
                        
                        # Exclude certain runtime parameters
                        keys_to_delete = []
                        exclusions_key_words = ['switch', 'load', 'select', 'max_workers', 'n_jobs']
                        for key, item in loaded_obj.__dict__.items():
                            for word in exclusions_key_words:
                                if word in key: keys_to_delete.append(key)
                        for key in keys_to_delete: del loaded_obj.__dict__[key]

                        self.update_nested_dict(self.__dict__, loaded_obj.__dict__)
                        if self.reconstructor_load_options['restore_environment']: 
                            self.logger.warning("Environment restoration not yet implemented. Skipping.")
                except Exception as e: 
                    # Some error handling:
                    # error = '[Errno 2] No such file or directory'
                    # if error in str(e): self.logger.warning(f"Reconstructor object not found at {dill_file_path}.")
                    # else: 
                    self.logger.warning(f"Error loading reconstructor object: {e}.")       
                    self.logger.warning(f"load_reconstructor option will be set to False for this run. Generating new reconstructor object.")
                    self.reconstructor_load_options['load_reconstructor'] = False
        return self

    def _restore_environment(self):
        # Write the requirements to a file
        requirements_file_path = os.path.join(self.reconstructor_dir, 'requirements.txt')
        with open(requirements_file_path, 'w') as f:
            f.write("\n".join(self.requirements))

        # Reinstall the requirements
        subprocess.run(['pip', 'install', '-r', requirements_file_path])

        # Check out the specific git version
        subprocess.run(['git', 'checkout', self.git_version])

        self.logger.info("Environment and git version restored")

    def bypass_to_templates(self):
        self.logger.info("Attempting to bypass pipeline steps to template extraction")
        try:
            assert self.reconstructor_load_options['load_reconstructor'], "Load reconstructor option is set to False. Cannot bypass to templates without reconstructor object."
            self.extract_templates(template_bypass=True)
            #self.extract_templates()
            self.logger.info("Bypass to templates successful")
            self.concatenate_switch = False
            self.sort_switch = False
            self.waveform_switch = False
            self.template_switch = False
        except Exception as e: self.logger.warning(f"Error bypassing to templates: {e}. Continuing with pipeline execution as normal.")
    
    def run_pipeline(self, **kwargs):
        # Pipeline switches
        self.concatenate_switch = kwargs.get('concatenate_switch', True)
        self.sort_switch = kwargs.get('sort_switch', True)
        self.waveform_switch = kwargs.get('waveform_switch', True)
        self.template_switch = kwargs.get('template_switch', True)
        self.recon_switch = kwargs.get('recon_switch', True)

        # Pipeline steps
        self.logger.info("Starting pipeline execution")        
        self.load_recordings()
        if self.reconstructor_load_options['load_reconstructor']: self.load_reconstructor() #TODO: Finish implementing environment load at some point. Not important.
        if self.reconstructor_load_options['load_templates_bypass']: self.bypass_to_templates() # Useful if sorting and waveform temp data have been dealt with but templates are still available
        if self.concatenate_switch: self.concatenate_recordings()
        if self.sort_switch: self.spikesort_recordings()
        if self.waveform_switch: self.extract_waveforms()
        if self.template_switch: self.extract_templates()        
        if self.recon_switch: self.analyze_and_reconstruct()
        if self.save_reconstructor_object: self.save_reconstructor()
        if self.run_lean: self.clean_up()
        self.logger.info("Pipeline execution completed")

    def clean_up(self):
        self.logger.info("Cleaning up intermediate files")
        if os.path.exists(self.waveforms_dir):
            shutil.rmtree(self.waveforms_dir)
            self.logger.debug("Removed waveforms directory: %s", self.waveforms_dir)
        self.logger.info("Clean-up completed")