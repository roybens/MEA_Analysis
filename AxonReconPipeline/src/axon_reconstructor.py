import os
import shutil
import logging
import spikeinterface.sorters as ss
import lib_helper_functions as helper
import lib_sorting_functions as sorter
import lib_waveform_functions as waveformer
import lib_template_functions as templater
import lib_axon_velocity_functions as axoner
import spikeinterface.full as si
import h5py
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from MEAProcessingLibrary import mea_processing_library as MPL

class AxonReconstructor:
    def __init__(self, h5_parent_dirs, allowed_scan_types=['AxonTracking'], stream_select=None, **kwargs):
        self.logger = self.setup_logger(kwargs.get('logger_level', 'INFO'))

        self.logger.info("Initializing AxonReconstructor")

        self.h5_parent_dirs = h5_parent_dirs
        self.allowed_scan_types = allowed_scan_types
        self.stream_select = stream_select
        self.unit_limit = kwargs.get('unit_limit', None)

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
        self.recon_dir = kwargs.get('recon_dir', './AxonReconPipeline/data/temp_data/reconstructions')

        self.n_jobs = kwargs.get('n_jobs', 4)
        self.max_workers = kwargs.get('max_workers', 32)
        assert self.max_workers % self.n_jobs == 0, "max_workers must be a multiple of n_jobs to avoid oversubscription."

        self.sorting_params = self.get_sorting_params(kwargs)
        self.te_params = self.get_te_params(kwargs)
        self.qc_params = self.get_qc_params(kwargs)
        self.av_params = self.get_av_params(kwargs)

        self.continuous_h5_file_info = None
        self.sorting_lists_by_scan = None
        self.waveforms = None
        self.templates = None
        self.reconstruction_results = None

        self.logger.debug("AxonReconstructor initialized with parameters: %s", self.__dict__)

    def setup_logger(self, logger_level='INFO'):
        logger = logging.getLogger('axon_reconstructor')
        if not logger.handlers:
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            logger.setLevel(logger_level.upper())        
            stream_handler = logging.StreamHandler()
            stream_handler.setLevel(logger_level.upper())
            stream_handler.setFormatter(formatter)        
            logger.addHandler(stream_handler)
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
        }
        default_params.update(kwargs.get('te_params', {}))
        for key, value in default_params.items():
            for k, v in kwargs.get('te_params', {}).items():
                if k == key: default_params[key] = v
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
            'detect_threshold': 0.01,
            'kurt_threshold': 0.1,
            'peak_std_threshold': 0.8,
            'upsample': 5,
            'neighbor_radius': 100,
            'r2_threshold': 0.8
        }
        av_params.update(kwargs.get('av_params', {}))
        for key, value in av_params.items():
            for k, v in kwargs.get('av_params', {}).items():
                if k == key: av_params[key] = v
        return av_params

    def load_recordings(self):
        self.logger.info("Loading continuous h5 files")
        
        self.recordings = {}        
        for dir in self.h5_parent_dirs:
            h5_file_paths = MPL.extract_raw_h5_filepaths(dir)
            for h5_path in h5_file_paths:
                h5_details = MPL.extract_recording_details(h5_path)
                date = h5_details[0]['date']
                chip_id = h5_details[0]['chipID']
                run_id = h5_details[0]['runID']
                scan_type = h5_details[0]['scanType']
                if scan_type in self.sorting_params['allowed_scan_types']:
                    device, recording_segments, stream_count, rec_counts = MPL.load_recordings(h5_path, stream_select=self.stream_select, logger=self.logger)
                    self.recordings[f"{date}_{chip_id}_{run_id}"] = {
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

        #self.logger.debug("Loaded continuous h5 files: %s", self.recordings)
        self.logger.info("Completed loading continuous h5 files")

    def concatenate_recordings(self):
        self.logger.info("Concatenating recordings")
        self.multirecordings = {}
        for recording in self.recordings.values():
            date = recording['date']
            chip_id = recording['chip_id']
            scanType = recording['scanType']
            run_id = recording['run_id']
            h5_path = recording['h5_path']
            h5 = h5py.File(h5_path)
            stream_ids = list(h5['wells'].keys())
            if self.stream_select is not None: 
                stream_ids = [stream_ids[self.stream_select]] if isinstance(stream_ids[self.stream_select], str) else stream_ids[self.stream_select]
            streams = {}
            for stream_id in stream_ids:
                recording_segments = recording['streams'][stream_id]['recording_segments']
                self.logger.info(f'Concatenating recording segments for {date}_{chip_id}_{run_id} stream {stream_id}')
                multirecording, common_el, multirec_save_path = sorter.concatenate_recording_segments(
                    h5_path, recording_segments, stream_id, save_dir=self.recordings_dir, logger=self.logger)  
                streams[stream_id] = {
                    'multirecording': multirecording,
                    'common_el': common_el,
                    'multirec_save_path': multirec_save_path
                }    
            self.multirecordings[f"{date}_{chip_id}_{run_id}"] = {
                'scanType': scanType,
                'date': date,
                'chip_id': chip_id,
                'run_id': run_id,
                'streams': streams
            }
        self.logger.info("Completed concatenation of recordings")
    
    def spikesort_recordings(self):
        self.logger.info("Starting spike sorting process")
        self.sortings = {} 
        for multirec in self.multirecordings.values():
            date = multirec['date']
            chip_id = multirec['chip_id']
            scanType = multirec['scanType']
            run_id = multirec['run_id']
            self.logger.info(f'Generating spike sorting objects for {date}_{chip_id}_{run_id}')
            streams = {}
            for stream_id, stream in multirec['streams'].items():
                mr = stream['multirecording']
                spikesorting_root = os.path.join(self.sorting_params['sortings_dir'], f'{date}/{chip_id}/{scanType}/{run_id}')
                sorting, stream_sort_path = sorter.sort_multirecording(mr, stream_id, save_root=spikesorting_root, sorting_params=self.sorting_params, logger=self.logger)
                streams[stream_id] = {
                    'sorting_path': stream_sort_path,
                    'sorting': sorting
                } 
            self.sortings[f"{date}_{chip_id}_{run_id}"] = {
                'scanType': scanType,
                'date': date,
                'chip_id': chip_id,
                'run_id': run_id,
                'spiksorting_root': spikesorting_root,
                'streams': streams
            }
        self.logger.debug("Completed spike sorting")

    def extract_waveforms(self):
        self.waveforms = {}
        self.logger.info("Extracting waveforms")
        for key, sorting in self.sortings.items():
            multirecs = self.multirecordings[key]
            streams = {}
            for stream_id, multirec in multirecs['streams'].items():
                mr = multirec['multirecording']
                sort = sorting['streams'][stream_id]['sorting']
                cleaned_sorting = waveformer.select_units(sort, **self.qc_params)
                cleaned_sorting = si.remove_excess_spikes(cleaned_sorting, mr) # Relevant if last spike time == recording_length
                cleaned_sorting.register_recording(mr)
                segment_sorting = si.SplitSegmentSorting(cleaned_sorting, mr)
                self.logger.info(f'Extracting waveforms from stream: {stream_id}')
                h5_path = self.recordings[key]['h5_path']
                waveforms = waveformer.extract_unit_waveforms(h5_path, stream_id, segment_sorting, save_root=self.waveforms_dir, logger=self.logger, te_params = self.te_params) # TODO: make distinct wf_params
                streams[stream_id] = waveforms
                if self.stream_select is not None:
                    formatted = 'well{:03}'.format(self.stream_select)
                    if stream_id == formatted:
                        break
            self.waveforms[key] = {
                'scanType': multirecs['scanType'],
                'date': multirecs['date'],
                'chip_id': multirecs['chip_id'],
                'run_id': multirecs['run_id'],
                'waveforms_dir': self.waveforms_dir,
                'streams': streams
            }
        self.logger.info("Completed waveform extraction")

    def extract_templates(self):
        self.logger.info("Extracting templates")
        self.templates = {}
        for key, waveforms in self.waveforms.items():
            streams = {}
            for stream_id, wfs in waveforms['streams'].items():
                if self.stream_select is not None: 
                    if stream_id != 'well{:03}'.format(self.stream_select): 
                        continue
                self.logger.info(f'Extracting templates from stream: {stream_id}')
                h5_path = self.recordings[key]['h5_path']
                sorting = self.sortings[key]['streams'][stream_id]['sorting']
                multirec = self.multirecordings[key]['streams'][stream_id]['multirecording']
                unit_templates = templater.extract_templates(
                    multirec, sorting, wfs, h5_path, stream_id, save_root=self.templates_dir, 
                    te_params=self.te_params, qc_params=self.qc_params, unit_limit=self.unit_limit, logger=self.logger)
                streams[stream_id] = {'units': unit_templates}
            self.templates[key] = {
                'scanType': waveforms['scanType'],
                'date': waveforms['date'],
                'chip_id': waveforms['chip_id'],
                'run_id': waveforms['run_id'],
                'templates_dir': self.templates_dir,
                'streams': streams
            }
        self.logger.debug("Completed template extraction")

    def analyze_and_reconstruct(self):
        self.logger.info("Reconstructing Axonal Morphology")
        
        self.reconstruction_results = axoner.analyze_and_reconstruct(
            self.templates,
            recon_dir=self.recon_dir,
            params=self.av_params,
            stream_select=self.stream_select,
            n_jobs=self.n_jobs,
            logger=self.logger
        )
        self.logger.debug("Completed analysis and reconstruction")
    
    def run_pipeline(self):
        self.logger.info("Starting pipeline execution")        
        self.load_recordings()
        self.concatenate_recordings()
        self.spikesort_recordings()
        self.extract_waveforms()
        self.extract_templates()        
        self.analyze_and_reconstruct()

        if self.run_lean:
            self.clean_up()

        self.logger.info("Pipeline execution completed")

    def clean_up(self):
        self.logger.info("Cleaning up intermediate files")
        if os.path.exists(self.sortings_dir):
            shutil.rmtree(self.sortings_dir)
            self.logger.debug("Removed sortings directory: %s", self.sortings_dir)
        if os.path.exists(self.waveforms_dir):
            shutil.rmtree(self.waveforms_dir)
            self.logger.debug("Removed waveforms directory: %s", self.waveforms_dir)
        if os.path.exists(self.templates_dir):
            shutil.rmtree(self.templates_dir)
            self.logger.debug("Removed templates directory: %s", self.templates_dir)
        self.logger.info("Clean-up completed")
