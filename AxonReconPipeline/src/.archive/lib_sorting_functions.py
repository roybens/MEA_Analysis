''' This file contains sorting-related functions that are used in the main pipeline functions. '''

'''imports'''
import sys
import os
import h5py
import spikeinterface.full as si
import spikeinterface.sorters as ss
import shutil
import numpy as np
from tqdm import tqdm
from pathlib import Path
import time

''' Local imports '''
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from MEAProcessingLibrary import mea_processing_library as MPL

'''logging setup'''
import logging
logger = logging.getLogger(__name__) #Create a logger
logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler()  # Create handlers, logs to console
stream_handler.setLevel(logging.DEBUG) # Set level of handlers
logger.addHandler(stream_handler) # Add handlers to the logger
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - %(module)s.%(funcName)s') # Create formatters and add it to handlers
stream_handler.setFormatter(formatter)

'''default settings'''
default_recording_dir = './AxonReconPipeline/data/temp_data/recordings'
default_sorting_dir = './AxonReconPipeline/data/temp_data/sortings'
default_n_jobs = 4

'''multirecording functions'''
def find_common_electrodes(rec_path, stream_id, scan_merge = False):
    """
    Function that returns the common electrodes of the successive axon scan recordings.
    
    Arguments
    ----------
    rec_path: str
        Path to the axon scan file.
    stream_id: str
        Well ID in the format "well***"; Well 1 would be "well001", Well 20 would be "well020"
        
    Returns
    ----------
    rec_names: list
        List of rec_names for the specified recording/well.
    common_el: list
        List of electrodes that are present in all axon scan recordings.
    """
    
    assert(os.path.exists(rec_path))
    
    h5 = h5py.File(rec_path)
    rec_names = list(h5['wells'][stream_id].keys())

    for i, rec_name in enumerate(rec_names):
        #rec_name = 'rec' + '%0*d' % (4, rec_id)
        try: 
            rec = si.MaxwellRecordingExtractor(rec_path, stream_id=stream_id, rec_name=rec_name)
            rec_el = rec.get_property("contact_vector")["electrode"]
            if i == 0:
                common_el = rec_el
            else:
                common_el = list(set(common_el).intersection(rec_el))
        except Exception as e:
            print(e)
            continue
            
    return rec_names, common_el
def concatenate_recording_slices(rec_path, stream_id, recording_dir = default_recording_dir, n_jobs = default_n_jobs):
    """
    Function that centers and concatenates the recordings of an axon scan for all common electrodes. 

    Arguments
    ----------
    rec_path: str
        Path to the axon scan file.
    stream_id: str
        Well ID in the format "well***"; Well 1 would be "well001", Well 20 would be "well020"
        
    Returns
    ----------
    multirecording: ConcatenatedRecordingSlice
        Concatenated recording across common electrodes (spikeinterface object)
    """
    
    '''Prepare the recording directory'''
    rec_names, common_el = find_common_electrodes(rec_path, stream_id)
    recording_details = MPL.extract_recording_details(rec_path)
    date = recording_details[0]['date']
    chip_id = recording_details[0]['chipID']
    scanType = recording_details[0]['scanType']
    run_id = recording_details[0]['runID']
    recording_dir_by_stream =recording_dir+f'/{date}/{chip_id}/{scanType}/{run_id}'
    recording_dir_by_stream = os.path.join(recording_dir_by_stream, stream_id)  

    '''Try to load multi-recording if it exists, otherwise create it'''
    multirec_name = f'{stream_id}_multirecording'
    multirec_save_path = recording_dir_by_stream+f'/{multirec_name}'
    # try: 
    #     assert False #Force the except block, remove later
    #     multirecording = si.load_extractor(multirec_save_path)
    #     logger.info(f"Loaded multirecording {multirec_name} from {multirec_save_path}")
    #     return multirecording, common_el
    # except:
    #     if os.path.exists(multirec_save_path):
    #         shutil.rmtree(multirec_save_path)
    #     os.makedirs(multirec_save_path, exist_ok=True)

    '''Center recordings and remove non-common electrodes'''
    rec_list = []          
    save_kwargs = dict(n_jobs=n_jobs)            
    for rec_name in rec_names:  
        try:
            rec_save_path = recording_dir_by_stream+f'/{rec_name}'
            rec = si.MaxwellRecordingExtractor(rec_path, stream_id=stream_id, rec_name=rec_name)                    
            try: 
                if not os.path.exists(rec_save_path):
                    raise FileNotFoundError
                rec_centered = si.load_extractor(rec_save_path)
                logger.info(f"Loaded centered recording {rec_name} from {rec_save_path}")
            except:
                if os.path.exists(rec_save_path):
                    shutil.rmtree(rec_save_path)
                
                
                chunk_size = np.min([10000, rec.get_num_samples()]) - 100 #Fallback for ultra short recordings (too little activity)
                logger.info(f"Centering recording {rec_name} with chunk size {chunk_size}")
                rec_centered = si.center(rec,chunk_size=chunk_size)                    
                #os.makedirs(rec_save_path, exist_ok=True) 
                #rec_centered.save(folder = rec_save_path, overwrite = True, verbose = True, **save_kwargs)
                #logger.info(f"Saved centered recording to {rec_save_path}")
            ch_id = rec.get_property("contact_vector")['device_channel_indices']
            rec_el = rec.get_property("contact_vector")["electrode"]                
            chan_idx = [np.where(rec_el == el)[0][0] for el in common_el]
            sel_channels = rec.get_channel_ids()[chan_idx]   
            rec_centered_sliced = rec_centered.channel_slice(sel_channels, renamed_channel_ids=list(range(len(chan_idx))))         
            rec_list.append(rec_centered_sliced)
        except Exception as e:
            print(e)
            continue 
        
    '''Concatenate the recordings across common electrodes and save the multirecording'''
    multirecording = si.concatenate_recordings(rec_list)
    logger.info(f"Concatenated recordings across common electrodes")
    #logger.info(f'Saving...')
    #multirecording.save(folder = multirec_save_path, overwrite = True, verbose = True, **save_kwargs)
    #logger.info(f"Saved multirecording to {multirec_save_path}")

    '''Optionally delete recordings to save space'''
    #TODO implement this

    return multirecording, common_el
def build_multirecording_objects(continuous_h5_dirs, allowed_scan_types=["AxonTracking"], recording_dir = default_recording_dir, stream_select = None, n_jobs = default_n_jobs):
    
    '''Generate multirecordings for all streams in all scans'''
    #multi_recs_by_scan = []
    #common_el_by_scan = []
    rec_list = [dir['h5_file_path'] for dir in continuous_h5_dirs.values()]
    for rec_path in rec_list:
        h5 = h5py.File(rec_path)
        stream_ids = list(h5['wells'].keys())
        if stream_select is not None:
            stream_ids = stream_ids[stream_select]
            if isinstance(stream_ids, str):
                stream_ids = [stream_ids]
        '''Iterate over all streams in the scan and build multirecordings for each stream'''
        multirecs_by_stream = []
        common_el_by_stream = []
        recording_details = MPL.extract_recording_details(rec_path)
        scanType = recording_details[0]['scanType']
        if scanType in allowed_scan_types:
            for stream_id in stream_ids:
                multirecording, common_el = concatenate_recording_slices(rec_path, stream_id,recording_dir = recording_dir, n_jobs = n_jobs)
                multirecs_by_stream.append(multirecording)
                common_el_by_stream.append(common_el)
           # multi_recs_by_scan.append(multirecs_by_stream)
            #common_el_by_scan.append(common_el_by_stream)
    #return multi_recs_by_scan, common_el_by_scan

'''sorting functions'''
def clean_sorting(rec, save_root, stream_id, sorter, sorter_params = dict(), clear_files=True, verbose=True):
    """
    Function that creates output folder if it does not exist, sorts the recording using the specified sorter
    and clears up large files afterwards. 

    Arguments
    ----------
    rec: MaxwellRecordingExtractor
        Recording to be sorted.
    save_root: str
        Root path where the sorted data will be stored. Stream name (i.e. well ID) will be appended.
    stream_id: str
        Well ID in the format "well***"; Well 1 would be "well001", Well 20 would be "well020"
    sorter: str
        Name of the sorter, as per the spikeinterface convention (e.g. 'kilosort2_5')
    sorter_params: dict
        Dictionary containing parameters for the spike sorter. If left empty, will default to default parameters.
    clear_files: bool
        Flag if large temporary files should be deleted after the sorting. Default: True
    verbose: bool
        Default: True

    Returns
    ----------
    sorting: Sorting object 
        Specific type depends on the sorter.
    """
    
    output_folder = Path(os.path.join(save_root, stream_id))
    sorter_output_file = os.path.join(output_folder, 'sorter_output', 'amplitudes.npy')
    sorting = []
    # Creates output folder if sorting has not yet been done
    if os.path.exists(sorter_output_file):
        return sorting
    elif (rec.get_total_duration() < 30):
        full_output_folder = Path(os.path.join(output_folder, 'sorter_output'))
        full_output_folder.mkdir(parents=True, exist_ok=True)
        np.save(sorter_output_file, np.empty(0)) #Empty file to indicate a failed sorting for future loops
        return sorting
    else:
        output_folder.mkdir(parents=True, exist_ok=True)
        raw_file = os.path.join(output_folder, 'sorter_output', 'recording.dat')
        wh_file = os.path.join(output_folder, 'sorter_output', 'temp_wh.dat')
        logger.debug(f"DURATION: {rec.get_num_frames() / rec.get_sampling_frequency()} s -- " f"NUM. CHANNELS: {rec.get_num_channels()}")

        # We use try/catch to not break loops when iterating over several sortings (e.g. when not all wells were recorded)
        try:
            t_start_sort = time.time()
            sorting = MPL.benshalom_kilosort2_docker_image(rec, output_folder=output_folder, verbose=verbose,)
            logger.debug(f"Spike sorting elapsed time {time.time() - t_start_sort} s")
            
            #Making sure we clean up the largest temporary files
            if clear_files & os.path.exists(wh_file):
                os.remove(wh_file)
            if clear_files & os.path.exists(raw_file):
                os.remove(raw_file)
        except Exception as e:
            print(e)
            if clear_files & os.path.exists(wh_file):
                os.remove(wh_file)
            if clear_files & os.path.exists(raw_file):
                os.remove(raw_file)
                
def sort_recording_list(path_list, save_root, sorter_params = dict(), clear_files=True, stream_select = None):
    """
    Function that iterates over a list of axon scans, finds common electrodes, concatenates and spike sorts the recording slices. 

    Arguments
    ----------
    path_list: list
        List of string referring to the paths of axon scan recording.
    save_path_changes: dict
        Dictionary containing keys 'pos' and 'vals' that indicate the changes to be made to the rec_path.
        Refer to the inidices after splitting the path by '/'.
    sorter: str
        Name of the sorter, as per the spikeinterface convention (e.g. 'kilosort2_5')
    sorter_params: dict
        Dictionary containing parameters for the spike sorter. If left empty, will default to default parameters.
    clear_files: bool
        Flag if large temporary files should be deleted after the sorting. Default: True
    verbose: bool
        Default: True

    Returns
    ----------
    sorting_list: list of sorting objects 
        Specific type depends on the sorter.
    
    """
    
    sorting_list = []
    
    if isinstance(path_list, str):
        path_list = [path_list]
    stream_ids = []
    
    for rec_path in tqdm(path_list, desc="Sorting recordings"):
        
        h5 = h5py.File(rec_path)
        #Check that all wells are recorded throughout all recordings (should not fail)
        stream_ids = list(h5['wells'].keys())
        if stream_select is not None:
            stream_ids = stream_ids[stream_select]
            if isinstance(stream_ids, str):
                stream_ids = [stream_ids]
        
        for stream_id in tqdm(stream_ids, desc="Sorting wells"):
            
            stream_id_path = os.path.join(save_root, stream_id)
            try:                    
                #sorting = ss.Kilosort2_5Sorter._get_result_from_folder(stream_id_path+'/sorter_output/')
                sorting = ss.Kilosort2Sorter._get_result_from_folder(stream_id_path+'/sorter_output/')
                sorting_list.append(sorting)
                continue                    
            except:
                pass
            
            sorter_output_file = Path(os.path.join(save_root, stream_id, 'sorter_output', 'amplitudes.npy'))
            if not os.path.exists(sorter_output_file):
                multirecording, common_el = concatenate_recording_slices(rec_path, stream_id)
                #sorting = 
                sorting = clean_sorting(multirecording, save_root, stream_id, sorter, sorter_params, clear_files=clear_files) #, verbose=verbose)
                sorting_list.append(sorting)            
    return sorting_list   
def spikesort_recordings(continuous_h5_dirs, sortings_dir = default_sorting_dir, allowed_scan_types=["AxonTracking"], stream_select = None):
    #Sorter params - dont really need to define this here, remove later
    sorter_params = si.get_default_sorter_params(si.Kilosort2_5Sorter)
    sorter_params['n_jobs'] = -1
    sorter_params['detect_threshold'] = 7
    sorter_params['minFR'] = 0.01
    sorter_params['minfr_goodchannels'] = 0.01
    sorter_params['keep_good_only'] = False
    sorter_params['do_correction'] = False
    verbose = True
    
    logger.info('Generating spike sorting objects and merging as needed.')    
    sorting_lists_by_scan = []
    for dir in continuous_h5_dirs.values():
        h5_file_path = dir['h5_file_path']
        scantype = dir['scanType']
        if scantype in allowed_scan_types:
            #logger.info(f"Processing {h5_file_path}...")
            #sorting_list_by_stream = 
            #extract_sortings_per_stream(h5_file_path, sortings_dir = sortings_dir, stream_select=stream_select)
            #sorting_lists_by_scan.append(sorting_list_by_stream)

            """ 
            Extract sorting objects per stream from sliced and cocatenated recordings (most/all fixed electrodes during axontracking scan)
            """
            recording_details = MPL.extract_recording_details(h5_file_path)
            date = recording_details[0]['date']
            chip_id = recording_details[0]['chipID']
            scanType = recording_details[0]['scanType']
            run_id = recording_details[0]['runID']            
            spikesorting_root = sortings_dir+f'/{date}/{chip_id}/{scanType}/{run_id}'

            #h5_file_paths = [dir['h5_file_path'] for dir in continuous_h5_dirs.values() if dir['scanType'] == "AxonTracking"]
            #sorting_list_by_stream = 
            sorting_list_by_stream = sort_recording_list(h5_file_path, save_root=spikesorting_root, sorter_params = sorter_params, stream_select = stream_select) 
            sorting_lists_by_scan.append(sorting_list_by_stream)
            #return sorting_list_by_stream
    return sorting_lists_by_scan