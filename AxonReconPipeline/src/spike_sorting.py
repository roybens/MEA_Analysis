# Source: https://github.com/hornauerp/axon_tracking/blob/main/axon_tracking/spike_sorting.py
import os, time, h5py
import spikeinterface.full as si
import numpy as np
from pathlib import Path
from tqdm import tqdm
from glob import glob

def sort_recording_list(path_list, save_path_changes, sorter, sorter_params = dict(), clear_files=True, verbose=True):
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
    
    for rec_path in tqdm(path_list, desc="Sorting recordings"):
        
        h5 = h5py.File(rec_path)
        #Check that all wells are recorded throughout all recordings (should not fail)
        stream_ids = list(h5['wells'].keys())

        save_root = convert_rec_path_to_save_path(rec_path, save_path_changes)
        
        for stream_id in tqdm(stream_ids, desc="Sorting wells"):
            sorter_output_file = Path(os.path.join(save_root, stream_id, 'sorter_output', 'amplitudes.npy'))
            if not os.path.exists(sorter_output_file):
                multirecording, common_el = concatenate_recording_slices(rec_path, stream_id)
                sorting = clean_sorting(multirecording, save_root, stream_id, sorter, sorter_params, clear_files=clear_files, verbose=verbose)
                sorting_list.append(sorting)
            
    return sorting_list



def convert_rec_path_to_save_path(rec_path, save_path_changes):
    """
    Function that converts a recording path to the corresponding save path for a spike sorting.

    Arguments
    ----------
    rec_path: str
        Path to the axon scan file.
    save_path_changes: dict
        Dictionary containing keys 'pos' and 'vals' that indicate the changes to be made to the rec_path.
        Refer to the inidices after splitting the path by '/'.
        
    Returns
    ----------
    save_path: str
        Root save path. Well ID will be appended during the sorting.
    """

    path_parts = rec_path.split('/')
    for x,y in zip(save_path_changes['pos'], save_path_changes['vals']):
        path_parts[x] = y
        
    save_path = os.path.join(*path_parts)
    
    return save_path

def find_common_electrodes(rec_path, stream_id):
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
        rec = si.MaxwellRecordingExtractor(rec_path, stream_id=stream_id, rec_name=rec_name)
        rec_el = rec.get_property("contact_vector")["electrode"]
        if i == 0:
            common_el = rec_el
        else:
            common_el = list(set(common_el).intersection(rec_el))
            
    return rec_names, common_el



def concatenate_recording_slices(rec_path, stream_id):
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

    rec_names, common_el = find_common_electrodes(rec_path, stream_id)
    if len(rec_names) == 1:
        rec = si.MaxwellRecordingExtractor(rec_path, stream_id=stream_id, rec_name=rec_names[0])
        return rec
    else:
        rec_list = []
        for rec_name in rec_names: 
            #rec_name = 'rec' + '%0*d' % (4, r)
            rec = si.MaxwellRecordingExtractor(rec_path, stream_id=stream_id, rec_name=rec_name)
                    
            ch_id = rec.get_property("contact_vector")['device_channel_indices']
            rec_el = rec.get_property("contact_vector")["electrode"]
            
            chan_idx = [np.where(rec_el == el)[0][0] for el in common_el]
            sel_channels = rec.get_channel_ids()[chan_idx]
            chunk_size = np.min([10000, rec.get_num_samples()]) - 100 #Fallback for ultra short recordings (too little activity)
            rec_centered = si.center(rec,chunk_size=chunk_size)
            rec_list.append(rec_centered.channel_slice(sel_channels, renamed_channel_ids=list(range(len(chan_idx)))))
        
        multirecording = si.concatenate_recordings(rec_list)
    
        return multirecording, common_el
    


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

        if verbose:
            print(f"DURATION: {rec.get_num_frames() / rec.get_sampling_frequency()} s -- "
                    f"NUM. CHANNELS: {rec.get_num_channels()}")

        # We use try/catch to not break loops when iterating over several sortings (e.g. when not all wells were recorded)
        try:
            t_start_sort = time.time()
            sorting = si.run_sorter(sorter, rec, output_folder=output_folder, verbose=verbose,
                                    **sorter_params)
            if verbose:
                print(f"\n\nSpike sorting elapsed time {time.time() - t_start_sort} s")
            
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
                
    return sorting

def generate_rec_list(path_parts):
    """
    Function that takes a list of strings (path parts) and finds all recordings matching the path pattern, and returns the stream ids for the first recordings.

    Arguments
    ----------
    path_parts: List of strings
        Parts of the path pattern to be concatenated to look for recordings matching the description (when using *wildcard)

    Returns
    ----------
    path_list: List of strings
        List of the recordings matching the pattern provided in path_parts.
    stream_ids: List of strings
        List of stream_ids (wells) recorded from the first recording.
    """
    path_pattern = os.path.join(*path_parts)
    path_list  = glob(path_pattern)
    h5 = h5py.File(path_list[0])
    stream_ids = list(h5['wells'].keys())
    path_list.sort()
    
    return path_list, stream_ids 

def concatenate_recording_list(path_list, stream_id):
    well_recording_list = []
    for rec_path in path_list: #Iterate over recordings to be concatenated
        try: # If not all wells were recorded, should be the only cause for an error
            rec = si.MaxwellRecordingExtractor(rec_path,stream_id=stream_id)
            well_recording_list.append(rec)
        except Exception:
            continue
               
    if len(well_recording_list) == len(path_list):
        multirecording = si.concatenate_recordings(well_recording_list)
    else:
        raise ValueError('Could not load all recordings!')
        
    saturated_count = find_saturated_channels(well_recording_list)
    clean_multirecording = multirecording.remove_channels(multirecording.get_channel_ids()[saturated_count>0])
    
    
    return multirecording

def find_saturated_channels(rec_list, threshold=0):
    """
    Function that creates output folder if it does not exist, sorts the recording using the specified sorter
    and clears up large files afterwards. 

    Arguments
    ----------
    rec_list: List of MaxwellRecordingExtractor objects.
        List of (potentially to be concatenated) recordings to be checked for saturated channels.
    threshold: float
        Maximum ratio of saturated signal for the channel to still be accepted as non-saturated. 

    Returns
    ----------
    saturated_count: np.array
        Number of recordings in which the saturation threshold was crossed (channel was considered to be saturated). Values go from 0 to len(rec_list). 
    """
    saturated_count = np.zeros((rec_list[0].get_num_channels()))
    
    for i in range(0, len(rec_list)):
        random_data = si.get_random_data_chunks(rec_list[i], num_chunks_per_segment = int((rec_list[i].get_total_duration()/60)))
        saturated = (np.sum((random_data == 0).astype("int16") + (random_data == 1023).astype("int16"),axis=0)) / random_data.shape[0]
        saturated_count += saturated > threshold
    return saturated_count