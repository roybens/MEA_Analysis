import os
import h5py
import spikeinterface.full as si
import spikeinterface.sorters as ss
import numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from MEAProcessingLibrary import mea_processing_library as MPL
import logging

default_recording_dir = './AxonReconPipeline/data/temp_data/recordings'
default_sorting_dir = './AxonReconPipeline/data/temp_data/sortings'

def find_common_electrodes(rec_path, stream_id, n_jobs=4, max_workers=24, logger=None):
    if logger:
        logger.info(f"Finding common electrodes in {rec_path} for stream {stream_id}")
    
    assert os.path.exists(rec_path)
    h5 = h5py.File(rec_path)
    rec_names = list(h5['wells'][stream_id].keys())

    def process_recording(rec_name):
        try:
            rec = si.MaxwellRecordingExtractor(rec_path, stream_id=stream_id, rec_name=rec_name)
            rec_el = rec.get_property("contact_vector")["electrode"]
            return rec_el
        except Exception as e:
            if logger:
                logger.warning(e)
            return None

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(process_recording, rec_names))
    
    common_el = set(results[0])
    for rec_el in results[1:]:
        try:
            rec_el = [el for el in rec_el if el is not None]
            common_el.intersection_update(rec_el)
        except:
            continue
    
    return rec_names, list(common_el)

def concatenate_recording_segments(h5_path, recording_segments, stream_id=None, save_dir=None, n_jobs=4, max_workers=24, logger=None):
    
    rec_names, common_el = find_common_electrodes(h5_path, stream_id, n_jobs, max_workers, logger=logger)
    multirec_name = f'{stream_id}_multirecording'
<<<<<<< HEAD
    h5_details = MPL.extract_recording_details(h5_path)
    date = h5_details[0]['date']
    chip_id = h5_details[0]['chipID']
    scanType = h5_details[0]['scanType']
    run_id = h5_details[0]['runID']
    multirec_save_path = os.path.join(save_dir, f'{date}/{chip_id}/{scanType}/{run_id}', multirec_name)
=======
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
>>>>>>> 57c3339 (updating memory handling)

    def process_rec_name(rec):
        try:
<<<<<<< HEAD
            chunk_size = min([10000, rec.get_num_samples()]) - 100
            rec_centered = si.center(rec, chunk_size=chunk_size)
            rec_el = rec.get_property("contact_vector")["electrode"]
=======
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
>>>>>>> 57c3339 (updating memory handling)
            chan_idx = [np.where(rec_el == el)[0][0] for el in common_el]
            sel_channels = rec.get_channel_ids()[chan_idx]
            rec_centered_sliced = rec_centered.channel_slice(sel_channels, renamed_channel_ids=list(range(len(chan_idx))))
            return rec_centered_sliced
        except Exception as e:
            if logger:
                logger.error(e)
            return None

    if logger: logger.info(f"Concatenating recording segments for {h5_path}, stream {stream_id}")
    else: print(f"Concatenating recording segments for {h5_path}, stream {stream_id}")
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        rec_list = list(executor.map(process_rec_name, recording_segments))
    
    rec_list = [rec for rec in rec_list if rec is not None]
    multirecording = si.concatenate_recordings(rec_list)
    return multirecording, common_el, multirec_save_path

def sort_multirecording(multirecording, stream_id, save_root, sorting_params=dict(), logger=None):
    if logger:
        logger.info(f"Sorting multi-recording for stream {stream_id}")

    stream_sort_path = os.path.join(save_root, stream_id)
    sorter_output_folder = os.path.join(stream_sort_path, 'sorter_output')
    sorter_output_file = os.path.join(sorter_output_folder, 'amplitudes.npy')

    if not sorting_params:
        sorting_params = ss.Kilosort2Sorter.default_params()
    verbose = sorting_params.get('verbose', True)

    # Try to load existing sorting if available
    try:
        if sorting_params.get('load_existing_sortings', False) and os.path.exists(sorter_output_folder):
            sorting = ss.Kilosort2Sorter._get_result_from_folder(sorter_output_folder)
            message = f"Loaded existing sorting from {sorter_output_folder}"
            if logger: logger.info(message)
            else: print(message)
            return sorting, stream_sort_path, message
    except Exception as e:
        if logger: logger.error(f"Error loading existing sorting for stream {stream_id}: {e}")
        else: print(f"Error loading existing sorting for stream {stream_id}: {e}")

    # If no existing sorting found or loading failed, perform sorting
    try:
        output_folder = Path(stream_sort_path)
        output_folder.mkdir(parents=True, exist_ok=True)
        sorting_params_filtered = {k: v for k, v in sorting_params.items() if k in ss.Kilosort2Sorter.default_params()}
        sorting = MPL.benshalom_kilosort2_docker_image(multirecording, sorting_params=sorting_params_filtered, output_folder=stream_sort_path, verbose=verbose)
        message = f"Completed sorting and saved results to {sorter_output_folder}"
        if logger: logger.info(message)
        else: print(message)
        return sorting, stream_sort_path, message
    except Exception as e:
        message = f"Error sorting multi-recording for stream {stream_id}: {e}"
        if logger: logger.error(message)
        else: print(message)        
        return None, stream_sort_path, message

