#general imports
import logging
import os
import glob
import shutil
import json
import re
import numpy as np
from multiprocessing import Pool, Manager
from tqdm import tqdm
import filecmp
from timeit import default_timer as timer
import pandas as pd
import h5py
#import docker

#spikeinterface imports
import spikeinterface
import spikeinterface.full as si
import spikeinterface.extractors as se
import spikeinterface.preprocessing as spre
import spikeinterface.sorters as ss
import spikeinterface.core as sc
import spikeinterface.postprocessing as sp
import spikeinterface.preprocessing as spre
import spikeinterface.qualitymetrics as qm

#Logger Setup
#Create a logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Create handlers
stream_handler = logging.StreamHandler()  # logs to console
#file_handler = logging.FileHandler('file.log')  # logs to a file

# Set level of handlers
stream_handler.setLevel(logging.DEBUG)
#file_handler.setLevel(logging.ERROR)

# Add handlers to the logger
logger.addHandler(stream_handler)
#logger.addHandler(file_handler)

# Create formatters and add it to handlers
#formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - %(module)s.%(funcName)s')
stream_handler.setFormatter(formatter)

# def extract_raw_h5_filepaths(directories):
#     """
#     This function finds all files named 'data.raw.h5' in the given directories and their subdirectories.

#     Parameters:
#     directories: The list of directories to search for 'data.raw.h5' files.

#     Returns:
#     h5_dirs: The list of paths to the found 'data.raw.h5' files.
#     """
#     logger.info(f"Extracting .h5 file paths from directories:")
#     h5_dirs = []
#     for directory in directories:
#         for root, dirs, files in os.walk(directory):
#             for file in files:
#                 if file == "data.raw.h5":
#                     h5_dirs.append(os.path.join(root, file))
#     return h5_dirs

def extract_raw_h5_filepaths(h5_dir):
    #walk through the directory and find all .h5 files
    h5_subdirs = []
    if h5_dir.endswith('.h5'): return [h5_dir]
    for root, dirs, files in os.walk(h5_dir): 
        for file in files:
            if file.endswith('.h5'):
                h5_subdirs.append(os.path.join(root, file))
    # assert len(h5_subdirs) > 0, "No .h5 files found in the directory."
    # assert len(h5_subdirs) == 1, "Ambiguous file selection. Multiple .h5 files found in the directory."
    return h5_subdirs

def extract_recording_details(h5_dirs):
    """
    This function extracts details about each recording from the given h5 directories.
    The details include the file path, run ID, scan type, chip ID, and recording date.

    Parameters:
    h5_dirs: The list of h5 directories.

    Returns:
    records: The list of dictionaries, where each dictionary contains details about a recording.
    """

    # If h5_dirs is a string, convert it to a list with a single element
    if isinstance(h5_dirs, str):
        h5_dirs = [h5_dirs]

    logger.info(f"Extracting recording details from h5 directories:")
    records = []
    for h5_dir in h5_dirs:
        try: assert '.h5' in h5_dir, "The input is not a list of h5 directories."
        except: 
            try: 
                h5_subdirs = extract_raw_h5_filepaths(h5_dir)
                assert len(h5_subdirs) > 0, "No .h5 files found in the directory."
                assert len(h5_subdirs) == 1, "Ambiguous file selection. Multiple .h5 files found in the directory."
                h5_dir = h5_subdirs[0]
            except: logger.error("Some error occurred during the extraction of .h5 file paths."); continue

        parent_dir = os.path.dirname(h5_dir)
        runID = os.path.basename(parent_dir)

        grandparent_dir = os.path.dirname(parent_dir)
        scan_type = os.path.basename(grandparent_dir)

        great_grandparent_dir = os.path.dirname(grandparent_dir)
        chipID = os.path.basename(great_grandparent_dir)

        ggg_dir = os.path.dirname(great_grandparent_dir)
        date = os.path.basename(ggg_dir)

        record = {'h5_file_path': h5_dir, 
                  'runID': runID, 
                  'scanType': scan_type, 
                  'chipID': chipID,
                  'date': date}
        records.append(record)

    return records

# def test_continuity(h5_file_path, verbose=False, stream_select = None, recordings = None, MaxID = None):
#     """
#     This function tests the continuity of a given h5 file by attempting to read recordings from the file until an error occurs.
#     It also counts the number of successful reads (recordings).

#     If the verbose flag is set to True, the function logs the number of recordings detected.

#     If a recording is an exception object, an error has occurred. The function logs the error and appends False to a list.
#     If the recording is successfully read, the function logs the success and appends True to the list.

#     After attempting to read all recordings, the function checks if all items in the list are True. If they are, all recordings are continuous, and the function logs this and returns True. If not all items are True, the data is not continuous, and the function logs this and returns False.

#     Parameters:
#     h5_file_path (str): The path to the h5 file to read from.
#     verbose (bool): If True, log detailed output. Default is False.

#     Returns:
#     bool: True if all recordings are continuous, False otherwise.
#     """
#     logger.info(f"Testing continuity of {h5_file_path}:")
#     stream_count, rec_count = count_wells_and_recs(h5_file_path, verbose = verbose, stream_select = stream_select)
#     if stream_count == 0:
#         logger.error("No recordings detected, none are continuous.")
#         return False  
#     #This part of the function might be entirely unnecessary:
#     TrueorFalse_list = []
#     recordings = []
#    # if MaxID is None and recordings is None: MaxID, recordings = load_recordings(h5_file_path, stream_select)

#     for stream_num in range(stream_count):
#         if stream_select is not None and stream_num != stream_select: 
#             logger.info(f"Skipping stream {stream_num}.")
#             continue #skip if stream_select is not None and stream_num is not stream_select, saves time on loading and verifying data
#         #for rec_num in range(rec_count[stream_num]):
#         for recording in recordings:
#             # try: 
#             #     recording, rec_name, stream_id = get_data_maxwell(h5_file_path, rec_num = rec_num, well_num = stream_num, verbose = verbose)
#             #     recordings = recordings.append(recording)
#             # except: recording, rec_name, stream_id = get_data_maxwell(h5_file_path, rec_num = rec_num, verbose = verbose)
#             if isinstance(recording, BaseException):
#                 e = recording
#                 if "Unable to open object (object 'routed' doesn't exist)" in str(e):
#                     logger.error("This error indicates that 'RecordSpikesOnly' was active during this recording. Data are not continuous.")
#                 else:
#                     logger.error("Unknown error")
#                     logger.error(e)
#                     logger.error("This error is unexpected. Please check the error message and try to resolve it.")
#                 TrueorFalse_list.append(False)
#             else:
#                 logger.info(f"Successfully read Stream ID: {stream_id}, Recording: {rec_name}, indicating continuity.")
#                 TrueorFalse_list.append(True)
#     #if all items in TrueorFalse_list are True, then the data is continuous
#     if all(TrueorFalse_list) and TrueorFalse_list != []:
#         logger.info("All recordings are continuous.")
#         return True, recordings
#     elif TrueorFalse_list == []:
#         logger.error("No recordings detected, none are continuous.")
#         return False, None
#     else:
#         logger.error("Data are not continuous.")
#         return False, None
#     #This part of the function might be entirely unnecessary:

def load_recordings(h5_file_path, stream_select=None, logger=None):
    """
    Loads all recordings from the specified file and identifies whether it's a MaxOne or MaxTwo file type.

    Parameters:
    h5_file_path (str): The path to the file to read from.
    stream_select (int, optional): Specific stream ID to select. Default is None.
    logger (Logger, optional): Logger for logging. Default is None.

    Returns:
    tuple: 
    (1) MaxID indicating file type (1 for MaxOne, 2 for MaxTwo) 
    (2) List of recordings loaded 
    (3) Total number of streams
    (4) List of recording counts per stream
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    h5_details = extract_recording_details(h5_file_path)
    h5_full_path = h5_details[0]['h5_file_path']
    scanType = h5_details[0]['scanType']
    logger.info(f"Scan Type: {scanType}")
    recordings = {}
    rec_counts = []

    try:
        recording = se.read_maxwell(h5_full_path, rec_name='rec0000', stream_id='well000')
        MaxID = 2
    except:
        try:
            recording = se.read_maxwell(h5_full_path, rec_name='rec0000')
            MaxID = 1
        except:
            logger.error("Error: Unable to read recording. Cannot identify as MaxOne or MaxTwo.")
            raise Exception("Error: Unable to read recording. Cannot identify as MaxOne or MaxTwo.")
            #return 0, [], 0, []

    if MaxID == 1:
        logger.info("MaxOne Detected.")
        expected_rec_count = len(h5py.File(h5_full_path)['recordings'])
        recs = []
        for rec_count in range(expected_rec_count):
            try:
                rec_name_str = f'rec{rec_count:04}'
                recording = se.read_maxwell(h5_full_path, rec_name=rec_name_str)
                recs.append(recording)
            except:
                continue
        recordings['well000'] = {'recording_segments': recs}
        return MaxID, recordings, 1, [expected_rec_count]

    elif MaxID == 2:
        logger.info("MaxTwo Detected.")
        with h5py.File(h5_full_path, 'r') as h5_file:
            expected_well_count = len(h5_file['wells'])
            streams_to_process = range(expected_well_count) if stream_select is None else [stream_select]
            for stream_id in streams_to_process:
                stream_id_str = f'well{stream_id:03}'
                expected_rec_names = list(h5_file['wells'][stream_id_str].keys())
                expected_rec_count = len(expected_rec_names)
                valid_rec_count = 0
                recs = []
                for rec_count in range(expected_rec_count):
                    try:
                        rec_name_str = f'rec{rec_count:04}'
                        if rec_name_str not in expected_rec_names:
                            continue
                        recording = se.read_maxwell(h5_full_path, rec_name=rec_name_str, stream_id=stream_id_str)
                        recs.append(recording)
                        valid_rec_count += 1
                        logger.debug(f"Stream ID: {stream_id_str}, Recording: {rec_name_str} loaded.")
                    except Exception as e:
                        logger.warning(e)
                        if "Unable to open object" in str(e):
                            logger.debug("This error may not be an issue. For some reason, some wells say they have segments when they don't.")
                rec_counts.append(valid_rec_count)
                recordings[stream_id_str] = {'recording_segments': recs}
                if stream_select is not None:
                    break
        return MaxID, recordings, expected_well_count, rec_counts


def count_wells_and_segs(h5_file_path, verbose=False, stream_select=None, recordings=None, MaxID=None):
    """
    This function counts the number of wells and recordings by utilizing the load_recordings function.

    Parameters:
    h5_file_path (str): The path to the file to read from.
    verbose (bool): If True, log detailed output. Default is False.
    stream_select (int, optional): Specific stream ID to select. Default is None.

    Returns:
    tuple: 
    (1) number of stream IDs detected 
    (2) and an array containing number of recording segments per stream ID.
    """
    if MaxID is None and recordings is None: MaxID, recordings = load_recordings(h5_file_path, stream_select)

    if MaxID == 0:
        return 0, 0

    stream_ids = []
    rec_counts = []

    if MaxID == 1:
        logger.info("MaxOne Detected.")    
        #Maxone
        while True:
            try:
                rec_name_str = f'rec{rec_count:04}'
                recording = se.read_maxwell(h5_file_path, rec_name=rec_name_str)
                #When counting recordings, Network scans dont require a rec_name and will loop endlessly if not handled
                if "Network" in h5_file_path:
                    assert rec_count == 0, "rec_count should be 0 before incrementing when 'network' is in h5_file_path"
                    rec_count = 1
                else:
                    rec_count += 1
            except:            
                rec_counts.append(rec_count)
                if rec_count == 0:
                    break            
                if verbose:
                    logger.info(f"{rec_count} recordings detected.")
        return 1, rec_counts if rec_counts else 0

    elif MaxID == 2:
        logger.info("MaxTwo Detected.")
        #MaxTwo
        while True:
            if stream_select is not None and stream_id != stream_select: stream_id += 1; continue #skip if stream_select is not None and stream_id is not stream_select, saves time on loading and verifying data
            try:
                stream_id_str = f'well{stream_id:03}'
                rec_name_str = f'rec{rec_count:04}'
                recording = se.read_maxwell(h5_file_path, rec_name=rec_name_str, stream_id=stream_id_str)
                #When counting recordings, Network scans dont require a rec_name and will loop endlessly if not handled
                if "Network" in h5_file_path:
                    assert rec_count == 0, "rec_count should be 0 before incrementing when 'network' is in h5_file_path"
                    rec_count = 1
                else: rec_count += 1
            except:            
                rec_counts.append(rec_count)
                if rec_count == 0: break            
                stream_ids.append(stream_id_str)
                if verbose: logger.info(f"Stream ID: {stream_id_str}, {rec_count} recordings detected.")
                rec_count = 0
                stream_id += 1
                if stream_select is not None: break #break if stream_select is not None
        if len(set(rec_counts)) > 1: logger.error(f"Warning: The number of recordings is not consistent across all stream IDs. Range: {min(rec_counts)}-{max(rec_counts)}")
        else: logger.info(f"Recordings per Stream ID: {rec_counts[0]}")
        logger.info(f"Stream IDs Detected: {len(stream_ids)}")    
        return len(stream_ids), rec_counts if rec_counts else 0

def count_wells_and_recs_old(h5_file_path, verbose=False, stream_select = None, recordings = None):
    """
    This function counts the number of wells (stream IDs) in a given file path by attempting to read recordings 
    from the file until an error occurs. It also counts the number of successful reads (recordings) for each stream ID.

    If the verbose flag is set to True, the function logs the stream ID and the number of recordings detected 
    for each stream ID.

    If the number of recordings is not consistent across all stream IDs, the function logs a warning and 
    the range of recordings detected. Otherwise, it logs the common number of recordings per stream ID.

    Finally, the function returns the total number of stream IDs detected and the number of recording segments 
    per stream ID.

    Parameters:
    h5_file_path (str): The path to the file to read from.
    verbose (bool): If True, log detailed output. Default is False.

    Returns:
    tuple: 
    (1) number of stream IDs detected 
    (2) and an array containing number of recording segments per stream ID.
    """
    logger.info(f"Counting wells and recordings in {h5_file_path}:")
    h5_details = extract_recording_details(h5_file_path)
    scanType = h5_details[0]['scanType']
    logger.info(f"Scan Type: {scanType}")
    stream_ids = []
    rec_counts = []
    stream_id = 0
    rec_count = 0

    #Check if MaxOne or MaxTwo
    stream_id_str = f'well{stream_id:03}'
    rec_name_str = f'rec{rec_count:04}'
    try: 
        recording = se.read_maxwell(h5_file_path, rec_name=rec_name_str, stream_id=stream_id_str)
        MaxID = 2
    except:
        try:
            recording = se.read_maxwell(h5_file_path, rec_name=rec_name_str)
            MaxID = 1
        except:
            logger.error("Error: Unable to read recording. Cannot identify as MaxOne or MaxTwo.")
            return 0, 0
    
    if MaxID == 1:
        logger.info("MaxOne Detected.")    
        #Maxone
        while True:
            try:
                rec_name_str = f'rec{rec_count:04}'
                recording = se.read_maxwell(h5_file_path, rec_name=rec_name_str)
                #When counting recordings, Network scans dont require a rec_name and will loop endlessly if not handled
                if "Network" in h5_file_path:
                    assert rec_count == 0, "rec_count should be 0 before incrementing when 'network' is in h5_file_path"
                    rec_count = 1
                else:
                    rec_count += 1
            except:            
                rec_counts.append(rec_count)
                if rec_count == 0:
                    break            
                if verbose:
                    logger.info(f"{rec_count} recordings detected.")
        return 1, rec_counts if rec_counts else 0

    elif MaxID == 2:
        logger.info("MaxTwo Detected.")
        #MaxTwo
        while True:
            if stream_select is not None and stream_id != stream_select: stream_id += 1; continue #skip if stream_select is not None and stream_id is not stream_select, saves time on loading and verifying data
            try:
                stream_id_str = f'well{stream_id:03}'
                rec_name_str = f'rec{rec_count:04}'
                recording = se.read_maxwell(h5_file_path, rec_name=rec_name_str, stream_id=stream_id_str)
                #When counting recordings, Network scans dont require a rec_name and will loop endlessly if not handled
                if "Network" in h5_file_path:
                    assert rec_count == 0, "rec_count should be 0 before incrementing when 'network' is in h5_file_path"
                    rec_count = 1
                else: rec_count += 1
            except:            
                rec_counts.append(rec_count)
                if rec_count == 0: break            
                stream_ids.append(stream_id_str)
                if verbose: logger.info(f"Stream ID: {stream_id_str}, {rec_count} recordings detected.")
                rec_count = 0
                stream_id += 1
                if stream_select is not None: break #break if stream_select is not None
        if len(set(rec_counts)) > 1: logger.error(f"Warning: The number of recordings is not consistent across all stream IDs. Range: {min(rec_counts)}-{max(rec_counts)}")
        else: logger.info(f"Recordings per Stream ID: {rec_counts[0]}")
        logger.info(f"Stream IDs Detected: {len(stream_ids)}")    
        return len(stream_ids), rec_counts if rec_counts else 0
    
def count_recording_segments(recording, verbose=True):
    """
    This function counts the number of recording segments in a given RecordingExtractor object by using the 
    get_num_segments method. 

    If the verbose flag is set to True, the function prints the number of recording segments detected.

    Finally, the function returns the total number of recording segments detected.

    Parameters:
    recording (RecordingExtractor): The RecordingExtractor object to count segments from.
    verbose (bool): If True, print detailed output. Default is True.

    Returns:
    int: The number of recording segments detected.
    """
    logger.info(f"Counting Segments in Recording:")
    num_segments = recording.get_num_segments()

    if verbose:
        logger.info(f"Number of segments in the recording: {num_segments}")

    return num_segments

def get_data_maxwell(h5_file_path, rec_num, well_num = None, verbose = False):
    """
    This function reads a recording from a given file path using the read_maxwell method from the RecordingExtractor object. 

    The function constructs the recording name and stream ID based on the provided recording number and well number. 
    If no well number is provided, the stream ID is set to None.

    If the verbose flag is set to True, the function prints detailed output.

    If an exception occurs while reading the recording, the function prints an error log and returns None for the recording.

    Parameters:
    file_path (str): The path to the file to read from.
    rec_num (int): The number of the recording to read from.
    well_num (int, optional): The number of the well to read from. Default is None.
    verbose (bool): If True, print detailed output. Default is False.

    Returns:
    tuple: A tuple containing the RecordingExtractor object, the recording name, and the stream ID. If an exception occurs, the RecordingExtractor object is None.
    """
    logger.info(f"Extracting Recording Object from h5_file:")
    rec_num =  str(rec_num).zfill(4)
    rec_name = 'rec' + rec_num
    stream_id='well' + str(well_num).zfill(3) if well_num is not None else None
    recording = None
    try:
        if well_num is not None:
            recording = se.read_maxwell(h5_file_path,rec_name=rec_name, stream_id=stream_id)
            if verbose: logger.info(f"Successfully read recording from well {well_num}.")
        else:
            recording = se.read_maxwell(h5_file_path,rec_name=rec_name)
            if verbose: logger.info("Successfully read recording.")
    except Exception as e:
        logger.error(f"Failed to read recording. Exception: {e}")
        return e, rec_name, stream_id
    return recording, rec_name, stream_id

def merge_mea_recording_segments(recordings, mode='concatenate'):
    """
    recordings: List of recording objects
    mode: Either 'append' or 'concatenate'. Determines how the recordings are merged.
    """
    logger.info(f"Merging Recording Objects ({mode}):")
    try:
        #(untested) for same set of channels in each recording, generate list of segments
        if mode == 'append':
            merged_recording = spikeinterface.append_recordings(recordings)
        #(untested) for same set of channels in each recording, generate single segment
        elif mode == 'concatenate':
            merged_recording = spikeinterface.concatenate_recordings(recordings)
        #for different sets of channels in each recording
        elif mode == 'aggregate':
            merged_recording = si.aggregate_channels(recordings)
        else:
            logger.error("Invalid mode. Must be either 'append' or 'concatenate'.")
    except Exception as e:
        merged_recording = None
        logger.error(f"Failed to merge recordings. Exception: {e}")

    return merged_recording

def merge_sortings(sortings, mode='aggregate'):
    """
    recordings: List of recording objects
    mode: Either 'append' or 'concatenate'. Determines how the recordings are merged.
    """
    logger.info(f"Merging Sorting Objects ({mode}):")
    try:
        #for different sets of channels in each recording
        if mode == 'aggregate':
            merged_sorting = si.aggregate_units(sortings)
        else:
            logger.error("Invalid mode. Must be either 'append' or 'concatenate'.")
    except Exception as e:
        merged_sorting = None
        logger.error(f"Failed to merge sortings. Exception: {e}")

    return merged_sorting

def get_channel_recording_stats(recording):
    """
    This function retrieves various statistics about a given recording.

    Parameters:
    recording: The recording object to get statistics for.

    Returns:
    fs: The sampling frequency of the recording.
    num_chan: The number of channels in the recording.
    channel_ids: The IDs of the channels in the recording.
    total_recording: The total duration of the recording in seconds.
    """
    logger.info(f"Getting Channel Recording Characteristics:")
    channel_ids = recording.get_channel_ids()
    fs = recording.get_sampling_frequency()
    num_chan = recording.get_num_channels()
    num_seg = recording.get_num_segments()
    total_recording = recording.get_total_duration()

    logger.info(f'Sampling frequency: {fs}')
    logger.info(f'Number of channels: {num_chan}')
    logger.info(f'Number of segments: {num_seg}')
    logger.info(f"Total recording duration: {total_recording} s")

    return fs, num_chan, channel_ids, total_recording

def preprocess_recording(recording): #AW25Jan24 - some hardcoded stuff here, discuss with Mandar later
    """
    This function performs preprocessing on a given recording. The preprocessing steps include:
    1. Bandpass filtering: This removes frequencies outside the range of 300 Hz to half the sampling frequency minus 1.
    2. Common median referencing (CMR): This is a technique used to reduce common noise sources. It works by subtracting the median of all channels from each individual channel.

    Parameters:
    recording: The recording object to preprocess.

    Returns:
    recording_cmr: The preprocessed recording object. It has been bandpass filtered and common median referenced.
    """
    recording_bp = spre.bandpass_filter(recording, freq_min=300, freq_max=(recording.sampling_frequency/2)-1)
    recording_cmr = spre.common_reference(recording_bp, reference='global', operator='median')
    recording_cmr.annotate(is_filtered=True)

    return recording_cmr

def prepare_recordings_for_merge(recording_list):
    """
    This function prepares a list of recordings for merging by getting a chunk of each recording and preprocessing it.
    It then performs quality checks to ensure that all recordings have the same sampling frequency, number of segments, data type, and number of samples.

    Parameters:
    recording_list: The list of recording objects to prepare for merging.

    Returns:
    recordings_to_merge: The list of preprocessed chunks of the recordings.
    """
    recordings_to_merge = []
    for recording in recording_list:
        # Get the recording statistics such as sampling frequency, number of channels, channel IDs, and total recording time
        fs, num_chan, channel_ids, total_rec_time = get_channel_recording_stats(recording)

        # Round the total recording time to the nearest whole number
        rounded_total_rec_time = round(total_rec_time)
        if total_rec_time > rounded_total_rec_time:
            time_in_s = rounded_total_rec_time
        else:
            time_in_s = total_rec_time

        # Define the start and end times for the recording chunk
        time_start = 0
        time_end = time_start + time_in_s

        # Get a chunk of the recording based on the start and end times and preprocess it
        recording_chunk = recording.frame_slice(start_frame=int(time_start * fs), end_frame=int(time_end * fs))
        recording_chunk = preprocess_recording(recording_chunk)
        recordings_to_merge.append(recording_chunk)

    # Quality checks
    fs = recordings_to_merge[0].get_sampling_frequency()
    num_segments = recordings_to_merge[0].get_num_segments()
    dtype = recordings_to_merge[0].get_dtype()

    ok1 = all(fs == rec.get_sampling_frequency() for rec in recordings_to_merge)
    ok2 = all(num_segments == rec.get_num_segments() for rec in recordings_to_merge)
    ok3 = all(dtype == rec.get_dtype() for rec in recordings_to_merge)
    ok4 = True
    for i_seg in range(num_segments):
        num_samples = recordings_to_merge[0].get_num_samples(i_seg)
        ok4 = all(num_samples == rec.get_num_samples(i_seg) for rec in recordings_to_merge)
        if not ok4:
            break

    if not (ok1 and ok2 and ok3 and ok4):
        raise ValueError("Recordings don't have the same sampling_frequency/num_segments/dtype/num samples")

    return recordings_to_merge

def preprocess_single_recording(recording):
    """
    This function prepares a single recording for merging by getting a chunk of the recording and preprocessing it.

    Parameters:
    recording: The recording object to prepare for merging.

    Returns:
    recording_chunk: The preprocessed chunk of the recording.
    """
    # Get the recording statistics such as sampling frequency, number of channels, channel IDs, and total recording time
    fs, num_chan, channel_ids, total_rec_time = get_channel_recording_stats(recording)

    # Round the total recording time to the nearest whole number
    rounded_total_rec_time = round(total_rec_time)
    if total_rec_time > rounded_total_rec_time:
        time_in_s = rounded_total_rec_time
    else:
        time_in_s = total_rec_time

    # Define the start and end times for the recording chunk
    time_start = 0
    time_end = time_start + time_in_s

    # Get a chunk of the recording based on the start and end times and preprocess it
    recording_chunk = recording.frame_slice(start_frame=int(time_start * fs), end_frame=int(time_end * fs))
    recording_chunk = preprocess_recording(recording_chunk)

    return recording_chunk

def run_kilosort2_docker_image(recording, output_folder, docker_image="spikeinterface/kilosort2-compiled-base:latest", verbose=False, logger=None):
    """
    Runs Kilosort2 sorter on the provided recording in chunks using a Docker image.

    Parameters:
    recording (RecordingExtractor): The recording extractor object.
    chunk_duration (float): Duration of each chunk in seconds.
    output_folder (str): The folder to save the sorting output.
    docker_image (str, optional): Docker image to use for Kilosort2. Default is "spikeinterface/kilosort2-compiled-base:latest".
    verbose (bool, optional): If True, enables verbose logging. Default is False.
    logger (Logger, optional): Logger for logging. Default is None.

    Returns:
    SortingExtractor: The sorting extractor object with the spike sorting results.
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    try:
        default_KS2_params = ss.Kilosort2Sorter.default_params()
        sorting = ss.run_kilosort2(recording, output_folder=output_folder, docker_image=docker_image, verbose=verbose, **default_KS2_params)
        logger.info("Kilosort2 processing complete")
        return sorting

    except Exception as e:
        logger.error(f"Error running Kilosort2 on recording: {e}")
        return None

# import os
# import shutil
# import spikeinterface.sorters as ss

def benshalom_kilosort2_docker_image(recording, output_folder, sorting_params=None, verbose = False):
    """
    Run Kilosort2 spike sorting using Docker.

    Parameters:
        recording: Recording object from SpikeInterface.
        output_folder: Path to the folder where the output will be saved.
        sorter_params: Dictionary of sorter parameters.

    Returns:
        sorting: Sorting object with the results.
    """
    docker_image = "rohanmalige/benshalom:v3"
    #verbose = sorting_params.pop('verbose', False)

    if sorting_params is None:
        sorting_params = ss.Kilosort2Sorter.default_params()
        sorting_params.update({
            'n_jobs': -1,
            'detect_threshold': 7,
            'minFR': 0.01,
            'minfr_goodchannels': 0.01,
            'keep_good_only': False,
            'do_correction': False
        })
    else:
        default_params = ss.Kilosort2Sorter.default_params()
        default_params.update(sorting_params)
        sorting_params = default_params

    #verbose = sorter_params.pop('verbose', False)

    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

    logger.info(f"Running Kilosort2 spike sorting using Docker Images: {docker_image}.")
    try:
        sorting = ss.run_sorter(
            sorter_name="kilosort2",
            recording=recording,
            output_folder=output_folder,
            docker_image=docker_image,
            verbose=verbose,
            **sorting_params
        )
    except Exception as e:
        logger.error(f"Error running Kilosort2: {e}")
        return None

    return sorting



def run_kilosort2_5_docker_image_GPUs(recording, output_folder, docker_image="spikeinterface/kilosort2_5-compiled-base:latest", verbose=False, num_gpus=1):
    # Update sorter parameters as needed
    sorter_params = si.get_default_sorter_params(si.Kilosort2_5Sorter)
    sorter_params['n_jobs'] = -1
    sorter_params['detect_threshold'] = 7
    sorter_params['minFR'] = 0.01
    sorter_params['minfr_goodchannels'] = 0.01
    sorter_params['keep_good_only'] = False
    sorter_params['do_correction'] = False

    # Add extra kwargs for GPU allocation
    extra_kwargs = {
        'docker_args': {
            'device_requests': [docker.types.DeviceRequest(count=num_gpus, capabilities=[['gpu']])]
        }
    }

    # Run Kilosort2.5
    sorting_KS2_5 = ss.run_kilosort2_5(
        recording,
        output_folder=output_folder,
        docker_image=docker_image,
        verbose=verbose,
        **sorter_params,
        #**extra_kwargs
    )

    return sorting_KS2_5

def extract_waveforms(recording,sorting,folder, load_if_exists = False, n_jobs = 4, sparse = True):
    job_kwargs = dict(n_jobs=n_jobs, chunk_duration="1s", progress_bar=True)
    #waveforms = si.extract_waveforms(recording,sorting_KS3,folder=folder,overwrite=True,**job_kwargs)
    waveforms = si.extract_waveforms(recording,sorting,folder=folder,overwrite=True, load_if_exists=load_if_exists, sparse = sparse, ms_before=1., ms_after=2.,allow_unfiltered=True,**job_kwargs)
    return waveforms

#This function supports generate_waveform_extractor_unit_by_unit
def copy_and_compare_files(args):
    old_name, new_name, file_pattern, unit_id, unit_folder, waveforms_folder, overwrite, log_messages = args
    do_very_shallow = True
    do_shallow = False
    do_deep = False
    if os.path.exists(new_name) and overwrite == False:
        try:                     

            #Very Shallow comparison
            # Compare file sizes
            if do_very_shallow:
                old_size = os.path.getsize(old_name)
                new_size = os.path.getsize(new_name)
                if old_size == new_size and do_shallow == False and do_deep == False:
                    log_messages.append(f"{new_name} already exists - skipping based on file size")
                    return
                elif old_size != new_size:
                    raise Exception(f"File sizes do not match for {old_name} and {new_name} based on file size alone")
                
            # Shallow comparisons
            #use filecmp.cmp to compare files, shallow = True                
            if do_shallow:
                #logger.debug(f"Shallow comparison of {old_name} and {new_name}")  
                cmp = filecmp.cmp(old_name, new_name, shallow=True)
                if cmp and do_deep == False and do_very_shallow == True:
                    log_messages.append(f"{new_name} already exists - skipping based on shallow comparison and file size")
                    return
                if cmp and do_deep == False and do_very_shallow == False:
                    log_messages.append(f"{new_name} already exists - skipping based on shallow comparison alone")
                    return
                elif not cmp:
                    raise Exception(f"Files do not match for {old_name} and {new_name} based on shallow comparison")
            
            #Actually, if full comparisons are needed, may as well just copy and overwrite the file, leaving this off for now
            # Full comparisons if all else fails.
            # Load the files as numpy arrays
            if do_deep:
                old_array = np.load(old_name)
                new_array = np.load(new_name)
                
                # Check if the shapes of the arrays are the same
                if old_array.shape != new_array.shape:
                    raise Exception(f"Shapes do not match for {old_name} and {new_name} based on deep comparison")
                
                # Check if the contents of the arrays are the same
                if np.array_equal(old_array, new_array):
                    log_messages.append(f"{new_name} already exists - skipping based on full comparison")
                    return
                else:
                    raise Exception(f"Contents do not match for {old_name} and {new_name} based on deep comparison")
        except:
            pass
    shutil.copy(old_name, new_name)
    #logger.info(f"{file_pattern[:-5]} for unit {unit_id} copied from {unit_folder} to {waveforms_folder}")
    log_messages.append(f"{file_pattern[:-5]} for unit {unit_id} copied from {unit_folder} to {waveforms_folder}")

#@profile
def generate_waveform_extractor_unit_by_unit(recording,sorting,folder, n_jobs = 4, units_per_extraction = 20, sparse = True, fresh_extractor = False, load_if_exists = False, job_kwargs = None):

    """
    This function generates a waveform extractor by processing units of data in chunks. 

    Parameters:
    recording: The recording from which to extract waveforms.
    sorting: The sorting of the recording.
    folder: The directory where the extracted waveforms will be stored.
    n_jobs: The number of jobs to run in parallel.
    units_per_extraction: The number of units to extract at a time.
    sparse: A boolean indicating whether to use sparse extraction.
    fresh_extractor: A boolean indicating whether to generate a new waveform extractor.
    load_if_exists: A boolean indicating whether to load the waveform extractor if it already exists.
    job_kwargs: Additional keyword arguments for the job.

    The function operates as follows:

    - If `load_if_exists` is `True`, it attempts to load an existing waveform extractor from the specified folder. If successful, it returns the extractor immediately.
    - If `fresh_extractor` is `True` or the specified folder does not exist, a new waveform extractor is generated from the full recording and sorting.
    - The function then processes units of data in chunks, defined by `units_per_extraction`, and extracts waveforms for these units.
    - Extracted waveforms from all unit folders are merged into a single folder.
    - The function verifies that the waveforms and sampled index files are contiguous. If they are not, it raises an exception.
    - The function attempts to load the full waveform extractor. If it fails to load, it raises an exception.
    - Temporary unit folders are deleted.
    - The waveform extractor is returned.
    """
    
    def try_load_waveform_extractor(folder, sorting):       

        try:
            # Load the waveform extractor
            #waveform_extractor = spikeinterface.load_extractor(folder)
            waveform_extractor = si.load_waveforms(folder)

            # Get the unit IDs from the waveform extractor and the sorting
            we_unit_ids = waveform_extractor.unit_ids
            sorting_unit_ids = sorting.get_unit_ids()

            # If the number of unit IDs is the same, return
            if all(a == b for a, b in zip(we_unit_ids, sorting_unit_ids)):
                logger.info(f"Waveform extractor loaded")
                return waveform_extractor
        except:
            logger.info(f"Waveform extractor not found")
            return None
        
    def check_extraction_status(unit_range_folder, unit_waveforms=None):
        extraction_info_file = unit_range_folder + "/extraction_info.json"
        # Check if the extraction info file exists
        if os.path.exists(extraction_info_file):
            with open(extraction_info_file, "r") as json_file:
                extraction_info = json.load(json_file)
                if extraction_info["status"] == "extraction successful":
                    match = re.search(r"units?_(\d+)(_to_(\d+))?", unit_range_folder)
                    if match:
                        start_unit = match.group(1)
                        end_unit = match.group(3) if match.group(3) else start_unit
                        logger.info(f"Previous extraction for units {start_unit} to {end_unit} was successful")
                    else:
                        logger.warning(f"Could not extract unit range from folder name: {unit_range_folder}")
                else:
                    raise Exception("Extraction was not successful")
        else:
            raise Exception("Extraction info file does not exist")
        if unit_waveforms is not None and unit_waveforms.get_num_segments() == 0:
            raise Exception("Waveforms are empty")
        return
    
    def flex_load_waveforms(unit_ids, folder, tot_num_units):      
        
        unit_id = unit_ids[0]
        while unit_id in unit_ids:
            # Try to load the widest range of units first
            for j in range(tot_num_units, unit_id, -1):
                unit_range_folder = folder + f"_unit_by_unit_temp/units_{unit_id}_to_{j}"
                if os.path.exists(unit_range_folder):
                    try:
                        unit_waveforms = si.load_waveforms(unit_range_folder)
                        # Extraction quality testing
                        check_extraction_status(unit_range_folder, unit_waveforms)
                        logger.info(f"Waveforms for units {unit_id} to {j} loaded from {unit_range_folder}")
                        unit_id = j + 1  # Update the next unit in sequence                        
                        return unit_waveforms, unit_id
                    except:
                        pass  # If loading fails, try the next range
            else:
                # If no range works, try to load a single unit
                unit_folder = folder + f"_unit_by_unit_temp/unit_{unit_id}"
                if os.path.exists(unit_folder):
                    try:
                        unit_waveforms = si.load_waveforms(unit_folder)
                       # Extraction quality testing
                        check_extraction_status(unit_range_folder, unit_waveforms)      
                        logger.info(f"Waveform for unit {unit_id} loaded from {unit_folder}")
                        unit_id += 1  # Update the next unit in sequence
                        return unit_waveforms, unit_id
                    except:
                        raise Exception(f"Failed to load waveforms for unit {unit_id}")
                else:
                    raise Exception(f"No folder found for unit {unit_id}")
        return unit_waveforms, unit_id
    
    def extract_waveforms_from_unit_range(unit_range, recording, sorting, folder, sparse, job_kwargs):
        unit_range_folder = folder + f"_unit_by_unit_temp/units_{unit_range[0]}_to_{unit_range[-1]}"
        logger.info(f"Extracting waveforms for units {unit_range[0]} to {unit_range[-1]} to {unit_range_folder}")
        unit_sorting = sorting.select_units(unit_range)
        unit_waveforms = si.extract_waveforms(recording, 
                                              unit_sorting, 
                                              folder=unit_range_folder,
                                              precompute_template=None, 
                                              overwrite=True, 
                                              load_if_exists=False, 
                                              sparse=sparse, 
                                              ms_before=1., 
                                              ms_after=2., 
                                              allow_unfiltered=True, 
                                              **job_kwargs)        
        # Create a JSON file to confirm successful extraction
        extraction_info = {
            "units": unit_range.tolist(),  # convert numpy array to list
            "folder": unit_range_folder,
            "status": "extraction successful"
        }
        with open(unit_range_folder + "/extraction_info.json", "w") as json_file:
            json.dump(extraction_info, json_file)
        logger.info(f"Waveforms extracted")
        return unit_waveforms
        
    def extract_waveforms_for_unit_ranges_to_temp(sorting, folder, n_jobs, units_per_extraction, load_if_exists, recording, sparse, job_kwargs):
        
        #Define default job_kwargs
        if job_kwargs is None:
            job_kwargs = dict(
                n_jobs=n_jobs, 
                chunk_duration="1s", 
                progress_bar=True
                )
        
        # Extract waveforms for each unit
        unit_folders = []        
        unit_ids = sorting.get_unit_ids()
        # Process the units in chunks
        i = 0
        tot_num_units = len(unit_ids)
        logger.info(f"Total number of units: {tot_num_units}")
        while i < tot_num_units:
            # Select up to units_per_extraction units
            unit_range = unit_ids[i:i+units_per_extraction]           
            if load_if_exists:
                try:
                    # Flexibly try to load the waveforms from the folder                
                    unit_waveforms, unit_id = flex_load_waveforms(unit_range, folder, tot_num_units)
                    #update i to the next unit in sequence
                    i = unit_id
                except:
                    unit_waveforms = extract_waveforms_from_unit_range(unit_range, recording, sorting, folder, sparse, job_kwargs)
                    i += units_per_extraction
            else:
                unit_waveforms = extract_waveforms_from_unit_range(unit_range, recording, sorting, folder, sparse, job_kwargs)
                i += units_per_extraction
        return unit_waveforms

    def create_waveform_extractor(recording, sorting, folder):
        # Remove the folder if it already exists
        if os.path.exists(folder):
            shutil.rmtree(folder)

        # Create a WaveformExtractor 
        waveform_extractor = sc.WaveformExtractor.create(
            recording, sorting, folder=folder, allow_unfiltered=True
        )

        # Set the parameters for the waveform extractor
        waveform_extractor.set_params(
            ms_before=1.,
            ms_after=2.,
        )
        
        return

    def merge_waveforms(folder, temp_units_folder, n_jobs, overwrite=True):
        # Create a waveforms folder at the same level as the unit folders
        waveforms_folder = os.path.join(folder, "waveforms")
        if not os.path.exists(waveforms_folder):
            os.makedirs(waveforms_folder)

        # From all unit_folder in unit_folders, move and merge the waveforms folder with higher waveforms folder we just created
        unit_folders = [os.path.join(temp_units_folder, name) for name in os.listdir(temp_units_folder) if os.path.isdir(os.path.join(temp_units_folder, name))]
        #sort the unit_folders
        unit_folders = sorted(unit_folders, key=lambda name: int(''.join(filter(str.isdigit, name))))

        # Create a manager object
        manager = Manager()

        # Create a list that all processes can write to
        log_messages = manager.list()
        
        # Create a pool of worker processes
        logger.info(f"Copying and merging waveforms folders from unit folders to template folder.")
        n_processes = min(n_jobs, len(unit_folders))
        n_processes = 1
        with Pool(processes=n_processes) as pool:
            # Prepare the arguments for each task
            tasks = []
            for unit_folder in unit_folders:
                unit_waveforms_folder = os.path.join(unit_folder, "waveforms")
                if os.path.exists(unit_waveforms_folder):
                    try: 
                        check_extraction_status(unit_folder)
                    except: 
                        logger.error(f"Invalid or incomplete waveforms folder found in {unit_folder} - skipping")
                        continue

                    # Copy the params.json file if it doesn't exist in the destination folder
                    if not os.path.exists(os.path.join(folder, "params.json")):
                        shutil.copy(os.path.join(unit_folder, "params.json"), folder)
                        logger.info(f"params.json file copied to {folder}")                

                    # Prepare the tasks for the sampled_index_i.npy and waveforms_i.npy files
                    for file_pattern in ["sampled_index_*.npy", "waveforms_*.npy"]:
                        for old_name in glob.glob(os.path.join(unit_waveforms_folder, file_pattern)):
                            # Get the unit ID from the unit folder name
                            unit_id = os.path.basename(old_name).split('_')[1]
                            new_name = os.path.join(waveforms_folder, os.path.basename(old_name))
                            tasks.append((old_name, new_name, file_pattern, unit_id, unit_folder, waveforms_folder, overwrite, log_messages))

            # Perform the tasks in parallel
            #debug
            #tasks = tasks[:10]
            #debug
            results = []
            for _ in tqdm(pool.imap_unordered(copy_and_compare_files, tasks), total=len(tasks)):
                results.append(_)

        # Log all the messages
        for message in log_messages:
            logger.info(message)

        logger.info(f"Unit waveforms copied. Full waveform extractor generated.")

    def verify_waveforms_contiguity(sorting, folder):
        unit_ids = sorting.get_unit_ids()
        waveforms_folder = os.path.join(folder, "waveforms")
        if not os.path.exists(waveforms_folder):
            raise Exception(f"Waveforms folder not found in {folder}, cannot verify contiguity")

        waveforms_files = [os.path.basename(file) for file in glob.glob(waveforms_folder + "/waveforms_*.npy")]
        sampled_index_files = [os.path.basename(file) for file in glob.glob(waveforms_folder + "/sampled_index_*.npy")]
        #waveforms_files.sort()
        #sampled_index_files.sort()
        if len(waveforms_files) == len(sampled_index_files) == len(unit_ids):
            missing_waveforms = [i for i in unit_ids if f"waveforms_{i}.npy" not in waveforms_files]
            missing_sampled_index = [i for i in unit_ids if f"sampled_index_{i}.npy" not in sampled_index_files]

            if not missing_waveforms and not missing_sampled_index:
                logger.info(f"Waveforms and sampled index files are contiguous.")
                return True
            else:
                if missing_waveforms:
                    logger.info(f"Missing waveforms files for unit IDs: {missing_waveforms}")
                if missing_sampled_index:
                    logger.info(f"Missing sampled index files for unit IDs: {missing_sampled_index}")
                return False
        else:
            logger.error(f"Waveforms and sampled index files are not contiguous.")
            return False

    #Begin the function
    logger.info(f"Extracting waveforms unit by unit:")
    
    # Check if the waveform extractor already exists
    logger.info(f"load_if_exists set to {load_if_exists}")
    if load_if_exists:        
        logger.info(f"Trying to load waveform extractor from {folder}")
        # Create a WaveformExtractor
        waveform_extractor = try_load_waveform_extractor(folder, sorting)
        if waveform_extractor is None:
            # Handle the case where the waveform extractor could not be loaded
            pass
        else: 
            return waveform_extractor

    # Extract waveforms for ranges of units
    # Choose the number of units to extract at a time
    units_per_extraction = int(units_per_extraction)
    logger.info(f"Extracting waveforms for up to {units_per_extraction} units at a time")
    extract_waveforms_for_unit_ranges_to_temp(sorting, folder, n_jobs, units_per_extraction, load_if_exists, recording, sparse, job_kwargs)
        
    #Generate a waveform extractor from full recording and sorting    
    if fresh_extractor:
        logger.info(f"fresh_extractor set to {fresh_extractor}, clearing folder {folder} and generating new waveform extractor.")
        logger.info(f"Generating new waveform extractor template from full recording and sorting.")
        create_waveform_extractor(recording, sorting, folder)
    elif os.path.exists(folder) == False:
        logger.info(f"Generating new waveform extractor template from full recording and sorting.")
        create_waveform_extractor(recording, sorting, folder)
    
    # from all unit_folder in unit_folders, move and merge the waveforms folder with higher waveforms folder we just created
    logger.info(f"Copying and merging waveforms folders from unit folders to template folder.")
    temp_units_folder = folder + f"_unit_by_unit_temp"
    overwrite = False
    merge_waveforms(folder, temp_units_folder, n_jobs, overwrite)
    
    #verify that waveforms_* and sampled_index_* files are present and contiguous
    # Create a waveforms folder at the same level as the unit folders
    logger.info(f"Verifying that waveforms and sampled index files are contiguous.")
    contiguity_status = verify_waveforms_contiguity(sorting, folder)
    if not contiguity_status:
        raise Exception("Waveforms and sampled index files are not contiguous")
    logger.info(f"Contiguity verified.")
    
    #Attempting to load the full waveform extractor
    logger.info(f"Attempting to load the full waveform extractor.") 
    # Create a WaveformExtractor
    waveform_extractor = try_load_waveform_extractor(folder, sorting)
    if waveform_extractor is None:
        raise Exception(f"Failed to load waveform extractor from {folder}") 
    
    logger.info(f"Deleting temporary unit folders.")
    # Remove unit folders and all their contents
    temp_units_folder = folder + f"_unit_by_unit_temp"
    shutil.rmtree(temp_units_folder)
    logger.info(f"Temporary unit folder {temp_units_folder} deleted.")
    logger.info(f"Waveform extractor successfully generated in chunks!")
    return waveform_extractor

#Waveform post-processing:
def get_quality_metrics(waveforms): 
    
    #Define default job_kwargs
    job_kwargs = dict(
        n_jobs=4, 
        chunk_duration="1s", 
        progress_bar=True
        )
    
    #TO DO: to fix the SNR showing as INF by looking at the extensions of waveforms
    # This is because the noise levels are 0.0
    # Similar issue with https://github.com/SpikeInterface/spikeinterface/issues/1467 
    # Best solution proposed was to stick with applying a gain to the entire recording before all preprocessing
    sp.compute_spike_amplitudes(waveforms,**job_kwargs)
    metrics = qm.compute_quality_metrics(waveforms, metric_names=['num_spikes','firing_rate', 'presence_ratio', 'snr',
                                                       'isi_violation', 'amplitude_cutoff','amplitude_median'], **job_kwargs)

    return metrics

def remove_similar_templates(waveforms):

    matrix = sp.compute_template_similarity(waveforms,load_if_exists=True)
    metrics = qm.compute_quality_metrics(waveforms,load_if_exists=True)
    temp_metrics = sp.compute_template_metrics(waveforms,load_if_exists=True)
    n = matrix.shape[0]

    # Find indices where values are greater than 0.5 and not on the diagonal
    indices = [(i,j) for i in range(1,n) for j in range(i-1) if matrix[i,j]>0.7]

    removables =[]
    # Print the indices
    for ind in indices:
        print(temp_metrics.index[ind[0]],temp_metrics.index[ind[1]])
        if metrics['amplitude_median'].loc[temp_metrics.index[ind[0]]] < metrics['amplitude_median'].loc[temp_metrics.index[ind[1]]]:
            smaller_index = temp_metrics.index[ind[0]]
        else:
            smaller_index = temp_metrics.index[ind[1]]

        removables.append(smaller_index)
    return removables

def remove_violated_units(metrics):

    """
    Removing based on Refractory violations, Firing_rate , snr_ratio
    amplitude_cutoff_thresh = 0.1
    isi_violations_ratio_thresh = 1
    presence_ratio_thresh = 0.9
    firing_rate = 0.1
    num_spikes = 200

    Returns an updated metrics dataframe
    
    """
    amplitude_cutoff_thresh = 0.1
    isi_violations_ratio_thresh = 1
    presence_ratio_thresh = 0.9
    firing_rate = 0.1
    num_spikes = 200
    our_query = f"(num_spikes > {num_spikes})&(amplitude_cutoff < {amplitude_cutoff_thresh}) & (isi_violations_ratio < {isi_violations_ratio_thresh}) & (presence_ratio > {presence_ratio_thresh}) & (firing_rate > {firing_rate})"

    metrics = metrics.query(our_query)

    return metrics


def postprocess_waveforms(waveforms, folder, get_quality=True, remove_violations=True, remove_similar=True):
    try: waveforms_good = si.load_waveforms(folder)
    except:
        if get_quality:
            # Start timer for quality metrics
            start_quality = timer()
            qual_metrics = get_quality_metrics(waveforms)
            # End timer for quality metrics
            end_quality = timer()
            logger.debug("Quality metrics computation takes", end_quality - start_quality)

        if remove_violations:
            # Start timer for removing violated units
            start_violations = timer()
            update_qual_metrics = remove_violated_units(qual_metrics)
            non_violated_units  = update_qual_metrics.index.values
            # End timer for removing violated units
            end_violations = timer()
            logger.debug("Removing violated units takes", end_violations - start_violations)

        if remove_similar:
            # Start timer for removing similar templates
            start_similar = timer()
            redundant_units = remove_similar_templates(waveforms)
            logger.info(f"redundant-units : {redundant_units}")
            non_violated_units = [item for item in non_violated_units if item not in redundant_units]
            for unit in redundant_units:
                try:
                    update_qual_metrics = update_qual_metrics.drop(unit)
                except Exception as e:
                    continue
            # End timer for removing similar templates
            end_similar = timer()
            logger.debug("Removing similar templates takes", end_similar - start_similar)

        # Start timer for computing template metrics
        start_template = timer()
        template_metrics = sp.compute_template_metrics(waveforms)
        template_metrics = template_metrics.loc[update_qual_metrics.index.values]
        # End timer for computing template metrics
        end_template = timer()
        logger.debug("Computing template metrics takes", end_template - start_template)

        # Create a new folder to store good waveforms
        if os.path.exists(folder):
            # If the new folder already exists, remove it
            shutil.rmtree(folder)
        # Create the directory
        os.makedirs(folder)        
        
        #Export the metrics to the sorted_unit_metrics folder
        sorted_unit_metrics_folder = folder + "/sorted_unit_metrics"
        os.makedirs(sorted_unit_metrics_folder)
        update_qual_metrics.to_excel(f"{sorted_unit_metrics_folder}/quality_metrics.xlsx")
        template_metrics.to_excel(f"{sorted_unit_metrics_folder}/template_metrics.xlsx")

        # Start the timer to measure the time taken for the following operations
        start = timer()
        # Get the extremum channel for each unit. This is the channel where the waveform reaches its minimum or maximum value.
        unit_extremum_channel = spikeinterface.full.get_template_extremum_channel(waveforms, peak_sign='both')
        # Filter the unit_extremum_channel dictionary to keep only the units that are in the non_violated_units list
        unit_extremum_channel = {key:value for key,value in unit_extremum_channel.items() if key in non_violated_units}
        # Select the good units and save them to the new folder
        waveforms_good = waveforms.select_units(non_violated_units,new_folder=folder)

        # End the timer and print the time taken for the operations
        end = timer()
        print("Removing redundant items takes", end - start)
    return waveforms_good