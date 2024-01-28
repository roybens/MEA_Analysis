#general imports
import logging

#spikeinterface imports
import spikeinterface
import spikeinterface.full as si
import spikeinterface.extractors as se
import spikeinterface.preprocessing as spre

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
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
stream_handler.setFormatter(formatter)

def count_wells_and_recs(h5_file_path, verbose=False):
    """
    This function counts the number of wells (stream IDs) in a given file path by attempting to read recordings 
    from the file until an error occurs. It also counts the number of successful reads (recordings) for each stream ID.

    If the verbose flag is set to True, the function prints the stream ID and the number of recordings detected 
    for each stream ID.

    If the number of recordings is not consistent across all stream IDs, the function prints a warning and 
    the range of recordings detected. Otherwise, it prints the common number of recordings per stream ID.

    Finally, the function returns the total number of stream IDs detected and the number of recording segments 
    per stream ID.

    Parameters:
    file_path (str): The path to the file to read from.
    verbose (bool): If True, print detailed output. Default is False.

    Returns:
    tuple: 
    (1) number of stream IDs detected 
    (2) and an array containing number of recording segments per stream ID.
    """
    logger.info(f"Counting wells and recordings in {h5_file_path}:")
    stream_ids = []
    rec_counts = []
    stream_id = 0
    rec_count = 0
    while True:
        try:
            stream_id_str = f'well{stream_id:03}'
            rec_name_str = f'rec{rec_count:04}'
            recording = se.read_maxwell(h5_file_path, rec_name=rec_name_str, stream_id=stream_id_str)
            rec_count += 1
        except:
            if rec_count == 0:
                break
            stream_ids.append(stream_id_str)
            rec_counts.append(rec_count)
            if verbose:
                logger.info(f"Stream ID: {stream_id_str}, {rec_count} recordings detected.")
            rec_count = 0
            stream_id += 1
    if len(set(rec_counts)) > 1: 
        logger.error(f"Warning: The number of recordings is not consistent across all stream IDs. Range: {min(rec_counts)}-{max(rec_counts)}")
    else: 
        logger.info(f"Recordings per Stream ID: {rec_count}")
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
    return recording, rec_name, stream_id

def merge_mea_recording_segments(recordings, mode='concatenate'):
    """
    recordings: List of recording objects
    mode: Either 'append' or 'concatenate'. Determines how the recordings are merged.
    """
    logger.info(f"Merging Recording Objects ({mode}):")
    try:
        if mode == 'append':
            merged_recording = spikeinterface.append_recordings(recordings)
        elif mode == 'concatenate':
            merged_recording = spikeinterface.concatenate_recordings(recordings)
        elif mode == 'aggregate':
            merged_recording = si.aggregate_channels(recordings)
        else:
            logger.error("Invalid mode. Must be either 'append' or 'concatenate'.")
    except Exception as e:
        merged_recording = None
        logger.error(f"Failed to merge recordings. Exception: {e}")

    return merged_recording

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