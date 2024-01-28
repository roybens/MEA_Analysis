#general imports
import logging

#spikeinterface imports

#local imports
import mea_processing_library as MPL

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

verbose = True
def main(h5_file_path, MaxWell_ID):
    # 1. Count the number of wells and recordings in the file
    stream_count, rec_count = MPL.count_wells_and_recs(h5_file_path, verbose = verbose)

    #2. Collect recording objects in each well.
    # - Pre-process and quality check each set of recordings. 
    # - Merge them into a single recording object. Collected merged recording objects in a list.    
    logger.info(f"Collecting, pre-processing, and merging recordings in each well:")
    merged_recordings_by_stream = []
    for stream_num in range(stream_count):
        recordings_to_merge = []
        for rec_num in range(rec_count[stream_num]):
            logger.info(f"stream_num: {stream_num}, rec_num: {rec_num}")
            recording, rec_name, stream_id = MPL.get_data_maxwell(h5_file_path, rec_num = rec_num, well_num = stream_num, verbose = True)                           
            recordings_to_merge.append(recording)
        #Pre-process and quality check recordings before merging
        recordings_to_merge = MPL.prepare_recordings_for_merge(recordings_to_merge)
        #Merge using aggergate method (allowing different channels in each recording)
        merged_recording = MPL.merge_mea_recording_segments(recordings_to_merge, mode = 'aggregate')
        MPL.count_recording_segments(merged_recording)
        merged_recordings_by_stream.append(merged_recording)

    
        

    print("boop")
    #MPL.get_data_maxwell(file_path, recnumber, well_num=well_number)