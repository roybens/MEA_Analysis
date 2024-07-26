#general imports
import logging
import os

#spikeinterface imports
import spikeinterface
import spikeinterface.extractors as se
import spikeinterface.sorters as ss
import spikeinterface.full as si

#local imports
from MEAProcessingLibrary.mea_processing_library import 
import lib_helper_functions as helper
import spike_sorting as hp #hornauerp, https://github.com/hornauerp/axon_tracking/blob/main/axon_tracking/spike_sorting.py

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

def process_and_merge_recordings(h5_file_path, stream_count, rec_count, recordings_dir, date, chip_id, scanType, runID, verbose):
    recs_to_sort_by_stream = []
    merged_recordings_by_stream = []
    for stream_num in range(stream_count):
        #debug
        if stream_num > 0: break
        #debug
        #Create a list of recording objects to merge
        recordings_to_merge = []
        recordings_to_sort = []
        for rec_num in range(rec_count[stream_num]):
            logger.info(f"stream_num: {stream_num}, rec_num: {rec_num}")
            #allow Maxone and Maxtwo files to be processed:
            try: recording, rec_name, stream_id = MPL.get_data_maxwell(h5_file_path, rec_num = rec_num, well_num = stream_num, verbose = verbose)
            except: recording, rec_name, stream_id = MPL.get_data_maxwell(h5_file_path, rec_num = rec_num, verbose = verbose)                                          
            recordings_to_merge.append(recording)
            #pre-process single recordings
            recording_chunk = MPL.preprocess_single_recording(recording)
            recordings_to_sort.append(recording_chunk)
        #Add recordings to list of recordings to sort
        recs_to_sort_by_stream.append(recordings_to_sort)
        stream_dir = f"Well#{str(stream_num+1)}"
        relative_path = recordings_dir + f"/{date}/{chip_id}/{scanType}/{runID}/{stream_dir}/aggregate_recording"
        #Check if recording object already exists
        if os.path.exists(relative_path):
            logger.info(f"Recording object already exists for {relative_path}")
            try:
                logger.info(f"Attempting to load recording object for {relative_path}")
                merged_recording = spikeinterface.load_extractor(relative_path)
                merged_recordings_by_stream.append(merged_recording)
                logger.info(f"Recording object loaded for {relative_path}")
                continue
            except:
                logger.info(f"Recording object could not be loaded for {relative_path}, will be re-created")            
        #Pre-process and quality check recordings before merging
        recordings_to_merge = MPL.prepare_recordings_for_merge(recordings_to_merge)
        #Merge using aggergate method (allowing different channels in each recording)
        merged_recording = MPL.merge_mea_recording_segments(recordings_to_merge, mode = 'aggregate')
        #Save recording objects
        logger.info(f"Saving recording objects")
        merged_recordings_by_stream.append(merged_recording)
        #merged_recording.save(folder = relative_path, overwrite = True, verbose = verbose)
    return recs_to_sort_by_stream, merged_recordings_by_stream

def sort_recordings(merged_recordings_by_stream, stream_count, rec_count, spikesorting_dir, date, chip_id, scanType, runID, recs_to_sort_by_stream, verbose):
    sortings_by_stream = []
    for stream_num in range(stream_count):
        #debug
        if stream_num > 0: break
        #debug
        merged_recording = merged_recordings_by_stream[stream_num]
        stream_dir = f"Well#{str(stream_num+1)}"
        sortings_to_merge = []
        merged_path = spikesorting_dir + f"/{date}/{chip_id}/{scanType}/{runID}/{stream_dir}/aggregate_sorting"
        spike_sort_merged = True
        if spike_sort_merged:
            #for rec_num in range(rec_count[stream_num]):            
                #rec_dir = f"Rec#{str(rec_num+1)}"
                #relative_path = spikesorting_dir + f"/{date}/{chip_id}/{scanType}/{runID}/{stream_dir}/sorting_{rec_dir}"
                #merged_recording = recs_to_sort_by_stream[stream_num][rec_num]
                #Check if recording object already exists
            if os.path.exists(merged_path):
                logger.info(f"Sorting object already exists for {merged_path}")
                try:
                    logger.info(f"Attempting to load sorting object for {merged_path}")
                    sorting = ss.Kilosort2Sorter._get_result_from_folder(merged_path+'/sorter_output/')
                    sortings_to_merge.append(sorting)
                    logger.info(f"Sorting object loaded for {merged_path}")
                    continue
                except:
                    logger.info(f"Sorting object could not be loaded for {merged_path}, will be re-created")
                    # Check if kilosort_output_folder exists
                    if os.path.exists(merged_path):
                        # Clear the kilosort_output_folder folder
                        print("Clearing kilosort_output_folder folder")
                        helper.empty_directory(merged_path)                 
                        print("Deleting recording folder")
                        os.rmdir(merged_path)        
            logger.info(f"Running kilosort2.5 on {merged_path}")
            merged_sorting = MPL.run_kilosort2_5_docker_image_GPUs(merged_recording, 
                                                                #chunk_duration = 60, 
                                                                output_folder=merged_path, 
                                                                verbose=verbose,
                                                                num_gpus=2)
            #sortings_to_merge.append(sorting)
            #Merge using aggergate method (allowing different channels in each recording)
            logger.info(f"Merging recording sorted")
            #merged_sorting = MPL.merge_sortings(sortings_to_merge, mode = 'aggregate')   
        else: 
            for rec_num in range(rec_count[stream_num]):            
                rec_dir = f"Rec#{str(rec_num+1)}"
                relative_path = spikesorting_dir + f"/{date}/{chip_id}/{scanType}/{runID}/{stream_dir}/sorting_{rec_dir}"
                recording = recs_to_sort_by_stream[stream_num][rec_num]
                #Check if recording object already exists
                if os.path.exists(relative_path):
                    logger.info(f"Sorting object already exists for {relative_path}")
                    try:
                        logger.info(f"Attempting to load sorting object for {relative_path}")
                        sorting = ss.Kilosort2Sorter._get_result_from_folder(relative_path+'/sorter_output/')
                        sortings_to_merge.append(sorting)
                        logger.info(f"Sorting object loaded for {relative_path}")
                        continue
                    except:
                        logger.info(f"Sorting object could not be loaded for {relative_path}, will be re-created")
                        # Check if kilosort_output_folder exists
                        if os.path.exists(relative_path):
                            # Clear the kilosort_output_folder folder
                            print("Clearing kilosort_output_folder folder")
                            helper.empty_directory(relative_path)                 
                            print("Deleting recording folder")
                            os.rmdir(relative_path)        
                logger.info(f"Running kilosort2 on {relative_path}")
                sorting = MPL.run_kilosort2_docker_image(recording, chunk_duration = 60, output_folder=relative_path, verbose=verbose)
                sortings_to_merge.append(sorting)        
            #Merge using aggergate method (allowing different channels in each recording)
            logger.info(f"Merging sorting objects")
            merged_sorting = MPL.merge_sortings(sortings_to_merge, mode = 'aggregate')        
        logger.info(f"Post-sorting processing")        
        merged_sorting = merged_sorting.remove_empty_units()
        merged_sorting = spikeinterface.curation.remove_excess_spikes(merged_sorting, merged_recording)
        #Save recording objects
        logger.info(f"Saving sorting objects")        
        sortings_by_stream.append(merged_sorting)
        merged_sorting.save(folder = merged_path, overwrite = True, verbose = verbose)
        #logger.info(f"Sorting object does not exist for {relative_path}")
    return sortings_by_stream


verbose = True
def main(h5_file_paths, h5_file_path, recordings_dir = './AxonReconPipeline/data/temp_data/recordings', spikesorting_dir = './AxonReconPipeline/data/temp_data/sortings'):
    
    #Extract recording details from the .h5 file path
    recording_details = MPL.extract_recording_details(h5_file_path)
    date = recording_details[0]['date']
    chip_id = recording_details[0]['chipID']
    scanType = recording_details[0]['scanType']
    runID = recording_details[0]['runID']     

    # 1. Count the number of wells and recordings in the file
    stream_count, rec_count = MPL.count_wells_and_recs(h5_file_path, verbose = verbose)
    
    #Just Testing Somethings
    #This funciton mimics https://github.com/hornauerp/axon_tracking
    #si.Kilosort2_5Sorter.set_kilosort2_5_path('/home/phornauer/Git/Kilosort_2020b') #Change
    hornauerp = False
    if hornauerp:
        sorter_params = si.get_default_sorter_params(si.Kilosort2_5Sorter)
        sorter_params['n_jobs'] = -1
        sorter_params['detect_threshold'] = 7
        sorter_params['minFR'] = 0.01
        sorter_params['minfr_goodchannels'] = 0.01
        sorter_params['keep_good_only'] = False
        sorter_params['do_correction'] = False
        spikesorting_dir = './AxonReconPipeline/data/temp_data/sortings'
        spikesorting_root = spikesorting_dir+f'/{date}/{chip_id}/{scanType}/{runID}'
        sorting_list = hp.sort_recording_list(h5_file_paths,
                                            save_root=spikesorting_root,
                                            #save_path_changes= 5,
                                            sorter = 'kilosort2_5',
                                            sorter_params = sorter_params,
                                                verbose = verbose)
    else:
        pass
        
    #2. Collect recording objects in each well. 
    # - Pre-process and quality check each set of recordings. 
    # - Merge them into a single recording object. 
    # - Save pre-processed recording objects.   
    logger.info(f"Collecting, pre-processing, and merging recordings in each well:")
    # Check if recordings folder exists in ./data/temp_data, if not create it.    
    if not os.path.exists(recordings_dir):
        os.makedirs(recordings_dir)
    #Process each well
    recs_to_sort_by_stream, merged_recordings_by_stream  = process_and_merge_recordings(
        h5_file_path, 
        stream_count, 
        rec_count, 
        recordings_dir, date, chip_id, scanType, runID, verbose)
   
    #3. Peform spike sorting on each merged recording object.
    # - Save spike sorting objects.
    logger.info(f"Generating spikesorting objects for each well:")
    # Check if recordings folder exists in ./data/temp_data, if not create it.
    spikesorting_dir = './AxonReconPipeline/data/temp_data/sortings'
    if not os.path.exists(spikesorting_dir):
        os.makedirs(spikesorting_dir)
    #Process each well
    sortings_by_stream = sort_recordings(
        merged_recordings_by_stream, 
        stream_count, 
        rec_count, 
        spikesorting_dir, 
        date, chip_id, scanType, runID, recs_to_sort_by_stream, verbose)

    return merged_recordings_by_stream, sortings_by_stream