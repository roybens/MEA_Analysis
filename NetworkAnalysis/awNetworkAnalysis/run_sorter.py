import os
import glob
import spikeinterface.sorters as ss
from MEA_Analysis.MEAProcessingLibrary import mea_processing_library as mea
# pip install -e **/MEA_Analysis

import sys

def run_sorter(
    raw_data_path,
    sorted_output_dir,
    waveform_output_dir,
    use_docker=True,
    try_load=True,
    ):
    
    #init
    print(f'raw_data_path: {raw_data_path}')
    #print(f'sorted_output_dir: {sorted_output_dir}')
    print(f'sorted_output_dir: {sorted_output_dir}')
    print(f'waveform_output_dir: {waveform_output_dir}')
    
    # load recording 
    _, recordings, _, _ = mea.load_recordings(raw_data_path) # expects a string path to a .h5 file, gets dict of recordings - 1 for each well
    recording_details = mea.extract_recording_details(raw_data_path)[0]     # NOTE: this works for a list of dirs or a single dir - but treats single dir 
                                                                            # as a list of a single dir
    h5_file_path = recording_details['h5_file_path']
    runID = recording_details['runID']
    scanType = recording_details['scanType']
    chipID = recording_details['chipID']
    date = recording_details['date']
    projectName = recording_details['projectName'] 

    # spike sort recordings
    kilosort2_params = ss.Kilosort2Sorter.default_params()
    
    #modify params as needed
    kilosort2_params['keep_good_only'] = True # false by default
    kilosort2_params['n_jobs'] = 32
    # change chunk size from '1s' to half a second to help with memory issues
    kilosort2_params['chunk_duration'] = '500ms'
    #kilosort2_params['NT'] = 100 # number of timepoints in waveforms
    #kilosort2_params['wave_length'] = 81 # increase from default of 61 to 100 to help capture more information for long refractory periods typical of excitatory neurons
    #kilosort2_params['n_jobs'] = 64
    
    
    for wellid, recording in recordings.items():
        # print
        print('====================================================================================================')
        print(f'Processing well: {wellid}')
        
        
        for segment in recording:  #   NOTE: generally, for this project we're working with network recordings, so there's only one segment per recording
                                            #       still, in the even we use activity scans or axon tracking scans - there will be multiple segments that need to be 
                                            #       handled after spike sorting
            # define output folder
            output_folder = os.path.join(sorted_output_dir, projectName, date, chipID, scanType, runID, wellid) #same filepath structure as the raw data + wellid
            
            # try to load the kilosort2 results if they already exist
            load_success = False
            if try_load:
                print(f'try_load is true. Attempting to load kilosort2 results for {wellid}...')
                
                print(f'output_folder: {output_folder}') 
                
                print(f'output_folder exists?: {os.path.exists(output_folder)}')
                
                # import sys
                # sys.exit() #for debugging          
                
                # check if the output folder already exists
                if os.path.exists(output_folder):
                    try:
                        
                        sorted_data = mea.load_kilosort2_results(output_folder) # FIXME: this function may not exist right now. Need to check.
                        # print(f'Kilosort2 results already exist for {wellid}. Skipping...')
                        # continue
                        print(f'Kilosort2 results loaded for {wellid}.')
                        load_success = True
                    except:
                        print(f'Kilosort2 results failed to load for {wellid}. Running...')
                else:
                    print(f'Kilosort2 results do not exist for {wellid}. Running...')
                    
                # import sys
                # sys.exit() #for debugging  
            
            # run kilosort2
            if not load_success: #  # if the kilosort2 results don't already exist, run kilosort2
                try:
                    if use_docker:
                        sorted_data = mea.run_kilosort2_docker_image(
                            segment, # each segment should be a maxwellrecordingextractor object
                            output_folder=output_folder,
                            docker_image="spikeinterface/kilosort2-compiled-base:latest", #FIXME: this is the wrong image... pretty sure
                            verbose=True, 
                            logger=None, 
                            params=kilosort2_params,
                            
                        )
                    else:
                        print('running kilosort2 without docker...')
                        #sys.exit()
                        sorted_data = mea.kilosort2_wrapper(
                            segment, # each segment should be a maxwellrecordingextractor object
                            output_folder=output_folder,
                            sorting_params=kilosort2_params,
                            verbose=True,
                        )
                except Exception as e:
                    print(f'Kilosort2 failed for {wellid} with error: {e}')
                    continue
                
            # assert  number of units is greater than 0
            try:
                unit_ids = sorted_data.get_unit_ids()
                assert len(unit_ids) > 0, 'No units found in sorted data.'
            except Exception as e:
                print(f'failed to get unit ids for {wellid} with error: {e}')
                print(f'skipping {wellid}...')
                continue                       
            
            # extract waveforms
            try:
                # TODO: decide if this is needed. I think it's probably not since these are continuous 1-segment recordings.
                # import spikeinterface.full as si
                # seg_index = 0 #this should always be zero since we're working with network recordings right?
                # chunk_size = min(10000, segment.get_num_samples()) - 100
                # rec_centered = si.center(segment, chunk_size=chunk_size)
                # seg_sort = si.SelectSegmentSorting(sorted_data, seg_index)
                # seg_sort = si.remove_excess_spikes(seg_sort, rec_centered)
                # seg_sort.register_recording(rec_centered)
                waveform_output = os.path.join(waveform_output_dir, projectName, date, chipID, scanType, runID, wellid)
                mea.extract_waveforms(segment, sorted_data, waveform_output, n_jobs=16, ms_before=1, ms_after=5) # aw 2025-02-14 13:26:25 increates ms_before and after to try and better capture excitatory neurons
            except Exception as e:
                print(f'Waveform extraction failed for {wellid} with error: {e}')
                continue
            
            # export to phy
            # import spikeinterface as si # core module only
            # from spikeinterface.postprocessing import compute_spike_amplitudes, compute_principal_components
            # from spikeinterface.exporters import export_to_phy
            # sorting_analyzer = si.create_sorting_analyzer(sorting=sorted_data, recording=segment)

            # # some computations are done before to control all options
            # sorting_analyzer.compute(['random_spikes', 'waveforms', 'templates', 'noise_levels'])
            # _ = sorting_analyzer.compute('spike_amplitudes')
            # _ = sorting_analyzer.compute('principal_components', n_components = 5, mode="by_channel_local")

            # # the export process is fast because everything is pre-computed
            # phy_folder = os.path.join(output_folder, 'phy')
            # export_to_phy(sorting_analyzer=sorting_analyzer, output_folder=phy_folder)
                
        # print big separator so its easy to see where each well starts and ends
        print('====================================================================================================')
    
    # attempt to find the log file and rename it for the current well
    def capture_logs(sorted_output_dir, projectName, date, chipID, scanType, runID, wellid):
        log_dir = os.path.join(sorted_output_dir, 'logs')
        log_files = glob.glob(os.path.join(log_dir, 'test_spikesort.log'))
        if len(log_files) > 0:
            for log_file in log_files:
                from datetime import datetime
                now = datetime.now()
                #make folder named for yyyymmdd
                log_date = now.strftime('%Y%m%d')
                log_dir = os.path.join(log_dir, log_date)
                if not os.path.exists(log_dir):
                    os.makedirs(log_dir)
                #rename log file 
                #new_log_file = os.path.join(log_dir, f'{projectName}_{date}_{chipID}_{scanType}_{runID}_{wellid}_spikesort.log')
                new_log_file = os.path.join(log_dir, f'{projectName}_{date}_{chipID}_{scanType}_{runID}_spikesort.log')
                os.rename(log_file, new_log_file)
    try: capture_logs(sorted_output_dir, projectName, date, chipID, scanType, runID, wellid)
    except: pass
    
    # print done
    print('done')
    
    # # return sorted_data
    # return sorted_data                          