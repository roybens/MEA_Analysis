import numpy as np
import matplotlib.pyplot as plt
from tsmoothie.smoother import GaussianSmoother
import spikeinterface
import spikeinterface.full as si
import spikeinterface.extractors as se
import spikeinterface.sorters as ss
import spikeinterface.comparison as sc
import spikeinterface.widgets as sw
import spikeinterface.postprocessing as sp
import spikeinterface.preprocessing as spre
import spikeinterface.qualitymetrics as qm
import lib_helper_functions as helper
from pathlib import Path
from timeit import default_timer as timer
import multiprocessing
import os
from datetime import datetime
import logging
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
logger = logging.getLogger()
logger.setLevel(logging.INFO)



job_kwargs = dict(n_jobs=4, chunk_duration="1s", progress_bar=True)
#Set the job_kwargs for the extract_waveforms function
#for block-wise axon track test
#job_kwargs_we = dict(n_jobs=4, chunk_duration="1s", progress_bar=True, method = "radius", radius_um = 300)
#for full axon track test
job_kwargs_we = dict(n_jobs=4, chunk_duration="1s", progress_bar=True, method = "radius", radius_um = 2000)
#note: default method is "radius" and default radius_um is 100S


def save_to_zarr(filepath,op_folder):
    """
    saves and compresses the file into the zarr format. 
    IP: file and the corresponding output folder.
    """
    from flac_numcodecs import Flac
    compressor = Flac(level =8)
    rec_original = se.read_maxwell(filepath)
    rec_int32 = spre.scale(rec_original, dtype="int32")
    # remove the 2^15 offset
    rec_rm_offset = spre.scale(rec_int32, offset=-2 ** 15)
    # now we can safely cast to int16
    rec_int16 = spre.scale(rec_rm_offset, dtype="int16")
    recording_zarr = rec_int16.save(format = "zarr",folder=op_folder,compressor=compressor,
                                    channel_chunk_size =2,n_jobs=64,chunk_duration="1s")
    return recording_zarr

def export_to_phy_datatype(we):

    from spikeinterface.exporters import export_to_phy
    export_to_phy(we,output_folder='/home/mmpatil/Documents/spikesorting/MEA_Analysis/Python/phy_folder',**job_kwargs)

def get_channel_recording_stats(recording):
    
    channel_ids = recording.get_channel_ids()
    fs = recording.get_sampling_frequency()
    num_chan = recording.get_num_channels()
    num_seg = recording.get_num_segments()
    total_recording = recording.get_total_duration()

    #print('Channel ids:', channel_ids)
    print('Sampling frequency:', fs)
    print('Number of channels:', num_chan)
    print('Number of segments:', num_seg)
    print(f"total_recording: {total_recording} s")
    return fs,num_chan,channel_ids, total_recording

def preprocess(recording):  ## some hardcoded stuff.
    """
    Does the bandpass filtering and the common average referencing of the signal
    
    """
    recording_bp = spre.bandpass_filter(recording, freq_min=300, freq_max=(recording.sampling_frequency/2)-1)
    

    recording_cmr = spre.common_reference(recording_bp, reference='global', operator='median')

    recording_cmr.annotate(is_filtered=True)

    return recording_cmr

def get_kilosort_result(folder):

    sorter = ss.Kilosort3Sorter._get_result_from_folder(folder)
    return sorter

def get_waveforms_result(folder,with_recording= True,sorter = None):

    waveforms = si.load_waveforms(folder,with_recording=with_recording,sorting=sorter)

    return waveforms

def run_kilosort(recording,output_folder):
    print("hello I am here")
    #docker_image = 'rohanmalige/benshalom:v3'
    docker_image = 'spikeinterface/kilosort2-compiled-base:latest'
    #docker_image = "kilosort2-maxwellcomplib:latest"
    default_KS2_params = ss.get_default_sorter_params('kilosort2')
    #default_KS3_params = ss.get_default_sorter_params('kilosort3')
    #default_KS2_params['keep_good_only'] = True
    #default_KS2_params['detect_threshold'] = 12
    #default_KS2_params['projection_threshold']=[18, 10]
    #default_KS2_params['preclust_threshold'] = 14
    #run_sorter = ss.run_kilosort2(recording, output_folder=output_folder, docker_image= "kilosort2-maxwellcomplib:latest",verbose=True, **default_KS2_params)
    try: run_sorter = ss.run_kilosort2(recording, output_folder=output_folder, docker_image= docker_image,verbose=True, **default_KS2_params)
    except Exception as e:
        print(e)
        helper.empty_directory(output_folder)
        os.rmdir(output_folder) 
        try: run_sorter = ss.run_kilosort2(recording, output_folder=output_folder, docker_image= docker_image,verbose=True, **default_KS2_params)
        except Exception as e:
            print(e)
    #run_sorter = ss.run_kilosort3(recording, output_folder=output_folder, docker_image= "spikeinterface/kilosort-compiled-base:latest",verbose=True, **default_KS3_params)

    sorting_KS3 = ss.Kilosort3Sorter._get_result_from_folder(output_folder+'/sorter_output/')
    return sorting_KS3

def extract_waveforms(recording,sorting_KS3,folder):
    folder = Path(folder)

    #waveforms = si.extract_waveforms(recording,sorting_KS3,folder=folder,overwrite=True,**job_kwargs)
    waveforms = si.extract_waveforms(recording,sorting_KS3,folder=folder,overwrite=True, sparse = True, ms_before=1., ms_after=2.,allow_unfiltered=True,**job_kwargs_we)
    return waveforms


def get_template_metrics(waveforms):
    return sp.compute_template_metrics(waveforms,load_if_exists=True)

def get_temp_similarity_matrix(waveforms):
    return sp.compute_template_similarity(waveforms,load_if_exists=True)


def get_quality_metrics(waveforms): 
    
    #TO DO: to fix the SNR showing as INF by looking at the extensions of waveforms
    # This is because the noise levels are 0.0
    # Similar issue with https://github.com/SpikeInterface/spikeinterface/issues/1467 
    # Best solution proposed was to stick with applying a gain to the entire recording before all preprocessing
    sp.compute_spike_amplitudes(waveforms,**job_kwargs)
    metrics = qm.compute_quality_metrics(waveforms, metric_names=['num_spikes','firing_rate', 'presence_ratio', 'snr',
                                                       'isi_violation', 'amplitude_cutoff','amplitude_median'], **job_kwargs)

    return metrics

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

def analyse_waveforms_sigui(waveforms) :
    import spikeinterface_gui
    # This creates a Qt app
    waveforms.run_extract_waveforms(**job_kwargs)
    app = spikeinterface_gui.mkQApp() 

    # create the mainwindow and show
    win = spikeinterface_gui.MainWindow(waveforms)
    win.show()
    # run the main Qt6 loop
    app.exec_()
    return


def get_unique_templates_channels(good_units, waveform):
    """
    Analyses all the units and their corresponding extremum channels.. Removes units which correspond to same channel keeping the highest amplitude one.
    
    ToDO : have to do some kind of peak matching to remove units which have same extremum electrode.
    """
    
    #get_extremum_channels.
    unit_extremum_channel =spikeinterface.full.get_template_extremum_channel(waveform, peak_sign='neg')
    #Step 1: keep only units that are in good_units 
    unit_extremum_channel = {key:value for key,value in unit_extremum_channel.items() if key in good_units}
    print(f"extremum channel : {unit_extremum_channel}")
    #Step3: get units that correspond to same electrodes.
    output_units = [[key for key, value in unit_extremum_channel.items() if value == v] for v in set(unit_extremum_channel.values()) if list(unit_extremum_channel.values()).count(v) > 1]
    print(f"Units that correspond to same electrode: {output_units}")
    #Step 3: get the metrics
    

    #Step4: select best units with same electrodes ( based on amp values)
    output=[]
    if output_units :
        for sublist in output_units :
            amp_max = 0 
            for unit in sublist:
                if metrics['amplitude_median'][int(unit)] > amp_max :
                    amp_max = metrics['amplitude_median'][int(unit)]
                    reqd_unit = unit
            output.append(reqd_unit)
    print(f"Best unit among the same electrodes {output}")
    #Step 5 --> unit_extremum_channel - output_units + output
    output_units = [element for sublist in output_units for element in sublist]
    new_list = [ item for item in output_units if item not in output]

    required_templates = {key:value for key,value in unit_extremum_channel.items() if key not in new_list}

    return required_templates

def get_channel_locations_mapping(recording):

    channel_locations = recording.get_channel_locations()
    channel_ids = recording.get_channel_ids()
    channel_locations_mappings= {channel_id: location for location, channel_id in zip(channel_locations, channel_ids)}
    return channel_locations_mappings

def get_data_maxwell(file_path,rec_num, well_num = None):

    rec_num =  str(rec_num).zfill(4)
    rec_name = 'rec' + rec_num
    stream_id='well' + str(well_num).zfill(3) if well_num is not None else None
    if well_num is not None:
        try:
            recording = se.read_maxwell(file_path,rec_name=rec_name, stream_id=stream_id)
        except Exception as e:
            print(f"Failed to read recording with well: {e}")
    else:
        try:
            recording = se.read_maxwell(file_path,rec_name=rec_name)
        except Exception as e:
            print(f"Failed to read recording without well: {e}")
    return recording,rec_name, stream_id


def routine_sequential(file_path,number_of_configurations,time_in_s):

    for rec_number in range(number_of_configurations):
        process_block(rec_number,file_path,time_in_s)
    return

 
def routine_parallel(file_path,number_of_configurations,time_in_s):

    inputs = [(x, file_path, time_in_s) for x in range(number_of_configurations)]
    pool = multiprocessing.Pool(processes=4)
    results = pool.starmap(process_block, inputs)
    # close the pool of processes
    pool.close()
    
    # wait for all the processes to finish
    pool.join()

    return results

def preprocess_waveforms(waveforms, recording_dir, rec_name, scan_type, run_id, get_quality=True, remove_violations=True, remove_similar=True):
    new_folder=recording_dir+f'/{scan_type}{run_id}_waveforms_good_'+ rec_name
    try: waveforms_good = si.load_waveforms(new_folder)
    except:
        if get_quality:
            # Start timer for quality metrics
            start_quality = timer()
            qual_metrics = get_quality_metrics(waveforms)
            # End timer for quality metrics
            end_quality = timer()
            logging.debug("Quality metrics computation takes", end_quality - start_quality)

        if remove_violations:
            # Start timer for removing violated units
            start_violations = timer()
            update_qual_metrics = remove_violated_units(qual_metrics)
            non_violated_units  = update_qual_metrics.index.values
            # End timer for removing violated units
            end_violations = timer()
            logging.debug("Removing violated units takes", end_violations - start_violations)

        if remove_similar:
            # Start timer for removing similar templates
            start_similar = timer()
            redundant_units = remove_similar_templates(waveforms)
            logging.info(f"redundant-units : {redundant_units}")
            non_violated_units = [item for item in non_violated_units if item not in redundant_units]
            for unit in redundant_units:
                try:
                    update_qual_metrics = update_qual_metrics.drop(unit)
                except Exception as e:
                    continue
            # End timer for removing similar templates
            end_similar = timer()
            logging.debug("Removing similar templates takes", end_similar - start_similar)

        # Start timer for computing template metrics
        start_template = timer()
        template_metrics = sp.compute_template_metrics(waveforms)
        template_metrics = template_metrics.loc[update_qual_metrics.index.values]
        # End timer for computing template metrics
        end_template = timer()
        logging.debug("Computing template metrics takes", end_template - start_template)

        # Uncomment these lines if you want to save the quality metrics and template metrics to Excel files
        # update_qual_metrics.to_excel(f"/home/mmp/disktb/mmpatil/MEA_Analysis/Python/sorted_unit_metrics/quality_metrics_{run_id}.xlsx")
        # template_metrics.to_excel(f"/home/mmp/disktb/mmpatil/MEA_Analysis/Python/sorted_unit_metrics/template_metrics_{run_id}.xlsx")

        # Uncomment this line if you want to get the unique templates channels
        # template_channel_dict = get_unique_templates_channels(non_violated_units,waveforms)

        # If you want to get the non-redundant templates, uncomment the following line
        # non_redundant_templates = list(template_channel_dict.keys())

        # Start the timer to measure the time taken for the following operations
        start = timer()

        # Get the extremum channel for each unit. This is the channel where the waveform reaches its minimum or maximum value.
        unit_extremum_channel = spikeinterface.full.get_template_extremum_channel(waveforms, peak_sign='both')

        # Filter the unit_extremum_channel dictionary to keep only the units that are in the non_violated_units list
        unit_extremum_channel = {key:value for key,value in unit_extremum_channel.items() if key in non_violated_units}

        # If you want to select the good units and save them to a new folder, uncomment the following lines
        new_folder=recording_dir+f'/{scan_type}{run_id}_waveforms_good_'+ rec_name
        if os.path.exists(new_folder):
            # If the new folder already exists, clear it
            helper.empty_directory(new_folder)
            # Remove the directory
            os.rmdir(new_folder)
        # Select the good units and save them to the new folder
        waveforms_good = waveforms.select_units(non_violated_units,new_folder=recording_dir+f'/{scan_type}{run_id}_waveforms_good_'+rec_name)

        # End the timer and print the time taken for the operations
        end = timer()
        print("Removing redundant items takes", end - start)
    return waveforms_good

def process_block(recnumber, scan_type, chip_id, date, file_path, sorting_folder = "./data/temp_data/sorting/", clear_temp_files=False, well_number=None):
    
    def clear_temp_folders(sorting_folder, clear_temp_files=False):
        # Check if the sorting_folder exists and is not empty
        if helper.isexists_folder_not_empty(sorting_folder):
            # If the clear_temp_files flag is set, clear the sorting folder
            if clear_temp_files:
                print("Clearing sorting folder")
                helper.empty_directory(sorting_folder)                

    def attempt_get_sorting_records(recording_dir, kilosort_output_folder):
        # Attempt to load existing sorting recordings
        print("Attempting to load extant sorting recordings")
        try:
            # Try to load the sorting recordings from the specified path
            sorting_KS3 = ss.Kilosort2Sorter._get_result_from_folder(kilosort_output_folder+'/sorter_output/')
            print("Loaded sorting output")
            # If the sorting recordings are loaded successfully, return them
            return sorting_KS3
        except:
            # If an error occurs while loading the sorting recordings, print an error message
            print("Failed to load sorting recordings")

            # Check if kilosort_output_folder exists
            if os.path.exists(kilosort_output_folder):
                # Clear the kilosort_output_folder folder
                print("Clearing kilosort_output_folder folder")
                helper.empty_directory(kilosort_output_folder)                 
                print("Deleting recording folder")
                os.rmdir(kilosort_output_folder)
            else:
                print("Recording directory does not exist")
            return None


    

    # If the sorting folder does not exist or is empty, run the sorting routine
    print("Running sorting routine")

    # Get the recording data from the Maxwell file
    recording, rec_name, well_name = get_data_maxwell(file_path, recnumber, well_num=well_number)

    # import json
    # # Define a JSON encoder for numpy arrays
    # class NumpyEncoder(json.JSONEncoder):
    #     def default(self, obj):
    #         if isinstance(obj, np.ndarray):
    #             return obj.tolist()
    #         return json.JSONEncoder.default(self, obj)

    # # Define the path to the JSON file
    # json_file_path = f"{sorting_folder}{chip_id}_{date}_Max2_{well_name}/recordings_{rec_name}.json"

    # # Write the recording data to the JSON file
    # with open(json_file_path, 'w') as json_file:
    #     json.dump(recording, json_file, cls=NumpyEncoder)

    # Log the name of the recording being processed
    logging.info(f"Processing recording: {rec_name}")

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

    # Get a chunk of the recording based on the start and end times
    recording_chunk = recording.frame_slice(start_frame=int(time_start * fs), end_frame=int(time_end * fs))

    # Preprocess the recording chunk
    recording_chunk = preprocess(recording_chunk)

    # Get the current working directory
    current_directory = os.getcwd()

    try:
        # Extract the run ID from the file path using a regular expression
        pattern = r"/(\d+)/data.raw.h5"
        run_id = int(re.search(pattern, file_path).group(1))

        # If the well name is not None, append it to the directory and sort ID names
        if well_name is not None:
            recording_dir = f"{sorting_folder}{chip_id}_{date}_Max2_{well_name}"
            sort_id = f"{chip_id}_{date}_Max2_{well_name}"
        else:
            # If the well name is None, do not append it to the directory and sort ID names
            recording_dir = f"{sorting_folder}{chip_id}_{date}_Max1"
            sort_id = f"{chip_id}_{date}_Max1"

        # Define the path for the Kilosort output folder
        kilosort_output_folder = recording_dir + f'/{scan_type}{run_id}_kilosort2_{rec_name}'

        # Clear the temporary sorting folder before starting the new sorting process
        clear_temp_folders(sorting_folder, clear_temp_files=clear_temp_files)
        sortingKS3 = attempt_get_sorting_records(recording_dir, kilosort_output_folder)
        
        # Create a new directory for the current recording. If the directory already exists, this will throw an error.
        try: os.mkdir(recording_dir, 0o777)
        except: 
            print("Directory already exists")
            pass
        
        # Change the current working directory to the new directory. All subsequent file operations will take place in this directory.
        #os.chdir(recording_dir)

        if sortingKS3 is None:     
            # Start timer for sorting
            start_sorting = timer()
            sortingKS3 = run_kilosort(recording_chunk,output_folder=kilosort_output_folder)
            sortingKS3 = sortingKS3.remove_empty_units()
            sortingKS3 = spikeinterface.curation.remove_excess_spikes(sortingKS3,recording_chunk)
            # End timer for sorting
            end_sorting = timer()
            logging.debug("Sorting takes", end_sorting - start_sorting)
        else:
            logging.info("Pre-loaded sortingKS3 found. Time saved on sorting.")

        # Start timer for waveform extraction
        start_waveform_extraction = timer()
        waveform_folder = recording_dir+f'/{scan_type}{run_id}_waveforms_'+ rec_name
        logging.info("Extracting waveforms")
        #try to delete waveforms folder if it exists
        try: waveforms = si.load_waveforms(waveform_folder)
        except: 
            if os.path.exists(waveform_folder):
                helper.empty_directory(waveform_folder)
                os.rmdir(waveform_folder)
            waveforms = extract_waveforms(recording_chunk,sortingKS3,folder = waveform_folder)
        # End timer for waveform extraction
        end_waveform_extraction = timer()
        logging.debug("Waveform extraction takes", end_waveform_extraction - start_waveform_extraction)
        waveforms_good = preprocess_waveforms(waveforms, recording_dir, rec_name, scan_type, run_id, get_quality=True, remove_violations=True, remove_similar=False)

        fig, ax = plt.subplots()
        si.plot_probe_map(recording_chunk, with_channel_ids=False, ax=ax)
        ax.set_xlim(0, 3850)
        ax.set_ylim(0, 2100) #microns

        #save figure
        waveform_good_folder=recording_dir+f'/{scan_type}{run_id}_waveforms_good_'+ rec_name
        fig.savefig(waveform_folder+f'/probe_map.png', dpi=300)
        fig.savefig(waveform_good_folder+f'/probe_map.png', dpi=300)


        
        #si.widgets.plot_unit_waveforms_density_map(waveforms_good)

       

        #channel_location_dict = get_channel_locations_mapping(recording_chunk)
        channel_locations = recording.get_channel_locations()
        # New dictionary with combined information
        # new_dict = {}

        # # Iterate over each template in the template_channel_dict dictionary
        # for template, channel in unit_extremum_channel.items():

        #     # If this channel is not already in the new dictionary, add it
        #     if template not in new_dict:
        #         new_dict[template] = {}

        #     # Add an entry for this template and its corresponding location to the new dictionary
        #     new_dict[template][channel] = [int(channel_location_dict[channel][0]/17.5),int(channel_location_dict[channel][1]/17.5)]
        
        # electrodes = []
        # # Iterate over each template in the template_channel_dict dictionary
        # for template, channel in waveforms.items():
        #     # Add an entry for this template and its corresponding location to the new dictionary
        #     electrodes.append(220* int(channel_location_dict[channel][1]/17.5)+int(channel_location_dict[channel][0]/17.5))
        # if clear_temp_files:
        #     helper.empty_directory(sorting_folder)
        # os.chdir(current_directory)
        return recording_dir
        # # file_name = '/mnt/disk15tb/mmpatil/MEA_Analysis/Python/Electrodes/Electrodes_'+rec_name
        # # helper.dumpdicttofile(new_dict,file_name)
        # 

       
        # return electrodes, len(update_qual_metrics)
    except Exception as e:
        logging.info(e)
        logging.info(f"Error in {rec_name} processing. Continuing to the next block")
        if clear_temp_files:
            helper.empty_directory(sorting_folder+'/block_'+rec_name)
        e = "The sorting failed."
        return e,e
    #helper.dumpdicttofile(failed_sorting_rec,'./failed_recording_sorting')

if __name__ =="__main__" :

    #routine_sequential()  #arguments
    #development code
    continuous_dir = ["/mnt/ben-shalom_nas/irc/media/harddrive8tb/Adam/MEA_AxonTraceDevScans/230828/16862/ActivityScan/000014/data.raw.h5",
                    "/mnt/ben-shalom_nas/irc/media/harddrive8tb/Adam/MEA_AxonTraceDevScans/230828/16862/Network/000015/data.raw.h5", 
                    "/mnt/disk20tb/PrimaryNeuronData/Maxtwo/FolicAcid/FolicAcid/230921/M05506/Network/000008/data.raw.h5"]    
    process_block
