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
from spikeinterface.sorters import run_sorter_local
import helper_functions as helper
from pathlib import Path
from timeit import default_timer as timer
import multiprocessing
import os
import sys
from datetime import datetime
import logging
import re
import pickle
import scipy.io as sio
import argparse
import pandas as pd
from spikeinterface.curation import CurationSorting
import pdb
from collections import defaultdict
#os.environ['HDF5_PLUGIN_PATH']='/home/mmp/Documents/Maxlab/so/'
# Configure the logger
# Manually clear the log file
with open('./application.log', 'w'):
    pass
logging.basicConfig(
    filename='./application.log',  # Log file name
    level=logging.DEBUG,  # Log level
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # Log format
    datefmt='%Y-%m-%d %H:%M:%S' # Timestamp format
)

logger = logging.getLogger("mea_pipeline")





job_kwargs = dict(n_jobs=32, chunk_duration="1s", progress_bar=False)


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
    recording_bp = spre.bandpass_filter(recording, freq_min=300, freq_max=6000)
    

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
   
    default_KS2_params = ss.get_default_sorter_params('kilosort2')
    default_KS2_params['keep_good_only'] = True
    # default_KS2_params['detect_threshold'] = 12
    # default_KS2_params['projection_threshold']=[18, 10]
    # default_KS2_params['preclust_threshold'] = 14
    run_sorter=run_sorter_local("kilosort2",recording, output_folder=output_folder, delete_output_folder=False,verbose=True,with_output=True,**default_KS2_params)
    #run_sorter = ss.run_kilosort2(recording, output_folder=output_folder, docker_image= "si-98-ks2-maxwell",verbose=True, **default_KS2_params)
    sorting_KS3 = ss.Kilosort3Sorter._get_result_from_folder(output_folder+'/sorter_output/')
    return sorting_KS3

def extract_waveforms(recording,sorting_KS3,folder):
    folder = Path(folder)

    waveforms = si.extract_waveforms(recording,sorting_KS3,folder=folder,overwrite=True,**job_kwargs)
    #waveforms = si.extract_waveforms(recording,sorting_KS3,folder=folder,overwrite=True, sparse = True, ms_before=1., ms_after=2.,allow_unfiltered=True,**job_kwargs)
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
def find_connected_components(pairs):
    # Create a graph using a dictionary
    graph = defaultdict(list)
    
    # Populate the graph with edges
    for x, y in pairs:
        graph[x].append(y)
        graph[y].append(x)  # Assuming undirected graph for transitive relationships

    # Function to perform DFS and find connected components
    def dfs(node, visited, component):
        visited.add(node)
        component.append(node)
        for neighbour in graph[node]:
            if neighbour not in visited:
                dfs(neighbour, visited, component)

    # List to store all components
    components = []
    visited = set()

    # Iterate through all nodes in the graph
    for node in graph:
        if node not in visited:
            component = []
            dfs(node, visited, component)
            components.append(component)

    return components

def merge_similar_templates(sorting,waveforms,sim_score =0.7):

    matrix = sp.compute_template_similarity(waveforms,load_if_exists=True)
    metrics = qm.compute_quality_metrics(waveforms,load_if_exists=True)
    temp_metrics = sp.compute_template_metrics(waveforms,load_if_exists=True)
    n = matrix.shape[0]

    # Find indices where values are greater than 0.5 and not on the diagonal
    indices = [(i,j) for i in range(1,n) for j in range(i-1) if matrix[i,j]>sim_score]

    pairs = []
    # Print the indices
    for ind in indices:
        logging.debug(f"{temp_metrics.index[ind[0]],temp_metrics.index[ind[1]]}")
        pairs.append((temp_metrics.index[ind[0]],temp_metrics.index[ind[1]]))
    connected_components = find_connected_components(pairs)
    cs = CurationSorting(sorting)
    for elements in connected_components:
        cs.merge(elements)
    return cs.sorting

  


def remove_similar_templates(waveforms,sim_score =0.7):

    matrix = sp.compute_template_similarity(waveforms,load_if_exists=True)
    metrics = qm.compute_quality_metrics(waveforms,load_if_exists=True)
    temp_metrics = sp.compute_template_metrics(waveforms,load_if_exists=True)
    n = matrix.shape[0]

    # Find indices where values are greater than 0.5 and not on the diagonal
    indices = [(i,j) for i in range(1,n) for j in range(i-1) if matrix[i,j]>sim_score]

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

def analyse_waveforms_sigui(waveforms_folder) :
    import spikeinterface_gui
    waveforms = si.load_waveforms(waveforms_folder)
    # This creates a Qt app
   # waveforms.run_extract_waveforms(**job_kwargs)
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

def get_data_maxwell(file_path,rec_num):

    rec_num =  str(rec_num).zfill(4)
    rec_name = 'rec' + rec_num
    recording = se.read_maxwell(file_path,rec_name=rec_name)
    return recording,rec_name




def process_block(file_path,time_in_s= 300,recnumber=0, sorting_folder = "./Sorting_Intermediate_files",clear_temp_files=True):
    
    
    #check if sorting_folder exists and empty
    if helper.isexists_folder_not_empty(sorting_folder):
        logging.info("clearing sorting folder")
        helper.empty_directory(sorting_folder)
    
    recording,rec_name = get_data_maxwell(file_path,recnumber)
    logging.info(f"Processing recording: {rec_name}")
    fs, num_chan, channel_ids, total_rec_time= get_channel_recording_stats(recording)
    if total_rec_time < time_in_s:
        time_in_s = total_rec_time
    time_start = 0
    time_end = time_start+time_in_s
    recording_chunk = recording.frame_slice(start_frame= int(time_start*fs),end_frame=int(time_end*fs))
    recording_chunk = preprocess(recording_chunk)
    
    current_directory =  os.path.dirname(os.path.abspath(__file__))
    try:
        pattern = r"/(\d+)/data.raw.h5"
        run_id = int(re.search(pattern, file_path).group(1))
        file_pattern = os.path.dirname(file_path)

        # Split the pattern into parts using '/' as the separator
        parts = file_pattern.split('/')

        # Extract the desired pattern
        desired_pattern = '/'.join(parts[-6:])  # Assuming you want the last 6 parts

        dir_name = sorting_folder
        #os.mkdir(dir_name,0o777,)
        os.chdir(dir_name)
        kilosort_output_folder = f"{current_directory}/../AnalyzedData/{desired_pattern}/kilosort2_{rec_name}"
        start = timer()
        sortingKS3 = run_kilosort(recording_chunk,output_folder='./Kilosort_tmp')
        logging.debug("Sorting complete")
        sortingKS3 = sortingKS3.remove_empty_units()
        sortingKS3 = spikeinterface.curation.remove_excess_spikes(sortingKS3,recording_chunk) #Sometimes KS returns spikes outside the number of samples. < https://github.com/SpikeInterface/spikeinterface/pull/1378>
        
        sortingKS3= sortingKS3.save(folder = kilosort_output_folder)
        waveform_folder =dir_name+'/waveforms_'+ rec_name
        logging.info("Extracting waveforms")
        waveforms = extract_waveforms(recording_chunk,sortingKS3,folder = waveform_folder)
        end = timer()
        logging.debug(f"Sort and extract waveforms takes: { end - start}")
        
        start = timer()
        
        qual_metrics = get_quality_metrics(waveforms)  
        
        update_qual_metrics = remove_violated_units(qual_metrics)
        non_violated_units  = update_qual_metrics.index.values
        #check for template similarity.
        curation_sortin_obj = merge_similar_templates(sortingKS3,waveforms)
        curation_sortin_obj.save(folder=f"{current_directory}/../AnalyzedData/{desired_pattern}/sorting")
        curatedmerged = curation_sortin_obj.get_unit_ids()
        end = timer()
         #todo: need to extract metrics here.
        logging.debug(f"Removing redundant items takes{end - start}")                                            #todo: need to extract metrics here.
        non_violated_units_new = [item for item in curatedmerged if item not in non_violated_units]
        
        
        waveforms= extract_waveforms(recording_chunk,curation_sortin_obj,folder = waveform_folder)
        waveform_good = waveforms.select_units(non_violated_units_new,new_folder=f"{current_directory}/../AnalyzedData/{desired_pattern}/waveforms_good")
        template_metrics = sp.compute_template_metrics(waveform_good)
        #template_metrics = template_metrics.loc[update_qual_metrics.index.values]
        qual_metrics = qm.compute_quality_metrics(waveform_good ,metric_names=['num_spikes','firing_rate', 'presence_ratio', 'snr',
                                                       'isi_violation', 'amplitude_cutoff','amplitude_median'])  ## to do : have to deal with NAN values
        
        template_metrics.to_excel(f"{current_directory}/../AnalyzedData/{desired_pattern}/template_metrics.xlsx")
        locations = sp.compute_unit_locations(waveform_good)
        qual_metrics['location_X'] = locations[:,0]
        qual_metrics['location_Y'] = locations[:,1]
        qual_metrics.to_excel(f"{current_directory}/../AnalyzedData/{desired_pattern}/quality_metrics.xlsx")
        ax = plt.subplot(111)

        sw.plot_probe_map(recording_chunk,ax=ax,with_channel_ids=False)
        for x,y in locations:
            ax.scatter(x,y, s=1)
        plt.savefig(f"{current_directory}/../AnalyzedData/{desired_pattern}/locations.pdf")
        plt.clf()
        #template_channel_dict = get_unique_templates_channels(non_violated_units,waveforms)
        #non_redundant_templates = list(template_channel_dict.keys())
        # extremum_channel_dict = 
        # ToDo Until the peak matching routine is implemented. We use this.
        unit_extremum_channel =spikeinterface.full.get_template_extremum_channel(waveforms, peak_sign='both')
        #Step 1: keep only units that are in good_units 
        unit_extremum_channel = {key:value for key,value in unit_extremum_channel.items() if key in non_violated_units}
        #waveform_good = waveforms.select_units(non_violated_units,new_folder=dir_name+'/waveforms_good_'+rec_name)
        
        #get the spike trains
        os.makedirs(f"{current_directory}/../AnalyzedData/{desired_pattern}/Spike_trains/",mode=0o777, exist_ok=True)
        for idx, unit_id in enumerate(non_violated_units):
            #print(unit_id)
            spike_train = sortingKS3.get_unit_spike_train(unit_id,start_frame=0*fs,end_frame=time_end*fs)
            #print(spike_train)
            if len(spike_train) > 0:
                spike_times = spike_train / float(fs)
                mat_filename = f"{current_directory}/../AnalyzedData/{desired_pattern}/Spike_trains/{unit_id}.mat"
                sio.savemat(mat_filename,{'spike_times':spike_times,'units':[unit_id]*len(spike_times)})
        
        
        channel_location_dict = get_channel_locations_mapping(recording_chunk)
        # file_name = '/mnt/disk15tb/mmpatil/MEA_Analysis/Python/Electrodes/Electrodes_'+rec_name
        # helper.dumpdicttofile(new_dict,file_name)
        electrodes = []

        # Iterate over each template in the template_channel_dict dictionary
        for template, channel in unit_extremum_channel.items():
            # Add an entry for this template and its corresponding location to the new dictionary
            electrodes.append(220* int(channel_location_dict[channel][1]/17.5)+int(channel_location_dict[channel][0]/17.5))
        electrode_data = {'electrodes':electrodes}
        #helper.dumpdicttofile(electrode_data,f"{current_directory}/../AnalyzedData/{desired_pattern}/associated_electrodes.json")
        sio.savemat(f"{current_directory}/../AnalyzedData/{desired_pattern}/associated_electrodes.mat",electrode_data)
        os.chdir(current_directory)
        if clear_temp_files:
            helper.empty_directory(sorting_folder)
        return electrodes, len(update_qual_metrics)
    except Exception as e:
        logging.info(e)
        logging.info(f"Error in {rec_name} processing. Continuing to the next block")
        if clear_temp_files:
            helper.empty_directory(sorting_folder)
        e = "The sorting failed."
        return e,e
    #helper.dumpdicttofile(failed_sorting_rec,'./failed_recording_sorting')


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
    


    



def main():

    """
    This direct main function run would be silent run..on server.

    """
    
    parser = argparse.ArgumentParser(description="Process inpout filepaths")
    # Check if the correct number of command-line arguments is provided
   
     # Positional argument (mandatory)
    parser.add_argument('file_dir_path', type=str, help='Path to the file or directory', nargs='?')

    parser.add_argument('-r', '--reference', type=str, help='Path to the reference file (optional)')

    parser.add_argument('-t','--type',nargs='+',type=str,default=['network today', 'network today/best'], help='Array types (Optional)')
   


    args = parser.parse_args()

    # Check if the mandatory file path/dir path is provided
    if args.file_dir_path is None:
        logger.info("Usage: python script.py <file_or_folder_path> ")
        sys.exit(1)

    # Get the user-provided path from the command-line argument
    path = args.file_dir_path
    
    # Check if the path exists
    if not os.path.exists(path):
        logger.info(f"The specified path '{path}' does not exist.")
        sys.exit(1)

    # Check if the path is a file
    if os.path.isfile(path):
        logger.debug(f"'{path}' is a file.")
        # Perform actions for a file here
        process_block(path,clear_temp_files=False)


    # Check if the path is a folder
    elif os.path.isdir(path):
        
        logger.debug(f"'{path}' is a folder.")
        # Perform actions for a folder here
        file_name_pattern = "data.raw.h5"
        subfolder_name = "Network"
        result = helper.find_files_with_subfolder(path, file_name_pattern, subfolder_name)
        
        if args.reference:
            
            data_f = args.reference
        else:
            logger.error("Reference file not provided")
            sys.exit(1)
        array_types = args.type
        data_df = pd.read_excel(data_f)
        network_today_runs = data_df[data_df['Assay'].str.lower().isin(array_types)]

        # Extract the 'Run #' values from the filtered DataFrame
        network_today_run_numbers = network_today_runs['Run #'].to_list()
        for path in result:   ##TO DO: check if the run number is in ref file.
            try:
                parts = path.split("/")
                run_id = int(parts[-2])
                if run_id in network_today_run_numbers:
                    _ = process_block(path) 
                else:
                    logger.info(f"{run_id} not a network assay")   
                    continue
            except Exception as e:
                logger.info(e)
                continue
                
    # If it's neither a file nor a folder, display an error message
    else:
        logger.info(f"'{path}' is neither a file nor a folder.")
        sys.exit(1)
    # You can add your code to perform specific actions for files and folders here
    
if __name__ =="__main__" :
    main()