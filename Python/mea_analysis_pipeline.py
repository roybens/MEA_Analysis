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
import helper_functions as helper
from pathlib import Path
from timeit import default_timer as timer
import multiprocessing
import os

FOLDER = ''

job_kwargs = dict(n_jobs=64, chunk_duration="1s", progress_bar=True)


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
    return fs,num_chan,channel_ids

def preprocess(recording):  ## some hardcoded stuff.
    """
    Does the bandpass filtering and the common average referencing of the signal
    
    """
    recording_bp = spre.bandpass_filter(recording, freq_min=300, freq_max=6000)
    recording_bp.annotate(is_filtered=True)

    recording_cmr = spre.common_reference(recording_bp, reference='global', operator='median')

    return recording_cmr

def get_kilosort_result(folder):

    sorter = ss.Kilosort3Sorter._get_result_from_folder(folder)
    return sorter

def get_waveforms_result(sorter,folder):

    waveforms = si.load_waveforms(folder,sorting=sorter)

    return waveforms

def run_kilosort(recording,output_folder):
   
    default_KS3_params = ss.get_default_sorter_params('kilosort3')
    default_KS3_params['keep_good_only'] = True
    default_KS3_params['detect_threshold'] = 24
    default_KS3_params['projection_threshold']=[30, 30]
    default_KS3_params['preclust_threshold'] = 26
    run_sorter = ss.run_kilosort3(recording, output_folder=output_folder, docker_image= "kilosort3-maxwellcomplib:latest",verbose=True, **default_KS3_params)
    sorting_KS3 = ss.Kilosort3Sorter._get_result_from_folder(output_folder+'/sorter_output/')
    return sorting_KS3

def extract_waveforms(recording,sorting_KS3,folder):
    folder = Path(folder)

    waveforms = si.extract_waveforms(recording,sorting_KS3,folder=folder,overwrite=True, ms_before=1., ms_after=2.,allow_unfiltered=True,**job_kwargs)

    return waveforms


def remove_violated_units(sortingKS3, waveforms):

    """
    Removing based on Refractory violations, Firing_rate , snr_ratio
    
    """
    all_units = sortingKS3.get_unit_ids()
    print(f"Total units counts {len(all_units)}")
    print(f"Units : {all_units}")
    rp_contamination,rp_violation = qm.compute_refrac_period_violations(waveforms)
    print(f"rp_contaminaiton : {rp_contamination}")

    refrct_violated_units = [unit for unit,ratio in rp_contamination.items() if ratio >0.0]

    snr_ratio = qm.compute_snrs(waveforms)
    print(f"SNR ratio : {snr_ratio}")
    snr_violated_units = [unit for unit,ratio in snr_ratio.items() if ratio >5]

    firing_rate = qm.compute_firing_rates(waveforms)
    print(f"Firing rate : {firing_rate}")
    fr_violated_units = [unit for unit,ratio in firing_rate.items() if ratio <=0.1]

    violated_units = list(set(refrct_violated_units+fr_violated_units))
    print(f"violated units {violated_units} removed")

    
    good_units = [units for units in all_units if units not in violated_units ]

    return good_units




def analyse_waveforms_sigui(waveforms) :
    import spikeinterface_gui
    # This creates a Qt app
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
    metrics = qm.compute_quality_metrics(waveform,**job_kwargs)

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


def process_block(recnumber,file_path,time_in_s):
    time_start = 0
    time_end = time_start+time_in_s

    #failed_sorting_rec = {}
 
        
    recording,rec_name = get_data_maxwell(file_path,recnumber)
    print(f"Processing recording: {rec_name}")
    fs, num_chan, channel_ids = get_channel_recording_stats(recording)
    recording_chunk = recording.frame_slice(start_frame= time_start*fs,end_frame=time_end*fs)
    recording_chunk = preprocess(recording_chunk)
    start = timer()
    current_directory = os.getcwd()
    try:
        dir_name = '/home/mmpatil/Documents/spikesorting/MEA_Analysis/Python/sorter_workspace/block_'+rec_name
        os.mkdir(dir_name,0o777)
        os.chdir(dir_name)
        kilosort_output_folder = dir_name+'/kilosort3_'+rec_name
        sortingKS3 = run_kilosort(recording_chunk,output_folder=kilosort_output_folder)
        waveform_folder =dir_name+'/waveforms_'+ rec_name
        waveforms = extract_waveforms(recording,sortingKS3,folder = waveform_folder)
    except Exception as e:
        print(e)
        print(f"Error in {rec_name} processing. Continuing to the next block")
        #failed_sorting_rec[rec_name] = e
        x = 2
        return x
    end = timer()
    print("Sort and extract waveforms takes", end - start)
    start = timer()
    non_violated_units = remove_violated_units(sortingKS3,waveforms)
    template_channel_dict = get_unique_templates_channels(non_violated_units,waveforms)
    non_redundant_templates = list(template_channel_dict.keys())
    end = timer()
    print("Removing redundant items takes", end - start)
    start = timer()
    waveform_good = waveforms.select_units(non_redundant_templates,new_folder='waveforms_good_'+rec_name)
    waveform_good.run_extract_waveforms()
    channel_location_dict = get_channel_locations_mapping(recording_chunk)
    # New dictionary with combined information
    new_dict = {}

    # Iterate over each template in the template_channel_dict dictionary
    for template, channel in template_channel_dict.items():

        # If this channel is not already in the new dictionary, add it
        if channel not in new_dict:
            new_dict[channel] = {}

        # Add an entry for this template and its corresponding location to the new dictionary
        new_dict[channel][template] = [int(channel_location_dict[channel][0]/17.5),int(channel_location_dict[channel][1]/17.5)]
    

    os.chdir(current_directory)
    file_name = 'Electrodes_'+rec_name
    helper.dumpdicttofile(new_dict,file_name)
    end = timer()
    print("Saving electrode configuraiton takes: ",end -start)
    x =1
    return 1

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
    


    



if __name__ =="__main__" :

    routine_sequential()  #arguments