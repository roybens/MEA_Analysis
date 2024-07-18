import os

import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt

def read_file(fn):
    df = pd.read_csv(fn)
    tmp_dict = df.to_dict(orient='list')
    for key in tmp_dict.keys():
        tmp_dict[key] = tmp_dict[key]
    return tmp_dict
def make_mea_dicts(fn, wt_chips):
    raw_dict = read_file(fn)
    #print(raw_dict)
    wt_inds = []
    cond_inds = []
    for i,chip_id in enumerate(raw_dict['Wellplate ID']):
        #print(chip_id)
        if chip_id in wt_chips:
            wt_inds.append(i)
        else:
            cond_inds.append(i)
    print(f'wt_inds is:{wt_inds} and cond_inds is:{cond_inds}')
    processed_dict = {}
    for curr_key,curr_val in raw_dict.items():
        wt_vals = []
        cond_vals = []
        for i,elem in enumerate(curr_val):
            if i in wt_inds:
                wt_vals.append(float(elem))
            else:
                cond_vals.append(float(elem))
        processed_dict[curr_key] = [wt_vals,cond_vals]
    return processed_dict
def plot_from_dict(curr_dict, key, axs = None, fig = None, plt_prefix = './Plots/'):
    if not axs:
        fig,axs = plt.subplots(1,1)
    [wt_vals,cond_vals] = curr_dict[key]
    wt_x = range(len(wt_vals))
    cond_x = range(len(wt_vals),len(wt_vals) + len(cond_vals))
    if 'Mean' in key:
        key_name = key[5:]
        std_key = f'{key_name} std'
        cv_key = f'{key_name} CV'
        [wt_err, cond_err] = curr_dict[std_key]
        axs.bar(wt_x, wt_vals,yerr = wt_err, color='black')
        axs.bar(cond_x, cond_vals, yerr = cond_err, color='red')
        [wt_cv, cond_cv] = curr_dict[cv_key]
        all_x = list(wt_x) + list(cond_x)
        all_y = wt_vals + cond_vals
        all_cvs = wt_cv + cond_cv
        print(f'all_x {all_x}, all_y {all_y}, all_cvs {all_cvs}')
        for i in all_x:
            axs.text(all_x[i], all_y[i]*0.2, f'CV={all_cvs[i]}', rotation=90,color = 'white')
    else:
        axs.bar(wt_x, wt_vals, color='black')
        axs.bar(cond_x, cond_vals, color='red')
    axs.set_title(key)
    fn = f'{plt_prefix}{key}.pdf'
    fig.savefig(fn)

def plot_activity_analysis(fn,wt_chip,plt_fldr = './Plots/'):
    os.makedirs(plt_fldr,exist_ok=True)
    act_dict = make_mea_dicts(fn, wt_chips)
    plot_from_dict(act_dict,'Active Area',plt_prefix=plt_fldr)
    plot_from_dict(act_dict,'Mean Firing Rate',plt_prefix=plt_fldr)
    plot_from_dict(act_dict, 'Mean Spike Amplitude',plt_prefix=plt_fldr)
    plot_from_dict(act_dict, 'Mean ISI',plt_prefix=plt_fldr)

def plot_network_analysis(fn,wt_chips,plt_fldr = './Plots/'):
    os.makedirs(plt_fldr, exist_ok=True)
    network_dict = make_mea_dicts(fn, wt_chips)
    plot_from_dict(network_dict, 'Burst Frequency',plt_prefix=plt_fldr)
    plot_from_dict(network_dict, 'Spikes within Bursts',plt_prefix=plt_fldr)
    plot_from_dict(network_dict, 'Mean Spikes per Burst',plt_prefix=plt_fldr)
    plot_from_dict(network_dict, 'Mean Spikes per Burst per Electrode',plt_prefix=plt_fldr)
    plot_from_dict(network_dict, 'Mean Burst Duration',plt_prefix=plt_fldr)
    plot_from_dict(network_dict, 'Mean Burst Peak Firing Rate',plt_prefix=plt_fldr)
    plot_from_dict(network_dict, 'Mean IBI',plt_prefix=plt_fldr)
    plot_from_dict(network_dict, 'Mean ISI within Burst',plt_prefix=plt_fldr)
    plot_from_dict(network_dict, 'Mean ISI outside Burst',plt_prefix=plt_fldr)
def find_all_files(folder,file_name):
    ans = glob.glob(f'{folder}/**/{file_name}', recursive = True)
    return ans

wt_chips = [16364,16461,16384]

fn = '/Users/rbenshalom/Library/CloudStorage/Box-Box/MEA_Data/2022_06_10_AS_Rat/20220613_171016/csv/activity_summary_metrics.csv'
fn = '/Users/rbenshalom/Library/CloudStorage/Box-Box/MEA_Data/2022_06_10_AS_Rat/20220614_131237/csv/network_summary_metrics.csv'
#plot_activity_analysis(fn,wt_chips,plt_fldr='./Plots/Div19_activity/')
plot_network_analysis(fn,wt_chips,plt_fldr='./Plots/Div19_network/')










# # Save the spike array to a npz file
        # np.savez(f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}/spike_array_uncompressed.npz", spike_array=spike_array)

        # # Assuming spike_array is your original array and original_fs is the original sampling frequency

        # # Define the target sampling frequency
        # target_fs = 100

        # # Calculate the resampling factor
        # resampling_factor = int(fs / target_fs)

        # # Calculate the number of bins needed
        # num_bins = spike_array.shape[1] // resampling_factor

        # # Reshape the array to the number of bins and check for any spike in each bin
        # downsampled_spike_array = np.any(spike_array[:, :num_bins * resampling_factor].reshape(spike_array.shape[0], num_bins, resampling_factor), axis=2).astype(int)

        # # Save the downsampled spike array to a npz file
        # np.savez(f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}/spike_array_downsampled.npz", spike_array=downsampled_spike_array)



    # os.makedirs(f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/Spike_trains/",mode=0o777, exist_ok=True)   #need to think of a better optimize d way to save. compare with export to phy
        # for idx, unit_id in enumerate(non_violated_units):
        #     #print(unit_id)
        #     spike_train = sortingKS3.get_unit_spike_train(unit_id,start_frame=0*fs,end_frame=time_end*fs)
        #     #print(spike_train)
        #     if len(spike_train) > 0:
        #         spike_times = spike_train / float(fs)
        #         #rohan made change here
        #         mat_filename = f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/Spike_trains/{unit_id}.mat"
        #         sio.savemat(mat_filename,{'spike_times':spike_times,'units':[unit_id]*len(spike_times)})
        
        
        # channel_location_dict = get_channel_locations_mapping(recording_chunk)
        # # file_name = '/mnt/disk15tb/mmpatil/MEA_Analysis/Python/Electrodes/Electrodes_'+rec_name
        # # helper.dumpdicttofile(new_dict,file_name)
        # electrodes = []

        # # Iterate over each template in the template_channel_dict dictionary
        # for template, channel in unit_extremum_channel.items():
        #     # Add an entry for this template and its corresponding location to the new dictionary
        #     electrodes.append(220* int(channel_location_dict[channel][1]/17.5)+int(channel_location_dict[channel][0]/17.5))
        # electrode_data = {'electrodes':electrodes}
        # #helper.dumpdicttofile(electrode_data,f"{current_directory}/../AnalyzedData/{desired_pattern}/associated_electrodes.json")
        # #rohan made change here
        # sio.savemat(f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/associated_electrodes.mat",electrode_data)