from cmath import exp
import os
from sqlite3 import DateFromTicks

import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import shutil


#NETWORK Per-electrode analysis

exportedDataFolder = 'C:/Users/Tim/Documents/Dev/MEA Data/ExampleRawData/' #This is the folder with raw data straight from MXWbio computer

electrodeData_dir = "1000_electrode_data_test/" #This is the name of the folder which will be created to contain all 'electrod_metrics.csv' files as well as averaged data csv


def compile_electrodeData(exportedDataFolder):
    outputfolderpath = os.path.join(exportedDataFolder,electrodeData_dir) #New folder where compiled csvs will be stored
    outputfolder = os.mkdir(outputfolderpath) #make the new folder that we just made the path for
    i=0
    for dirs, subdirs, files, in os.walk(exportedDataFolder): #read all files and folders 
        for f in files:
            f_name, f_ext = os.path.splitext(f) #split names so we can rename
            if f_name.endswith('electrode_metrics'): #only choose electrode_metrics sheets
                new_name = '{}_{}{}'.format(f_name,i,f_ext) #string which will be new name
                oldpath = os.path.join(exportedDataFolder, dirs,f) #where the old file was
                newpath = os.path.join(exportedDataFolder, dirs,new_name) #current path of new file with new file name
                newfilename = os.rename(oldpath,newpath)
                destpath = os.path.join(outputfolderpath,new_name) # final new path of new file/filename
                shutil.copy(newpath, destpath) #copy the renamed file to new destination
        i=i+1
    return exportedDataFolder


def read_electrode_data(exportedDataFolder):
    compile_electrodeData(exportedDataFolder)
    working_folder = os.chdir(os.path.join(exportedDataFolder,electrodeData_dir)) #folder with electrode data .csvs
    files = os.listdir(working_folder) #get list of files in directory
    compiled_df = pd.DataFrame() #empty df to concat to
    colsToAverage = ['Wellplate ID','Electrode Spike Count', 'Electrode Firing Rate', 'Electrode Spike Amplitude', 'Mean Electrode ISI'] #which columns to average
    date = ['Date'] #date column to append to averages.

    for i in range (0,len(files)):
        edf = pd.read_csv(files[i], skiprows=1, parse_dates=['Date']).iloc[:,:-4] #skip top line, drop the last 4 columns, read date
        edf_means = edf[colsToAverage].mean() #average specified columns
        edf_means = pd.concat([edf.iloc[2][date], edf_means]) #add date to top of averaged columns
        compiled_df = pd.concat([compiled_df, edf_means], axis=1) #add each new chip averages to one compiled dataframe

    compiled_df2 = compiled_df.transpose().reset_index(drop=True)
    compiled_df2.sort_values(['Date','Wellplate ID'], inplace=True, ascending=True)
    print(compiled_df2)

    compiled_df2.to_csv('compiledMeanData.csv')

    return compiled_df2

def plot_perElectrode_activity(compiled_df2):

    #plot stuff
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    # ax1.bar(edf['Instance'], edf['Electrode Spike Count'], color='red')
    # ax2.bar(edf['Instance'], edf['Electrode Spike Amplitude'], color='blue')
    # ax3.bar(edf['Instance'], edf['Electrode Firing Rate'], color='green')
    # ax4.bar(edf['Instance'], edf['Mean Electrode ISI'], color='black')
    # plt.show()
    return


#compile_electrodeData(exportedDataFolder)
read_electrode_data(exportedDataFolder)













############################################################################
######Roy's code for graphing Activity and Network
def read_file(fn):
    df = pd.read_csv(fn, skiprows=1) #get rid of top row and use row 2 as index
    df = df.iloc[:,:-4] #drop the last 4 columns of useless data
    df = df.drop(['Date'], axis=1) #drop the 'Date' column because it can't convert date to float later
    tmp_dict = df.to_dict(orient='list')
    for key in tmp_dict.keys():
        tmp_dict[key] = tmp_dict[key]
    return tmp_dict
def make_mea_dicts(fn, wt_chips):
    raw_dict = read_file(fn)
    print(raw_dict)
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

def plot_activity_analysis(fn,wt_chips,plt_fldr = './Plots/'):
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

wt_chips = [16657,16744,16665]

#fn = 'C:/Users/Tim/Documents/Dev/MEA Data/SYNGAP MEA data for scripts 062422/20220624_120702/csv/activity_summary_metrics.csv'
#fn = 'C:/Users/Tim/Documents/Dev/MEA Data/SYNGAP MEA data for scripts 062422/20220624_120937/csv/network_summary_metrics.csv'



#plot_activity_analysis(fn,wt_chips,plt_fldr='./Plots/syngapActivity_TF/')
#plot_network_analysis(fn,wt_chips,plt_fldr='./Plots/syngapNetwork_TF/')

