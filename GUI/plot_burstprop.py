import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import stats
import json

# This script reads the compiled csv made by compileNetworkFiles_JL.m and a reference note
# to plot the wt vs. het burst properties overdays

# setting starts here
###############################################################################################################################################



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import stats
# This script reads the compiled csv made by compileNetworkFiles_JL.m and a reference note
# to plot the wt vs. het burst properties overdays

# Read the JSON file
with open('net_plt_settings.json', 'r') as file:
    settings_data = json.load(file)
# set plot saving dir
opDir = settings_data['opDir']
# set data and reference note dir
data_f = opDir + 'Network_outputs/Compiled_Networks.csv'
#data_f = '/home/jonathan/Documents/Scripts/Matlab/scripts_output/CDKL5/ActivityScan_outputs/Compiled_ActivityScan.csv'
reference_f = settings_data['refDir']
#reference_f = '/home/jonathan/Documents/Scripts/Python/CDKL5_Notes.xlsx'

#opDir = '/home/jonathan/Documents/Scripts/Matlab/scripts_output/CDKL5/'

# set exclude lists
if settings_data['excludeChips'] != '':
    chip_exclude = settings_data['excludeChips']
else :
    chip_exclude =[]
if settings_data['excludeRunIDs']!='':
    run_exclude = settings_data['excludeRunIDs']
else:
    run_exclude =[]
#print(chip_exclude)
#print(run_exclude)
#if chip_exclude[0] =='':
    #  chip_exclude = []
#if run_exclude[0] =='':
    #  run_exclude =[]


# set the keywords for assay type
assay_type_keywords = settings_data['assayTypes']
#print(assay_type_keywords)
#print(type(assay_type_keywords))
# setting ends here
###############################################################################################################################################

# make output dir if not existed
opt_dir = opDir + 'Network_outputs/BurstProperty_graphs/'
if not os.path.exists(opt_dir):
    os.makedirs(opt_dir)

# read files
data_df = pd.read_csv(data_f)
ref_df = pd.read_excel(reference_f)

# sort df based on Run IDs and reindex
data_df = data_df.sort_values(by=['Run_ID'],ascending=True,ignore_index=True)
ref_df = ref_df.sort_values(by=['Run #'],ascending=True,ignore_index=True)

#To do : need to check why NAN
#data_df = data_df['DataFrame Column'].fillna(0)
data_df = data_df.replace(np.NaN,0.0)
print(data_df)

#combine data df with reference note
assay_l = []
genotype_l = []
for i in data_df.index:
    if data_df.loc[i]['Run_ID'] in list(ref_df['Run #']):
        temp_df = ref_df[ref_df['Run #'] == data_df.loc[i]['Run_ID']]
        assay = str(temp_df['Assay'].unique()[0])
        assay_l.append(assay)
        genotype = str(temp_df['Neuron Source'].unique()[0])
        genotype_l.append(genotype)
assay_l = [x.lower() for x in assay_l]
genotype_l = [x.lower() for x in genotype_l]

df = data_df.assign(Assay=assay_l, Genotype=genotype_l)

def plot_network_graph(working_df,output_type, assay_type):
    #extract data based on assay_type
    #print(working_df)
    #print(assay_type)
    df = working_df[working_df['Assay'].str.lower().str.contains(assay_type.lower())]
    #print(df)
    #create assay title
    assay_title = 'Network ' + assay_type.title()
    #assay_type = df['Assay'].unique()[0].title()
    #define input based on output
    if output_type == 'IBI':
        title = 'IBI'
    if output_type == 'Burst_Peak':
        title = 'Burst Peak'
    if output_type == 'Number_Bursts':
        title = 'Number Bursts'
    if output_type == 'Spike_per_Burst':
        title = 'Spike per Burst'

    #div array
    div = df['DIV'].unique()
    total_div = len(div)
    #output array
    wt = []
    het = []
    for i in div: #fill data from data frame
        temp_df = df.loc[(df['DIV'] == i) & (df['Genotype'] == 'wt cortex')]
        wt.append(np.array(temp_df[output_type]))
        temp_df = df.loc[(df['DIV'] == i) & (df['Genotype'] == 'het cortex')]
        het.append(np.array(temp_df[output_type]))

    ##plot
    # bar width
    w = total_div/32
    #create x-coordinates of bars
    x_day = [] #creat x-axis values
    for i in range(0,len(div)):
        x_day.append(int(div[i]))
    x_wt = []
    x_het = []   
    x_d = list(range(1,total_div+1)) #create x-axis bar centers
    for i in x_d:
        x_wt.append(i-w*.75)
        x_het.append(i+w*.75)
    #data series
    y_wt = wt
    y_het = het
    file_path = title + '_y_wt.txt'
    y_wt_list = [list(row) for row in y_wt]

    with open(file_path, 'w') as f:
        for row in y_wt_list:
            f.write('\t'.join(map(str, row)) + '\n')
    
    file_path = title + '_y_het.txt'

    y_het_list = [list(row) for row in y_het]

    with open(file_path, 'w') as f:
        for row in y_het_list:
            f.write('\t'.join(map(str, row)) + '\n')
    #do some ttest
   
    #plotting
    fig, ax = plt.subplots()
    #plot WT bar
    mean_wt = [np.mean([n for n in yi if np.isfinite(n)]) for yi in y_wt]
    yerr_wt = [np.std([n for n in yi if np.isfinite(n)], ddof=1)/np.sqrt(np.size(yi)) for yi in y_wt]
    n_wt = [len(yi) for yi in y_wt]
    mean_het = [np.mean([n for n in yi if np.isfinite(n)]) for yi in y_het]
    yerr_het = [np.std([n for n in yi if np.isfinite(n)], ddof=1)/np.sqrt(np.size(yi)) for yi in y_het]
    n_het = [len(yi) for yi in y_het]
    output_file = title + '_wt_statistics.txt'
    with open(output_file, 'w') as file:
        file.write("WT Statistics\n")
        file.write("Mean: " + ", ".join([str(m) for m in mean_wt]) + "\n")
        file.write("SEM: " + ", ".join([str(sem) for sem in yerr_wt]) + "\n")
        file.write("Sample Size (n): " + ", ".join([str(n) for n in n_wt]) + "\n")
    output_file = title + '_het_statistics.txt'
    with open(output_file, 'w') as file:
        file.write("HeT Statistics\n")
        file.write("Mean: " + ", ".join([str(m) for m in mean_het]) + "\n")
        file.write("SEM: " + ", ".join([str(sem) for sem in yerr_het]) + "\n")
        file.write("Sample Size (n): " + ", ".join([str(n) for n in n_het]) + "\n")
    
    p_values =[]

    for i in range(len(mean_wt)):

        t_stat,p_value = stats.ttest_ind_from_stats(mean_wt[i],yerr_wt[i],n_wt[i],
                                                    
                                                    mean_het[i],yerr_het[i],n_het[i])
        p_values.append(p_value)
    
    ax.bar(x_wt, 
            height=mean_wt,
            yerr= yerr_wt,    # error bars
            capsize=3, # error bar cap width in points
            width=w,    # bar width
            color=(0,0,0,0),  # face color transparent
            edgecolor='black',
            ecolor='black')
    #plot HET bar
    ax.bar(x_het, 
            height=mean_het,
            yerr=yerr_het,    # error bars
            capsize=3, # error bar cap width in points
            width=w,    # bar width
            color=(0,0,0,0),  # face color transparent
            edgecolor='red',
            ecolor='red')
    #plot wt and het scatters
    for i in range(len(x_wt)):
        wt_scatter = ax.scatter(x_wt[i]+np.zeros(y_wt[i].size), y_wt[i], color='black', label='WT', s=20)
    for i in range(len(x_het)):
        het_scatter = ax.scatter(x_het[i]+np.zeros(y_het[i].size), y_het[i], color='red', label='HET', s=20)
    for i in range(len(x_wt)):
        # wt_data = [n for n in y_wt[i] if np.isfinite(n)]
        # het_data = [n for n in y_het[i] if np.isfinite(n)]
        # t_stat, p_value = stats.ttest_ind(wt_data, het_data)

        maxim = max(np.max(y_wt[i]), np.max(y_het[i]))
        p_value = p_values[i]
        if p_value > 0.05:
            ax.plot([x_wt[i], x_het[i]], [maxim + 0.1*maxim] * 2, 'k', linewidth=1.5)
            ax.text((x_wt[i] + x_het[i]) / 2, maxim + 0.15*maxim, "ns", ha='center', va='bottom', fontsize=10)
            continue
        elif p_value <= 0.001 :
            ax.plot([x_wt[i], x_het[i]], [maxim+ 0.1*maxim] * 2, 'k', linewidth=1.5)
            ax.text((x_wt[i] + x_het[i]) / 2, maxim + 0.15*maxim, "***", ha='center', va='bottom', fontsize=10)
            continue
        elif p_value <= 0.01 :
            ax.plot([x_wt[i], x_het[i]], [maxim+ 0.1*maxim] * 2, 'k', linewidth=1.5)
            ax.text((x_wt[i] + x_het[i]) / 2, maxim + 0.15*maxim, "**", ha='center', va='bottom', fontsize=10)
            continue
        elif p_value <= 0.05 :
            ax.plot([x_wt[i], x_het[i]], [maxim+ 0.1*maxim] * 2, 'k', linewidth=1.5)
            ax.text((x_wt[i] + x_het[i]) / 2, maxim + 0.15*maxim, "*", ha='center', va='bottom', fontsize=10)
            continue
    df.replace(np.inf,)
    #axis scaling
    xmin = 0
    xmax = (max(df['DIV']) - xmin)*1.25
    ymin = 0
    ymax = (max(df[output_type]) - ymin)*1.25
    #labelings
    ax.legend(handles=[wt_scatter,het_scatter])
    plt.title(assay_title + ' ' + title)
    plt.xlabel('DIV')
    plt.ylabel(title)
    plt.xticks(x_d, x_day)
    plt.axis([xmin, total_div + 1, ymin, ymax])

    #save plot
    plt.savefig(opt_dir + '/' + assay_title + ' ' + title +'.pdf', dpi=300)

#exclude chip ids and runs that are in the exclude list
exclude_l = []
for i in df.index:
    if df['Chip_ID'][i] in chip_exclude or df['Run_ID'][i] in run_exclude:
        df = df.drop(index = i)

#Run grapher
#assay_types = list(df['Assay'].unique())
for i in assay_type_keywords:
    #working_df = df[df['Assay'].str.lower().str.contains(i.lower())]
    plot_network_graph(df, 'IBI',i)
    plot_network_graph(df, 'Burst_Peak',i)
    plot_network_graph(df, 'Number_Bursts',i)
    plot_network_graph(df, 'Spike_per_Burst',i)


