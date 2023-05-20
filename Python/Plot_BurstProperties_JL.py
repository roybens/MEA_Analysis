import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# This script reads the compiled csv made by compileNetworkFiles_JL.m and a reference note
# to plot the wt vs. het burst properties overdays

# setting starts here
###############################################################################################################################################

# set data and reference note dir
data_f = '/home/jonathan/Documents/Scripts/Matlab/scripts_output/ADNP/Network_outputs/Compiled_Networks.csv'
reference_f = '/home/jonathan/Documents/Scripts/Python/ADNP_Notes.xlsx'
# set plot saving dir
opDir = '/home/jonathan/Documents/Scripts/Matlab/scripts_output/ADNP/'

# set exclude lists
chip_exclude = []
run_exclude = []

# set the keywords for assay type
#assay_type_keywords = ['today','last','best']
assay_type_keywords = ['today','last']

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
    df = working_df[working_df['Assay'].str.lower().str.contains(assay_type.lower())]
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

    #plotting
    fig, ax = plt.subplots()
    #plot WT bar
    ax.bar(x_wt, 
            height=[np.mean([n for n in yi if np.isfinite(n)]) for yi in y_wt],
            yerr=[np.std([n for n in yi if np.isfinite(n)], ddof=1)/np.sqrt(np.size(yi)) for yi in y_wt],    # error bars
            capsize=3, # error bar cap width in points
            width=w,    # bar width
            color=(0,0,0,0),  # face color transparent
            edgecolor='black',
            ecolor='black')
    #plot HET bar
    ax.bar(x_het, 
            height=[np.mean([n for n in yi if np.isfinite(n)]) for yi in y_het],
            yerr=[np.std([n for n in yi if np.isfinite(n)], ddof=1)/np.sqrt(np.size(yi)) for yi in y_het],    # error bars
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
    plt.savefig(opt_dir + '/' + assay_title + ' ' + title +'.png', dpi=300)

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