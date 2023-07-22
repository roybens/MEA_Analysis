import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from datetime import datetime
import shutil

##################################################################################
#################################### Actions #####################################
#Set actions here
#set this True to copy and rename the raw data files into a single folder
rename_copies_and_export_compiled_data = False

#set this True to make graphs (Network graphs, Active Electrodes bar and scatters graphs)
make_graphs = True
##################################################################################
##################################################################################


##################################################################################
############################## set custom variables ##############################

##### custom directories #####
#set a directory of where the Exported_Data was stored and where the new files will be added at
prj_dir = 'C:/Users/jonat/Documents/RoyBensLab/MEA/graphs_test'

#set a folder name to store all the renamed copies of electrode_data.csv
electrodeData_dir = '1000_electrode_data'

#set a folder name to store all the graphs
graph_dir = 'Graphs'

#set a file name for the compiled mean data
electrodeData_file_name = 'compiledMeanData.csv'

##### custom chip info #####
#set a dictionary that for Wellplate ID information
#format: {Wellplate ID:[Chip number, Genotype]}
prj_dic = {16657:[7,'WT'],
         16757:[8,'HET'],
         16744:[9,'WT'],
         16742:[10,'HET'],
         16665:[11,'WT'],
         16963:[12,'HET'],
         16850:[20,'WT'],
         16930:[24,'HET'],
         }

##### custom start date #####
#set the date for Day in Vitro 0
div0_date = '2022-08-04' #format 'yyyy-mm-dd'


##### custom average columns #####
#set the desired colomns in electrode_data.csv to take averages
colsToAverage = ['Electrode Spike Count', 'Electrode Firing Rate', 'Electrode Spike Amplitude', 'Mean Electrode ISI', 'Active Electrode'] #which columns to average

##### columns to graphs #####
#set the desired graphs outputs
#only works for ['Electrode Spike Count', 'Electrode Firing Rate', 'Electrode Spike Amplitude', 'Mean Electrode ISI']
colsToGraph = ['Electrode Spike Count', 'Electrode Firing Rate', 'Electrode Spike Amplitude', 'Mean Electrode ISI']
##################################################################################
##################################################################################







##################################################################################
##################################################################################
##################################################################################
#################################### functions ####################################
#make copy of electode_metrics.csv and complie means
def make_copies_compile_edata():
    #set output folder directory and create a output folder if not already exsited
    os.chdir(prj_dir)
    if electrodeData_dir not in os.listdir():
        os.mkdir(electrodeData_dir)
    export_dir = os.path.join(prj_dir,electrodeData_dir)
    #create exported df
    e_df_dic = {'Date':[],
                'Div':[],
                'Assay':[],
                'Chip Number':[],
                'Wellplate ID':[],
                'Genotype':[]}
    for i in range(len(colsToAverage)):
        e_df_dic.update({colsToAverage[i]:[]})
    e_df = pd.DataFrame(data = e_df_dic)
    #walk through all files
    for dirs, subdirs, files, in os.walk(prj_dir):
        for f in files:
            f_name, f_ext = os.path.splitext(f)
            if f_name.endswith('electrode_metrics'):
                #read electrode_metrics.csv and save date, div, wellplate id, chip number, genotype
                os.chdir(dirs)
                df = pd.read_csv(f, skiprows=1).iloc[:,:-4]
                date = df['Date'].unique()[0]
                div = (datetime.strptime(date, "%Y-%m-%d") - datetime.strptime(div0_date, "%Y-%m-%d")).days
                wellid = df['Wellplate ID'][1]
                chipnumb = prj_dic.get(wellid)[0]
                genotype = prj_dic.get(wellid)[1]
                #determine assay type
                assay_type = 'Unknown'
                if 'analysis_params.csv' in os.listdir(dirs):
                    t_df = pd.read_csv('analysis_params.csv', skiprows = 1)
                    assay_type = t_df['Analysis Type'][0]
                #complie means from electrode_metrics.csv
                edf_means = df[colsToAverage].mean()
                #rename and check if file is duplitcated 
                n_date = date[:4] + date[5:7] + date[8:10]
                new_name = f'{n_date}_{f_name}_{assay_type}_{str(chipnumb).zfill(2)}_{genotype}{f_ext}'
                file_duplicate_count = len(e_df.loc[(e_df['Div'] == div) & (e_df['Wellplate ID'] == wellid)])
                if file_duplicate_count > 0: #if file is duplicated, rename according to duplicated number
                    new_name = f'{n_date}_{f_name}_{assay_type}_{str(chipnumb).zfill(2)}_{genotype} ({file_duplicate_count+1}){f_ext}'
                #update it to exported df
                edf_list = [date, div, assay_type, chipnumb, wellid, genotype] + list (edf_means)
                e_df.loc[len(e_df.index)] = edf_list
                #make a copy of renamed electrode_metrics.csv files into a single folder
                olddir = os.path.join(prj_dir, dirs,f) # where the old file was
                destdir = os.path.join(f'{prj_dir}/{electrodeData_dir}',new_name) # where the copy of file will be
                shutil.copy(olddir, destdir) # make copy
    #convert active electrodes to in scale of 100%
    e_df['Active Electrode'] = e_df['Active Electrode']*100
    #Organize df based on Div and Chip Number
    e_df = e_df.sort_values(by=['Div', 'Chip Number']).reset_index(drop=True)
    #Export .csv file
    e_df.to_csv(prj_dir+'/'+electrodeData_file_name)

#plot network graphs
def plot_network_graph(inp):
    #set output folder directory and create a output folder if not already exsited
    os.chdir(prj_dir)
    if graph_dir not in os.listdir():
        os.mkdir(graph_dir)
    export_dir = os.path.join(prj_dir, graph_dir)
    #read csv
    os.chdir(prj_dir)
    df = pd.read_csv(electrodeData_file_name)
    #define input based on output
    if inp == 'Electrode Spike Count':
        out = 'Network Spike Count'
        inp_title = 'Electrode Spike Count'
    if inp == 'Electrode Firing Rate':
        out = 'Network Electrode Firing Rate'
        inp_title = 'Electrode Firing Rate (Hz)'
    if inp == 'Electrode Spike Amplitude':
        out = 'Network Spike Amplitude'
        inp_title = 'Electrode Spike Amplitude (uV)'
    if inp == 'Mean Electrode ISI':
        out = 'Network Electrode ISI'
        inp_title = 'Mean Electrode ISI'
    #create x and y axis
    x_div = df['Div'].unique()
    wt_data = []
    het_data = []
    #fill in x and y axis data
    x = np.arange(len(x_div))
    for i in x_div:
        temp_df = df.loc[(df['Div'] == i) & (df['Genotype'] == 'WT')]
        wt_data.append(np.array(temp_df[inp]))
        temp_df = df.loc[(df['Div'] == i) & (df['Genotype'] == 'HET')]
        het_data.append(np.array(temp_df[inp]))
    #plotting
    fig, ax = plt.subplots()
    #plot WT bar
    ax.bar(x-0.15, 
        height=[np.mean(yi) for yi in wt_data],
        yerr=[np.std(yi) for yi in wt_data],    # error bars
        capsize=3, # error bar cap width in points
        width=0.25,    # bar width
        color=(0,0,0,0),  # face color transparent
        edgecolor='black',
        ecolor='black')
    #plot HET bar
    ax.bar(x+0.15, 
            height=[np.mean(yi) for yi in het_data],
        yerr=[np.std(yi) for yi in het_data],    # error bars
        capsize=3, # error bar cap width in points
        width=0.25,    # bar width
        color=(0,0,0,0),  # face color transparent
        edgecolor='red',
        ecolor='red')
    #plot wt and het scatters
    for i in range(len(wt_data)):
        wt = ax.scatter(x[i]-0.15+np.zeros(wt_data[i].size), wt_data[i], color='black', label='WT', s=20)
    for i in range(len(het_data)):
        het = ax.scatter(x[i]+0.15+np.zeros(het_data[i].size), het_data[i], color='red', label='HET', s=20)
    #axis scaling
    xmin = -1
    xmax = len(x_div)
    ymin = 0
    ymax = (max(df[inp]) - ymin)*1.25
    #labelings
    ax.legend(handles=[wt,het])
    plt.title(out)
    plt.xlabel('Day in Vitro')
    plt.ylabel(inp_title)
    plt.xticks(x, x_div)
    plt.axis([xmin, xmax, ymin, ymax])
    #save plot
    plt.savefig(export_dir +'/' + out +'.png', dpi=300)

def plot_active_electrode_bars():
    #set output folder directory and create a output folder if not already exsited
    os.chdir(prj_dir)
    if graph_dir not in os.listdir():
        os.mkdir(graph_dir)
    export_dir = os.path.join(prj_dir, graph_dir)
    #read csv
    os.chdir(prj_dir)
    df = pd.read_csv('compiledMeanData.csv')
    #create x and y axis
    x_div = df['Div'].unique()
    bar_data = []
    x = np.arange(len(x_div))
    for w_id in prj_dic.keys():
        well_data = []
        for day in x_div:
            ae_df = df.loc[(df['Div'] == day) & (df['Wellplate ID'] == w_id)]
            ae = 0
            if len(ae_df) > 0:
                ae = ae_df['Active Electrode'].unique()[-1]
            well_data.append(ae)
        bar_data.append(well_data)
    #plot bars
    fig, ax = plt.subplots()
    well_ids = list(prj_dic.keys())
    colors = plt.cm.viridis(np.linspace(0, 1, len(well_ids)))
    for i in (range(len(well_ids))):
            plt.bar(x-0.35+i*0.8/len(well_ids),
            bar_data[i],
            width = 0.6/len(well_ids),
            color = colors[i],
            label = well_ids[i])
    #axis scaling
    xmin = -1
    xmax = len(x_div)
    ymin = 0
    ymax = 100
    #labelings
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.175, box.width, box.height * 0.85])
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.125),
            fancybox=True, shadow=True, ncol=5)
    plt.title('Active Electrodes')
    plt.xlabel('Day in Vitro')
    plt.ylabel('Active Electrodes (%)')
    plt.xticks(x, x_div)
    plt.axis([xmin, xmax, ymin, ymax])
    #save plot
    plt.savefig(export_dir +'/' + 'Active Electrodes Bar Chart' + '.png', dpi=300)

def plot_active_electrode_scatter():
    #set output folder directory and create a output folder if not already exsited
    os.chdir(prj_dir)
    if graph_dir not in os.listdir():
        os.mkdir(graph_dir)
    export_dir = os.path.join(prj_dir, graph_dir)
    #read csv
    os.chdir(prj_dir)
    df = pd.read_csv('compiledMeanData.csv')
    #create x and y axis
    x_div = df['Div'].unique()
    scatter_data = []
    x = np.arange(len(x_div))
    for w_id in prj_dic.keys():
        well_data = []
        for day in x_div:
            ae_df = df.loc[(df['Div'] == day) & (df['Wellplate ID'] == w_id)]
            ae = 0
            if len(ae_df) > 0:
                ae = ae_df['Active Electrode'].unique()[-1]
            well_data.append(ae)
        scatter_data.append(well_data)
    #plot scatters
    fig, ax = plt.subplots()
    well_ids = list(prj_dic.keys())
    colors = plt.cm.viridis(np.linspace(0, 1, len(well_ids)))
    for i in (range(len(well_ids))):
        plt.plot(x, scatter_data[i],color = colors[i])
        plt.scatter(x, scatter_data[i],color = colors[i])
    #axis scaling
    xmin = -1
    xmax = len(x_div)
    ymin = -0.05
    ymax = 100
    #labelings
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.175, box.width, box.height * 0.85])
    #plt.legend(well_ids)
    plt.legend(well_ids, loc='upper center', bbox_to_anchor=(0.5, -0.125),
            fancybox=True, shadow=True, ncol=5)
    plt.title('Active Electrodes')
    plt.xlabel('Day in Vitro')
    plt.ylabel('Active Electrodes (%)')
    plt.xticks(x, x_div)
    plt.axis([xmin, xmax, ymin, ymax])
    #save plot
    plt.savefig(export_dir + '/' + 'Active Electrodes Scatter Plot' + '.png', dpi=300)
##################################################################################



##################################################################################
#call functions
if rename_copies_and_export_compiled_data: 
    make_copies_compile_edata()
if make_graphs:
    for i in colsToGraph:
        plot_network_graph(i)
    plot_active_electrode_bars()
    plot_active_electrode_scatter()
