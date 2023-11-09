import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import stats
import pdb
import json
from math import sqrt
# This script reads the compiled csv made by compileNetworkFiles_JL.m and a reference note
# to plot the wt vs. het burst properties overdays

# setting starts here
###############################################################################################################################################
# Read the JSON file
with open('act_plt_settings.json', 'r') as file:
    settings_data = json.load(file)
## set plot saving dir
opDir = settings_data['opDir']
# set data and reference note dir
data_f = opDir + 'ActivityScan_outputs/Compiled_ActivityScan.csv'
#data_f = '/home/jonathan/Documents/Scripts/Matlab/scripts_output/CDKL5/ActivityScan_outputs/Compiled_ActivityScan.csv'
reference_f = settings_data['refDir']
#reference_f = '/home/jonathan/Documents/Scripts/Python/CDKL5_Notes.xlsx'

#opDir = '/home/jonathan/Documents/Scripts/Matlab/scripts_output/CDKL5/'

input_string1 = input("Enter the assay type that need to be plotted")
# set exclude lists
if not input_string1:
    print("no  assay value inputted")
    exit(0)
elif len(input_string1.split(','))> 1:
    assay_type_keywords = [x.lower().strip() for x in input_string1.split(',')]
else:
    assay_type_keywords=[input_string1.lower().strip()]

input_string2 = input("Enter comma-separated chip ids to exclude (hit enter if none )")
# set exclude lists
if not input_string2:
    chip_exclude =[]
else :
    try:
        chip_exclude =[x for x in input_string2.split(',')]
        print(f"chips being excluded are {chip_exclude}")
    except Exception:
        print("Invalid input")
input_string3 = input("Enter comma-separated run_ids to exclude (hit enter if none )")
# set exclude lists
if not input_string3:
    run_exclude =[]
else :
    try:
        chip_exclude =[x for x in input_string3.split(',')]
        print(f"chips being excluded are {chip_exclude}")
    except Exception:
        print("Invalid input")

# setting ends here
###############################################################################################################################################

# make output dir if not existed
opt_dir = opDir + 'ActivityScan_outputs/meanActivityProperty_graphs/'
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

df = data_df.assign(Assay=assay_l)
#pdb.set_trace()
print(df.columns)
def plot_network_graph(working_df,output_type, assay_type):
    #extract data based on assay_type
    df = working_df[working_df['Assay'].str.lower().str.contains(assay_type.lower())]
    #create assay title
    assay_title = assay_type.title()
    #define input based on output
    if output_type == 'Mean_FiringRate':
        title = 'Firing Rate'
    if output_type == 'Mean_SpikeAmplitude':
        title = 'Amplitude'

    #div array
    div = df['DIV'].unique()
    total_div = len(div)
    # Find unique values in 'Genotype' column and count them
    unique_genotypes = df['NeuronType'].unique()
    total_genotypes = len(unique_genotypes)

    # Print the number of unique genotypes
    print(f"Number of unique Genotypes: {total_genotypes}")

    # Initialize output arrays for each unique genotype
    output_arrays = {genotype: [] for genotype in unique_genotypes}
    print(unique_genotypes)
    # Fill data from data frame
    for i in div:
        for genotype in unique_genotypes:
            temp_df = df.loc[(df['DIV'] == i) & (df['NeuronType'] == genotype)]
            output_arrays[genotype].append(np.array(temp_df[output_type]))

    ##plot
    # bar width
    w = total_div/32
    #create x-coordinates of bars
   # Create x-coordinates of bars
    x_day = [int(d) for d in div]
    x_genotype = {genotype: [] for genotype in unique_genotypes}
    x_d = list(range(1, total_div + 1))
    
    # Assign x-coordinates for each genotype
    for i, genotype in enumerate(unique_genotypes):
        for x in x_d:
            x_genotype[genotype].append(x + w * (i - len(unique_genotypes) / 2))

    #plotting
    fig, ax = plt.subplots()
    # Generate a list of distinct colors based on the number of genotypes
    colors = [plt.colormaps['Set1'](i) for i in np.linspace(0, 1, len(unique_genotypes))]# Using a colormap to generate colors
    colors2 = [plt.colormaps['Set2'](i) for i in np.linspace(0, 1, len(unique_genotypes))]#
    # Plot data for each genotype
    mean_data_all ={}
    yerr_data_all = {}
    n_data_all={}
    #plot WT bar
    for i,genotype in enumerate(unique_genotypes):
            y_data = output_arrays[genotype]
            #print("type: ",type(genotype))
            # Calculate statistics
            mean_data = [np.mean([n for n in yi if np.isfinite(n)]) for yi in y_data]
            yerr_data = [np.std([n for n in yi if np.isfinite(n)], ddof=1)/np.sqrt(np.size(yi)) for yi in y_data]
            n_data = [len(yi) for yi in y_data]
            # Store statistics in dictionaries
            mean_data_all[genotype] = mean_data
            yerr_data_all[genotype] = yerr_data
            n_data_all[genotype] = n_data
            # Save statistics to file
            output_file = f"intermediate_files/{title}_{genotype}_statistics.txt"
            with open(output_file, 'w') as file:
                file.write(f"{genotype} Statistics\n")
                file.write("Mean: " + ", ".join([str(m) for m in mean_data]) + "\n")
                file.write("SEM: " + ", ".join([str(sem) for sem in yerr_data]) + "\n")
                file.write("Sample Size (n): " + ", ".join([str(n_data)]) + "\n")

            # Plot bars
            ax.bar(x_genotype[genotype],
                height=mean_data,
                yerr=yerr_data,
                capsize=3,
                width=w,
                color=colors[i],
                edgecolor='black',
                ecolor='black',label=genotype)

            # Plot scatter points
            for j in range(len(x_genotype[genotype])):
                ax.scatter(x_genotype[genotype][j] + np.zeros(y_data[j].size), y_data[j], s=20,color=colors2[i])

    #perform ttest
    for i in range(len(x_d)):
        maxim = max([max( output_arrays[genotype][i] )for genotype in unique_genotypes])
        count = 1
        p_values = []
        for j, genotype1 in enumerate(unique_genotypes):
            for k, genotype2 in enumerate(unique_genotypes):
                if j < k:
                    #pdb.set_trace()
                    #print("mean_data_all",mean_data_all[genotype1])
                    #print("type:",type(genotype1))
                    mean1, sem1, n1 = mean_data_all[genotype1][i], yerr_data_all[genotype1][i], n_data_all[genotype1][i]
                    mean2, sem2, n2 = mean_data_all[genotype2][i], yerr_data_all[genotype2][i], n_data_all[genotype2][i]
                    #t_stat, p_value = stats.ttest_ind_from_stats(mean1, sem1, n1, mean2, sem2, n2)
                    sed = sqrt(sem1**2.0 + sem2**2.0)
                    t_stat = (mean1 - mean2) / sed
                    # degrees of freedom
                    degreef = n1+n2 - 2
                    alpha=0.05
                    # calculate the critical value
                    cv = stats.t.ppf(1.0 - alpha, degreef)
                    # calculate the p-value
                    p_value = (1.0 - stats.t.cdf(abs(t_stat), degreef)) * 2.0
                    p_values.append([mean1,sem1,mean2,sem2,p_value])

                    # Plot significance
                    #maxim = max(np.max(output_arrays[genotype1][i]), np.max(output_arrays[genotype2][i]))
                    x1, x2 = x_genotype[genotype1][i], x_genotype[genotype2][i]
                    ax.plot([x1, x2], [maxim + 0.1*maxim*(count)] * 2, 'k', linewidth=1.5)
                    sign = "***" if p_value <= 0.001 else "**" if p_value <= 0.01 else "*" if p_value <= 0.05 else "ns"
                    
                    ax.text((x1 + x2) / 2, maxim +0.1*maxim*(count), sign, ha='center', va='bottom', fontsize=10)
                    count = count +1
    
                    with open(output_file, 'a') as file:
                                file.write(f"P values:{p_values} \n")

    # Axis scaling and labeling
    xmin = 0
    xmax = (max(df['DIV']) - xmin)*1.25
    ymin = 0
    ymax = (max(df[output_type]) - ymin)*1.4

    plt.title(assay_title + ' ' + title)
    plt.xlabel('DIV')
    plt.ylabel(title)
    plt.xticks(x_d, x_day)
    plt.axis([xmin, total_div + 1, ymin, ymax])
    plt.legend(title='type',loc='upper right')


    svg_dir = os.path.join(opt_dir, 'svg')
    jpg_dir = os.path.join(opt_dir, 'jpg')

    # Create the directories if they do not exist
    os.makedirs(svg_dir, exist_ok=True)
    os.makedirs(jpg_dir, exist_ok=True)

    # Now, save the figures to the specified format and directory
    plt.savefig(os.path.join(svg_dir, f"{assay_title} {title}.pdf"), dpi=300, format='svg')
    plt.savefig(os.path.join(jpg_dir, f"{assay_title} {title}.jpg"), dpi=300, format='jpg')
    return fig

#exclude chip ids and runs that are in the exclude list
exclude_l = []
for i in df.index:
    if df['Chip_ID'][i] in chip_exclude or df['Run_ID'][i] in run_exclude:
        df = df.drop(index = i)

#Run grapher
#assay_types = list(df['Assay'].unique())
for i in assay_type_keywords:
    #working_df = df[df['Assay'].str.lower().str.contains(i.lower())]
    plot_network_graph(df, 'Mean_FiringRate',i)
    plot_network_graph(df, 'Mean_SpikeAmplitude',i)
 