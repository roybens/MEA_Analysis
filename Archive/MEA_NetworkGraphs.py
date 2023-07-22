import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from datetime import datetime

#define file name
f = 'compiledMeanData.xlsx'

def plot_network_graph(file,out): #input: (file name, plot output) eg. ('compiledMeanData.xlsx', 'Network Spike Count')
    df = pd.read_excel(file)

    #define input based on output
    if out == 'Network Spike Count':
        inp = 'Electrode Spike Count'
        inp_title = 'Electrode Spike Count'
    if out == 'Network Electrode Firing Rate':
        inp = 'Electrode Firing Rate'
        inp_title = 'Electrode Firing Rate (Hz)'
    if out == 'Network Spike Amplitude':

        inp = 'Electrode Spike Amplitude'
        inp_title = 'Electrode Spike Amplitude (uV)'
    if out == 'Network Electrode ISI':
        inp = 'Mean Electrode ISI'
        inp_title = 'Mean Electrode ISI'
    
    #convert date to day(s) from start date
    df['start date'] = df['Date'].min() # start date as earliest date
    d = 4 #Assume the first recodring day is the d-th day
    df['Day'] = (df['Date'] - df['start date']) / np.timedelta64(1, 'D') + d # Day as delta-date

    ##form plot data structure
    #day array
    day = df['Day'].unique()
    total_days = len(day)
    #output array
    out_wt = [] #empty wt array
    out_het = [] #empty het array
    for i in day: #fill data from data frame
        temp_df = df.loc[(df['Day'] == i) & (df['Column1'] == 'WT')]
        out_wt.append(np.array(temp_df[inp]))
        temp_df = df.loc[(df['Day'] == i) & (df['Column1'] == 'HET')]
        out_het.append(np.array(temp_df[inp]))

    ##plot
    # bar width
    w = total_days/32
    #create x-coordinates of bars
    x_day = [] #creat x-axis values
    for i in range(0,len(day)):
        x_day.append(int(day[i]))
    x_wt = []
    x_het = []   
    x_d = list(range(1,total_days+1)) #create x-axis bar centers
    for i in x_d:
        x_wt.append(i-w*.75)
        x_het.append(i+w*.75)
    #data series
    y_wt = out_wt
    y_het = out_het

    #plotting
    fig, ax = plt.subplots()
    #plot WT bar
    ax.bar(x_wt, 
           height=[np.mean(yi) for yi in y_wt],
           yerr=[np.std(yi) for yi in y_wt],    # error bars
           capsize=3, # error bar cap width in points
           width=w,    # bar width
           color=(0,0,0,0),  # face color transparent
           edgecolor='black',
           ecolor='black')
    #plot HET bar
    ax.bar(x_het, 
           height=[np.mean(yi) for yi in y_het],
           yerr=[np.std(yi) for yi in y_het],    # error bars
           capsize=3, # error bar cap width in points
           width=w,    # bar width
           color=(0,0,0,0),  # face color transparent
           edgecolor='red',
           ecolor='red')
    #plot wt and het scatters
    for i in range(len(x_wt)):
        wt = ax.scatter(x_wt[i]+np.zeros(y_wt[i].size), y_wt[i], color='black', label='WT', s=20)
    for i in range(len(x_het)):
        het = ax.scatter(x_het[i]+np.zeros(y_het[i].size), y_het[i], color='red', label='HET', s=20)

    #axis scaling
    xmin = 0
    xmax = (max(df['Day']) - xmin)*1.25
    ymin = 0
    ymax = (max(df[inp]) - ymin)*1.25
    #labelings
    ax.legend(handles=[wt,het])
    plt.title(out)
    plt.xlabel('Day in Vitro')
    plt.ylabel(inp_title)
    plt.xticks(x_d, x_day)
    plt.axis([xmin, total_days + 1, ymin, ymax])
    #save plot
    plt.savefig(out +'.png', dpi=300)


#Call the functions based on output types
plot_network_graph(f,'Network Spike Count')
plot_network_graph(f,'Network Electrode Firing Rate')
plot_network_graph(f,'Network Spike Amplitude')
plot_network_graph(f,'Network Electrode ISI')
