import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multicomp as multi

'''logging setup'''
import logging
logger = logging.getLogger(__name__) #Create a logger
logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler()  # Create handlers, logs to console
stream_handler.setLevel(logging.DEBUG) # Set level of handlers
logger.addHandler(stream_handler) # Add handlers to the logger
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - %(module)s.%(funcName)s') # Create formatters and add it to handlers
stream_handler.setFormatter(formatter)

# Function to extract lengths and categorize wells
def test_oneway_ANOVA_and_annotate_plot(data):
    # Perform one-way ANOVA
    f_val, p_val = stats.f_oneway(data[data['Category'] == 'Homo']['Length'],
                                  data[data['Category'] == 'HET']['Length'],
                                  data[data['Category'] == 'WT']['Length'])
    print('One-way ANOVA result: F =', f_val, ', p =', p_val)

    # If p_val < 0.05, perform post-hoc test
    if p_val < 0.05:
        print('Performing Tukey HSD test...')
        print('Null hypothesis: All group means are equal.')
        mc = multi.MultiComparison(data['Length'], data['Category'])
        result = mc.tukeyhsd()
        print(result)

        # Annotate plot with significance
        categories = ['Homo', 'HET', 'WT']
        category_coordinates = {category: i for i, category in enumerate(categories)}
        y_max = data['Length'].max()

        sig_height = None
        y_diff = None
        for i, reject in enumerate(result.reject):
            if reject:
                # Get the coordinates of the categories being compared
                group1, group2 = mc.groupsunique[np.triu_indices(len(mc.groupsunique), 1)[0][i]], mc.groupsunique[np.triu_indices(len(mc.groupsunique), 1)[1][i]]
                coord1, coord2 = category_coordinates[group1], category_coordinates[group2]

                # Draw a horizontal line between the pairs
                if y_diff is None: ymin, ymax = plt.ylim() #avoid scaling difference between sigbars through multiple loop iters
                y_diff = ymax-ymin
                if sig_height is None: sig_height = y_max + (i+1)*10
                else: sig_height = sig_height + y_diff*.10 #scale by plot height
                if sig_height >= ymax: plt.ylim(ymax*1.10, ymin)
                plt.plot([coord1, coord2], [sig_height, sig_height], color='k')
                # Place the asterisk at the center of the line
                plt.text((coord1 + coord2) / 2, sig_height, '*', ha='center', va='bottom', color='k', fontsize=20)
    # Call the function after plotting
    #perform_stat_tests_and_annotate_plot(data)
def extract_velocities(dfs, min_branches=1, dense_percentile=0):
    
    lengths = []
    categories = []
    unit_ids = []     
    for i, df in enumerate(dfs):
        # Determine category based on well index
        if i % 3 == 0:
            category = 'Homo'
        elif i % 3 == 1:
            category = 'HET'
        else:
            category = 'WT'

        #get max channel_density
        densities = []
        for j, density in enumerate(df['channel_density']):
            densities.append(eval(density)[0])
        max_channel_density = max(densities)       
        
        # Extract lengths and append to the lists
        for j, vel_list in enumerate(df['velocity']):
            if len(eval(vel_list))<min_branches: continue # exclude units with only one branch
            if eval(df['channel_density'][j])[0] < max_channel_density*dense_percentile: continue #exclude low density reconstructions
            unit_ids.extend(eval(df['unit_ids'][j])*len(eval(vel_list)))
            lengths.extend(eval(vel_list))  # Convert string representation of list to actual list
            categories.extend([category] * len(eval(vel_list)))
            
            # Check for negative values in vel_list
            for vel in eval(vel_list):
                if vel < 0:
                    print(f"Negative value found in unit_id: {eval(df['unit_ids'][j])}, stream_id: {i}")
            
    return lengths, categories, unit_ids
def extract_lengths(dfs, min_branches=1, dense_percentile=0):
    
    lengths = []
    categories = []
    unit_ids = []     
    for i, df in enumerate(dfs):
        # Determine category based on well index
        if i % 3 == 0:
            category = 'Homo'
        elif i % 3 == 1:
            category = 'HET'
        else:
            category = 'WT'

        #get max channel_density
        densities = []
        for j, density in enumerate(df['channel_density']):
            densities.append(eval(density)[0])
        max_channel_density = max(densities)       
        
        # Extract lengths and append to the lists
        for j, length_list in enumerate(df['length']):
            if len(eval(length_list))<min_branches: continue # exclude units with only one branch
            if eval(df['channel_density'][j])[0] < max_channel_density*dense_percentile: continue #exclude low density reconstructions
            unit_ids.extend(eval(df['unit_ids'][j])*len(eval(length_list)))
            lengths.extend(eval(length_list))  # Convert string representation of list to actual list
            categories.extend([category] * len(eval(length_list)))
            
    return lengths, categories, unit_ids
def plot_branch_lengths_by_gene(file_paths, min_branches=0, dense_percentile=0, ylim=(0, 1000)):
    # Extract lengths and categories
    lengths, categories, unit_ids = extract_lengths(file_paths, min_branches=min_branches, dense_percentile=dense_percentile)

    # Create a DataFrame for plotting
    data = pd.DataFrame({
        'Length': lengths,
        'Category': categories,
        'unit_ids': unit_ids
    })

    # Plotting
    plt.figure(figsize=(12, 8))
    # Set custom pastel colors for seaborn plots
    pastel_colors = ['#FFB6C1', '#ADD8E6', '#90EE90']

    plt.figure(figsize=(12, 8))
    # Boxplot with custom pastel colors
    sns.boxplot(x='Category', y='Length', data=data, showmeans=False, meanline=False, showfliers=False, palette=pastel_colors)

    # Stripplot with jitter and color coding based on unit_id
    unique_unit_ids = data['unit_ids'].unique()
    palette = sns.color_palette('husl', len(unique_unit_ids))  # Create a color palette based on the number of unique unit_ids
    unit_id_color_map = dict(zip(unique_unit_ids, palette))

    sns.stripplot(x='Category', y='Length', data=data, jitter=True, alpha=0.5,
                  palette=unit_id_color_map, hue='unit_ids', dodge=True, legend=False)

    plt.title('Branch Length Comparison by Well Category')
    plt.xlabel('Category')
    plt.ylabel('Branch Length')
    #plt.ylim(50, 1600)
    plt.ylim(ylim)
    plt.grid(False)
    #plt.legend(title='unit_ids', bbox_to_anchor=(1.05, 1), loc='upper left')  # Adjust legend position
    test_oneway_ANOVA_and_annotate_plot(data)
    plt.show()
def plot_velocities_by_gene(file_paths, min_branches=0, dense_percentile=0, ylim=(0, 1000)):
    # Extract lengths and categories
    lengths, categories, unit_ids = extract_velocities(file_paths, min_branches=min_branches, dense_percentile=dense_percentile)

    # Create a DataFrame for plotting
    data = pd.DataFrame({
        'Length': lengths,
        'Category': categories,
        'unit_ids': unit_ids
    })

    # Plotting
    plt.figure(figsize=(12, 8))
    # Set custom pastel colors for seaborn plots
    pastel_colors = ['#FFB6C1', '#ADD8E6', '#90EE90']

    plt.figure(figsize=(12, 8))
    # Boxplot with custom pastel colors
    sns.boxplot(x='Category', y='Length', data=data, showmeans=False, meanline=False, showfliers=False, palette=pastel_colors)

    # Stripplot with jitter and color coding based on unit_id
    unique_unit_ids = data['unit_ids'].unique()
    palette = sns.color_palette('husl', len(unique_unit_ids))  # Create a color palette based on the number of unique unit_ids
    unit_id_color_map = dict(zip(unique_unit_ids, palette))

    sns.stripplot(x='Category', y='Length', data=data, jitter=True, alpha=0.5,
                  palette=unit_id_color_map, hue='unit_ids', dodge=True, legend=False)

    plt.title('Axon Velocity Comparison by Well Category')
    plt.xlabel('Category')
    plt.ylabel('Axon Velocity')
    #plt.ylim(-200, 1000)
    plt.ylim(ylim)
    plt.grid(False)
    test_oneway_ANOVA_and_annotate_plot(data)
    #plt.legend(title='unit_ids', bbox_to_anchor=(1.05, 1), loc='upper left')  # Adjust legend position
    plt.show()