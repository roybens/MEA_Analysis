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

## Old fuctions

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

## Data Extraction Functions
def load_data(file_paths, categories, divs):
    data = []
    for file_path, category, div in zip(file_paths, categories, divs):
        df = pd.read_csv(file_path)
        df['Category'] = category
        df['DIV'] = div
        data.append(df)
    return data

def extract_branch_lengths(data):
    lengths = []
    categories = []
    DIVs = []
    for df in data:
        category = df['Category'].iloc[0]
        div = df['DIV'].iloc[0]
        density_threshold = df['channel_density'].apply(eval).apply(lambda x: x[0]).quantile(0.95)
        for length_list, velocity_list, density in zip(df['length'], df['velocity'], df['channel_density']):
            if eval(density)[0] >= density_threshold:
                filtered_lengths = [l for l, v in zip(eval(length_list), eval(velocity_list)) if v > 0]  # Filter out lengths with non-positive velocities
                lengths.extend(filtered_lengths)
                categories.extend([category] * len(filtered_lengths))
                DIVs.extend([div] * len(filtered_lengths))
    return pd.DataFrame({'Length': lengths, 'Category': categories, 'DIV': DIVs})

def extract_branch_velocities(data):
    velocities = []
    categories = []
    DIVs = []
    for df in data:
        category = df['Category'].iloc[0]
        div = df['DIV'].iloc[0]
        density_threshold = df['channel_density'].apply(eval).apply(lambda x: x[0]).quantile(0.95)
        for velocity_list, density in zip(df['velocity'], df['channel_density']):
            if eval(density)[0] >= density_threshold:
                filtered_velocities = [v for v in eval(velocity_list) if v > 0]  # Filter out non-positive velocities
                velocities.extend(filtered_velocities)
                categories.extend([category] * len(filtered_velocities))
                DIVs.extend([div] * len(filtered_velocities))
    return pd.DataFrame({'Velocity': velocities, 'Category': categories, 'DIV': DIVs})

def extract_number_of_branches(data):
    branches = []
    categories = []
    DIVs = []
    for df in data:
        category = df['Category'].iloc[0]
        div = df['DIV'].iloc[0]
        density_threshold = df['channel_density'].apply(eval).apply(lambda x: x[0]).quantile(0.95)
        for branches_list, velocity_list, density in zip(df['branch_id'], df['velocity'], df['channel_density']):
            if eval(density)[0] >= density_threshold:
                num_valid_branches = len([v for v in eval(velocity_list) if v > 0])  # Count branches with positive velocities
                branches.append(num_valid_branches)
                categories.append(category)
                DIVs.append(div)
    return pd.DataFrame({'Branches': branches, 'Category': categories, 'DIV': DIVs})



import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import sem, ttest_ind
import statsmodels.stats.multitest as mt

def add_significance_annotations(ax, data, y, group_by):
    pairs = []
    p_values = []
    
    # Compare WT vs HET within each DIV
    divs = data['DIV'].unique()
    for div in divs:
        data_div = data[data['DIV'] == div]
        categories = data_div[group_by].unique()
        if len(categories) == 2:  # Ensure there are exactly two categories to compare
            cat1, cat2 = categories
            data1 = data_div[data_div[group_by] == cat1][y]
            data2 = data_div[data_div[group_by] == cat2][y]
            t_stat, p_val = ttest_ind(data1, data2)
            pairs.append((div, cat1, cat2))
            p_values.append(p_val)
            print(f'Comparing {cat1} vs {cat2} at DIV {div}: t-statistic={t_stat}, p-value={p_val}')
    
    # Adjust p-values for multiple comparisons
    reject, p_vals_corrected, _, _ = mt.multipletests(p_values, alpha=0.05, method='bonferroni')
    
    # Annotate plot with significance lines and stars
    for (div, cat1, cat2), p_val, reject_h0 in zip(pairs, p_vals_corrected, reject):
        if reject_h0:
            x1 = divs.tolist().index(div)
            y_max = data[data['DIV'] == div][y].max()
            y, h, col = y_max + 1, 1, 'k'
            ax.plot([x1 - 0.2, x1 + 0.2], [y + h, y + h], lw=1.5, c=col)
            ax.text(x1, y + h, "*", ha='center', va='bottom', color=col)
            print(f'Significant difference between {cat1} and {cat2} at DIV {div} after Bonferroni correction: adjusted p-value={p_val}')

def plot_with_significance(data, y, title, ylabel):
    plt.figure(figsize=(12, 8))
    ax = sns.barplot(x='DIV', y=y, hue='Category', data=data, palette=['#FF0000', '#0000FF'], alpha=0.6, errorbar='sd')
    
    sns.stripplot(x='DIV', y=y, hue='Category', data=data, jitter=True, dodge=True, marker='o', alpha=0.6, palette=['#FF0000', '#0000FF'], ax=ax, edgecolor='gray')
    
    # means = data.groupby(['DIV', 'Category'])[y].mean().unstack()
    # errors = data.groupby(['DIV', 'Category'])[y].apply(sem).unstack()
    
    # for i, div in enumerate(means.index):
    #     for j, category in enumerate(means.columns):
    #         ax.errorbar(x=i + (j - 0.5) * 0.25, y=means.loc[div, category], yerr=errors.loc[div, category], fmt='none', c='black', capsize=5)
    
    add_significance_annotations(ax, data, y, 'Category')
    
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel('DIV')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[0:2], labels[0:2], title='Category')
    plt.show()

    return ax

def count_units_and_branches(data):
    unit_counts = data.groupby(['DIV', 'Category']).size().unstack()
    branch_counts = data.groupby(['DIV', 'Category']).sum().unstack()
    return unit_counts, branch_counts

# # Count units and branches included in each DIV for WT and HET
# unit_counts_lengths, branch_counts_lengths = count_units_and_branches(branch_lengths_data)
# unit_counts_velocities, branch_counts_velocities = count_units_and_branches(branch_velocities_data)

# print("Unit counts (Lengths):")
# print(unit_counts_lengths)

# print("\nBranch counts (Lengths):")
# print(branch_counts_lengths)

# print("\nUnit counts (Velocities):")
# print(unit_counts_velocities)

# print("\nBranch counts (Velocities):")
# print(branch_counts_velocities)

# print("\nUnit counts (Number of Branches):")
# print(unit_counts_branches)

# print("\nBranch counts (Number of Branches):")
# print(branch_counts_branches)


