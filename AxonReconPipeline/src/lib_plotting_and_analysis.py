import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import sem, ttest_ind
import statsmodels.stats.multitest as mt
import logging
import ast

# Logging setup
logger = logging.getLogger(__name__)  # Create a logger
logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler()  # Create handlers, logs to console
stream_handler.setLevel(logging.DEBUG)  # Set level of handlers
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - %(module)s.%(funcName)s')
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)  # Add handlers to the logger

#treat data
def count_all_units_and_branches(data):
    """Count all units and branches in the dataset."""
    
    # Ensure `data` is a list of DataFrames
    if isinstance(data, pd.DataFrame):
        #data = [data]
        pass
    
    total_units = 0
    total_branches = 0
    for df in data:
        df['branch_id'] = df['branch_id'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
        units_per_well = len(df['unit_ids'])
        total_units += units_per_well
        for i in df.index:  # use df.index to handle data indices changing after filtering
            total_branches += len(df['branch_id'][i])
    return total_units, total_branches

def validate_template_densities(data, threshold=0.95):
    """Validate and filter density data based on given criteria."""
    total_units, total_branches = count_all_units_and_branches(data)
    print('Unvalidated:')
    print(f'Total units: {total_units}, Total branches: {total_branches}')
    
    valid_data = []
    for df in data:
        df['channel_density'] = df['channel_density'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
        density_threshold = df['channel_density'].apply(lambda x: x[0]).quantile(threshold)
        valid_df = df[df['channel_density'].apply(lambda x: x[0]) >= density_threshold]
        valid_df.index = range(len(valid_df))  # Reset index
        valid_data.append(valid_df)
    
    valid_units, valid_branches = count_all_units_and_branches(valid_data)
    print('After density validation:')
    print(f'Total units: {valid_units}, Total branches: {valid_branches}')
    
    return valid_data

def validate_branch_lengths_data(data):
    """Validate and filter branch lengths data based on given criteria."""
    #total_units, total_branches = count_all_units_and_branches(data)
    total_units, total_branches = len(data['unit_ids']), len(data['branch_id'])
    print('Before branch length validation:')
    print(f'Total units: {total_units}, Total branches: {total_branches}')
    
    data['length'] = data['length'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    valid_data = data[data['length'].apply(lambda x: all(l > 0 for l in x))]
    
    valid_units, valid_branches = count_all_units_and_branches([valid_data])
    print('After branch length validation:')
    print(f'Total units: {valid_units}, Total branches: {valid_branches}')
    
    return valid_data

def validate_velocity_data(data):
    """Validate and filter velocity data based on given criteria."""
    total_units, total_branches = count_all_units_and_branches(data)
    print('Before velocity validation:')
    print(f'Total units: {total_units}, Total branches: {total_branches}')
    
    data['velocity'] = data['velocity'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    valid_data = data[data['velocity'].apply(lambda x: all(v > 0 for v in x))]
    
    valid_units, valid_branches = count_all_units_and_branches([valid_data])
    print('After velocity validation:')
    print(f'Total units: {valid_units}, Total branches: {valid_branches}')
    
    return valid_data

def validate_number_of_branches_data(data):
    """Validate and filter number of branches data based on given criteria."""
    total_units, total_branches = count_all_units_and_branches(data)
    print('Before number of branches validation:')
    print(f'Total units: {total_units}, Total branches: {total_branches}')
    
    data['branch_id'] = data['branch_id'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    valid_data = data[data['branch_id'].apply(lambda x: len(x) > 0)]
    
    valid_units, valid_branches = count_all_units_and_branches([valid_data])
    print('After number of branches validation:')
    print(f'Total units: {valid_units}, Total branches: {valid_branches}')
    
    return valid_data

# Colors
color_map = {'WT': '#0000FF', 'Het': '#FF0000', 'Homo': '#FFA500'}

def load_data(file_paths, categories, divs):
    """Load data from CSV files and append category and DIV information."""
    data = []
    for file_path, category, div in zip(file_paths, categories, divs):
        df = pd.read_csv(file_path)
        df['Category'] = category
        df['DIV'] = div
        data.append(df)
    return data

def extract_branch_lengths(data):
    """Extract branch lengths for high-density regions with positive velocities."""
    unit_ids, lengths, categories, DIVs = [], [], [], []
    for df in data:
        category, div = df['Category'].iloc[0], df['DIV'].iloc[0]
        #density_threshold = df['channel_density'].apply(eval).apply(lambda x: x[0]).quantile(0.95)
        for unit_id, length_list, velocity_list, density in zip(df['unit_ids'], df['length'], df['velocity'], df['channel_density']):
            #if eval(density)[0] >= density_threshold:
            #filtered_lengths = [l for l, v in zip(eval(length_list), eval(velocity_list)) if v > 0]
            unit_ids.extend([unit_id] * len(length_list))
            lengths.extend(length_list)
            categories.extend([category] * len(length_list))
            DIVs.extend([div] * len(length_list))
    return pd.DataFrame({'unit_ids': unit_ids, 'Length': lengths, 'Category': categories, 'DIV': DIVs})

def extract_branch_velocities(data):
    """Extract branch velocities for high-density regions with positive velocities."""
    velocities, categories, DIVs = [], [], []
    for df in data:
        category, div = df['Category'].iloc[0], df['DIV'].iloc[0]
        # density_threshold = df['channel_density'].apply(eval).apply(lambda x: x[0]).quantile(0.95)
        for velocity_list, density in zip(df['velocity'], df['channel_density']):
            #if eval(density)[0] >= density_threshold:
            #filtered_velocities = [v for v in eval(velocity_list) if v > 0]
            velocities.extend(velocity_list)
            categories.extend([category] * len(velocity_list))
            DIVs.extend([div] * len(velocity_list))
    return pd.DataFrame({'Velocity': velocities, 'Category': categories, 'DIV': DIVs})

def extract_number_of_branches(data):
    """Extract the number of branches with positive velocities for high-density regions."""
    branches, categories, DIVs = [], [], []
    for df in data:
        category, div = df['Category'].iloc[0], df['DIV'].iloc[0]
        #density_threshold = df['channel_density'].apply(eval).apply(lambda x: x[0]).quantile(0.95)
        for branches_list, velocity_list, density in zip(df['branch_id'], df['velocity'], df['channel_density']):
            #if eval(density)[0] >= density_threshold:
            #num_valid_branches = len([v for v in eval(velocity_list) if v > 0])
            branches.append(velocity_list)
            categories.append(category)
            DIVs.append(div)
    return pd.DataFrame({'Branches': branches, 'Category': categories, 'DIV': DIVs})

def identify_outliers(series):
    """Identify outliers using the IQR method."""
    Q1 = series.quantile(0.25)
    Q3 = series.quantile(0.75)
    IQR = Q3 - Q1
    return (series < (Q1 - 1.5 * IQR)) | (series > (Q3 + 1.5 * IQR))

def add_significance_annotations(ax, data, y, group_by):
    """Add significance annotations to the plot."""
    pairs, p_values = [], []

    divs = data['DIV'].unique()
    for div in divs:
        data_div = data[data['DIV'] == div]
        categories = data_div[group_by].unique()
        if len(categories) == 2:
            cat1, cat2 = categories
            data1, data2 = data_div[data_div[group_by] == cat1][y], data_div[data_div[group_by] == cat2][y]
            t_stat, p_val = ttest_ind(data1, data2)
            pairs.append((div, cat1, cat2))
            p_values.append(p_val)

    reject, p_vals_corrected, _, _ = mt.multipletests(p_values, alpha=0.05, method='bonferroni')
    
    for (div, cat1, cat2), p_val, reject_h0 in zip(pairs, p_vals_corrected, reject):
        if reject_h0:
            x1 = list(divs).index(div)
            y_max = data[data['DIV'] == div][y].max()
            y, h, col = y_max + 1, 1, 'k'
            ax.plot([x1 - 0.2, x1 + 0.2], [y + h, y + h], lw=1.5, c=col)
            ax.text(x1, y + h, "*", ha='center', va='bottom', color=col)

def plot_with_significance(data, y, title, ylabel):
    """Plot data with significance annotations and outliers as diamonds."""
    plt.figure(figsize=(12, 8))
    sorted_data = data[data['Category'].isin(['WT', 'Het', 'Homo'])].copy()
    sorted_data['Category'] = pd.Categorical(sorted_data['Category'], categories=['WT', 'Het', 'Homo'], ordered=True)
    
    palette = [color_map[cat] for cat in sorted_data['Category'].cat.categories]

    
    ax = sns.barplot(x='DIV', y=y, hue='Category', data=sorted_data, palette=palette, alpha=0.6, errorbar='se')
    # kwargs = {
    #     'outlier_pop': 2
    # }
    sns.stripplot(x='DIV', y=y, hue='Category', data=sorted_data, jitter=True, 
                  dodge=True, marker='D', alpha=0.6, palette=palette, ax=ax, edgecolor='gray', 
                  #**kwargs
                  )  # marker='D' for diamonds
    
    add_significance_annotations(ax, sorted_data, y, 'Category')
    
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel('DIV')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[0:3], labels[0:3], title='Category')
    plt.show()

    return ax

def count_units_and_branches(data):
    """Count units and branches grouped by DIV and Category."""
    unit_counts = data.groupby(['DIV', 'Category']).size().unstack()
    branch_counts = data.groupby(['DIV', 'Category']).sum().unstack()
    return unit_counts, branch_counts


