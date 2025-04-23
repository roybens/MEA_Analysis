import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import sem, ttest_ind
import statsmodels.stats.multitest as mt
import logging
import ast
import pingouin as pg

# Logging setup
#reset logger
import logging
logger = logging.getLogger(__name__)  # Create a logger
logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler()  # Create handlers, logs to console
#stream_handler.setLevel(logging.DEBUG)  # Set level of handlers
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - %(module)s.%(funcName)s')
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)  # Add handlers to the logger

'''Analyze and Treat Data'''
def load_data(file_paths, categories, divs):
    """Load data from CSV files and append category and DIV information."""
    data = []
    for file_path, category, div in zip(file_paths, categories, divs):
        df = pd.read_csv(file_path)
        df.insert(0, 'category', category) #Insert category column at position 0
        df.insert(1, 'DIV', div)
        df = df.to_dict()
        #eval all lists in the dataframe
        for key, value in df.items(): 
            #if isinstance(value, str):
                #df[key] = [eval(val) for val in value]
            #if isinstance(value, dict):
            check = True
            try: eval(df[key][0]) #check if eval works
            except: check = False
            if check: df[key] = [eval(val) for val in value.values() if isinstance(val, str) and val.startswith('[')]
        data.append(df)
    return data

def count_units_and_branches(data, verbose=False):
    """Count units and branches grouped by DIV and Category."""
    #Old method:
    # unit_count = 0
    # branch_count = 0
    # for df in data:
    #     required_columns = ['DIV', 'Category', 'unit_ids', 'branch_id']
    #     if not all(column in df.columns for column in required_columns):
    #         raise ValueError(f"Data must contain the following columns: {required_columns}")
    
    #     unit_count += len(df['unit_ids'].index)
    #     branch_counts = [len(eval(branch_list)) for branch_list in df['branch_id']]
    #     branch_count += sum(branch_counts)

    #New method:
    tot_unit_count = 0
    tot_branch_count = 0
    for well in data:
        unit_count = len(well['unit_ids'])
        branch_count = [len(branch_list) for branch_list in well['branch_id']]
        branch_count = sum(branch_count)
        if verbose: print(f'DIV: {well["DIV"][0]}, Category: {well["category"][0]}, Units: {unit_count}, Branches: {branch_count}, Mean Branches per unit: {branch_count / unit_count}')

        # Update total counts
        tot_unit_count += unit_count
        tot_branch_count += branch_count    
    #print(f'Total units: {tot_unit_count}, Total branches: {tot_branch_count}')
    return tot_unit_count, tot_branch_count

def validate_template_densities(data, std=1, verbose=False, method='by_well'):
    """Validate and filter density data based on given criteria."""
    print('Validating templates...')
    print(f'Method: {method}')
    
    def get_density_threshold(density_values):
        density_values = pd.Series(density_values)
        
        # Calculate the mean and standard deviation
        mean_value = density_values.mean()
        std_dev = density_values.std()

        # Define the threshold as one standard deviation below the mean
        density_threshold = mean_value - std_dev*std
        return density_threshold
    
    def filter_data(data, density_values, threshold):
        """Filter data based on density threshold."""
        before_count = len(data["unit_ids"])
        valid_indices = [idx for idx, val in enumerate(density_values) if val >= threshold]
        valid_return = {key: [val[idx] for idx in valid_indices] for key, val in data.items()}
        after_count = len(valid_return["unit_ids"])
        if verbose: print(f'Before: {before_count} units, After: {after_count} units')
        return valid_return
    
    valid_data = data.copy()

    if method == 'by_well': #TODO: Update and implement if needed.
        raise NotImplementedError("Method 'by_well' is not implemented yet.")
        for i, well in enumerate(data):
            # Extract and evaluate 'channel_density'
            density_values = [d[0] for d in well['channel_density']]
            density_threshold = pd.Series(density_values).quantile(threshold)
            valid_data[i] = filter_data(well, density_values, density_threshold)
            if verbose: print(f'Well {i}, Density threshold: {density_threshold}')

    elif method == 'by_DIV':
        divs = set(well['DIV'][0] for well in data)
        for div in divs:
            div_data = [well for well in data if well['DIV'][0] == div]
            pre_unit_count, _ = count_units_and_branches(div_data, verbose=False)
            try: assert len(div_data) > 0, f"No data found for DIV {div}"
            except AssertionError as e: print(e); continue
            combined_density_values = [d[0] for well in div_data for d in well['channel_density']]
            combined_density_values.sort()
            density_threshold = get_density_threshold(combined_density_values)
            for i, well in enumerate(div_data):
                density_values = [d[0] for d in well['channel_density']]
                density_values.sort()
                valid_data[valid_data.index(well)] = filter_data(well, density_values, density_threshold)
            if verbose: 
                print(f'DIV {div}, Density threshold: {density_threshold}')
            div_data = [well for well in valid_data if well['DIV'][0] == div]
            post_unit_count, _ = count_units_and_branches(div_data, verbose=False)
            print(f'DIV {div}, Units before: {pre_unit_count}, Units after: {post_unit_count}')

    elif method == 'by_category': #TODO: Update and implement if needed.
        raise NotImplementedError("Method 'by_category' is not implemented yet.")
        categories = set(well['category'][0] for well in data)
        for category in categories:
            category_data = [well for well in data if well['category'][0] == category]
            try: assert len(category_data) > 0, f"No data found for category {category}"
            except AssertionError as e: print(e); continue
            combined_density_values = [d[0] for well in category_data for d in well['channel_density']]
            density_threshold = pd.Series(combined_density_values).quantile(threshold)
            for i, well in enumerate(category_data):
                density_values = [d[0] for d in well['channel_density']]
                valid_data[valid_data.index(well)] = filter_data(well, density_values, density_threshold)
            if verbose:
                print(f'Category {category}, Density threshold: {density_threshold}')

    elif method == 'global':
        combined_density_values = [d[0] for well in data for d in well['channel_density']]
        density_threshold = get_density_threshold(combined_density_values)
        for i, well in enumerate(data):
            density_values = [d[0] for d in well['channel_density']]
            valid_data[i] = filter_data(well, density_values, density_threshold)
        if verbose:
            print(f'Global density threshold: {density_threshold}')

    valid_units, valid_branches = count_units_and_branches(valid_data)
    #print('After density validation:')
    print(f'Total units: {valid_units}, Total branches: {valid_branches}')
    
    return valid_data

def validate_branch_lengths_data(data, min_length=0, stdevs=3, outlier_method='by_DIV', remove_outlier_data=False):
    """Validate and filter branch lengths data based on given criteria."""
    #Old method:
    # data_copy = data.copy()
    # total_units, total_branches = len(set(data_copy['unit_ids'])), len(data_copy['branch_id'])
    # print('Before branch length validation:')
    # print(f'Total units: {total_units}, Total branches: {total_branches}')
    
    # data_copy['length'] = data_copy['length'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    # valid_data = data_copy[data_copy['length'].apply(lambda x: x > 0)]
    # valid_data = valid_data[np.abs(valid_data['length'] - valid_data['length'].mean()) <= (stdevs * valid_data['length'].std())]

    # valid_units, valid_branches = len(set(valid_data['unit_ids'])), len(valid_data['branch_id'])
    # print('After branch length validation:')
    # print(f'Total units: {valid_units}, Total branches: {valid_branches}')

    #New method:
    print('Validating branch lengths...')
    def remove_short_branches(data, min_length):
        valid_data = data.copy()
        for i, well in enumerate(data):
            for unit_idx, branch_lengths in enumerate(well['length']):
                valid_lengths = []
                valid_branch_ids = []
                valid_velocities = []
                for j, branch_length in enumerate(branch_lengths):
                    if branch_length > min_length:
                        valid_lengths.append(branch_length)
                        valid_branch_ids.append(well['branch_id'][unit_idx][j])
                        valid_velocities.append(well['velocity'][unit_idx][j])
                valid_data[i]['length'][unit_idx] = valid_lengths
                valid_data[i]['branch_id'][unit_idx] = valid_branch_ids
                valid_data[i]['velocity'][unit_idx] = valid_velocities
        return valid_data

    def remove_outliers(data, stdevs, method='by_DIV'):
        def filter_data(well, branch_lengths, mean_length, std_length):
            valid_lengths = []
            valid_branch_ids = []
            valid_velocities = []
            for j, length in enumerate(branch_lengths):
                if np.abs(length - mean_length) <= stdevs * std_length:
                    valid_lengths.append(length)
                    valid_branch_ids.append(well['branch_id'][unit_idx][j])
                    valid_velocities.append(well['velocity'][unit_idx][j])
            return valid_lengths, valid_branch_ids, valid_velocities

        valid_data = data.copy()

        if method == 'by_DIV':
            divs = set(well['DIV'][0] for well in data)
            for div in divs:
                div_data = [well for well in data if well['DIV'][0] == div]
                _, pre_branch_count = count_units_and_branches(div_data)
                #print(f'DIV {div}: Before outlier removal: {pre_branch_count} branches')
                for well in div_data:
                    for unit_idx, branch_lengths in enumerate(well['length']):
                        mean_length = np.mean(branch_lengths)
                        std_length = np.std(branch_lengths)
                        valid_lengths, valid_branch_ids, valid_velocities = filter_data(well, branch_lengths, mean_length, std_length)
                        valid_data[valid_data.index(well)]['length'][unit_idx] = valid_lengths
                        valid_data[valid_data.index(well)]['branch_id'][unit_idx] = valid_branch_ids
                        valid_data[valid_data.index(well)]['velocity'][unit_idx] = valid_velocities
                div_data = [well for well in valid_data if well['DIV'][0] == div]
                _, post_branch_count = count_units_and_branches(div_data)
                #print(f'DIV {div}: After outlier removal: {post_branch_count} branches')
                print(f'DIV {div}: Before: {pre_branch_count} branches, After: {post_branch_count} branches')

        elif method == 'by_well':
            for i, well in enumerate(data):
                for unit_idx, branch_lengths in enumerate(well['length']):
                    mean_length = np.mean(branch_lengths)
                    std_length = np.std(branch_lengths)
                    valid_lengths, valid_branch_ids, valid_velocities = filter_data(well, branch_lengths, mean_length, std_length)
                    valid_data[i]['length'][unit_idx] = valid_lengths
                    valid_data[i]['branch_id'][unit_idx] = valid_branch_ids
                    valid_data[i]['velocity'][unit_idx] = valid_velocities

        return valid_data
    
    valid_data = data.copy()
    _, pre_branch_count = count_units_and_branches(valid_data)
    valid_data = remove_short_branches(valid_data, min_length)
    if remove_outlier_data: valid_data = remove_outliers(valid_data, stdevs, method=outlier_method)
    _, post_branch_count = count_units_and_branches(valid_data)
    print(f'Before: {pre_branch_count} branches, After: {post_branch_count} branches')   
    return valid_data

def validate_velocity_data(data, stdevs=3, outlier_method='by_DIV', remove_outlier_data=False):
    """Validate and filter velocity data based on given criteria."""
    print('Validating velocity data...')
    
    def remove_negative_velocities(data):
        valid_data = data.copy()
        for i, well in enumerate(data):
            for unit_idx, velocities in enumerate(well['velocity']):
                valid_velocities = []
                valid_branch_ids = []
                valid_lengths = []
                for j, velocity in enumerate(velocities):
                    if velocity >= 0:
                        valid_velocities.append(velocity)
                        valid_branch_ids.append(well['branch_id'][unit_idx][j])
                        valid_lengths.append(well['length'][unit_idx][j])
                valid_data[i]['velocity'][unit_idx] = valid_velocities
                valid_data[i]['branch_id'][unit_idx] = valid_branch_ids
                valid_data[i]['length'][unit_idx] = valid_lengths
        return valid_data

    def remove_outliers(data, stdevs, method='by_DIV'):
        def filter_data(well, velocities, mean_velocity, std_velocity):
            valid_velocities = []
            valid_branch_ids = []
            valid_lengths = []
            for j, velocity in enumerate(velocities):
                if np.abs(velocity - mean_velocity) <= stdevs * std_velocity:
                    valid_velocities.append(velocity)
                    valid_branch_ids.append(well['branch_id'][unit_idx][j])
                    valid_lengths.append(well['length'][unit_idx][j])
            return valid_velocities, valid_branch_ids, valid_lengths

        valid_data = data.copy()

        if method == 'by_DIV':
            divs = set(well['DIV'][0] for well in data)
            for div in divs:
                div_data = [well for well in data if well['DIV'][0] == div]
                _, pre_branch_count = count_units_and_branches(div_data)
                for well in div_data:
                    for unit_idx, velocities in enumerate(well['velocity']):
                        mean_velocity = np.mean(velocities)
                        std_velocity = np.std(velocities)
                        valid_velocities, valid_branch_ids, valid_lengths = filter_data(well, velocities, mean_velocity, std_velocity)
                        valid_data[valid_data.index(well)]['velocity'][unit_idx] = valid_velocities
                        valid_data[valid_data.index(well)]['branch_id'][unit_idx] = valid_branch_ids
                        valid_data[valid_data.index(well)]['length'][unit_idx] = valid_lengths
                div_data = [well for well in valid_data if well['DIV'][0] == div]
                _, post_branch_count = count_units_and_branches(div_data)
                print(f'DIV {div}: Before: {pre_branch_count} branches, After: {post_branch_count} branches')

        elif method == 'by_well':
            for i, well in enumerate(data):
                for unit_idx, velocities in enumerate(well['velocity']):
                    mean_velocity = np.mean(velocities)
                    std_velocity = np.std(velocities)
                    valid_velocities, valid_branch_ids, valid_lengths = filter_data(well, velocities, mean_velocity, std_velocity)
                    valid_data[i]['velocity'][unit_idx] = valid_velocities
                    valid_data[i]['branch_id'][unit_idx] = valid_branch_ids
                    valid_data[i]['length'][unit_idx] = valid_lengths

        return valid_data
    
    valid_data = data.copy()
    _, pre_branch_count = count_units_and_branches(valid_data)
    valid_data = remove_negative_velocities(valid_data)
    if remove_outlier_data:
        valid_data = remove_outliers(valid_data, stdevs, method=outlier_method)
    _, post_branch_count = count_units_and_branches(valid_data)
    print(f'Before: {pre_branch_count} branches, After: {post_branch_count} branches')   
    return valid_data

def validate_and_extract_branch_counts(data, verbose=False):
    """Validate and filter number of branches data based on given criteria."""
    print('Validating and extracting branch counts...')
    valid_data = data.copy()
    initial_unit_count, initial_branch_count = count_units_and_branches(valid_data)

    filtered_data = []
    for well in valid_data:
        unit_count = len(well['unit_ids'])
        branch_counts = [len(branch_list) for branch_list in well['branch_id']]
        
        if any(branch_counts):
            well['branch_counts'] = branch_counts
            filtered_data.append(well)

    final_unit_count, final_branch_count = count_units_and_branches(filtered_data)

    print(f'Before: {initial_unit_count} units, After: {final_unit_count} units')
    #print(f'Before: {initial_branch_count} branches, After: {final_branch_count} branches')
    
    return filtered_data

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

def test_significance(data, dep_var, group_by, print_stats=False, method='auto'):
    """
    Test significance between groups using t-tests or ANOVA.
    
    Parameters:
    - data: DataFrame containing the data.
    - dep_var: The dependent variable column name.
    - group_by: The column name to group by.
    - print_stats: Whether to print statistical results.
    - method: 'ttest' for t-tests, 'anova' for ANOVA, 'auto' to decide based on the number of categories.
    
    Returns:
    - pairs: List of compared pairs (for t-tests or Tukey's HSD after ANOVA).
    - p_values: List of p-values for each comparison.
    """
    divs = data['DIV'].unique()
    pairs, p_values = [], []

    for div in divs:
        data_div = data[data['DIV'] == div]
        categories = data_div[group_by].unique()
        
        if method == 'ttest' or (method == 'auto' and len(categories) == 2):
            pairs_div = [(cat1, cat2) for i, cat1 in enumerate(categories) for cat2 in categories[i+1:]]
            
            for (cat1, cat2) in pairs_div:
                data1 = data_div[data_div[group_by] == cat1][dep_var].dropna().astype(float)
                data2 = data_div[data_div[group_by] == cat2][dep_var].dropna().astype(float)
                t_stat, p_val = stats.ttest_ind(data1, data2)
                pairs.append((div, cat1, cat2))
                p_values.append(p_val)
                if print_stats: 
                    print(f'Comparing {cat1} vs {cat2} at DIV {div}: t-statistic={t_stat}, p-value={p_val}')
        
        elif method == 'anova' or (method == 'auto' and len(categories) > 2):
            data_groups = [data_div[data_div[group_by] == cat][dep_var].dropna().astype(float) for cat in categories]
            f_stat, p_val = stats.f_oneway(*data_groups)
            if print_stats: print(f'ANOVA for DIV {div}: F-statistic={f_stat}, p-value={p_val}')
            
            if p_val < 0.05:  # Perform post-hoc test if ANOVA is significant
                # Ensure the data is clean and properly typed
                cleaned_data = data_div.dropna(subset=[dep_var, group_by])
                endog = cleaned_data[dep_var].astype(float)
                groups = cleaned_data[group_by].astype(str)
                
                tukey_result = pairwise_tukeyhsd(endog=endog, groups=groups, alpha=0.05)
                for result in tukey_result.summary().data[1:]:
                    cat1, cat2, _, p_val, _, _, _ = result
                    pairs.append((div, cat1, cat2))
                    p_values.append(p_val)
                    if print_stats:
                        print(f'Comparing {cat1} vs {cat2} at DIV {div}: p-value={p_val}')
    
    return pairs, p_values

def identify_outliers(series):
    """Identify outliers using Tukey's method."""
    q1 = np.percentile(series, 25)
    q3 = np.percentile(series, 75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    return (series < lower_bound) | (series > upper_bound)

def add_significance_annotations(ax, data, dep_var, group_by, print_stats=False, bar_width=0.8/3, stagger=5):
    """Add significance annotations to the plot."""
    print("Disclaimer: Significance levels are indicated as follows: * for p <= 0.05, ** for p <= 0.01, *** for p <= 0.001, **** for p <= 0.0001.")
    
    pairs, p_values = test_significance(data, dep_var, group_by, print_stats=print_stats)
    divs = data['DIV'].unique()
    annotations_per_div = {div: 0 for div in divs}

    # Get the y range of the plot
    y_max = ax.get_ylim()[1]

    for (div, cat1, cat2), p_val in zip(pairs, p_values):
        data_div = data[data['DIV'] == div]
        data1 = data_div[data_div[group_by] == cat1][dep_var].dropna().astype(float)
        data2 = data_div[data_div[group_by] == cat2][dep_var].dropna().astype(float)
        if p_val <= 0.05:
            stars = ""
            if p_val <= 0.0001:
                stars = "****"
            elif p_val <= 0.001:
                stars = "***"
            elif p_val <= 0.01:
                stars = "**"
            elif p_val <= 0.05:
                stars = "*"
            
            if stars:
                #divs = divs.sort()
                sorted_divs = sorted(divs)
                x_ind = list(sorted_divs).index(div)
                categories = data[data['DIV'] == div][group_by].unique()
                cat1_index = list(categories).index(cat1)
                cat2_index = list(categories).index(cat2)

                div_pos = x_ind
                div_width = len(categories) * bar_width
                div_start = div_pos - div_width / 2
                bar1_pos = div_start + (cat1_index) * bar_width + bar_width / 2
                bar2_pos = div_start + (cat2_index) * bar_width + bar_width / 2
                x1 = bar1_pos
                x2 = bar2_pos
                
                y_max_sem1 = stats.sem(data1)   #get the y_min values for vertical lines
                y_max_sem2 = stats.sem(data2)
                assert y_max_sem1 > 0, f"y_max_sem1={y_max_sem1}"
                assert y_max_sem2 > 0, f"y_max_sem2={y_max_sem2}"
                
                y_mean1 = data1.mean()
                y_mean2 = data2.mean()
                y_min1 = y_mean1 + y_max_sem1/2 + y_max*0.05
                y_min2 = y_mean2 + y_max_sem2/2 + y_max*0.05
                
                stagger = y_max * 0.05 #add stagger to y_min values for vertical coords
                stagger_offset = annotations_per_div[div] * stagger
                y = y_max + 1 + stagger_offset
                h, col = 1, 'k'
                assert y_min1<y+h, f"y_min1={y_min1}, y+h={y+h}"
                assert y_min2<y+h, f"y_min2={y_min2}, y+h={y+h}"
                
                ax.plot([x1, x2], [y + h, y + h], lw=1.5, c=col) #horizontal line
                ax.plot([x1, x1], [y_min1, y + h], lw=1.5, ls='dotted', c=col) #vertical line
                ax.plot([x2, x2], [y_min2, y + h], lw=1.5, ls='dotted', c=col) #vertical line
                ax.text((x1 + x2) / 2, y + h, stars, ha='center', va='bottom', color=col) #stars

                annotations_per_div[div] += 1

                if print_stats: print(f'Significant difference between {cat1} and {cat2} at DIV {div}: p-value={p_val}, stars={stars}')

                #break

def plot_with_significance(data, y, title, ylabel, scatter=True, cat_order=None, print_stats=False, div_width=.8, exclude_outliers=False):
    """Plot data with significance annotations and outliers as diamonds."""
    df = data.copy()
    
    plt.figure(figsize=(12, 8))
    cats = df['category'].unique()
    cat_list = []
    if cat_order is not None: #cats = [cat for cat in cat_order if cat in cats]
        for cat in cat_order:
            if cat not in cats: logger.warning(f"Category '{cat}' not found in data.")
            else: cat_list.append(cat)
        cats = cat_list
    df = df.sort_values(by='category', key=lambda x: x.map({cat: i for i, cat in enumerate(cats)}))
    
    sorted_data = df[df['category'].isin(cats)].copy()
    sorted_data['category'] = pd.Categorical(sorted_data['category'], categories=cats, ordered=True)    
    palette = sns.color_palette("husl", len(cats))
    
    df['outlier'] = df.groupby(['DIV', 'category'])[y].transform(identify_outliers)

    if exclude_outliers: 
        sorted_data = sorted_data[~df['outlier']]
        title += ' (Excluding Outliers)'
    
    if scatter:
        ax = sns.stripplot(x='DIV', y=y, hue='category', data=sorted_data[~df['outlier']], jitter=True, dodge=True, marker='o', alpha=0.6, palette=palette, edgecolor='gray')
        # Highlight outliers
        ax = sns.stripplot(x='DIV', y=y, hue='category', data=sorted_data[df['outlier']], jitter=True, dodge=True, marker='D', palette=palette, edgecolor='red', size=6, ax=ax, legend=False)
        ax = sns.barplot(x='DIV', y=y, hue='category', data=sorted_data, palette=palette, alpha=0.6, errorbar='se', 
                         capsize=0.1, 
                         errwidth=1, 
                         dodge=True, ax=ax, width=div_width)
    else:
        ax = sns.barplot(x='DIV', y=y, hue='category', data=sorted_data, palette=palette, 
                        capsize=0.1,
                        errwidth=1, 
                        alpha=0.6, errorbar='se', width=div_width)
    
    bar_width = div_width / len(cats)
    add_significance_annotations(ax, sorted_data, y, 'category', print_stats=print_stats, bar_width=bar_width)
    
    # bar_width = div_width / len(cats)
    # add_significance_annotations(ax, sorted_data, y, 'category', print_stats=print_stats, bar_width=bar_width)
    
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel('DIV')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[0:3], labels[0:3], title='Category')
    plt.show()

    return ax

# Example usage
# data = pd.read_csv('your_data.csv')
# plot_with_significance(data, 'dependent_variable', 'Plot Title', 'Y-axis Label', scatter=True, cat_order=['WT', 'KO', 'TG'])