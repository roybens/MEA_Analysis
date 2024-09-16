import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multicomp as multi
import logging

# Logging setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
logger.addHandler(stream_handler)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - %(module)s.%(funcName)s')
stream_handler.setFormatter(formatter)

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

# Colors
color_map = {'Homo': '#FFA500', 'HET': '#FF0000', 'WT': '#0000FF'}

def test_oneway_ANOVA_and_annotate_plot(data):
    f_val, p_val = stats.f_oneway(data[data['Category'] == 'Homo']['Length'],
                                  data[data['Category'] == 'HET']['Length'],
                                  data[data['Category'] == 'WT']['Length'])
    print('One-way ANOVA result: F =', f_val, ', p =', p_val)

    if p_val < 0.05:
        print('Performing Tukey HSD test...')
        print('Null hypothesis: All group means are equal.')
        mc = multi.MultiComparison(data['Length'], data['Category'])
        result = mc.tukeyhsd()
        print(result)

        categories = ['Homo', 'HET', 'WT']
        category_coordinates = {category: i for i, category in enumerate(categories)}
        y_max = data['Length'].max()

        sig_height = None
        y_diff = None
        for i, reject in enumerate(result.reject):
            if reject:
                group1, group2 = mc.groupsunique[np.triu_indices(len(mc.groupsunique), 1)[0][i]], mc.groupsunique[np.triu_indices(len(mc.groupsunique), 1)[1][i]]
                coord1, coord2 = category_coordinates[group1], category_coordinates[group2]

                if y_diff is None: ymin, ymax = plt.ylim()
                y_diff = ymax - ymin
                if sig_height is None: sig_height = y_max + (i + 1) * 10
                else: sig_height = sig_height + y_diff * 0.10
                if sig_height >= ymax: plt.ylim(ymax * 1.10, ymin)
                plt.plot([coord1, coord2], [sig_height, sig_height], color='k')
                plt.text((coord1 + coord2) / 2, sig_height, '*', ha='center', va='bottom', color='k', fontsize=20)

def extract_lengths(dfs, min_branches=1, dense_percentile=0):
    lengths, categories, unit_ids = [], [], []
    for i, df in enumerate(dfs):
        category = ['Homo', 'HET', 'WT'][i % 3]

        densities = [eval(density)[0] for density in df['channel_density']]
        max_channel_density = max(densities)

        for j, length_list in enumerate(df['length']):
            if len(eval(length_list)) < min_branches: continue
            if eval(df['channel_density'][j])[0] < max_channel_density * dense_percentile: continue
            unit_ids.extend(eval(df['unit_ids'][j]) * len(eval(length_list)))
            lengths.extend(eval(length_list))
            categories.extend([category] * len(eval(length_list)))

    return lengths, categories, unit_ids

def plot_branch_lengths_by_gene(file_paths, min_branches=0, dense_percentile=0, ylim=(0, 1000)):
    lengths, categories, unit_ids = extract_lengths(file_paths, min_branches=min_branches, dense_percentile=dense_percentile)
    data = pd.DataFrame({'Length': lengths, 'Category': categories, 'unit_ids': unit_ids})

    plt.figure(figsize=(12, 8))
    sns.boxplot(x='Category', y='Length', data=data, showmeans=False, meanline=False, showfliers=False, palette=color_map)

    unique_unit_ids = data['unit_ids'].unique()
    palette = sns.color_palette('husl', len(unique_unit_ids))
    unit_id_color_map = dict(zip(unique_unit_ids, palette))

    sns.stripplot(x='Category', y='Length', data=data, jitter=True, alpha=0.5, palette=unit_id_color_map, hue='unit_ids', dodge=True, legend=False)

    plt.title('Branch Length Comparison by Well Category')
    plt.xlabel('Category')
    plt.ylabel('Branch Length')
    plt.ylim(ylim)
    plt.grid(False)
    test_oneway_ANOVA_and_annotate_plot(data)
    plt.show()

def plot_velocities_by_gene(file_paths, min_branches=0, dense_percentile=0, ylim=(0, 1000)):
    lengths, categories, unit_ids = extract_branch_velocities(file_paths, min_branches=min_branches, dense_percentile=dense_percentile)
    data = pd.DataFrame({'Length': lengths, 'Category': categories, 'unit_ids': unit_ids})

    plt.figure(figsize=(12, 8))
    sns.boxplot(x='Category', y='Length', data=data, showmeans=False, meanline=False, showfliers=False, palette=color_map)

    unique_unit_ids = data['unit_ids'].unique()
    palette = sns.color_palette('husl', len(unique_unit_ids))
    unit_id_color_map = dict(zip(unique_unit_ids, palette))

    sns.stripplot(x='Category', y='Length', data=data, jitter=True, alpha=0.5, palette=unit_id_color_map, hue='unit_ids', dodge=True, legend=False)

    plt.title('Axon Velocity Comparison by Well Category')
    plt.xlabel('Category')
    plt.ylabel('Axon Velocity')
    plt.ylim(ylim)
    plt.grid(False)
    test_oneway_ANOVA_and_annotate_plot(data)
    plt.show()
