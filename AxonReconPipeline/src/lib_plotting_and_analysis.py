import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

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
def extract_lengths(dfs):
    lengths = []
    categories = []
    for i, df in enumerate(dfs):
        # Determine category based on well index
        if i % 3 == 0:
            category = 'Homo'
        elif i % 3 == 1:
            category = 'HET'
        else:
            category = 'WT'
        
        # Extract lengths and append to the lists
        for length_list in df['length']:
            lengths.extend(eval(length_list))  # Convert string representation of list to actual list
            categories.extend([category] * len(eval(length_list)))
    return lengths, categories

def plot_branch_lengths_by_gene(dfs):
        # Extract lengths and categories
        lengths, categories = extract_lengths(dfs)

        # Create a DataFrame for plotting
        data = pd.DataFrame({
            'Length': lengths,
            'Category': categories
        })

        # Plotting
        plt.figure(figsize=(12, 8))
        sns.boxplot(x='Category', y='Length', data=data, showmeans=True, meanline=True, showfliers=False)
        sns.stripplot(x='Category', y='Length', data=data, jitter=True, color='black', alpha=0.5)

        # Adding SEM bars
        category_means = data.groupby('Category')['Length'].mean()
        category_sems = data.groupby('Category')['Length'].sem()
        for i, category in enumerate(category_means.index):
            plt.errorbar(i, category_means[category], yerr=category_sems[category], fmt='none', capsize=5, color='red')

        plt.title('Branch Length Comparison by Well Category')
        plt.xlabel('Category')
        plt.ylabel('Branch Length')
        plt.grid(False)
        plt.show()
