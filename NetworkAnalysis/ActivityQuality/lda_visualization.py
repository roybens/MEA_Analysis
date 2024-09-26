# # Apply LDA to multiple factors according to columns
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# from sklearn.preprocessing import StandardScaler, LabelEncoder
# from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
# from scipy import stats
# from matplotlib.patches import Ellipse


# def preprocess_data(data, feature_names):
#     # Check for and handle infinite values
#     data.replace([np.inf, -np.inf], np.nan, inplace=True)
#     data.dropna(subset=feature_names, inplace=True)  # Ensure no NaNs in features

#     # Ensure all data is numeric and finite
#     data[feature_names] = data[feature_names].apply(pd.to_numeric, errors='coerce')
#     data.dropna(subset=feature_names, inplace=True)

#     # Standardize the features
#     scaler = StandardScaler()
#     X_scaled = scaler.fit_transform(data[feature_names])
#     return X_scaled

# def safe_covariance(sub_data):
#     if sub_data.shape[0] < 2 or np.var(sub_data) == 0:
#         return np.eye(sub_data.shape[1])  # Return identity if not enough data or no variance
#     covariance = np.cov(sub_data.T)
#     # Add a small jitter in the diagonal elements if necessary
#     epsilon = 1e-5
#     np.fill_diagonal(covariance, np.diag(covariance) + epsilon)
#     return covariance

# def plot_biplot(data, feature_names, group_column):
#     scaler = StandardScaler()
#     X_scaled = scaler.fit_transform(data[feature_names])
#     lda = LDA(n_components=2)
#     X_lda = lda.fit_transform(X_scaled, data[group_column])
#     lda_df = pd.DataFrame(X_lda, columns=['LD1', 'LD2'])
#     lda_df[group_column] = data[group_column]
    
#     # Decrease scaling factor to make ellipses smaller
#     scaling_factor = np.sqrt(5.991)  # Adjust this factor. 
#     # 5.991 is the critical value for 95% confidence level from the chi-square distribution (χ²) with 2 degrees of freedom.
#     # 2.991 is the critical value for 90% confidence level from the chi-square distribution (χ²) with 2 degrees of freedom.
    
#     explained_variance_ratio = lda.explained_variance_ratio_ * 100
#     plt.figure(figsize=(12, 9))  # You can adjust the figure size to better fit your needs
#     palette = sns.color_palette("Set1", n_colors=len(lda_df[group_column].unique()))
#     color_map = dict(zip(lda_df[group_column].unique(), palette))
    
#     sns.scatterplot(x='LD1', y='LD2', hue=group_column, data=lda_df, palette=color_map, s=50, alpha=0.7)
    
#     # Collect all ellipse data to determine plot limits
#     ellipse_extents = []
#     for group in lda_df[group_column].unique():
#         sub_data = lda_df[lda_df[group_column] == group]
#         covariance = np.cov(sub_data[['LD1', 'LD2']].T)
#         lambda_, v = np.linalg.eig(covariance)
#         lambda_ = np.sqrt(lambda_) * scaling_factor
#         ell = Ellipse(xy=(np.mean(sub_data['LD1']), np.mean(sub_data['LD2'])),
#                       width=lambda_[0]*2, height=lambda_[1]*2,
#                       angle=np.rad2deg(np.arccos(v[0, 0])), edgecolor=color_map[group], fc='None', lw=2, linestyle='--')
#         plt.gca().add_patch(ell)
#         ellipse_extents.append((np.mean(sub_data['LD1']) - lambda_[0], np.mean(sub_data['LD1']) + lambda_[0]))
#         ellipse_extents.append((np.mean(sub_data['LD2']) - lambda_[1], np.mean(sub_data['LD2']) + lambda_[1]))

#     # Determine new axis limits
#     ellipse_extents = np.array(ellipse_extents)
#     xlim = [ellipse_extents[:,0].min(), ellipse_extents[:,1].max()]
#     ylim = [ellipse_extents[:,0].min(), ellipse_extents[:,1].max()]
#     # Set axis limits based on the actual range of LDA components with some padding
#     plt.xlim([lda_df['LD1'].min() - 2, lda_df['LD1'].max() + 2])
#     plt.ylim([lda_df['LD2'].min() - 2, lda_df['LD2'].max() + 2])


#     plt.xlabel(f'LD1 ({explained_variance_ratio[0]:.2f}%)')
#     plt.ylabel(f'LD2 ({explained_variance_ratio[1]:.2f}%)')
#     plt.title(f'LDA Analysis of {group_column}')
#     plt.legend(title=group_column)
#     plt.grid(True)
#     plt.show()

# def plot_density(data, feature_names, group_column):
#     X = data[feature_names]
#     le = LabelEncoder()
#     y = le.fit_transform(data[group_column])
#     y_labels = le.inverse_transform(y)
#     lda = LDA()
#     X_lda = lda.fit_transform(X, y)
#     lda_df = pd.DataFrame(X_lda, columns=['LD1'])
#     lda_df[group_column] = y_labels
#     fig, ax = plt.subplots(figsize=(10, 6))
#     unique_labels = np.unique(y_labels)
#     palette = sns.color_palette("Set1", len(unique_labels))
#     color_map = dict(zip(unique_labels, palette))
#     medians = {}
#     for label in unique_labels:
#         class_data = lda_df[lda_df[group_column] == label]['LD1']
#         kde = stats.gaussian_kde(class_data)
#         x_range = np.linspace(lda_df['LD1'].min(), lda_df['LD1'].max(), 100)
#         y_values = kde(x_range)
#         y_normalized = y_values / y_values.max() * 0.4
#         ax.fill_between(x_range, y_normalized, alpha=0.3, color=color_map[label], label=f'{label} Density')
#         medians[label] = np.median(class_data)
#     for label in unique_labels:
#         scatter_data = lda_df[lda_df[group_column] == label]
#         ax.scatter(scatter_data['LD1'], np.zeros(len(scatter_data)), color=color_map[label], s=50, edgecolor='k', label=label, alpha=0.5)
#     for label in unique_labels:
#         median = medians[label]
#         ax.axvline(median, color=color_map[label], linestyle='--', linewidth=1)
#         ax.text(median, 0.45, f'{label}\nMedian: {median:.2f}', color=color_map[label], ha='center', va='top', fontsize=8, rotation=90)
#     ax.set_title(f'LDA Analysis of {group_column}')
#     ax.set_xlabel('Linear Discriminant 1')
#     ax.set_ylabel('Density')
#     ax.set_ylim([-0.05, 0.5])
#     handles, labels = ax.get_legend_handles_labels()
#     by_label = dict(zip(labels, handles))
#     ax.legend(by_label.values(), by_label.keys(), loc='upper right', title=group_column)
#     plt.show()

# def lda_plot(file_path, feature_names, group_column):
#     try:
#         data = pd.read_csv(file_path)
#         data = preprocess_data(data, feature_names)
#         if len(data[group_column].unique()) == 2:
#             plot_density(data, feature_names, group_column)
#         elif len(data[group_column].unique()) > 2:
#             plot_biplot(data, feature_names, group_column)
#     except Exception as e:
#         print(f"An error occurred: {e}")

# Apply LDA to multiple factors according to columns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from scipy import stats
from matplotlib.patches import Ellipse

def plot_biplot(data, feature_names, group_column):
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(data[feature_names])
    lda = LDA(n_components=2)
    X_lda = lda.fit_transform(X_scaled, data[group_column])
    lda_df = pd.DataFrame(X_lda, columns=['LD1', 'LD2'])
    lda_df[group_column] = data[group_column]
    
    # Decrease scaling factor to make ellipses smaller
    scaling_factor = np.sqrt(5.991)  # Adjust this factor. 
    # 5.991 is the critical value for 95% confidence level from the chi-square distribution (χ²) with 2 degrees of freedom.
    # 2.991 is the critical value for 90% confidence level from the chi-square distribution (χ²) with 2 degrees of freedom.
    
    explained_variance_ratio = lda.explained_variance_ratio_ * 100
    plt.figure(figsize=(12, 9))  # You can adjust the figure size to better fit your needs
    palette = sns.color_palette("Set1", n_colors=len(lda_df[group_column].unique()))
    color_map = dict(zip(lda_df[group_column].unique(), palette))
    
    sns.scatterplot(x='LD1', y='LD2', hue=group_column, data=lda_df, palette=color_map, s=50, alpha=0.7)
    
    # Collect all ellipse data to determine plot limits
    ellipse_extents = []
    for group in lda_df[group_column].unique():
        sub_data = lda_df[lda_df[group_column] == group]
        covariance = np.cov(sub_data[['LD1', 'LD2']].T)
        lambda_, v = np.linalg.eig(covariance)
        lambda_ = np.sqrt(lambda_) * scaling_factor
        ell = Ellipse(xy=(np.mean(sub_data['LD1']), np.mean(sub_data['LD2'])),
                      width=lambda_[0]*2, height=lambda_[1]*2,
                      angle=np.rad2deg(np.arccos(v[0, 0])), edgecolor=color_map[group], fc='None', lw=2, linestyle='--')
        plt.gca().add_patch(ell)
        ellipse_extents.append((np.mean(sub_data['LD1']) - lambda_[0], np.mean(sub_data['LD1']) + lambda_[0]))
        ellipse_extents.append((np.mean(sub_data['LD2']) - lambda_[1], np.mean(sub_data['LD2']) + lambda_[1]))

    # Determine new axis limits
    ellipse_extents = np.array(ellipse_extents)
    xlim = [ellipse_extents[:,0].min(), ellipse_extents[:,1].max()]
    ylim = [ellipse_extents[:,0].min(), ellipse_extents[:,1].max()]
    # Set axis limits based on the actual range of LDA components with some padding
    plt.xlim([lda_df['LD1'].min() - 2, lda_df['LD1'].max() + 2])
    plt.ylim([lda_df['LD2'].min() - 2, lda_df['LD2'].max() + 2])


    plt.xlabel(f'LD1 ({explained_variance_ratio[0]:.2f}%)')
    plt.ylabel(f'LD2 ({explained_variance_ratio[1]:.2f}%)')
    plt.title(f'LDA Analysis of {group_column}')
    plt.legend(title=group_column)
    plt.grid(True)
    plt.show()

def plot_density(data, feature_names, group_column):
    X = data[feature_names]
    le = LabelEncoder()
    y = le.fit_transform(data[group_column])
    y_labels = le.inverse_transform(y)
    lda = LDA()
    X_lda = lda.fit_transform(X, y)
    lda_df = pd.DataFrame(X_lda, columns=['LD1'])
    lda_df[group_column] = y_labels
    fig, ax = plt.subplots(figsize=(10, 6))
    unique_labels = np.unique(y_labels)
    palette = sns.color_palette("Set1", len(unique_labels))
    color_map = dict(zip(unique_labels, palette))
    medians = {}
    for label in unique_labels:
        class_data = lda_df[lda_df[group_column] == label]['LD1']
        kde = stats.gaussian_kde(class_data)
        x_range = np.linspace(lda_df['LD1'].min(), lda_df['LD1'].max(), 100)
        y_values = kde(x_range)
        y_normalized = y_values / y_values.max() * 0.4
        ax.fill_between(x_range, y_normalized, alpha=0.3, color=color_map[label], label=f'{label} Density')
        medians[label] = np.median(class_data)
    for label in unique_labels:
        scatter_data = lda_df[lda_df[group_column] == label]
        ax.scatter(scatter_data['LD1'], np.zeros(len(scatter_data)), color=color_map[label], s=50, edgecolor='k', label=label, alpha=0.5)
    for label in unique_labels:
        median = medians[label]
        ax.axvline(median, color=color_map[label], linestyle='--', linewidth=1)
        ax.text(median, 0.45, f'{label}\nMedian: {median:.2f}', color=color_map[label], ha='center', va='top', fontsize=8, rotation=90)
    ax.set_title(f'LDA Analysis of {group_column}')
    ax.set_xlabel('Linear Discriminant 1')
    ax.set_ylabel('Density')
    ax.set_ylim([-0.05, 0.5])
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right', title=group_column)
    plt.show()

def lda_plot(file_path, feature_names, group_column):
    """
    Parameters:
    file_path (str): Path to the CSV file containing the dataset.
    feature_names (list): List of column names to use as features for LDA.
    group_column (str): Column name to use for grouping the data in the LDA analysis.

    Returns:
    Depending on the number of groups, displays either a biplot or a density plot.
    """

    data = pd.read_csv(file_path)
    data = data.replace(np.nan, 0.0)

    if len(data[group_column].unique()) == 2:
        plot_density(data, feature_names, group_column)
    elif len(data[group_column].unique()) > 2:
        plot_biplot(data, feature_names, group_column)