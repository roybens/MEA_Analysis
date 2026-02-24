import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import pearsonr, zscore
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# Load the Excel files
quality_metrics = pd.read_excel("quality_metrics.xlsx", index_col=0)
template_metrics = pd.read_excel("template_metrics.xlsx", index_col=0)

# Join the data on the index (unit numbers)
merged_data = quality_metrics.join(template_metrics, how="inner")

# Open a PDF for saving all figures
with PdfPages("complete_analysis.pdf") as pdf:
    
    # 1. Amplitude Median Statistical Inference
    amplitudes = merged_data['amplitude_median'].dropna()
    mean_amp, median_amp = amplitudes.mean(), amplitudes.median()
    std_amp, min_amp, max_amp = amplitudes.std(), amplitudes.min(), amplitudes.max()
    z_scores = np.abs(zscore(amplitudes))
    outliers = amplitudes[z_scores > 3]

    # Amplitude statistics and outliers as text figure
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis("off")
    stats_text = (
        f"Amplitude Median Summary Statistics:\n"
        f"Mean: {mean_amp:.2f}\nMedian: {median_amp:.2f}\n"
        f"Std Dev: {std_amp:.2f}\nMin: {min_amp:.2f}\nMax: {max_amp:.2f}\n\n"
        f"Outliers Detected (Z-score > 3):"
    )
    for idx, outlier in outliers.items():
        stats_text += f"\nIndex {idx}: Amplitude {outlier:.2f}"
    ax.text(0.1, 0.5, stats_text, fontsize=12, va="top")
    pdf.savefig()
    plt.close()
    
    # Amplitude Median Histogram
    plt.figure(figsize=(8, 6))
    plt.hist(amplitudes, bins=30, color="blue", alpha=0.7)
    plt.title("Histogram of Amplitude Median")
    plt.xlabel("Amplitude")
    plt.ylabel("Frequency")
    pdf.savefig()
    plt.close()
    
    # Amplitude Median Boxplot
    plt.figure(figsize=(8, 6))
    plt.boxplot(amplitudes, vert=False)
    plt.title("Boxplot of Amplitude Median")
    plt.xlabel("Amplitude")
    pdf.savefig()
    plt.close()

    # 2. Firing Rate Statistical Inference
    firing_rates = merged_data['firing_rate'].dropna()
    mean_firing, median_firing = firing_rates.mean(), firing_rates.median()
    std_firing, min_firing, max_firing = firing_rates.std(), firing_rates.min(), firing_rates.max()
    z_scores_firing = np.abs(zscore(firing_rates))
    outliers_firing = firing_rates[z_scores_firing > 3]

    # Firing Rate statistics and outliers as text figure
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis("off")
    stats_text = (
        f"Firing Rate Summary Statistics:\n"
        f"Mean: {mean_firing:.2f}\nMedian: {median_firing:.2f}\n"
        f"Std Dev: {std_firing:.2f}\nMin: {min_firing:.2f}\nMax: {max_firing:.2f}\n\n"
        f"Outliers Detected (Z-score > 3):"
    )
    for idx, outlier in outliers_firing.items():
        stats_text += f"\nIndex {idx}: Firing Rate {outlier:.2f}"
    ax.text(0.1, 0.5, stats_text, fontsize=12, va="top")
    pdf.savefig()
    plt.close()

    # Firing Rate Histogram
    plt.figure(figsize=(8, 6))
    plt.hist(firing_rates, bins=30, color="green", alpha=0.7)
    plt.title("Histogram of Firing Rate")
    plt.xlabel("Firing Rate")
    plt.ylabel("Frequency")
    pdf.savefig()
    plt.close()
    
    # Firing Rate Boxplot
    plt.figure(figsize=(8, 6))
    plt.boxplot(firing_rates, vert=False)
    plt.title("Boxplot of Firing Rate")
    plt.xlabel("Firing Rate")
    pdf.savefig()
    plt.close()

    # 3. Relationship Between Firing Rate and Amplitude Median
    stat_data = merged_data[['firing_rate', 'amplitude_median']].dropna()
    correlation, p_value = pearsonr(stat_data['firing_rate'], stat_data['amplitude_median'])

    # Scatter Plot with Regression Line
    plt.figure(figsize=(8, 6))
    sns.regplot(x='firing_rate', y='amplitude_median', data=stat_data, scatter_kws={'alpha':0.5})
    plt.title(f"Scatter Plot of Firing Rate vs. Amplitude Median\nCorrelation: {correlation:.2f}, p-value: {p_value:.2e}")
    plt.xlabel("Firing Rate")
    plt.ylabel("Amplitude Median")
    pdf.savefig()
    plt.close()

    # Correlation statistics as text figure
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis("off")
    stats_text = (
        f"Relationship Between Firing Rate and Amplitude Median:\n"
        f"Correlation Coefficient: {correlation:.2f}\n"
        f"P-value: {p_value:.2e}\n\n"
    )
    if p_value < 0.05:
        stats_text += "There is a statistically significant association between firing rate and amplitude."
    else:
        stats_text += "There is no statistically significant association between firing rate and amplitude."
    ax.text(0.1, 0.5, stats_text, fontsize=12, va="top")
    pdf.savefig()
    plt.close()

    # 4. Clustering and Classification
    features_with_amplitude = merged_data[['firing_rate', 'half_width', 'num_negative_peaks', 'num_positive_peaks', 
                                           'peak_to_valley', 'peak_trough_ratio', 'recovery_slope', 
                                           'repolarization_slope', 'num_spikes', 'amplitude_median']].dropna()

    features_without_amplitude = merged_data[['firing_rate', 'half_width', 'num_negative_peaks', 'num_positive_peaks', 
                                              'peak_to_valley', 'peak_trough_ratio', 'recovery_slope', 
                                              'repolarization_slope', 'num_spikes']].dropna()

    # Standardize features
    scaler = StandardScaler()
    X_scaled_with_amplitude = scaler.fit_transform(features_with_amplitude)
    X_scaled_without_amplitude = scaler.fit_transform(features_without_amplitude)

    # PCA for visualization
    pca = PCA(n_components=2)

    # K-Means Clustering with 'amplitude_median'
    X_pca_with_amplitude = pca.fit_transform(X_scaled_with_amplitude)
    kmeans_with_amplitude = KMeans(n_clusters=2, random_state=42)
    kmeans_labels_with_amplitude = kmeans_with_amplitude.fit_predict(X_scaled_with_amplitude)

    # Plot K-Means clustering results with 'amplitude_median'
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=X_pca_with_amplitude[:, 0], y=X_pca_with_amplitude[:, 1], hue=kmeans_labels_with_amplitude, palette='viridis')
    plt.title("K-Means Clustering (With 'amplitude_median')")
    pdf.savefig()
    plt.close()

    # K-Means Clustering without 'amplitude_median'
    X_pca_without_amplitude = pca.fit_transform(X_scaled_without_amplitude)
    kmeans_without_amplitude = KMeans(n_clusters=2, random_state=42)
    kmeans_labels_without_amplitude = kmeans_without_amplitude.fit_predict(X_scaled_without_amplitude)

    # Plot K-Means clustering results without 'amplitude_median'
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=X_pca_without_amplitude[:, 0], y=X_pca_without_amplitude[:, 1], hue=kmeans_labels_without_amplitude, palette='viridis')
    plt.title("K-Means Clustering (Without 'amplitude_median')")
    pdf.savefig()
    plt.close()

    # Random Forest Classifier and Feature Importance
    rf_with_amplitude = RandomForestClassifier(random_state=42)
    rf_with_amplitude.fit(X_scaled_with_amplitude, kmeans_labels_with_amplitude)
    feature_importance_with = pd.DataFrame({'Feature': features_with_amplitude.columns, 'Importance': rf_with_amplitude.feature_importances_})
    feature_importance_with = feature_importance_with.sort_values(by='Importance', ascending=False)

    # Plot feature importance with 'amplitude_median'
    plt.figure(figsize=(8, 6))
    sns.barplot(x='Importance', y='Feature', data=feature_importance_with, palette="viridis")
    plt.title("Feature Importance from Random Forest (With 'amplitude_median')")
    pdf.savefig()
    plt.close()

    rf_without_amplitude = RandomForestClassifier(random_state=42)
    rf_without_amplitude.fit(X_scaled_without_amplitude, kmeans_labels_without_amplitude)
    feature_importance_without = pd.DataFrame({'Feature': features_without_amplitude.columns, 'Importance': rf_without_amplitude.feature_importances_})
    feature_importance_without = feature_importance_without.sort_values(by='Importance', ascending=False)

    # Plot feature importance without 'amplitude_median'
    plt.figure(figsize=(8, 6))
    sns.barplot(x='Importance', y='Feature', data=feature_importance_without, palette="viridis")
    plt.title("Feature Importance from Random Forest (Without 'amplitude_median')")
    pdf.savefig()
    plt.close()