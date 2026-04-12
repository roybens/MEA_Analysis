import matplotlib
matplotlib.use('Agg')  # Set the backend to 'Agg' before importing pyplot
import matplotlib.pyplot as plt
import io
import base64
import numpy as np
import os
import sys
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
from scipy.stats import sem
from dash import html

# Import MEAPlotter from IPNAnalysis
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from IPNAnalysis.meaplotter import MEAPlotter


def _fig_to_dash_outputs(fig):
    """Convert a matplotlib figure to (html.Img element, svg_bytes, png_bytes)."""
    png_buf = io.BytesIO()
    fig.savefig(png_buf, format='png', bbox_inches='tight')
    png_buf.seek(0)
    img_data = base64.b64encode(png_buf.read()).decode('utf-8')
    img_element = html.Img(src='data:image/png;base64,{}'.format(img_data))
    png_bytes = png_buf.getvalue()

    svg_buf = io.BytesIO()
    fig.savefig(svg_buf, format='svg', bbox_inches='tight')
    svg_buf.seek(0)
    svg_bytes = svg_buf.getvalue()

    return img_element, svg_bytes, png_bytes


def plot_metrics_graphs(data, divs, metrics, ordered_genotypes, selected_colors):
    """
    Plot grouped bar charts (mean ± SEM + scatter) per DIV using MEAPlotter.

    Uses t-test for two groups and one-way ANOVA + Tukey HSD for three or more
    groups to annotate significance brackets.

    Parameters
    ----------
    data : pd.DataFrame
        Must contain 'DIV' and 'NeuronType' columns alongside the metric columns.
    divs : list
        DIV values to include on the x-axis.
    metrics : list[str]
        Column names to plot, one figure per metric.
    ordered_genotypes : list[str]
        NeuronType labels in the desired display order.
    selected_colors : dict
        Mapping of genotype label -> colour string or hex code.

    Returns
    -------
    images : list[html.Img]
    svg_bytes_list : list[bytes]
    png_bytes_list : list[bytes]
    """
    images = []
    svg_bytes_list = []
    png_bytes_list = []

    if not ordered_genotypes:
        return [html.Div('No data to display.')], [], []

    div_order = [int(d) for d in divs]
    plotter = MEAPlotter(
        group_order=ordered_genotypes,
        palette=selected_colors,
        div_order=div_order,
    )

    for metric in metrics:
        if metric not in data.columns:
            images.append(html.Div(f'Metric {metric} not found in data.'))
            svg_bytes_list.append(b'')
            png_bytes_list.append(b'')
            continue

        fig, ax = plt.subplots()
        try:
            plotter.plot_bars_by_div(
                df=data,
                div_col='DIV',
                group_col='NeuronType',
                y=metric,
                div_order=div_order,
                group_order=ordered_genotypes,
                palette=selected_colors,
                title=metric,
                ax=ax,
                annotate=True,
                clip_upper=False,
                show_upper_outliers=False,
            )
        except Exception as e:
            print(f"Error plotting metric {metric} across DIVs {div_order}: {type(e).__name__}: {e}")
            plt.close(fig)
            continue

        img_element, svg_bytes, png_bytes = _fig_to_dash_outputs(fig)
        images.append(img_element)
        svg_bytes_list.append(svg_bytes)
        png_bytes_list.append(png_bytes)
        plt.close(fig)

    return images, svg_bytes_list, png_bytes_list


# Backward-compatible aliases retained so any external callers continue to work.
def plot_activity_graphs(data, divs, metrics, ordered_genotypes, selected_colors):
    """Alias for `plot_metrics_graphs`. Kept for backward compatibility."""
    return plot_metrics_graphs(data, divs, metrics, ordered_genotypes, selected_colors)


def plot_bar_with_p_values(data, divs, metrics, ordered_genotypes, selected_colors):
    """Alias for `plot_metrics_graphs`. Kept for backward compatibility."""
    return plot_metrics_graphs(data, divs, metrics, ordered_genotypes, selected_colors)

def plot_isi_graph(data_df, selected_divs, metric, selected_colors):
    all_svg_bytes_list = []
    filtered_df = data_df[data_df['DIV'].isin(selected_divs)]
    neuron_types = filtered_df['NeuronType'].unique()
    max_edges = None
    for item in filtered_df.to_dict('records'):
        if max_edges is None or len(item['networkAPFreqEdges']) > len(max_edges):
            max_edges = item['networkAPFreqEdges']

    plt.figure(figsize=(11, 8))
    aggregate_data = {nt: [] for nt in neuron_types}

    for item in filtered_df.to_dict('records'):
        bins = item[metric]
        edges = item[metric.replace('Bins', 'Edges')]
        neuron_type = item['NeuronType']

        mask = bins != 0
        filtered_bins = bins[mask]
        try:
            filtered_edges = edges[:-1][mask]
        except:
            print("Error: mask and edges length don't match")
            print("Edges length:", len(edges[:-1]), "Mask length:", len(mask))
            continue

        normalized_bins = min_max_normalize(filtered_bins)
        smoothed_bins = gaussian_filter1d(normalized_bins, sigma=1)
        interp = interp1d(filtered_edges, smoothed_bins, kind='previous', bounds_error=False, fill_value=0)
        standardized_bins = interp(max_edges[:-1])

        aggregate_data[neuron_type].append(standardized_bins)

    for neuron_type, bins_list in aggregate_data.items():
        if bins_list:
            bins_array = np.array(bins_list)
            mean_bins = np.mean(bins_array, axis=0)
            stderr_bins = sem(bins_array, axis=0)

            color = selected_colors.get(neuron_type, 'black')
            plt.plot(max_edges[:-1], mean_bins, color=color, label=neuron_type)
            plt.fill_between(max_edges[:-1], mean_bins - stderr_bins, mean_bins + stderr_bins, color=color, alpha=0.3)

    plt.xscale('log')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Percentage Count')
    plt.title(f"{metric.replace('APFreqBins', ' AP Frequencies')}")
    plt.legend(loc='upper right')
    plt.xlim(0.1, 500)
    plt.tight_layout()

    img_bytes = io.BytesIO()
    plt.savefig(img_bytes, format='png')
    img_bytes.seek(0)
    img_data = base64.b64encode(img_bytes.read()).decode('utf-8')
    img_element = html.Img(src='data:image/png;base64,{}'.format(img_data))

    svg_bytes = io.BytesIO()
    plt.savefig(svg_bytes, format='svg')
    svg_bytes.seek(0)
    svg_data = svg_bytes.getvalue()
    all_svg_bytes_list.append((f"{metric}_plot.svg", svg_data))

    plt.close()

    return img_element, svg_bytes, img_bytes.getvalue(), all_svg_bytes_list

def min_max_normalize(data):
    min_val = np.min(data)
    max_val = np.max(data)
    return (data - min_val) / (max_val - min_val)

