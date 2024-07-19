import matplotlib
matplotlib.use('Agg')  # Set the backend to 'Agg' before importing pyplot
import matplotlib.pyplot as plt
import io
import base64
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
from scipy.stats import sem
from dash import html
from scipy.stats import ttest_ind, t

def plot_activity_graphs(data, div, output_types, ordered_genotypes, selected_colors):
    markers = ['o', 's', 'D', '^', 'v', 'x', '*', 'H', '8', 'p']
    images = []
    svg_bytes_list = []
    png_bytes_list = []
    for output_type in output_types:
        total_genotypes = len(ordered_genotypes)
        if total_genotypes == 0:
            return [html.Div('No data to display.')], [], []

        output_arrays = {genotype: [] for genotype in ordered_genotypes}
        chip_arrays = {genotype: [] for genotype in ordered_genotypes}
        well_arrays = {genotype: [] for genotype in ordered_genotypes}

        for i in div:
            for genotype in ordered_genotypes:
                temp_df = data.loc[(data['DIV'] == i) & (data['NeuronType'] == genotype)]
                if output_type not in temp_df.columns:
                    return [html.Div(f'Output type {output_type} not found in data.')], [], []
                output_arrays[genotype].append(np.array(temp_df[output_type]))
                chip_arrays[genotype].append(np.array(temp_df['Chip_ID']))
                well_arrays[genotype].append(np.array(temp_df['Well']))

        bar_width = 0.25
        gap_between_bars = 0
        gap_between_groups = bar_width

        total_bar_group_width = total_genotypes * bar_width + (total_genotypes - 1) * gap_between_bars
        total_plot_width = len(div) * (total_bar_group_width + gap_between_groups) + gap_between_groups
        x_genotype = {genotype: [] for genotype in ordered_genotypes}
        base_x_coordinate = np.arange(len(div)) * (total_bar_group_width + gap_between_groups) + gap_between_groups + bar_width
        offset = (total_genotypes * bar_width + (total_genotypes - 1) * gap_between_bars) / 2
        centered_x = base_x_coordinate - offset + bar_width / 2
        for i, genotype in enumerate(ordered_genotypes):
            x_genotype[genotype] = centered_x + i * (bar_width + gap_between_bars)

        fig, ax = plt.subplots()

        mean_data_all = {}
        yerr_data_all = {}
        n_data_all = {}

        for i, genotype in enumerate(ordered_genotypes):
            y_data = output_arrays[genotype]
            chipy_data = chip_arrays[genotype]
            welly_data = well_arrays[genotype]

            mean_data = [np.nanmean(yi) if len(yi) > 0 else np.nan for yi in y_data]
            yerr_data = [np.nanstd(yi, ddof=1) / np.sqrt(len(yi)) if len(yi) > 1 else np.nan for yi in y_data]
            n_data = [len(yi) for yi in y_data]

            mean_data_all[genotype] = mean_data
            yerr_data_all[genotype] = yerr_data
            n_data_all[genotype] = n_data

            if len(x_genotype[genotype]) != len(mean_data):
                print(f"Length mismatch for genotype {genotype}: {len(x_genotype[genotype])} vs {len(mean_data)}")

            ax.bar(x_genotype[genotype], mean_data, yerr=yerr_data, capsize=3, width=bar_width, color=selected_colors.get(genotype, 'blue'), edgecolor='black', ecolor='black', label=genotype, alpha=0.6)

            jitter_amount = 0.05
            for j in range(len(x_genotype[genotype])):
                combined_data = [str(chip) + str(well) for chip, well in zip(chipy_data[j], welly_data[j])]
                for k in range(len(y_data[j])):
                    ax.scatter(x_genotype[genotype][j] + np.random.uniform(-jitter_amount, jitter_amount, 1), y_data[j][k], s=10, color=selected_colors.get(genotype, 'blue'), marker=markers[i % len(markers)])

        valid_arrays = [array for genotype_arrays in output_arrays.values() for array in genotype_arrays if len(array) > 0]
        if valid_arrays:
            maxim = max(max(array[np.isfinite(array)]) for array in valid_arrays)
        else:
            maxim = 0
            print("No valid arrays found.")

        for i in range(len(base_x_coordinate)):
            count = 1
            for j, genotype1 in enumerate(ordered_genotypes):
                for k, genotype2 in enumerate(ordered_genotypes):
                    if j < k:
                        mean1, sem1, n1 = mean_data_all[genotype1][i], yerr_data_all[genotype1][i], n_data_all[genotype1][i]
                        mean2, sem2, n2 = mean_data_all[genotype2][i], yerr_data_all[genotype2][i], n_data_all[genotype2][i]
                        if np.isnan(mean1) or np.isnan(mean2):
                            continue  # Skip t-test if any mean is NaN
                        sed = np.sqrt(sem1 ** 2.0 + sem2 ** 2.0)
                        t_stat = (mean1 - mean2) / sed
                        degreef = n1 + n2 - 2
                        alpha = 0.05
                        cv = t.ppf(1.0 - alpha, degreef)
                        p_value = (1.0 - t.cdf(abs(t_stat), degreef)) * 2.0
                        sign = "***" if p_value <= 0.001 else "**" if p_value <= 0.01 else "*" if p_value <= 0.05 else "ns"
                        if sign != 'ns':
                            ax.plot([x_genotype[genotype1][i], x_genotype[genotype2][i]], [maxim + 0.05 * maxim * count] * 2, 'k', linewidth=1.3)
                            ax.text((x_genotype[genotype1][i] + x_genotype[genotype2][i]) / 2, maxim + 0.05 * maxim * count, sign, ha='center', va='bottom', fontsize=7)
                            ax.axvline(x_genotype[genotype1][i], color='black', linestyle=':', linewidth=0.5)
                            #ax.axvline(x_genotype[genotype2][i], color='black', linestyle=':', linewidth=0.5)
                        count += 1

        plt.title(f"{output_type}", fontsize=14)
        plt.xlabel('DIV', fontsize=12)
        plt.ylabel(f"{output_type}", fontsize=12)
        plt.xticks(base_x_coordinate, div, fontsize=10)
        plt.legend(title='NeuronType', loc='upper left', bbox_to_anchor=(1, 1), fontsize='small')
        plt.tight_layout()

        img_bytes = io.BytesIO()
        plt.savefig(img_bytes, format='png')
        img_bytes.seek(0)
        img_data = base64.b64encode(img_bytes.read()).decode('utf-8')
        img_element = html.Img(src='data:image/png;base64,{}'.format(img_data))
        images.append(img_element)
        png_bytes_list.append(img_bytes.getvalue())

        svg_bytes = io.BytesIO()
        plt.savefig(svg_bytes, format='svg')
        svg_bytes.seek(0)
        svg_data = svg_bytes.getvalue()
        svg_bytes_list.append(svg_data)

        plt.close(fig)

    return images, svg_bytes_list, png_bytes_list
def plot_bar_with_p_values(data, divs, metrics, ordered_genotypes, selected_colors):
    markers = ['o', 's', 'D', '^', 'v', 'x', '*', 'H', '8', 'p']

    images = []
    svg_bytes_list = []
    png_bytes_list = []

    for output_type in metrics:
        total_genotypes = len(ordered_genotypes)
        #print(f"Number of unique Genotypes: {total_genotypes}")

        output_arrays = {genotype: [] for genotype in ordered_genotypes}
        chip_arrays = {genotype: [] for genotype in ordered_genotypes}
        well_arrays = {genotype: [] for genotype in ordered_genotypes}

        for i in divs:
            for genotype in ordered_genotypes:
                temp_df = data.loc[(data['DIV'] == i) & (data['NeuronType'].str.strip() == genotype)]
                output_arrays[genotype].append(np.array(temp_df[output_type]))
                chip_arrays[genotype].append(np.array(temp_df['Chip_ID']))
                well_arrays[genotype].append(np.array(temp_df['Well']))

        bar_width = 0.20
        gap_between_bars = 0.0
        gap_between_groups = bar_width

        total_bar_group_width = total_genotypes * bar_width + (total_genotypes - 1) * gap_between_bars
        total_plot_width = len(divs) * (total_bar_group_width + gap_between_groups) + gap_between_groups

        x_genotype = {genotype: [] for genotype in ordered_genotypes}
        base_x_coordinate = np.arange(len(divs)) * (total_bar_group_width + gap_between_groups) + gap_between_groups + bar_width
        offset = (total_genotypes * bar_width + (total_genotypes - 1) * gap_between_bars) / 2
        centered_x = base_x_coordinate - offset + bar_width / 2
        for i, genotype in enumerate(ordered_genotypes):
            x_genotype[genotype] = centered_x + (i) * (bar_width + gap_between_bars)

        fig, ax = plt.subplots()
        mean_data_all = {}
        yerr_data_all = {}
        n_data_all = {}

        for i, genotype in enumerate(ordered_genotypes):
            y_data = output_arrays[genotype]
            chipy_data = chip_arrays[genotype]
            welly_data = well_arrays[genotype]
            mean_data = [np.mean(yi[np.isfinite(yi)]) if len(yi[np.isfinite(yi)]) > 0 else np.nan for yi in y_data]
            yerr_data = [np.std(yi[np.isfinite(yi)], ddof=1) / np.sqrt(np.size(yi[np.isfinite(yi)])) if len(yi[np.isfinite(yi)]) > 1 else np.nan for yi in y_data]
            n_data = [len(yi[np.isfinite(yi)]) for yi in y_data]
            mean_data_all[genotype] = mean_data
            yerr_data_all[genotype] = yerr_data
            n_data_all[genotype] = n_data

            ax.bar(x_genotype[genotype], mean_data, yerr=yerr_data, capsize=3, width=bar_width, color=selected_colors.get(genotype, 'blue'), label=genotype, alpha=0.7)

            jitter_amount = 0.05
            for j in range(len(x_genotype[genotype])):
                combined_data = [str(chip) + str(well) for chip, well in zip(chipy_data[j], welly_data[j])]
                for k in range(len(y_data[j])):
                    ax.scatter(
                        x_genotype[genotype][j] + np.random.uniform(-jitter_amount, jitter_amount, 1),
                        y_data[j][k],
                        s=10,
                        color=selected_colors.get(genotype, 'blue'),
                        marker=markers[i % len(markers)]
                    )

        valid_arrays = [array for genotype_arrays in output_arrays.values() for array in genotype_arrays if len(array) > 0]
        if valid_arrays:
            maxim = max(max(array) for array in valid_arrays)
        else:
            maxim = 0

        for i in range(len(base_x_coordinate)):
            count = 1
            p_values = []
            for j, genotype1 in enumerate(ordered_genotypes):
                for k, genotype2 in enumerate(ordered_genotypes):
                    if j < k:
                        mean1, sem1, n1 = mean_data_all[genotype1][i], yerr_data_all[genotype1][i], n_data_all[genotype1][i]
                        mean2, sem2, n2 = mean_data_all[genotype2][i], yerr_data_all[genotype2][i], n_data_all[genotype2][i]
                        if n1 > 1 and n2 > 1:
                            sed = np.sqrt(sem1 ** 2.0 + sem2 ** 2.0)
                            t_stat = (mean1 - mean2) / sed
                            degreef = n1 + n2 - 2

                            p_value = (1.0 - t.cdf(abs(t_stat), degreef)) * 2.0
                            p_values.append([mean1, sem1, mean2, sem2, p_value])

                            x1, x2 = x_genotype[genotype1][i], x_genotype[genotype2][i]
                            sign = "***" if p_value <= 0.001 else "**" if p_value <= 0.01 else "*" if p_value <= 0.05 else "ns"
                            if sign != 'ns':
                                ax.plot([x1, x2], [maxim + 0.05 * maxim * (count)] * 2, 'k', linewidth=1.5)
                                ax.text((x1 + x2) / 2, maxim + 0.05 * maxim * (count), f"{sign}", ha='center', va='bottom', fontsize=7)
                                ax.axvline(x1, color='black', linestyle=':', linewidth=0.5)
                                ax.axvline(x2, color='black', linestyle=':', linewidth=0.5)
                            count += 1
            y_ticks = ax.get_yticks()
            y_ticks_selected = [y_ticks[0], y_ticks[len(y_ticks) // 2], y_ticks[-1]]
            ax.set_yticks(y_ticks_selected)

        plt.title(f"{output_type}", fontsize=14)
        plt.xlabel('DIV', fontsize=12)
        plt.ylabel(f"{output_type}", fontsize=12)
        plt.xticks(base_x_coordinate, divs, fontsize=10)
        plt.legend(title='NeuronType', loc='upper left', bbox_to_anchor=(1, 1), fontsize='small')
        plt.tight_layout()
        plt.xlim([0, total_plot_width])

        img_bytes = io.BytesIO()
        plt.savefig(img_bytes, format='png')
        img_bytes.seek(0)
        img_data = base64.b64encode(img_bytes.read()).decode('utf-8')
        img_element = html.Img(src='data:image/png;base64,{}'.format(img_data))
        images.append(img_element)
        png_bytes_list.append(img_bytes.getvalue())

        svg_bytes = io.BytesIO()
        plt.savefig(svg_bytes, format='svg')
        svg_bytes.seek(0)
        svg_data = svg_bytes.getvalue()
        svg_bytes_list.append(svg_data)

        plt.close(fig)

    return images, svg_bytes_list, png_bytes_list
def plot_isi_graph(data_df, selected_divs, metric):
    filtered_df = data_df[data_df['DIV'].isin(selected_divs)]
    neuron_types = filtered_df['NeuronType'].unique()
    colors = {neuron_type: plt.cm.tab10(i) for i, neuron_type in enumerate(neuron_types)}
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
                print("error mask and edges length don't match")
                print("edges here", len(edges[:-1]), "mask length", len(mask))
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

            plt.plot(max_edges[:-1], mean_bins, color=colors[neuron_type], label=neuron_type)
            plt.fill_between(max_edges[:-1], mean_bins - stderr_bins, mean_bins + stderr_bins, color=colors[neuron_type], alpha=0.3)

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

    plt.close()

    return img_element, svg_bytes, img_bytes.getvalue()

def min_max_normalize(data):
    min_val = np.min(data)
    max_val = np.max(data)
    return (data - min_val) / (max_val - min_val)

