import pandas as pd
import numpy as np
import scipy.io
import io
import base64

def parse_contents(contents, filename):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
            df.columns = df.columns.str.strip()  # Strip whitespace from column names
            df['Chip_Well'] = df['Chip_ID'].astype(str) + '_' + df['Well'].astype(str)
        else:
            return html.Div([
                'Unsupported file format'
            ])
    except Exception as e:
        print(f"Error processing file {filename}: {e}")
        return html.Div([
            'There was an error processing this file.'
        ])

    return df

def parse_mat(contents):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    mat_data = scipy.io.loadmat(io.BytesIO(decoded))

    all_records = np.squeeze(mat_data['allExtMetrics'])
    refined_results = []
    
    for idx, record in enumerate(all_records):
        try:
            if record.size > 0:
                for item in record:
                    if item.size > 0:
                        div = item[0][0][0]['DIV'][0][0]
                        well = item[0][0][0]['Well'][0][0]
                        run_id = item[0][0][0]['Run_ID'][0][0]
                        chip_id = item[0][0][0]['Chip_ID'][0]
                        neuron_type = item[0][0][0]['NeuronType'][0]
                        network_ap_freq_bins = item[0][0][0]['networkAPFreqBins'][0]
                        network_ap_freq_edges = item[0][0][0]['networkAPFreqEdges'][0]
                        burst_ap_freq_bins = item[0][0][0]['burstAPFreqBins'][0]
                        burst_ap_freq_edges = item[0][0][0]['burstAPFreqEdges'][0]
                        nonburst_ap_freq_bins = item[0][0][0]['nonburstAPFreqBins'][0]
                        nonburst_ap_freq_edges = item[0][0][0]['nonburstAPFreqEdges'][0]
                        refined_entry = {
                            'Run_ID': run_id,
                            'Chip_ID': chip_id,
                            'DIV': div,
                            'WellRecord': well,
                            'NeuronType': neuron_type,
                            'networkAPFreqBins': network_ap_freq_bins,
                            'networkAPFreqEdges': network_ap_freq_edges,
                            'burstAPFreqBins': burst_ap_freq_bins,
                            'burstAPFreqEdges': burst_ap_freq_edges,
                            'nonburstAPFreqBins': nonburst_ap_freq_bins,
                            'nonburstAPFreqEdges': nonburst_ap_freq_edges
                        }
                        refined_results.append(refined_entry)
        except Exception as e:
            print(f"An error occurred at entry {idx + 1}: {str(e)}")
            continue
    
    data_df = pd.DataFrame.from_dict(refined_results)
    available_metrics = ['networkAPFreqBins', 'burstAPFreqBins', 'nonburstAPFreqBins']
    return data_df, available_metrics
