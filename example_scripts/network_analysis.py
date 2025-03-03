import numpy as np
from scipy.stats import norm
from scipy.signal import convolve
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from MEA_Analysis.IPNAnalysis.helper_functions import detect_bursts_statistics # NOTE: Mandar code, not mine.
import os

'''derive optimization targets'''
# def build_network_metric_targets_dict(network_metrics):
    
#     #initialize network_metric_targets
#     bursting_data = network_metrics['bursting_data']
#     mega_bursting_data = network_metrics['mega_bursting_data']
#     duration_seconds = network_metrics['timeVector'][-1]
    
#     network_metric_targets = {
#         #General Data
#         'source': network_metrics['source'], # 'simulated' or 'experimental'
#         #'timeVector': network_metrics['timeVector'],
        
#         # Spiking Data
#         'spiking_data': {
#             #'spike_times': network_metrics['spiking_data']['spike_times'],
#             #'spiking_times_by_unit': network_metrics['spiking_data']['spiking_times_by_unit'],
#             #'spiking_data_by_unit': network_metrics['spiking_data']['spiking_data_by_unit'],
#             'spiking_summary_data': {
#                 'MeanFireRate': {
#                     'target': network_metrics['spiking_data']['spiking_summary_data']['MeanFireRate'],
#                     'min': get_min(network_metrics['spiking_data']['spiking_data_by_unit'], 'FireRate'),
#                     'max': get_max(network_metrics['spiking_data']['spiking_data_by_unit'], 'FireRate'),
#                     #'E_I_ratio_assumption': 0.7, #E_I_ratio = 5  # 1:5 ratio of E to I neurons #TODO: need scientific basis for this assumption
                    
#                     ## TODO: Get scientific basis for these assumptions # aw 2025-01-26 21:01:23
#                     #E
#                     'max_E_assumption': 5, # Hz  #outside of bursts 
#                     'min_E_assumption': 0.5, # Hz
#                     'max_E_assumption_inBurst': 50, # Hz #inside of bursts
#                     'min_E_assumption_inBurst': 10, # Hz
                    
#                     #I
#                     'max_I_assumption': 10, # Hz
#                     'min_I_assumption': 1, # Hz
#                     'max_I_assumption_inBurst': 100, # Hz
#                     'min_I_assumption_inBurst': 20, # Hz                    
                                        
#                     'weight': 1, # TODO: update these with Nfactors
#                 },
#                 'CoVFireRate': {
#                     'target': network_metrics['spiking_data']['spiking_summary_data']['CoVFireRate'],
#                     #'min': get_min(network_metrics['spiking_data']['spiking_data_by_unit'], 'fr_CoV'),
#                     #'max': get_max(network_metrics['spiking_data']['spiking_data_by_unit'], 'fr_CoV'),
#                     'weight': 1, # TODO: update these with Nfactors
#                 },
#                 'MeanISI': {
#                     'target': network_metrics['spiking_data']['spiking_summary_data']['MeanISI'],
#                     'min': get_min(network_metrics['spiking_data']['spiking_data_by_unit'], 'meanISI'),
#                     'max': get_max(network_metrics['spiking_data']['spiking_data_by_unit'], 'meanISI'),
#                     'weight': 1, # TODO: update these with Nfactors
#                 },
#                 'CoV_ISI': {
#                     'target': network_metrics['spiking_data']['spiking_summary_data']['CoV_ISI'],
#                     # 'min': get_min(network_metrics['spiking_data']['spiking_data_by_unit'], 'isi_CoV'),
#                     # 'max': get_max(network_metrics['spiking_data']['spiking_data_by_unit'], 'isi_CoV'),
#                     'weight': 1, # TODO: update these with Nfactors
#                 },
#             },
#         },
        
        
#         #Bursting Data
#         'bursting_data': {
#             'bursting_summary_data': {
#                 'MeanBurstRate': {
#                     'target': bursting_data['bursting_summary_data'].get('mean_Burst_Rate'),
#                     # 'min': get_min_burst(bursting_data['bursting_data_by_unit'], 'bursts', duration_seconds), # NOTE: Calculated as individual unit burst participation rate.
#                     # 'max': get_max_burst(bursting_data['bursting_data_by_unit'], 'bursts', duration_seconds), # NOTE: Calculated as individual unit burst participation rate.
#                     'min': 0,
#                     'max': None,
#                     'weight': 1,
#                 },                
#                 'MeanWithinBurstISI': {
#                     'target': bursting_data['bursting_summary_data'].get('MeanWithinBurstISI'),
#                     'min': get_min(bursting_data['bursting_data_by_unit'], 'mean_isi_within'),
#                     'max': get_max(bursting_data['bursting_data_by_unit'], 'mean_isi_within'),
#                     'weight': 1,
#                 },
#                 'CovWithinBurstISI': {
#                     'target': bursting_data['bursting_summary_data'].get('CoVWithinBurstISI'),
#                     # 'min': get_min(bursting_data['bursting_data_by_unit'], 'cov_isi_within'),
#                     # 'max': get_max(bursting_data['bursting_data_by_unit'], 'cov_isi_within'),
#                     'weight': 1,
#                 },
#                 'MeanOutsideBurstISI': {
#                     'target': bursting_data['bursting_summary_data'].get('MeanOutsideBurstISI'),
#                     'min': get_min(bursting_data['bursting_data_by_unit'], 'mean_isi_outside'),
#                     'max': get_max(bursting_data['bursting_data_by_unit'], 'mean_isi_outside'),
#                     'weight': 1,
#                 },
#                 'CoVOutsideBurstISI': {
#                     'target': bursting_data['bursting_summary_data'].get('CoVOutsideBurstISI'),
#                     # 'min': get_min(bursting_data['bursting_data_by_unit'], 'cov_isi_outside'),
#                     # 'max': get_max(bursting_data['bursting_data_by_unit'], 'cov_isi_outside'),
#                     'weight': 1,
#                 },
#                 'MeanNetworkISI': {
#                     'target': bursting_data['bursting_summary_data'].get('MeanNetworkISI'),
#                     'min': get_min(bursting_data['bursting_data_by_unit'], 'mean_isi_all'),
#                     'max': get_max(bursting_data['bursting_data_by_unit'], 'mean_isi_all'),
#                     'weight': 1,
#                 },
#                 'CoVNetworkISI': {
#                     'target': bursting_data['bursting_summary_data'].get('CoVNetworkISI'),
#                     # 'min': get_min(bursting_data['bursting_data_by_unit'], 'cov_isi_all'),
#                     # 'max': get_max(bursting_data['bursting_data_by_unit'], 'cov_isi_all'),
#                     'weight': 1,
#                 },
#                 'NumUnits': {
#                     'target': bursting_data['bursting_summary_data'].get('NumUnits'),
#                     # 'min': 1,
#                     # 'max': None,
#                     # 'weight': 1,
#                 },
#                 'Number_Bursts': {
#                     'target': bursting_data['bursting_summary_data'].get('Number_Bursts'),
#                     # 'min': 1,
#                     # 'max': None,
#                     # 'weight': 1, #TODO: update these with Nfactors
#                 },
#                 'mean_IBI': {
#                     'target': bursting_data['bursting_summary_data'].get('mean_IBI'),
#                     'min': None,
#                     'max': None,
#                     'weight': 1, #TODO: update these with Nfactors
#                 },
#                 'cov_IBI': {
#                     'target': bursting_data['bursting_summary_data'].get('cov_IBI'),
#                     'min': None,
#                     'max': None,
#                     'weight': 1, #TODO: update these with Nfactors
#                 },
#                 'mean_Burst_Peak': {
#                     'target': bursting_data['bursting_summary_data'].get('mean_Burst_Peak'),
#                     'min': None,
#                     'max': None,
#                     'weight': 1, #TODO: update these with Nfactors
#                 },
#                 'cov_Burst_Peak': {
#                     'target': bursting_data['bursting_summary_data'].get('cov_Burst_Peak'),
#                     'min': None,
#                     'max': None,
#                     'weight': 1, #TODO: update these with Nfactors
#                 },
#                 'fano_factor': {
#                     'target': bursting_data['bursting_summary_data'].get('fano_factor'),
#                     'min': None,
#                     'max': None,
#                     'weight': 1, #TODO: update these with Nfactors
#                 },
#                 'baseline': {
#                     'target': bursting_data['bursting_summary_data'].get('baseline'),
#                     'min': None,
#                     'max': None,
#                     'weight': 1,
#                 },
#             },
#         },
        
#         #Mega Bursting Data
#         'mega_bursting_data': {
#             'bursting_summary_data': {
#                 'MeanBurstRate': { # WHole network burst rate
#                     'target': mega_bursting_data['bursting_summary_data'].get('mean_Burst_Rate'),
#                     # 'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'bursts'),
#                     # 'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'bursts'),
#                     'min' : 0,
#                     'max' : None,
#                     'weight': 1,
                

#                     # 'min': get_min_burst(mega_bursting_data['bursting_data_by_unit'], 'bursts', duration_seconds), # NOTE: Calculated as individual unit burst participation rate.
#                     # 'max': get_max_burst(mega_bursting_data['bursting_data_by_unit'], 'bursts', duration_seconds), # NOTE: Calculated as individual unit burst participation rate.
#                     # 'weight': 1,
                    
#                     # aw 2025-01-25 21:12:43 - #TODO: need to maybe rethink how I"m getting min and maxes for a lot of these data.
#                     # like in this case, using min and max burst participation rate for individual units to set min and max for
#                     # whole network burst rate doesnt make sense. Unless at least one unit participates in al bursts - max for individual
#                     # units will always be less than burst rate for whole network.
#                 },                  
#                 'MeanWithinBurstISI': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('MeanWithinBurstISI'),
#                     'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'mean_isi_within'),
#                     'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'mean_isi_within'),
#                     # 'min' : 0,
#                     # 'max' : None,
#                     'weight': 1,
#                 },
#                 'CovWithinBurstISI': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('CoVWithinBurstISI'),
#                     # 'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'cov_isi_within'),
#                     # 'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'cov_isi_within'),
#                     'weight': 1,
#                 },
#                 'MeanOutsideBurstISI': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('MeanOutsideBurstISI'),
#                     'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'mean_isi_outside'),
#                     'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'mean_isi_outside'),
#                     'weight': 1,
#                 },
#                 'CoVOutsideBurstISI': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('CoVOutsideBurstISI'),
#                     # 'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'cov_isi_outside'),
#                     # 'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'cov_isi_outside'),
#                     'weight': 1,
#                 },
#                 'MeanNetworkISI': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('MeanNetworkISI'),
#                     'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'mean_isi_all'),
#                     'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'mean_isi_all'),
#                     'weight': 1,
#                 },
#                 'CoVNetworkISI': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('CoVNetworkISI'),
#                     # 'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'cov_isi_all'),
#                     # 'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'cov_isi_all'),
#                     'weight': 1,
#                 },
#                 'NumUnits': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('NumUnits'),
#                     # 'min': 1,
#                     # 'max': None,
#                     # 'weight': 1,
#                 },
#                 'Number_Bursts': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('Number_Bursts'),
#                     # 'min': 1,
#                     # 'max': None,
#                     # 'weight': 1, #TODO: update these with Nfactors
#                 },
#                 'mean_IBI': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('mean_IBI'),
#                     'min': None,
#                     'max': None,
#                     'weight': 1, #TODO: update these with Nfactors
#                 },
#                 'cov_IBI': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('cov_IBI'),
#                     'min': None,
#                     'max': None,
#                     'weight': 1, #TODO: update these with Nfactors
#                 },
#                 'mean_Burst_Peak': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('mean_Burst_Peak'),
#                     'min': None,
#                     'max': None,
#                     'weight': 1, #TODO: update these with Nfactors
#                 },
#                 'cov_Burst_Peak': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('cov_Burst_Peak'),
#                     'min': None,
#                     'max': None,
#                     'weight': 1, #TODO: update these with Nfactors
#                 },
#                 'fano_factor': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('fano_factor'),
#                     'min': None,
#                     'max': None,
#                     'weight': 1, #TODO: update these with Nfactors
#                 },
#                 'baseline': {
#                     'target': mega_bursting_data['bursting_summary_data'].get('baseline'),
#                     'min': None,
#                     'max': None,
#                     'weight': 1,
#                 },
#             },        
#         }
#     }
#     return network_metric_targets

def get_min(data, key):
    data_min = min(unit[key] for unit in data.values() if unit[key] is not None and unit[key] > 0 and not (
        unit[key] != unit[key] or unit[key] == float('inf') or unit[key] == float('-inf')
        ))
    
    #expand the code above into a more readable format
    # data_min = []
    # for unit in data.values():
    #     if unit[key] is not None and unit[key] > 0 and not (unit[key] != unit[key] or unit[key] == float('inf') or unit[key] == float('-inf')):
    #         data_min.append(unit[key])
    # data_min = min(data_min)
    
    
    
    #if data is an array with one value, return that value
    # if isinstance(data_min, list) and len(data_min) == 1:
    #     data_min = data_min[0]
    # try: 
    #     len(data_min)
    #     print('data_min is an array')
    # except TypeError: pass
    data_min = handle_numpy_float64(data_min)
    #print('data_min:', data_min)
    return data_min

def get_max(data, key):
    data_list = [unit[key] for unit in data.values() if not (unit[key] is None or unit[key] != unit[key] or unit[key] == float('inf') or unit[key] == float('-inf'))]
    data_max = max(unit[key] for unit in data.values() if not (unit[key] is None or unit[key] != unit[key] or unit[key] == float('inf') or unit[key] == float('-inf')))
    #if data is an array with one value, return that value
    # if type(data_max) == list:
    #     data_max = data_max[0]
    data_max = handle_numpy_float64(data_max)
    return data_max

def get_min_burst(data, key, duration_seconds):
    
    #start min_rate at positive infinity
    min_rate = float('inf')
    for unit in data.values():
        num_bursts = len(unit[key])
        burst_participation_rate = num_bursts / duration_seconds        
        try:
            if burst_participation_rate < min_rate:
                min_rate = burst_participation_rate
        except:
            continue
    
    return min_rate

def get_max_burst(data, key, duration_seconds):
    
    #start min_rate at positive infinity
    #min_rate = float('inf')
    max_rate = 0
    for unit in data.values():
        num_bursts = len(unit[key])
        burst_participation_rate = num_bursts / duration_seconds        
        try:
            if burst_participation_rate > max_rate:
                max_rate = burst_participation_rate
        except:
            continue
    
    return max_rate  

def save_network_metric_dict_with_timestamp(network_metrics, network_metrics_targets, output_dir):
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    date = datetime.now().strftime("%Y%m%d")
    #output_path = output_dir + '/fitness_args_' + timestamp + '.py'
    output_path = os.path.join(output_dir, f'{date}_features.py')
    
    output_path = output_path.replace('-', '_') #conver all '-' to '_' in the timestamp
    print('output_path:', output_path)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    def format_dict(d, indent=0):
        formatted_str = ''        
        if indent==0: formatted_str += 'import numpy as np\n' # add import statement to the top of the file
        for key, value in d.items():
            if isinstance(value, dict):
                formatted_str += ' ' * indent + f"'{key}': {{\n" + format_dict(value, indent + 4) + ' ' * indent + '},\n'
                
                # try: key = int(key)
                # except: pass
                
                # if isinstance(key, str):
                #     formatted_str += ' ' * indent + f"'{key}': {{\n" + format_dict(value, indent + 4) + ' ' * indent + '},\n'
                # elif isinstance(key, int):
                #     formatted_str += ' ' * indent + f"{key}: {{\n" + format_dict(value, indent + 4) + ' ' * indent + '},\n'
                # else:
                #     print('Warning: Unrecognized key type:', type(key), 'key:', key)
                #     formatted_str += ' ' * indent + f"'{key}': {{\n" + format_dict(value, indent + 4) + ' ' * indent + '},\n'
            elif isinstance(value, list):
                # convert list to string representation
                list_str = ', '.join([str(item) for item in value])
                formatted_str += ' ' * indent + f"'{key}': [{list_str}],\n"
                
                # formatted_str += ' ' * indent + f"'{key}': ["
                # for i, item in enumerate(value):
                #     formatted_str += ' ' + f"{item},"
                #     # if i < len(value) - 1:
                #     #     formatted_str += ','
                #     #formatted_str += '\n'
                # # remove last comma
                # formatted_str = formatted_str[:-1]
                # formatted_str += '],\n'
            elif isinstance(value, (int, float)):
                formatted_str += ' ' * indent + f"'{key}': {value},\n"
            elif isinstance(value, str):
                formatted_str += ' ' * indent + f"'{key}': '{value}',\n"
            elif isinstance(value, np.ndarray):
                #formatted_str += ' ' * indent + f"'{key}': np.array({value.tolist()}),\n"
                #convert to tuple
                formatted_str += ' ' * indent + f"'{key}': {tuple(value.tolist())},\n"
            elif isinstance(value, tuple):
                formatted_str += ' ' * indent + f"'{key}': {value},\n"
            elif isinstance(value, bool):
                formatted_str += ' ' * indent + f"'{key}': {value},\n"
            elif isinstance(value, None.__class__):
                formatted_str += ' ' * indent + f"'{key}': None,\n"
            else:
                print(type(value))
                print('Warning: Unrecognized data type for key:', key, 'value:', value)
                formatted_str += ' ' * indent + f"'{key}': {value},\n"
        return formatted_str

    # Convert network_metric_targets to a formatted string
    formatted_fitness_args = 'fitness_args = {\n' + format_dict(network_metrics_targets, 4) + '}'
    
    # append some additional info to the formatted string
    formatted_fitness_args += (
        f'\n\n#=== FITNESS FUNCTION ARGUMENTS ===\n'
        f'# fitness function\n'
        f'fitnessFuncArgs = {{}}\n'
        f'fitnessFuncArgs[\'targets\'] = fitness_args\n'
        f'number_of_units = fitness_args[\'bursting_data\'][\'bursting_summary_data\'][\'NumUnits\'][\'target\']\n'
        f'fitnessFuncArgs[\'features\'] = {{\n'
        f'  # tweaking the ratio a bit based on Roy\'s suggestion\n'
        f'  \'num_excite\': len(fitness_args[\'excit_units\']),\n'
        f'  \'num_inhib\': len(fitness_args[\'inhib_units\']),\n'
        f'}}\n'
        f'fitnessFuncArgs[\'maxFitness\'] = 1000\n'
        f'#fitnessFuncArgs[\'skip_existing\'] = False\n'
    )

    with open(output_path, 'w') as f:
        f.write(formatted_fitness_args)

    print('Updated fitness args saved to:', output_path)

def handle_numpy_float64(data):
    try: 
        if isinstance(data, np.ndarray):
            data = data.tolist()
    except: 
        pass
    return data

'''sim funcs'''
def get_simulated_network_activity_metrics(simData=None, popData=None, **kwargs):
    import time
    
    #candidate_label = kwargs.get('candidate_label', None)
    print('') #for formatting
    print('Calculating Network Activity Metrics for Simulated Data...')
    start_time = time.time()
    assert simData is not None, 'No simData provided'
    #this part should be useful for fitness during simulation
    # if simData is None:
    #     try: 
    #         simData = sim.simData
    #         print('Using simData from netpyne.sim.simData')
    #     except:
    #         print('No simData provided or found in netpyne.sim.simData')
    #         return None

    #initialize network_data
    #network_data = init_network_data_dict()
    # # HACK:
    # from RBS_network_models.network_analysis import init_network_data_dict, extract_bursting_activity_data
    # # HACK:
    network_data = init_network_data_dict()
    rasterData = simData.copy()
    
    #check if rasterData['spkt'] is empty
    if len(rasterData['spkt']) == 0:
        print('No spike times found in rasterData')
        return None
    
    #convert time to seconds - get initially available data
    spike_times = np.array(rasterData['spkt']) / 1000
    timeVector = np.array(rasterData['t']) / 1000
    spike_times_by_unit = {int(i): spike_times[rasterData['spkid'] == i] for i in np.unique(rasterData['spkid'])} #mea_analysis_pipeline.py expects spike_times as dictionary    
    
    #extract spiking metrics from simulated data
    try: 
        extract_metrics_from_simulated_data(spike_times, timeVector, spike_times_by_unit, rasterData, popData, **kwargs)
    except Exception as e:
        print(f'Error extracting metrics from simulated data: {e}')
        pass
    
    #extract bursting metrics from simulated data (but this one works for both simulated and experimental data)
    try: 
        extract_bursting_activity_data(spike_times, spike_times_by_unit, **kwargs)
    except Exception as e:
        print(f'Error calculating bursting activity: {e}')
        pass
    
    #return network_data
    
    def convert_single_element_arrays(data):
        if isinstance(data, dict):
            for key, value in data.items():
                data[key] = convert_single_element_arrays(value)
        elif isinstance(data, np.ndarray) and data.size == 1:
            return data.item()
        return data

    network_data = convert_single_element_arrays(network_data)    
    print('Network Activity Metrics Calculated and Extracted!')
    print(f'Elapsed time: {time.time() - start_time} seconds')
    print('') #for formatting    
    return network_data

def get_CoV_fr_simulated(spike_times, total_duration, window_size=1.0): 
    """
    Calculate the Coefficient of Variation of Firing Rate (CoV FR) over time windows.

    Parameters:
    - spike_times: List or array of spike times for a single unit (in seconds).
    - total_duration: Total duration of the simulation (in seconds).
    - window_size: Size of the time window for calculating firing rates (default: 1.0 second).

    Returns:
    - CoV_fr: Coefficient of Variation of Firing Rate (float).
    """
    #TODO import window from convolution_params?
    if len(spike_times) == 0:
        # If no spikes, CoV FR is undefined (return NaN)
        return np.nan

    # Divide the total duration into non-overlapping time windows
    #window_size = 0.1
    total_duration = total_duration[-1]
    num_windows = int(total_duration / window_size)
    window_edges = np.linspace(0, total_duration, num_windows + 1)

    # Count spikes in each window
    spike_counts = np.histogram(spike_times, bins=window_edges)[0]

    # Convert spike counts to firing rates (spikes per second)
    firing_rates = spike_counts / window_size

    # Compute CoV FR: standard deviation divided by mean
    mean_fr = np.mean(firing_rates)
    std_fr = np.std(firing_rates)

    # Avoid division by zero if mean_fr is 0
    if mean_fr == 0:
        return np.nan

    CoV_fr = std_fr / mean_fr
    return CoV_fr

def calculate_simulated_network_activity_metrics(spike_times_by_unit):
    # Initialize spiking data dictionary
    network_data['simulated_data']['spiking_data_by_unit'] = {}
    
    # Iterate over each unit to calculate individual metrics
    for unit, spike_times in spike_times_by_unit.items():
        if len(spike_times) < 2:
            meanISI = np.nan
            CoV_ISI = np.nan
        else:
            isi = np.diff(spike_times)
            meanISI = np.mean(isi)
            CoV_ISI = np.std(isi) / meanISI  # Correct CoV calculation for ISI
        
        fr = len(spike_times) / network_data['timeVector'][-1]  # Firing rate (spikes per second)
        
        # Calculate CoV FR using the get_CoV_fr function with a 1-second time window
        CoV_fr = get_CoV_fr_simulated(spike_times, network_data['timeVector'])

        # Store calculated metrics for each unit
        network_data['simulated_data']['spiking_data_by_unit'][unit] = {
            'FireRate': fr,
            'CoV_fr': CoV_fr,
            'meanISI': meanISI,
            'CoV_ISI': CoV_ISI,
            'spike_times': spike_times,
        }
        
        # Determine population type
        if unit in network_data['simulated_data']['E_Gids']:
            network_data['simulated_data']['spiking_data_by_unit'][unit]['pop'] = 'E'
        elif unit in network_data['simulated_data']['I_Gids']:
            network_data['simulated_data']['spiking_data_by_unit'][unit]['pop'] = 'I'
        else:
            raise ValueError(f'Unit {unit} not found in E_Gids or I_Gids')

    # Extract E and I gids
    E_Gids = network_data['simulated_data']['E_Gids']
    I_Gids = network_data['simulated_data']['I_Gids']

    # Calculate mean and CoV metrics for excitatory and inhibitory populations
    E_CoV_ISI = np.nanmean([network_data['simulated_data']['spiking_data_by_unit'][unit]['CoV_ISI'] 
                            for unit in E_Gids if unit in network_data['simulated_data']['spiking_data_by_unit']])
    I_CoV_ISI = np.nanmean([network_data['simulated_data']['spiking_data_by_unit'][unit]['CoV_ISI'] 
                            for unit in I_Gids if unit in network_data['simulated_data']['spiking_data_by_unit']])
    MeanISI_E = np.nanmean([network_data['simulated_data']['spiking_data_by_unit'][unit]['meanISI'] 
                            for unit in E_Gids if unit in network_data['simulated_data']['spiking_data_by_unit']])
    MeanISI_I = np.nanmean([network_data['simulated_data']['spiking_data_by_unit'][unit]['meanISI'] 
                            for unit in I_Gids if unit in network_data['simulated_data']['spiking_data_by_unit']])
    CoVFiringRate_E = np.nanmean([network_data['simulated_data']['spiking_data_by_unit'][unit]['CoV_fr'] 
                                  for unit in E_Gids if unit in network_data['simulated_data']['spiking_data_by_unit']])
    CoVFiringRate_I = np.nanmean([network_data['simulated_data']['spiking_data_by_unit'][unit]['CoV_fr'] 
                                  for unit in I_Gids if unit in network_data['simulated_data']['spiking_data_by_unit']])
    
    # Add population-level metrics to the network data
    network_data['simulated_data']['CoV_ISI_E'] = E_CoV_ISI
    network_data['simulated_data']['CoV_ISI_I'] = I_CoV_ISI
    network_data['simulated_data']['MeanISI_E'] = MeanISI_E
    network_data['simulated_data']['MeanISI_I'] = MeanISI_I
    network_data['simulated_data']['CoVFireRate_E'] = CoVFiringRate_E
    network_data['simulated_data']['CoVFireRate_I'] = CoVFiringRate_I
            
def extract_metrics_from_simulated_data(spike_times, timeVector, spike_times_by_unit, rasterData, popData, **kwargs):
    
    #add immediately available data to network_data
    network_data['source'] = 'simulated'
    network_data['timeVector'] = timeVector #seconds
    network_data['spiking_data']['spike_times'] = spike_times
    network_data['simulated_data']['soma_voltage'] = rasterData['soma_voltage']
    network_data['spiking_data']['spiking_times_by_unit'] = spike_times_by_unit

    #add uniquely simulated data to network_data
    network_data['simulated_data']['E_Gids'] = popData['E']['cellGids']
    network_data['simulated_data']['I_Gids'] = popData['I']['cellGids']
    network_data['simulated_data']['MeanFireRate_E'] = rasterData['popRates']['E']
    network_data['simulated_data']['MeanFireRate_I'] = rasterData['popRates']['I']
    
    #get simulated network activity metrics
    calculate_simulated_network_activity_metrics(spike_times_by_unit)
    
    #derive remaining spiking data from simulated data
    num_E = len(network_data['simulated_data']['E_Gids'])
    num_I = len(network_data['simulated_data']['I_Gids'])
    
    #overall mean firing rate
    E_popRates = rasterData['popRates']['E']
    I_popRates = rasterData['popRates']['I']
    network_data['spiking_data']['spiking_summary_data']['MeanFireRate'] = (
        (E_popRates * num_E + I_popRates * num_I) / (num_E + num_I)
        )
    
    #overall CoV firing rate
    E_CoV = network_data['simulated_data']['CoVFireRate_E']
    I_CoV = network_data['simulated_data']['CoVFireRate_I']
    network_data['spiking_data']['spiking_summary_data']['CoVFireRate'] = (
        (E_CoV * num_E + I_CoV * num_I) / (num_E + num_I)
        )
    
    #overall mean ISI
    E_meanISI = network_data['simulated_data']['MeanISI_E']
    I_meanISI = network_data['simulated_data']['MeanISI_I']
    network_data['spiking_data']['spiking_summary_data']['MeanISI'] = (
        (E_meanISI * num_E + I_meanISI * num_I) / (num_E + num_I)
        )
    
    #overall CoV ISI
    E_CoV_ISI = network_data['simulated_data']['CoV_ISI_E']
    I_CoV_ISI = network_data['simulated_data']['CoV_ISI_I']
    network_data['spiking_data']['spiking_summary_data']['CoV_ISI'] = (
        (E_CoV_ISI * num_E + I_CoV_ISI * num_I) / (num_E + num_I)
        )
    
    #spiking data by unit - just go through the simulated version of the data but de-identify pop
    simulated_spiking_data_by_unit = network_data['simulated_data']['spiking_data_by_unit']
    spiking_data_by_unit = {}
    for unit, unit_data in simulated_spiking_data_by_unit.items():
        spiking_data_by_unit[unit] = {
            'FireRate': unit_data['FireRate'],
            'CoV_fr': unit_data['CoV_fr'],
            'meanISI': unit_data['meanISI'],
            'CoV_ISI': unit_data['CoV_ISI'],
            'spike_times': unit_data['spike_times'],
        }
    network_data['spiking_data']['spiking_data_by_unit'] = spiking_data_by_unit

def get_simulated_network_activity_metrics(conv_params, mega_params, simData=None, popData=None, **kwargs):
    import time
    
    #initialize network_data
    network_data = init_network_data_dict()
    
    #add conv_params to kwargs
    kwargs['conv_params'] = conv_params
    kwargs['mega_params'] = mega_params
    
    #candidate_label = kwargs.get('candidate_label', None)
    print('') #for formatting
    print('Calculating Network Activity Metrics for Simulated Data...')
    start_time = time.time()
    assert simData is not None, 'No simData provided'
    
    # copy simData
    rasterData = simData.copy()
    
    #check if rasterData['spkt'] is empty
    if len(rasterData['spkt']) == 0:
        print('No spike times found in rasterData')
        return None
    
    #convert time to seconds - get initially available data
    spike_times = np.array(rasterData['spkt']) / 1000
    timeVector = np.array(rasterData['t']) / 1000
    spike_times_by_unit = {int(i): spike_times[rasterData['spkid'] == i] for i in np.unique(rasterData['spkid'])} #mea_analysis_pipeline.py expects spike_times as dictionary    
    
    #extract spiking metrics from simulated data
    try: 
        extract_metrics_from_simulated_data(spike_times, timeVector, spike_times_by_unit, rasterData, popData, **kwargs)
    except Exception as e:
        print(f'Error extracting metrics from simulated data: {e}')
        pass
    
    #extract bursting metrics from simulated data (but this one works for both simulated and experimental data)
    try: 
        extract_bursting_activity_data(spike_times, spike_times_by_unit, **kwargs)
    except Exception as e:
        print(f'Error calculating bursting activity: {e}')
        pass
    
    #return network_data
    
    def convert_single_element_arrays(data):
        if isinstance(data, dict):
            for key, value in data.items():
                data[key] = convert_single_element_arrays(value)
        elif isinstance(data, np.ndarray) and data.size == 1:
            return data.item()
        return data

    network_data = convert_single_element_arrays(network_data)    
    print('Network Activity Metrics Calculated and Extracted!')
    print(f'Elapsed time: {time.time() - start_time} seconds')
    print('') #for formatting    
    return network_data

'''Plotting Functions'''
def plot_network_summary(
    network_metrics,
    bursting_plot_path=None,
    bursting_fig_path=None,
    save_path=None,
    mode='2p',    
    ):
    # assertions
    assert save_path is not None, f"Error: save_path must be specified in plot_network_activity"
    
    #init paths
    if bursting_plot_path is None and bursting_fig_path is None:
        
        #HACK: dumb way to to allow for pdf or png save paths - fix later I guess.
        if '.pdf' in save_path: save_path = save_path.replace('.pdf', '.npy') #switch it back to .npy so that it gets handled correctly
        if '.png' in save_path: save_path = save_path.replace('.png', '.npy') 
        
        assert 'network_metrics' in save_path, f"Error: save_path must contain 'network_metrics'" # HACK: this is a hack to make sure we're saving in the right place
        split_path = save_path.split('network_metrics/')
        plot_path = split_path[0] + 'network_plots/' + split_path[1]
        raster_plot_path = plot_path.replace('.npy', '_raster_plot.pdf')
        # raster_fig_path = plot_path.replace('.npy', '_raster_fig.pkl')
        # raster_plot_path_png = plot_path.replace('.npy', '_raster_plot.png')
        bursting_plot_path = plot_path.replace('.npy', '_bursting_plot.pdf')
        bursting_fig_path = plot_path.replace('.npy', '_bursting_fig.pkl')
        bursting_plot_path_png = plot_path.replace('.npy', '_bursting_plot.png')        
    
    # import plotting functions
    from MEA_Analysis.NetworkAnalysis.awNetworkAnalysis.network_analysis import plot_raster_plot_experimental, plot_network_bursting_experimental
    
    # init metrics
    spiking_data_by_unit = network_metrics['spiking_data']['spiking_data_by_unit']
    bursting_ax = network_metrics['bursting_data']['bursting_summary_data']['ax']
    mega_bursting_ax = network_metrics['mega_bursting_data']['bursting_summary_data']['ax']
    
    # choose mode:
    # mode = '2p' - 2 panels
    # mode = '3p' - 3 panels
    if mode == '2p':
        
        # init plot
        fig, ax = plt.subplots(2, 1, figsize=(16, 9))
        
        # plot raster plot
        print("Generating raster plot...")
        #spiking_data_by_unit = network_metrics['spiking_data']['spiking_data_by_unit']
        ax[0] = plot_raster_plot_experimental(ax[0], spiking_data_by_unit)

        # plot network bursting plots
        print("Generating network bursting plot...")
        ax[1] = plot_network_bursting_experimental(ax[1], bursting_ax, mega_ax=mega_bursting_ax)
        #mega_bursting_ax = mega_metrics['bursting_data']['bursting_summary_data']['ax']
        #ax[1] = plot_network_bursting_experimental(ax[1], mega_bursting_ax)
        
        # adjust axes
        print("Adjusting axes...")
        # ensure x-axes are the same for both plots
        # get the narrowes x-axis range shared by both plots
        x_min = min([ax[0].get_xlim()[0], ax[1].get_xlim()[0]])
        x_max = max([ax[0].get_xlim()[1], ax[1].get_xlim()[1]])
        ax[0].set_xlim(x_min, x_max)
        ax[1].set_xlim(x_min, x_max)
        plt.tight_layout()
        
    elif mode == '3p':
        # init plot
        fig, ax = plt.subplots(3, 1, figsize=(16, 9))
        
        # plot raster plot
        print("Generating raster plot...")
        #spiking_data_by_unit = network_metrics['spiking_data']['spiking_data_by_unit']
        ax[0] = plot_raster_plot_experimental(ax[0], spiking_data_by_unit)
        
        # plot network bursting plots
        print("Generating network bursting plot...")
        ax[1] = plot_network_bursting_experimental(ax[1], bursting_ax)
        
        # plot network bursting plots
        print("Generating mega network bursting plot...")
        ax[2] = plot_network_bursting_experimental(ax[2], mega_bursting_ax, mode='mega')
        
        # adjust axes
        print("Adjusting axes...")
        # ensure x-axes are the same for both plots
        # get the narrowes x-axis range shared by both plots
        x_min = min([ax[0].get_xlim()[0], ax[1].get_xlim()[0], ax[2].get_xlim()[0]])
        x_max = max([ax[0].get_xlim()[1], ax[1].get_xlim()[1], ax[2].get_xlim()[1]])
        ax[0].set_xlim(x_min, x_max)
        ax[1].set_xlim(x_min, x_max)
        ax[2].set_xlim(x_min, x_max)
        plt.tight_layout()
        
    # #for debugging, limit x-axis to 35 seconds (in ms)
    try:
        ax[0].set_xlim(0, 35)
        ax[1].set_xlim(0, 35)
        ax[2].set_xlim(0, 35) # just using try except to catch the error if it's a 2 panel plot
    except:
        pass
            
    # make sure the plot path exists
    if not os.path.exists(os.path.dirname(raster_plot_path)):
        os.makedirs(os.path.dirname(raster_plot_path), exist_ok=True)
    
    # save plot as pdf
    fig.savefig(bursting_plot_path) #save as pdf
    print(f"Network summary plot saved to {bursting_plot_path}")
    
    # save plot as png
    fig.savefig(bursting_plot_path_png, dpi=600) #save as png
    print(f"Network summary plot saved to {bursting_plot_path_png}")   

def plot_raster(ax, spiking_data_by_unit, unit_ids=None, E_gids=None, I_gids=None, data_type='simulated'):
        """
        Plot a raster plot of spiking activity.

        Parameters:
            ax (matplotlib.axes.Axes): The matplotlib axis to draw the plot on.
            spiking_data_by_unit (dict): A dictionary where keys are neuron IDs (GIDs) 
                                        and values are spike time data for each unit.
            E_gids (list): List of excitatory neuron IDs.
            I_gids (list): List of inhibitory neuron IDs.

        Returns:
            matplotlib.axes.Axes: The axis with the raster plot.
        """
        
        #main logic
        if data_type == 'simulated':
            list_of_gids = [E_gids, I_gids]
            list_of_colors = ['y.', 'b.']
            list_of_labels = ['Excitatory', 'Inhibitory']
            def plot_mixed(list_of_gids, list_of_colors, list_of_labels):
                assert E_gids is not None, "E_gids must be provided"
                assert I_gids is not None, "I_gids must be provided"

                # Create a dict of all gids and their respective colors and labels
                gid_color_label_dict = {gid: (color, label) for gids, color, label in zip(list_of_gids, list_of_colors, list_of_labels) for gid in gids}
                
                #sort all gids by fr
                sorted_spike_data = {k: v for k, v in sorted(spiking_data_by_unit.items(), key=lambda item: item[1]['FireRate'])}
                sorted_key = 0
                for gid, data in sorted_spike_data.items():
                    spike_times = data['spike_times']
                    #sometimes spike_times is a single value, not a list - make it a list
                    if isinstance(spike_times, (int, float)):
                        spike_times = [spike_times]
                    color, label = gid_color_label_dict[gid]
                    ax.plot(spike_times, [sorted_key] * len(spike_times), f'{color}', markersize=2, label=label)
                    sorted_key += 1
            plot_mixed(list_of_gids, list_of_colors, list_of_labels)
            ax.set_title('Simulated Raster Plot')
            # ax.set_xlabel('Time (s)')
            ax.set_ylabel('Neuron GID')
        elif data_type == 'experimental':            
            color = 'k.'
            label = 'Experimental'
            def plot_exprimental(unit_ids, color, label):           
                #sort
                sorted_spike_data = {k: v for k, v in sorted(spiking_data_by_unit.items(), key=lambda item: item[1]['FireRate'])}
                sorted_key = 0
                for gid, data in sorted_spike_data.items():
                    spike_times = data['spike_times']
                    ax.plot(spike_times, [sorted_key] * len(spike_times), f'{color}', markersize=1, label=label)
                    sorted_key += 1
                
                # #plot
                # for gid in unit_ids:
                #     if gid in spiking_data_by_unit:
                #         spike_times = spiking_data_by_unit[gid]['spike_times']
                #         if isinstance(spike_times, (int, float)):
                #             spike_times = [spike_times]
                #         #ax.plot(spike_times, [gid] * len(spike_times), f'{color}', markersize=1, label=label if gid == gids[0] else None)
                #         ax.plot(spike_times, [gid] * len(spike_times), f'{color}', markersize=1, label=label)
            plot_exprimental(unit_ids, color, label)
            #ax.set_xlabel('Time (s)')
            ax.set_ylabel('Unit ID')
            ax.set_title('Reference Raster Plot')
        else:
            raise ValueError(f"Invalid data_type: {data_type}")
        
        ax.set_title('Raster Plot')
        ax.legend()
        
        # Set the marker size in the legend
        legend = ax.legend()
        for handle in legend.legendHandles:
            #handle._legmarker.set_markersize(5)
            handle._markersize = 5

        # only show unique legend items
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
           
        plt.tight_layout()
        
        # HACK: there must be an easier way to do this...
        # fix x and y lims
        true_max_x = 0
        true_min_x = 0
        for line in ax.lines:
            try:
                x_data = line.get_xdata()
                max_x = max(x_data)
                min_x = min(x_data)
                if max_x > true_max_x:
                    true_max_x = max_x
                if min_x < true_min_x:
                    true_min_x = min_x
            except:
                continue
        ax.set_xlim(true_min_x, true_max_x)
        
        true_max_y = 0
        true_min_y = 0
        for line in ax.lines:
            try:
                y_data = line.get_ydata()
                max_y = max(y_data)
                min_y = min(y_data)
                if max_y > true_max_y:
                    true_max_y = max_y
                if min_y < true_min_y:
                    true_min_y = min_y
            except:
                continue
            ax.set_ylim(true_min_y, true_max_y)       
        
        return ax

def plot_network_summary_slide(sim_data_path, 
                                #convolution_params_path, 
                                #reference_data_npy, 
                                average_fitness,
                                trim_start=0,
                                ):
        start_network_metrics = time()
        #from CDKL5_DIV21.src.conv_params import conv_params
        def get_network_metrics():
            #convolution_params = import_module_from_path(CONVOLUTION_PARAMS)
            #from DIV21.src.conv_params import conv_params
            from RBS_network_models.CDKL5.DIV21.src.conv_params import conv_params
            #from RBS_network_models.developing.CDKL5.DIV21.src.conv_params import conv_params
            #get network metrics
            kwargs = {
                'simData': sim.allSimData,
                'simConfig': sim.cfg,
                #'conv_params': convolution_params.conv_params,
                'conv_params': conv_params,
                'popData': sim.net.allPops,
                'cellData': sim.net.allCells, #not actually used in the fitness calculation, but whatever
            }
            #from DIV21.utils.fitness_helper import calculate_network_metrics
            #from RBS_network_models.developing.utils.fitness_helper import calculate_network_metrics
            from RBS_network_models.utils.fitness_helper import calculate_network_metrics
            error, kwargs = calculate_network_metrics(kwargs)
            #save network metrics
            network_metrics_path = sim_data_path.replace('_data', '_network_metrics')
            if '.pkl' in network_metrics_path:
                network_metrics_path = network_metrics_path.replace('.pkl', '.npy')
            elif '.json' in network_metrics_path:
                network_metrics_path = network_metrics_path.replace('.json', '.npy')
            np.save(network_metrics_path, kwargs)
            # with open(network_metrics_path, 'w') as f:
            #     json.dump(kwargs, f, indent=4)
            return error, kwargs
        error, kwargs = get_network_metrics()
        if error: return error
        # privileged_print("\tNetwork metrics calculated - kwargs dict created.")
        # privileged_print(f'\tTime taken: {time()-start_network_metrics} seconds')
        print("\tNetwork metrics calculated - kwargs dict created.")
        print(f'\tTime taken: {time()-start_network_metrics} seconds')

        #plot raster
        start_raster = time()
        def plot_simulated_raster_wrapper(ax = None, subplot=False, trim_start = None, **kwargs):
            network_metrics = kwargs['network_metrics']
            spiking_data_by_unit = network_metrics['spiking_data']['spiking_data_by_unit']
            popData = sim.net.allPops
            E_gids = popData['E']['cellGids']
            I_gids = popData['I']['cellGids']
            
            if ax is None:
                fig, ax_raster = plt.subplots(1, 1, figsize=(16, 4.5))
            else:
                ax_raster = ax
                subplot = True
                
            ax_raster = plot_raster(ax_raster, spiking_data_by_unit, E_gids=E_gids, I_gids=I_gids, data_type='simulated')
            
            # if trim_start, trim first x seconds from the start of the simulation
            if trim_start is not None and trim_start > 0 and trim_start < ax_raster.get_xlim()[1]:
                ax_raster.set_xlim(trim_start, ax_raster.get_xlim()[1])
            elif trim_start is not None and trim_start > 0 and trim_start > ax_raster.get_xlim()[1]:
                modified_trim = ax_raster.get_xlim()[1]*0.1
                ax_raster.set_xlim(modified_trim, ax_raster.get_xlim()[1])
                print('boop')                
            
            #break if subplot
            if subplot:
                #plt.close()
                return ax_raster

            if DEBUG_MODE:
                #save local for debugging
                #dev_dir = '/pscratch/sd/a/adammwea/workspace/RBS_neuronal_network_models/optimization_projects/CDKL5_DIV21/scripts'
                dev_dir = os.path.dirname(os.path.realpath(__file__))
                plt.savefig(os.path.join(dev_dir, '_raster_plot.png'), dpi=300)
                #save local for debugging
                
            #save wherever data is saved
            #sim_data_path = SIMULATION_RUN_PATH
            raster_plot_path = sim_data_path.replace('_data', '_raster_plot')
            #remove file type and replace with png
            if '.json' in raster_plot_path:
                raster_plot_path = raster_plot_path.replace('.json', '.png')
            elif '.pkl' in raster_plot_path:
                raster_plot_path = raster_plot_path.replace('.pkl', '.png')
            #raster_plot_path = raster_plot_path.replace('.json', '.png')
            plt.savefig(raster_plot_path, dpi=300)
            print(f"Raster plot saved to {raster_plot_path}")
            
            #save as pdf
            raster_plot_path = sim_data_path.replace('_data', '_raster_plot')
            #raster_plot_path = raster_plot_path.replace('.json', '.pdf')
            if '.json' in raster_plot_path:
                raster_plot_path = raster_plot_path.replace('.json', '.pdf')
            elif '.pkl' in raster_plot_path:
                raster_plot_path = raster_plot_path.replace('.pkl', '.pdf')
            plt.savefig(raster_plot_path)
            print(f"Raster plot saved to {raster_plot_path}")
            plt.close()
            
            raster_plots_paths = [
                raster_plot_path,
                raster_plot_path.replace('.pdf', '.png'),
            ]
            
            return raster_plots_paths
        raster_plot_paths = plot_simulated_raster_wrapper(trim_start = trim_start, **kwargs)
        privileged_print("\tIndividual raster plots saved.")
        privileged_print(f'\tTime taken: {time()-start_raster} seconds')

        #plot bursting summary
        start_bursting = time()
        def plot_simulated_bursting_wrapper(ax = None, subplot=False, trim_start = None, **kwargs):
            if ax is None:
                fig, new_ax = plt.subplots(1, 1)
                fig.set_size_inches(16, 4.5)
            else:
                new_ax = ax
                subplot = True #if ax is passed in, then we are plotting on a subplot
            
            #
            conv_params = kwargs['conv_params']
            SpikeTimes = kwargs['network_metrics']['spiking_data']['spiking_times_by_unit']
            #from DIV21.utils.fitness_helper import plot_network_activity_aw
            #from RBS_network_models.developing.utils.analysis_helper import plot_network_activity_aw
            from MEA_Analysis.NetworkAnalysis.awNetworkAnalysis.network_analysis import plot_network_activity_aw
            #bursting_ax, _ = plot_network_activity_aw(new_ax, SpikeTimes, **conv_params) #TODO need to make sure this function agrees with mandar. would be best if we shared a function here.
            bursting_ax = kwargs['network_metrics']['bursting_data']['bursting_summary_data']['ax']
            mega_ax = kwargs['network_metrics']['mega_bursting_data']['bursting_summary_data']['ax']
            
            # # HACK
            # mega_conv_params = kwargs['conv_params'].copy()
            # mega_conv_params['binSize'] *= 5
            # mega_conv_params['gaussianSigma'] *= 15
            #mega_ax, _ = plot_network_activity_aw(new_ax, SpikeTimes, **mega_conv_params) #TODO need to make sure this function agrees with mandar. would be best if we shared a function here.
            
            from MEA_Analysis.NetworkAnalysis.awNetworkAnalysis.network_analysis import plot_network_bursting_experimental
            new_ax = plot_network_bursting_experimental(new_ax, bursting_ax, mega_ax=mega_ax)            
            
            # if trim_start, trim first x seconds from the start of the simulation
            if trim_start is not None and trim_start > 0 and trim_start < new_ax.get_xlim()[1]:
                new_ax.set_xlim(trim_start, new_ax.get_xlim()[1])
            elif trim_start is not None and trim_start > 0 and trim_start > new_ax.get_xlim()[1]:
                modified_trim = new_ax.get_xlim()[1]*0.1
                new_ax.set_xlim(modified_trim, new_ax.get_xlim()[1])    
            
            new_ax.set_title('Bursting summary')
            new_ax.set_xlabel('Time (s)')
            new_ax.set_ylabel('Fire rate (Hz)')
            
            #break if subplot
            if subplot:
                #plt.close()
                return new_ax
            
            plt.tight_layout()
            
            if DEBUG_MODE:
                # Save local for debugging
                dev_dir = os.path.dirname(os.path.realpath(__file__))
                plt.savefig(os.path.join(dev_dir, '_bursting_plot.png'), dpi=300)
                # Save local for debugging
            
            # save wherever data is saved
            #sim_data_path = SIMULATION_RUN_PATH
            bursting_plot_path = sim_data_path.replace('_data', '_bursting_plot')
            #remove file type and replace with png
            #bursting_plot_path = bursting_plot_path.replace('.json', '.png')
            if '.json' in bursting_plot_path:
                bursting_plot_path = bursting_plot_path.replace('.json', '.png')
            elif '.pkl' in bursting_plot_path:
                bursting_plot_path = bursting_plot_path.replace('.pkl', '.png')
            plt.savefig(bursting_plot_path, dpi=300)
            print(f"Bursting plot saved to {bursting_plot_path}")
            
            #save as pdf
            bursting_plot_path = sim_data_path.replace('_data', '_bursting_plot')
            #bursting_plot_path = bursting_plot_path.replace('.json', '.pdf')
            if '.json' in bursting_plot_path:
                bursting_plot_path = bursting_plot_path.replace('.json', '.pdf')
            elif '.pkl' in bursting_plot_path:
                bursting_plot_path = bursting_plot_path.replace('.pkl', '.pdf')
            plt.savefig(bursting_plot_path)
            print(f"Bursting plot saved to {bursting_plot_path}")
            plt.close()
            
            bursting_plot_paths = [
                bursting_plot_path,
                bursting_plot_path.replace('.pdf', '.png'),
            ]
            
            return bursting_plot_paths    
        bursting_plot_paths = plot_simulated_bursting_wrapper(trim_start = trim_start, **kwargs)
        privileged_print("\tIndividual bursting plots saved.")
        privileged_print(f'\tTime taken: {time()-start_bursting} seconds')

        # combine plots into a single summary plot
        start_summary = time()
        def plot_simulation_summary(trim_start = None, **kwargs):
            fig, ax = plt.subplots(2, 1, figsize=(16, 9))
            raster_plot_ax, bursting_plot_ax = ax
            subplot = True
            raster_plot_ax = plot_simulated_raster_wrapper(ax=raster_plot_ax, subplot=subplot, trim_start = trim_start, **kwargs)
            bursting_plot_ax = plot_simulated_bursting_wrapper(ax=bursting_plot_ax, subplot=subplot, trim_start = trim_start, **kwargs)
            
            #make both plots share the same x-axis
            raster_plot_ax.get_shared_x_axes().join(raster_plot_ax, bursting_plot_ax)    
            plt.tight_layout()
            
            if DEBUG_MODE:
                # Save local for debugging
                dev_dir = os.path.dirname(os.path.realpath(__file__))
                plt.savefig(os.path.join(dev_dir, '_summary_plot.png'), dpi=300)
                # Save local for debugging
                
            # save wherever data is saved
            #sim_data_path = SIMULATION_RUN_PATH
            summary_plot_path = sim_data_path.replace('_data', '_summary_plot')
            #remove file type and replace with png
            #summary_plot_path = summary_plot_path.replace('.json', '.png')
            if '.json' in summary_plot_path:
                summary_plot_path = summary_plot_path.replace('.json', '.png')
            elif '.pkl' in summary_plot_path:
                summary_plot_path = summary_plot_path.replace('.pkl', '.png')
            plt.savefig(summary_plot_path, dpi=300)
            print(f"Summary plot saved to {summary_plot_path}")
            
            #save as pdf
            summary_plot_path = sim_data_path.replace('_data', '_summary_plot')
            #summary_plot_path = summary_plot_path.replace('.json', '.pdf')
            if '.json' in summary_plot_path:
                summary_plot_path = summary_plot_path.replace('.json', '.pdf')
            elif '.pkl' in summary_plot_path:
                summary_plot_path = summary_plot_path.replace('.pkl', '.pdf')
            plt.savefig(summary_plot_path)  
            print(f"Summary plot saved to {summary_plot_path}")
            plt.close()
            
            summary_plot_paths = [
                summary_plot_path,
                summary_plot_path.replace('.pdf', '.png'),
            ]    
            return summary_plot_paths
        plot_simulation_summary(trim_start=trim_start, **kwargs)
        privileged_print("\tSimulation summary plot saved.")
        privileged_print(f'\tTime taken: {time()-start_summary} seconds')

        # comparison plot against reference data snippet
        start_comparison = time()
        activate_print()
        def plot_comparision_plot(ax_list=None, subplot=False, trim_start = None, **kwargs):
            #activate_print() #for debugging
            #fig, ax = plt.subplots(4, 1, figsize=(16, 9))
            if ax_list is None:
                fig, ax_list = plt.subplots(4, 1, figsize=(16, 9))
                sim_raster_ax, sim_bursting_ax, ref_raster_ax, ref_bursting_ax = ax_list
            else:
                #assert that there be 4 axes in the ax list
                assert len(ax_list) == 4, "There must be 4 axes in the ax_list."
                sim_raster_ax, sim_bursting_ax, ref_raster_ax, ref_bursting_ax = ax_list
                subplot = True
            
            #plot simulated raster
            sim_raster_ax = plot_simulated_raster_wrapper(ax=sim_raster_ax, subplot=True, trim_start=trim_start, **kwargs)
            
            #plot simulated bursting
            sim_bursting_ax = plot_simulated_bursting_wrapper(ax=sim_bursting_ax, subplot=True, trim_start=trim_start, **kwargs)
            
            # plot reference raster
            def plot_reference_raster_wrapper(ax=None, subplot=False, sim_data_length=None, **kwargs):
                #load npy ref data
                print(ax)
                print(subplot)
                #print(kwargs)
                
                #load npy ref data
                ref_data = np.load(
                    REFERENCE_DATA_NPY, 
                    allow_pickle=True
                    ).item()
                network_metrics = ref_data
                spiking_data_by_unit = network_metrics['spiking_data']['spiking_data_by_unit'].copy()
                unit_ids = [unit_id for unit_id in spiking_data_by_unit.keys()]
                
                # Trim reference data to match the length of simulated data
                if sim_data_length is not None:
                    for unit_id in unit_ids:
                        #spiking_data_by_unit[unit_id] = spiking_data_by_unit[unit_id][:sim_data_length]
                        #max_spike_time_index = np.argmax(spiking_data_by_unit[unit_id])
                        spike_times = spiking_data_by_unit[unit_id]['spike_times']
                        spike_times = spike_times[spike_times < sim_data_length]
                        spiking_data_by_unit[unit_id]['spike_times'] = spike_times
                
                if ax is None:
                    fig, ax_raster = plt.subplots(1, 1, figsize=(16, 4.5))
                else:
                    ax_raster = ax
                    subplot = True
                    
                ax_raster = plot_raster(ax_raster, spiking_data_by_unit, unit_ids=unit_ids, data_type='experimental')
                
                if subplot:
                    #plt.close()
                    print("subplot")
                    return ax_raster
                
                if DEBUG_MODE:
                    #save local for debugging
                    #dev_dir = '/pscratch/sd/a/adammwea/workspace/RBS_neuronal_network_models/optimization_projects/CDKL5_DIV21/scripts'
                    dev_dir = os.path.dirname(os.path.realpath(__file__))
                    plt.savefig(os.path.join(dev_dir, '_ref_raster_plot.png'), dpi=300)
                    #save as pdf
                    plt.savefig(os.path.join(dev_dir, '_ref_raster_plot.pdf'))
                    #save local for debugging
                
                
                #TODO: semi unfinished compared to similar function for simulated raster
            #sim_data_length = len(kwargs['network_metrics']['spiking_data']['spiking_data_by_unit'][0])
            sim_data_length = kwargs['simConfig']['duration']/1000 #in seconds
            ref_raster_ax = plot_reference_raster_wrapper(ax=ref_raster_ax, subplot=True, sim_data_length = sim_data_length, **kwargs)
            
            #plot reference bursting
            def plot_reference_bursting_wrapper(ax=None, subplot=False, sim_data_length=None, trim_start = None, **kwargs):
                #ref_data = np.load(REFERENCE_DATA_NPY, allow_pickle=True).item()
                ref_data = np.load(REFERENCE_DATA_NPY, allow_pickle=True).item()
                network_metrics = ref_data
                conv_params = kwargs['conv_params']
                SpikeTimes = network_metrics['spiking_data']['spiking_times_by_unit']
                if ax is None:
                    fig, new_ax = plt.subplots(1, 1)
                    fig.set_size_inches(16, 4.5)
                else:
                    new_ax = ax
                    subplot = True
                    
                if sim_data_length is not None:
                    for unit_id in SpikeTimes.keys():
                        spike_times = SpikeTimes[unit_id]
                        spike_times = spike_times[spike_times < sim_data_length]
                        SpikeTimes[unit_id] = spike_times
                
                #from DIV21.utils.fitness_helper import plot_network_activity_aw
                #from RBS_network_models.developing.utils.analysis_helper import plot_network_activity_aw
                from MEA_Analysis.NetworkAnalysis.awNetworkAnalysis.network_analysis import plot_network_activity_aw
                new_ax, _ = plot_network_activity_aw(new_ax, SpikeTimes, **conv_params) #TODO need to make sure this function agrees with mandar. would be best if we shared a function here.
                
                new_ax.set_title('Bursting summary')
                new_ax.set_xlabel('Time (s)')
                new_ax.set_ylabel('Fire rate (Hz)')
                
                #break if subplot
                if subplot:
                    #plt.close()
                    return new_ax
                
                plt.tight_layout()
                
                if DEBUG_MODE:
                    # Save local for debugging
                    dev_dir = os.path.dirname(os.path.realpath(__file__))
                    plt.savefig(os.path.join(dev_dir, '_ref_bursting_plot.png'), dpi=300)
                    # plot as pdf
                    plt.savefig(os.path.join(dev_dir, '_ref_bursting_plot.pdf'))
                    # Save local for debugging
                
                # TODO: semi unfinished compared to similar function for simulated bursting
            sim_data_length = kwargs['simConfig']['duration']/1000 #in seconds
            ref_bursting_ax = plot_reference_bursting_wrapper(ax=ref_bursting_ax, subplot=True, sim_data_length=sim_data_length, trim_start = trim_start, **kwargs)
            
            # # remove the first second of x-axis for all plots
            # def set_xlim_to_first_value_over_one(ax):
            #     x_data = ax.lines[0].get_xdata()
            #     first_value_over_one = next((x for x in x_data if x > 1), None)
            #     if first_value_over_one is not None:
            #         ax.set_xlim(first_value_over_one, ax.get_xlim()[1])
            
            # set_xlim_to_first_value_over_one(sim_raster_ax)
            # set_xlim_to_first_value_over_one(ref_raster_ax)
            # set_xlim_to_first_value_over_one(sim_bursting_ax)
            # set_xlim_to_first_value_over_one(ref_bursting_ax)
            
            # #print axes limits
            # print(sim_raster_ax.get_xlim())
            # print(ref_raster_ax.get_xlim())
            # print(sim_bursting_ax.get_xlim())
            # print(ref_bursting_ax.get_xlim())
            #ref_bursting_ax.append(ref_bursting_ax.get_xlim()[0] + 1)
            
            # # Ensure y-axis is the same for bursting plots
            # sim_bursting_ylim = sim_bursting_ax.get_ylim()
            # ref_bursting_ylim = ref_bursting_ax.get_ylim()
            
            # ensure xaxis of refernce plots matches simulated raster
            ref_raster_ax.set_xlim(sim_raster_ax.get_xlim())
            ref_bursting_ax.set_xlim(sim_raster_ax.get_xlim())
            sim_bursting_ax.set_xlim(ref_bursting_ax.get_xlim())
            
            # Ensure y-axis is the same for bursting plots
            def adjust_ylim_based_on_xlim(ax):
                x_data = ax.lines[0].get_xdata()
                y_data = ax.lines[0].get_ydata()
                xlim = ax.get_xlim()
                filtered_y_data = [y for x, y in zip(x_data, y_data) if xlim[0] <= x <= xlim[1]]
                if filtered_y_data:
                    ax.set_ylim(min(filtered_y_data) * 0.8, max(filtered_y_data) * 1.2)

            adjust_ylim_based_on_xlim(sim_bursting_ax)
            adjust_ylim_based_on_xlim(ref_bursting_ax)  
            
            
            # Calculate the combined y-axis limits
            sim_bursting_ylim = sim_bursting_ax.get_ylim()
            ref_bursting_ylim = ref_bursting_ax.get_ylim()
            # Calculate the combined y-axis limits
            combined_ylim = [
                min(sim_bursting_ylim[0], ref_bursting_ylim[0]) * 0.8,
                max(sim_bursting_ylim[1], ref_bursting_ylim[1]) * 1.2
            ]
            
            # Set the same y-axis limits for both axes
            sim_bursting_ax.set_ylim(combined_ylim)
            ref_bursting_ax.set_ylim(combined_ylim)
            
            # Make all four plots share the same x-axis as simulated raster
            sim_raster_ax.get_shared_x_axes().join(sim_raster_ax, sim_bursting_ax)
            sim_raster_ax.get_shared_x_axes().join(sim_raster_ax, ref_raster_ax)
            sim_raster_ax.get_shared_x_axes().join(sim_raster_ax, ref_bursting_ax)  

            #
            if subplot:
                #plt.tight_layout()
                #plt.close()
                return [sim_raster_ax, sim_bursting_ax, ref_raster_ax, ref_bursting_ax] 
            
            plt.tight_layout()
            
            if DEBUG_MODE:
                # Save local for debugging
                dev_dir = os.path.dirname(os.path.realpath(__file__))
                fig.savefig(os.path.join(dev_dir, '_comparison_plot.png'), dpi=300)
                # Save local for debugging
            
            # save wherever data is saved
            #sim_data_path = SIMULATION_RUN_PATH
            comparison_plot_path = sim_data_path.replace('_data', '_comparison_plot')
            #remove file type and replace with png
            #comparison_plot_path = comparison_plot_path.replace('.json', '.png')
            if '.json' in comparison_plot_path:
                comparison_plot_path = comparison_plot_path.replace('.json', '.png')
            elif '.pkl' in comparison_plot_path:
                comparison_plot_path = comparison_plot_path.replace('.pkl', '.png')
            fig.savefig(comparison_plot_path, dpi=300)
            print(f"Comparison plot saved to {comparison_plot_path}")
            
            #save as pdf
            comparison_plot_path = sim_data_path.replace('_data', '_comparison_plot')
            #comparison_plot_path = comparison_plot_path.replace('.json', '.pdf')
            if '.json' in comparison_plot_path:
                comparison_plot_path = comparison_plot_path.replace('.json', '.pdf')
            elif '.pkl' in comparison_plot_path:
                comparison_plot_path = comparison_plot_path.replace('.pkl', '.pdf')
            fig.savefig(comparison_plot_path)
            print(f"Comparison plot saved to {comparison_plot_path}")
            plt.close()
        plot_comparision_plot(trim_start=trim_start, **kwargs)
        privileged_print("\tComparison plot saved.")
        privileged_print(f'\tTime taken: {time()-start_comparison} seconds')

        # build comparison summary slide
        start_summary_slide = time()
        def build_comparision_summary_slide(sim_data_path, trim_start = None,):
            fig, ax_list = plt.subplots(4, 1, figsize=(16, 9))
            
            #plot comparison plot
            plot_comparision_plot(ax_list=ax_list, subplot=True, trim_start=trim_start, **kwargs)
            
            #remove x-axis labels from all plots
            for ax in ax_list:
                ax.set_xlabel('')
                
            #remove titles from all plots
            for ax in ax_list:
                ax.set_title('')
                
            #remove x-axis ticks from all plots except the bottom one
            for ax in ax_list[:-1]:
                ax.set_xticks([])
                
            # add 'simulated' and 'reference' labels above each raster plot
            ax_list[0].text(0.5, 1.1, 'Simulated', ha='center', va='center', transform=ax_list[0].transAxes, fontsize=12)
            ax_list[2].text(0.5, 1.1, 'Reference', ha='center', va='center', transform=ax_list[2].transAxes, fontsize=12)
            
            # add time(s) label to the bottom plot
            ax_list[-1].set_xlabel('Time (s)')
                    
            #create space to the right of plots for text
            fig.subplots_adjust(right=0.6)
            
            #create some space at the bottom for one line of text
            fig.subplots_adjust(bottom=0.1)
            
            #add text
            spiking_summary = kwargs['network_metrics']['spiking_data']['spiking_summary_data']
            bursting_summary = kwargs['network_metrics']['bursting_data']['bursting_summary_data']
            simulated_spiking_data = kwargs['network_metrics']['simulated_data']
            #simulated_bursting_data = kwargs['network_metrics']['simulated_data']
            
            #TODO: if meanBurstrate is not in the spiking summary, then the following will fail, add nan in that case for now
            if 'mean_Burst_Rate' not in bursting_summary:
                bursting_summary['mean_Burst_Rate'] = np.nan                
            
            text = (
                # f"Summary of comparison between simulated and reference data:\n"
                # f"Simulated raster plot is shown in the top left, reference raster plot is shown in the bottom left.\n"
                # f"Simulated bursting plot is shown in the top right, reference bursting plot is shown in the bottom right.\n"
                f"Simulated spiking summary:\n"
                f"  - Mean firing rate: {spiking_summary['MeanFireRate']} Hz\n"
                f"  - Coefficient of variation: {spiking_summary['CoVFireRate']}\n"
                f"  - Mean ISI: {spiking_summary['MeanISI']} s\n"
                f"  - Coefficient of variation of ISI: {spiking_summary['CoV_ISI']}\n"
                f"  - Mean firing rate E: {simulated_spiking_data['MeanFireRate_E']} Hz\n"
                f"  - Mean CoV FR E: {simulated_spiking_data['CoVFireRate_E']}\n"
                f"  - Mean firing rate I: {simulated_spiking_data['MeanFireRate_I']} Hz\n"
                f"  - Mean CoV FR I: {simulated_spiking_data['CoVFireRate_I']}\n"
                f"  - Mean ISI E: {simulated_spiking_data['MeanISI_E']} s\n"
                f"  - CoV ISI E: {simulated_spiking_data['CoV_ISI_E']}\n"
                f"  - Mean ISI I: {simulated_spiking_data['MeanISI_I']} s\n"
                f"  - CoV ISI I: {simulated_spiking_data['CoV_ISI_I']}\n"
                
                f"Simulated bursting summary:\n"        
                #mean burst rate
                #number of units
                f"  - Number of units: {bursting_summary['NumUnits']}\n"
                #baseline
                f"  - Baseline: {bursting_summary['baseline']} Hz\n"
                #fanofactor
                f"  - Fanofactor: {bursting_summary['fano_factor']}\n"
                #num bursts
                f"  - Number of bursts: {bursting_summary['Number_Bursts']}\n"
                #mean burst rate
                f"  - Mean burst rate: {bursting_summary['mean_Burst_Rate']} bursts/second\n"     
                #mean burst peak
                f"  - Mean burst peak: {bursting_summary['mean_Burst_Peak']} Hz\n"
                #burst peak CoV
                f"  - CoV burst peak: {bursting_summary['cov_Burst_Peak']}\n"
                #mean IBI
                f"  - Mean IBI: {bursting_summary['mean_IBI']} s\n"
                #CoV IBI
                f"  - CoV IBI: {bursting_summary['cov_IBI']}\n"
                #mean within burst isi
                f"  - Mean within burst ISI: {bursting_summary['MeanWithinBurstISI']} s\n"
                #CoV within burst isi
                f"  - CoV within burst ISI: {bursting_summary['CoVWithinBurstISI']}\n"
                #mean outside burst isi
                f"  - Mean outside burst ISI: {bursting_summary['MeanOutsideBurstISI']} s\n"
                #CoV outside burst isi
                f"  - CoV outside burst ISI: {bursting_summary['CoVOutsideBurstISI']}\n"
                #mean whole network isi
                f"  - Mean network ISI: {bursting_summary['MeanNetworkISI']} s\n"
                #CoV whole network isi
                f"  - CoV network ISI: {bursting_summary['CoVNetworkISI']}\n"
                )
            
            #append reference metrics to text
            ref_data = np.load(REFERENCE_DATA_NPY, allow_pickle=True).item()
            ref_spiking_summary = ref_data['spiking_data']['spiking_summary_data']
            ref_bursting_summary = ref_data['bursting_data']['bursting_summary_data']    
            text += (
                f"\n"
                f"\n"
                f"Reference spiking summary:\n"
                f"  - Mean firing rate: {ref_spiking_summary['MeanFireRate']} Hz\n"
                f"  - Coefficient of variation: {ref_spiking_summary['CoVFireRate']}\n"
                f"  - Mean ISI: {ref_spiking_summary['MeanISI']} s\n"
                f"  - Coefficient of variation of ISI: {ref_spiking_summary['CoV_ISI']}\n"
                
                f"Reference bursting summary:\n"
                #number of units
                f"  - Number of units: {ref_bursting_summary['NumUnits']}\n"
                #baseline
                f"  - Baseline: {ref_bursting_summary['baseline']} Hz\n"
                #fanofactor
                f"  - Fanofactor: {ref_bursting_summary['fano_factor']}\n"
                #num bursts
                f"  - Number of bursts: {ref_bursting_summary['Number_Bursts']}\n"
                #mean burst rate
                #f"  - Mean burst rate: {ref_bursting_summary['mean_Burst_Rate']} Hz\n" #TODO: this is not in the data        
                #mean burst peak
                f"  - Mean burst peak: {ref_bursting_summary['mean_Burst_Peak']} Hz\n"
                #burst peak CoV
                f"  - CoV burst peak: {ref_bursting_summary['cov_Burst_Peak']}\n"
                #mean IBI
                f"  - Mean IBI: {ref_bursting_summary['mean_IBI']} s\n"
                #CoV IBI
                f"  - CoV IBI: {ref_bursting_summary['cov_IBI']}\n"
                #mean within burst isi
                f"  - Mean within burst ISI: {ref_bursting_summary['MeanWithinBurstISI']} s\n"
                #CoV within burst isi
                f"  - CoV within burst ISI: {ref_bursting_summary['CoVWithinBurstISI']}\n"
                #mean outside burst isi
                f"  - Mean outside burst ISI: {ref_bursting_summary['MeanOutsideBurstISI']} s\n"
                #CoV outside burst isi
                f"  - CoV outside burst ISI: {ref_bursting_summary['CoVOutsideBurstISI']}\n"
                #mean whole network isi
                f"  - Mean network ISI: {ref_bursting_summary['MeanNetworkISI']} s\n"
                #CoV whole network isi
                f"  - CoV network ISI: {ref_bursting_summary['CoVNetworkISI']}\n"
                )
            
            #add average fitness to text
            text += (
                f"\n"
                f"\n"
                f"\nAverage fitness: {average_fitness}"
                ) 
            
            #add text to the right of the plots
            fig.text(0.65, 0.5, text, ha='left', va='center', fontsize=9)
            
            #add data path at the bottom of the slide
            #fig.text(0.5, 0.05, f"Data path: {SIMULATION_RUN_PATH}", ha='center', va='center', fontsize=10)
            #go a little lower and add data path
            fig.text(0.5, 0.02, f"Data path: {sim_data_path}", ha='center', va='center', fontsize=9) 
            
            if DEBUG_MODE:
                #save local for debugging
                dev_dir = os.path.dirname(os.path.realpath(__file__))
                fig.savefig(os.path.join(dev_dir, '_comparison_summary_slide.png'), dpi=300)
                #save as pdf
            
            # save wherever data is saved
            sim_data_path = sim_data_path
            comparison_summary_slide_path = sim_data_path.replace('_data', '_comparison_summary_slide')
            #remove file type and replace with png
            #comparison_summary_slide_path = comparison_summary_slide_path.replace('.json', '.png')
            if '.json' in comparison_summary_slide_path:
                comparison_summary_slide_path = comparison_summary_slide_path.replace('.json', '.png')
            elif '.pkl' in comparison_summary_slide_path:
                comparison_summary_slide_path = comparison_summary_slide_path.replace('.pkl', '.png')
            fig.savefig(comparison_summary_slide_path, dpi=300)
            privileged_print(f"Comparison summary slide saved to {comparison_summary_slide_path}")
            
            #save as pdf
            comparison_summary_slide_path = sim_data_path.replace('_data', '_comparison_summary_slide')
            #comparison_summary_slide_path = comparison_summary_slide_path.replace('.json', '.pdf')
            if '.json' in comparison_summary_slide_path:
                comparison_summary_slide_path = comparison_summary_slide_path.replace('.json', '.pdf')
            elif '.pkl' in comparison_summary_slide_path:
                comparison_summary_slide_path = comparison_summary_slide_path.replace('.pkl', '.pdf')
            fig.savefig(comparison_summary_slide_path)
            privileged_print(f"Comparison summary slide saved to {comparison_summary_slide_path}")
            plt.close()
        build_comparision_summary_slide(sim_data_path, trim_start = trim_start)
        privileged_print("\tComparison summary slide saved.")
        privileged_print(f'\tTime taken: {time()-start_summary_slide} seconds')

def plot_raster_plot_experimental(ax, spiking_data_by_unit):
    """Plot a raster plot for spiking data."""
    
    # Calculate the average firing rate for each unit
    firing_rates = {}
    for gid in spiking_data_by_unit:
        spike_times = spiking_data_by_unit[gid]['spike_times']
        spike_times = [spike_times] if isinstance(spike_times, (int, float)) else spike_times
        firing_rate = len(spike_times) / (max(spike_times) - min(spike_times)) if len(spike_times) > 1 else 0
        firing_rates[gid] = firing_rate
    
    # Sort the units based on their average firing rates
    sorted_units = sorted(firing_rates, key=firing_rates.get)
    
    # Create a mapping from original gid to new y-axis position
    gid_to_ypos = {gid: pos for pos, gid in enumerate(sorted_units)}
    
    # Plot the units in the sorted order
    for gid in sorted_units:
        spike_times = spiking_data_by_unit[gid]['spike_times']
        spike_times = [spike_times] if isinstance(spike_times, (int, float)) else spike_times
        ax.plot(spike_times, [gid_to_ypos[gid]] * len(spike_times), 'b.', markersize=2)

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Unit ID (sorted by firing rate)')
    ax.set_title('Raster Plot')
    plt.tight_layout()
    return ax

def plot_network_bursting_experimental(ax, bursting_ax, mega_ax=None, mode='normal'):
    """Plot network bursting activity with clear differentiation between lines."""
    
    if mode == 'normal':
        # Copy ax features to the new ax
        ax.set_xlim(bursting_ax.get_xlim())
        ax.set_ylim(bursting_ax.get_ylim())
        ax.set_ylabel('Firing Rate (Hz)')
        ax.set_xlabel('Time (s)')
        ax.set_title(bursting_ax.get_title())

        # Plot Line 1 (blue)
        ax.plot(
            bursting_ax.get_lines()[0].get_xdata(),
            bursting_ax.get_lines()[0].get_ydata(),
            color='blue',
            #label='Bursting'
        )
        
        # plot bursting peaks
        ax.plot(
            bursting_ax.get_lines()[1].get_xdata(),
            bursting_ax.get_lines()[1].get_ydata(),
            'o', color='red', 
            label='Bursts'
        )

        # # Mark peaks with vertical dotted lines
        # for x_peak in bursting_ax.get_lines()[1].get_xdata():
        #     ax.axvline(x=x_peak, linestyle='dotted', color='red', alpha=0.8)

        # If mega_ax is not None, plot filled area and peaks
        if mega_ax is not None:
            # Plot the filled area behind other lines
            ax.fill_between(
                mega_ax.get_lines()[0].get_xdata(),
                mega_ax.get_lines()[0].get_ydata(),
                color='red',
                alpha=0.2,  # Transparency for the filled area
                label='Mega Bursts'
            )
            
            #outline the filled aread, dotted line
            ax.plot(
                mega_ax.get_lines()[0].get_xdata(),
                mega_ax.get_lines()[0].get_ydata(),
                color='red',
                linestyle='dotted'
            )

            # Mark mega peaks with vertical dotted lines
            for x_peak in mega_ax.get_lines()[1].get_xdata():
                ax.axvline(x=x_peak, linestyle='dotted', 
                        color='red', 
                        #label='Mega Bursts'
                        #alpha=0.8
                        )
                
            # remove redundant legend items
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys())

        # Add a legend for clarity
        ax.legend()

        return ax
    elif mode == 'mega':
        # this is for the case where we're plotting mega bursts only
        # Copy ax features to the new ax
        ax.set_xlim(bursting_ax.get_xlim())
        ax.set_ylim(bursting_ax.get_ylim())
        ax.set_ylabel('Firing Rate (Hz)')
        ax.set_xlabel('Time (s)')
        ax.set_title(bursting_ax.get_title())
        
        # Plot Line 1 (red)
        ax.plot(
            bursting_ax.get_lines()[0].get_xdata(),
            bursting_ax.get_lines()[0].get_ydata(),
            color='red',
            #label='Mega Bursts'
        )
        
        # plot bursting peaks
        ax.plot(
            bursting_ax.get_lines()[1].get_xdata(),
            bursting_ax.get_lines()[1].get_ydata(),
            'o', 
            color='orange',
            label='Mega Bursts'
        )
        
        # Mark peaks with vertical dotted lines
        for x_peak in bursting_ax.get_lines()[1].get_xdata():
            ax.axvline(x=x_peak, linestyle='dotted', color='red', alpha=0.8)
            
        # Add a legend for clarity
        ax.legend()
        
        return ax

'''Experimental Data Functions'''
def analyze_bursting_activity(spike_times, spike_times_by_unit, isi_threshold):
    '''Slightly modified version of the code from mea_analysis_pipeline.py'''
    #spike_times = {int(i): rasterData['spkt'][rasterData['spkid'] == i] for i in np.unique(rasterData['spkid'])}
    
    # convert spike_times to dictionary
    #spike_times_dict = {int(i): spike_times[spike_times['spkid'] == i] for i in np.unique(spike_times['spkid'])}
    #spike_times = spike_times_dict
    spike_times = spike_times_by_unit
    
    
    #burst_statistics = helper.detect_bursts_statistics(spike_times, isi_threshold=isi_threshold)
    burst_statistics = detect_bursts_statistics(spike_times, isi_threshold=isi_threshold)
    bursts_by_unit = [unit_stats['bursts'] for unit_stats in burst_statistics.values()]
    
    all_isis_within_bursts = np.concatenate([stats['isis_within_bursts'] for stats in burst_statistics.values() if stats['isis_within_bursts'].size > 0])
    all_isis_outside_bursts = np.concatenate([stats['isis_outside_bursts'] for stats in burst_statistics.values() if stats['isis_outside_bursts'].size > 0])
    all_isis = np.concatenate([stats['isis_all'] for stats in burst_statistics.values() if stats['isis_all'].size > 0])

    mean_isi_within_combined = np.mean(all_isis_within_bursts) if all_isis_within_bursts.size > 0 else np.nan
    cov_isi_within_combined = np.cov(all_isis_within_bursts) if all_isis_within_bursts.size > 0 else np.nan
    #fano_factor_within_combined = np.concatenate([[stats['fano_factor_within']] for stats in burst_statistics.values() if stats['fano_factor_within'] is not None])

    mean_isi_outside_combined = np.mean(all_isis_outside_bursts) if all_isis_outside_bursts.size > 0 else np.nan
    cov_isi_outside_combined = np.cov(all_isis_outside_bursts) if all_isis_outside_bursts.size > 0 else np.nan
    #fano_factor_outside_combined = np.concatenate([[stats['fano_factor_outside']] for stats in burst_statistics.values() if stats['fano_factor_outside'] is not None])

    mean_isi_all_combined = np.mean(all_isis) if all_isis.size > 0 else np.nan
    cov_isi_all_combined = np.cov(all_isis) if all_isis.size > 0 else np.nan    
    #fano_factor_all_combined = np.concatenate([[stats['fano_factor_all']] for stats in burst_statistics.values() if stats['fano_factor'] is not None])

    bursting_summary_data = {}
    bursting_summary_data['MeanWithinBurstISI'] = mean_isi_within_combined
    bursting_summary_data['CoVWithinBurstISI'] = cov_isi_within_combined
    bursting_summary_data['MeanOutsideBurstISI'] = mean_isi_outside_combined   
    bursting_summary_data['CoVOutsideBurstISI'] = cov_isi_outside_combined
    bursting_summary_data['MeanNetworkISI'] = mean_isi_all_combined
    bursting_summary_data['CoVNetworkISI'] = cov_isi_all_combined
    #bursting_summary_data['fano_factor_within'] = fano_factor_within_combined
    #bursting_summary_data['fano_factor_outside'] = fano_factor_outside_combined
    #bursting_summary_data['fano_factor_all'] = fano_factor_all_combined
    bursting_summary_data['NumUnits'] = len(spike_times)
    bursting_summary_data['Number_Bursts'] = sum(len(unit_stats['bursts']) for unit_stats in burst_statistics.values())
    bursting_summary_data['mean_IBI'] = np.mean(all_isis) if all_isis.size > 0 else np.nan
    bursting_summary_data['cov_IBI'] = np.cov(all_isis) if all_isis.size > 0 else np.nan
    
    bursting_data = {
        'bursting_summary_data': bursting_summary_data,
        'bursting_data_by_unit': burst_statistics,
    }

    return bursting_data

def plot_network_activity_aw(ax,SpikeTimes, min_peak_distance=1.0, 
                             binSize=0.1, 
                             gaussianSigma=0.16, 
                             thresholdBurst=1.2, 
                             prominence=1, 
                             figSize=(10, 6),
                             title='Network Activity'
                             ):
    # TODO: stole this code from MEA_pipline ipn analysis. I'm going to make some changes for myself. 
            # I need to put them back into the MEA pipeline at somepoint
            
    # # aw 2025-01-16 13:26:25 - its been a while since I made teh note above. I'm pretty much ready to put changes
    # back into the MEA pipeline.
    
    relativeSpikeTimes = []
    units = 0
    for unit_id, spike_times in SpikeTimes.items():
        temp_spike_times = spike_times
        #sometimes the spike times are in a arrays, somtimes they are floats       
        if isinstance(temp_spike_times, np.ndarray): relativeSpikeTimes.extend(temp_spike_times) 
        elif isinstance(temp_spike_times, float): relativeSpikeTimes.append(temp_spike_times)
        elif isinstance(temp_spike_times, int): relativeSpikeTimes.append(temp_spike_times)
        else:
            Warning(f'unit {unit_id} has spike times that are not float or array')
            continue
        units += 1 # Set the first spike time to 0
    assert all(isinstance(x, (int, float)) for x in relativeSpikeTimes), 'All elements in relativeSpikeTimes must be ints or floats'
    
    relativeSpikeTimes = np.array(relativeSpikeTimes)
    relativeSpikeTimes.sort() # Convert to NumPy array

    # Step 1: Bin all spike times into small time windows
    timeVector = np.arange(min(relativeSpikeTimes), max(relativeSpikeTimes), binSize)  # Time vector for binning
    binnedTimes, _ = np.histogram(relativeSpikeTimes, bins=timeVector)  # Bin spike times
    binnedTimes = np.append(binnedTimes, 0)  # Append 0 to match MATLAB's binnedTimes length

    # Step 2: Smooth the binned spike times with a Gaussian kernel
    kernelRange = np.arange(-3*gaussianSigma, 3*gaussianSigma + binSize, binSize)  # Range for Gaussian kernel
    kernel = norm.pdf(kernelRange, 0, gaussianSigma)  # Gaussian kernel
    kernel *= binSize  # Normalize kernel by bin size
    firingRate = convolve(binnedTimes, kernel, mode='same') / binSize  # Convolve and normalize by bin size
    firingRate = firingRate / units  # Convert to Hz

    # Create a new figure with a specified size (width, height)
    #fig, ax = plt.subplots(figsize=figSize)

    # Plot the smoothed network activity
    ax.plot(timeVector, firingRate, color='royalblue')
    # Restrict the plot to the first and last 100 ms
    ax.set_xlim([min(relativeSpikeTimes), max(relativeSpikeTimes)])
    ax.set_ylim([min(firingRate)*0.8, max(firingRate)*1.2])  # Set y-axis limits to min and max of firingRate
    ax.set_ylabel('Firing Rate [Hz]')
    ax.set_xlabel('Time [ms]')
    #ax.set_title('Network Activity', fontsize=11)
    ax.set_title(title, fontsize=11)

    # Step 3: Peak detection on the smoothed firing rate curve
    rmsFiringRate = np.sqrt(np.mean(firingRate**2))  # Calculate RMS of the firing rate

    peaks, properties = find_peaks(firingRate, prominence=prominence, distance=min_peak_distance)  # Find peaks above the threshold
    print(properties.keys())
    burstPeakTimes = timeVector[peaks]  # Convert peak indices to times
    #burstPeakValues = properties['prominences']  # Get the peak values
    burstPeakValues = firingRate[peaks]  # Get the peak values


    #Calculate the ISIs between spiketimes
    # Calculate the intervals between consecutive peaks
    intervals = np.diff(burstPeakTimes)

    # Calculate the mean interval between peaks
    mean_interburstinterval = np.mean(intervals)
    
    covariance_interburstinterval = np.cov(intervals)
    # Count the number of peaks
    num_peaks = len(peaks)

    # Calculate the mean peak height
    mean_peak_height = np.mean(burstPeakValues)
    cov_peak_height = np.cov(burstPeakValues)

    # # Print the results
    # print("Mean Interval between Peaks:", mean_interburstinterval)
    # print("Number of Peaks:", num_peaks)
    # print("Mean Peak Height:", mean_peak_height)

    # #calculate the mean isi and the variance of the isi
    # mean_isi =np.mean([SpikeTimes.values()])

    #spikecounts= [len(x) for x in SpikeTimes.values()]
    spikecounts = []
    for x in SpikeTimes.values():
        try: count = len(x)
        except:
            if isinstance(x, float): count = 1
            else: count = 0
        spikecounts.append(count)
    
    var_spikecounts = np.var(spikecounts)
    mean_spikecounts = np.mean(spikecounts)
    fanofact = var_spikecounts/mean_spikecounts
    
    #mean burst rate
    time_in_seconds = timeVector[-1] - timeVector[0]
    mean_burst_rate = num_peaks/time_in_seconds

    network_data = {
        "Number_Bursts": num_peaks,
        'mean_Burst_Rate': mean_burst_rate,
        "mean_IBI": mean_interburstinterval,
        "cov_IBI": covariance_interburstinterval,
        "mean_Burst_Peak": mean_peak_height, # TODO: add fitting func for this later
        "cov_Burst_Peak": cov_peak_height,
        "fano_factor": fanofact
    }   


    # Plot the threshold line and burst peaks
    #ax.plot(np.arange(timeVector[-1]), thresholdBurst * rmsFiringRate * np.ones(np.ceil(timeVector[-1]).astype(int)), color='gray')
    ax.plot(burstPeakTimes, burstPeakValues, 'or')  # Plot burst peaks as red circles    

    return ax,network_data

def analyze_convolved_spiking_signal(spike_times, spike_times_by_unit, min_peak_distance, 
                                     binSize, gaussianSigma, 
                                     thresholdBurst, prominence, title = 'Network Activity'):
    fig, ax = plt.subplots()
    #spike_times = {i: rasterData['spkt'][rasterData['spkid'] == i] for i in np.unique(rasterData['spkid'])} #plot_network_activity expects spike_times as dictionary
    #spike_times_dict = {int(i): spike_times[spike_times['spkid'] == i] for i in np.unique(spike_times['spkid'])}
    try:
        #ax, convolved_signal_metrics = helper.plot_network_activity(ax, spike_times_by_unit, min_peak_distance=min_peak_distance, binSize=binSize, gaussianSigma=gaussianSigma, thresholdBurst=thresholdBurst)
        ax, convolved_signal_metrics = plot_network_activity_aw(ax, spike_times_by_unit, min_peak_distance=min_peak_distance, 
                                                                binSize=binSize, gaussianSigma=gaussianSigma, 
                                                                thresholdBurst=thresholdBurst,
                                                                prominence=prominence,
                                                                title=title,
                                                                )
        
        #pull y data out of ax
        y_data = ax.lines[0].get_ydata()
        signal = y_data
        baseline = np.mean(signal)
        
        convolved_data = {
            'mean_Burst_Rate': convolved_signal_metrics['mean_Burst_Rate'],
            'mean_Burst_Peak': convolved_signal_metrics['mean_Burst_Peak'],
            'cov_Burst_Peak': convolved_signal_metrics['cov_Burst_Peak'],
            'fano_factor': convolved_signal_metrics['fano_factor'],
            'mean_IBI': convolved_signal_metrics['mean_IBI'],
            'cov_IBI': convolved_signal_metrics['cov_IBI'],
            'Number_Bursts': convolved_signal_metrics['Number_Bursts'],
            #'baseline': convolve_signal_get_baseline(spike_times, binSize=binSize, gaussianSigma=gaussianSigma),
            'baseline': baseline,
            #'fig': fig,
            'ax': ax,
        }
        return convolved_data
    except Exception as e:
        #set all convolved data to nan
        print(f'Error analyzing convolved spiking signal: {e}')
        print(f'might be cause by all spike times being in the same bin')
        convolved_data = {
            'mean_Burst_Rate': np.nan,
            'mean_Burst_Peak': np.nan,
            'cov_Burst_Peak': np.nan,
            'fano_factor': np.nan,
            'mean_IBI': np.nan,
            'cov_IBI': np.nan,
            'Number_Bursts': np.nan,
            'baseline': np.nan,
        }
        return convolved_data

def init_burst_analysis_params(isi_threshold=None):
    # try:
    #     #from simulate._config_files.burst_analysis_params import burst_analysis_params
    #     burst_analysis_params_temp = burst_analysis_params
    #     burst_analysis_params_temp['isi_threshold'] = isi_threshold if isi_threshold is not None else burst_analysis_params['isi_threshold']
    # except:
    #     burst_analysis_params_temp = {
    #         'isi_threshold': 0.1 if isi_threshold is None else isi_threshold,
    #     }
    
    burst_analysis_params_temp = {
        'isi_threshold': 0.1 if isi_threshold is None else isi_threshold,
    }    
    burst_analysis_params = burst_analysis_params_temp
    return burst_analysis_params

def extract_bursting_activity_data(spike_times, spike_times_by_unit, **kwargs):

    #conv_params = init_convolution_params()
    conv_params = kwargs.get('conv_params', None)
    try: conv_params = conv_params.conv_params if conv_params is not None else None
    except: conv_params = conv_params
    assert conv_params is not None, 'No convolution parameters provided'
    assert 'binSize' in conv_params, 'Convolution parameters must include "binSize"'
    assert 'gaussianSigma' in conv_params, 'Convolution parameters must include "gaussianSigma"'
    
    # # HACK
    # #make mega bursting params by copying the conv_params and multiplying binSize and gaussianSigma
    # mega_conv_params = conv_params.copy()
    # mega_conv_params['binSize'] = mega_conv_params['binSize'] * 5
    # mega_conv_params['gaussianSigma'] = mega_conv_params['gaussianSigma'] * 15
    
    mega_conv_params = kwargs.get('mega_params', None)
    try: mega_conv_params = mega_conv_params.mega_conv_params if mega_conv_params is not None else None
    except: mega_conv_params = mega_conv_params
    assert mega_conv_params is not None, 'No mega convolution parameters provided'
    assert 'binSize' in mega_conv_params, 'Mega convolution parameters must include "binSize"'
    assert 'gaussianSigma' in mega_conv_params, 'Mega convolution parameters must include "gaussianSigma"'
    
    #extract convolution parameters
    binSize = conv_params['binSize']
    gaussianSigma = conv_params['gaussianSigma']
    thresholdBurst = conv_params['thresholdBurst']
    min_peak_distance = conv_params['min_peak_distance']
    
    #unit-wise burst analysis
    burst_analysis_params = init_burst_analysis_params()
    isi_threshold = burst_analysis_params['isi_threshold']
    
    #convolve spiking signal and get baseline
    prominence = conv_params['prominence']
    convolved_data = analyze_convolved_spiking_signal(spike_times, spike_times_by_unit, min_peak_distance, 
                                                      binSize, gaussianSigma, 
                                                      thresholdBurst, prominence, 
                                                      #title='Network Activity (Convolved)'
                                                      )
    bursting_data = analyze_bursting_activity(spike_times, spike_times_by_unit, isi_threshold)
    
    #mega convovle and burst analysis
    prominence = mega_conv_params['prominence']
    mega_convolved_data = analyze_convolved_spiking_signal(spike_times, spike_times_by_unit,
                                                            mega_conv_params['min_peak_distance'], 
                                                            mega_conv_params['binSize'], 
                                                            mega_conv_params['gaussianSigma'], 
                                                            mega_conv_params['thresholdBurst'],
                                                            prominence,
                                                            title='Network Activity (Mega)'
                                                          )
    mega_bursting_data = analyze_bursting_activity(spike_times, spike_times_by_unit, isi_threshold)
    #mega_bursting_data = bursting_data.copy()

    #add convolved data to the bursting summary data
    bursting_summary_data = bursting_data['bursting_summary_data']
    for key in convolved_data.keys():
        bursting_summary_data[key] = convolved_data[key]
    bursting_data['bursting_summary_data'] = bursting_summary_data
    
    #mega_bursting_data = mega_bursting_data
    mega_bursting_summary_data = mega_bursting_data['bursting_summary_data']
    for key in mega_convolved_data.keys():
        mega_bursting_summary_data[key] = mega_convolved_data[key]
    
    # Verify that any single-element array is converted to a scalar
    for key in bursting_data['bursting_summary_data'].keys():
        value = bursting_data['bursting_summary_data'][key]
        if isinstance(value, np.ndarray) and value.size == 1:
            bursting_data['bursting_summary_data'][key] = value.item()  # .item() fetches the scalar value directly
            
    # verifying that any single-element array is converted to a scalar
    for key in mega_bursting_summary_data.keys():
        value = mega_bursting_summary_data[key]
        if isinstance(value, np.ndarray) and value.size == 1:
            mega_bursting_summary_data[key] = value.item()
    
    network_data['bursting_data'] = bursting_data
    network_data['mega_bursting_data'] = mega_bursting_data
    # network_data = {
    #     'bursting_data': bursting_data,
    # }

    #return network_data

def get_mean_fr(recording_object, sorting_object, sampling_rate=10000):
    total_fr = 0
    unit_count = 0
    units = sorting_object.get_unit_ids()
    for unit in units:
        spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
        assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
        if len(spike_train) > 1:
            duration = recording_object.get_total_duration()  # seconds
            fr = len(spike_train) / duration  # spikes per second
            if not np.isnan(fr) and not np.isinf(fr):
                total_fr += fr
                unit_count += 1
    return total_fr / unit_count if unit_count > 0 else None

def get_CoV_fr_experimental(recording_object, sorting_object, sampling_rate=10000, window_size=1.0):
    """
    Calculate the Coefficient of Variation of Firing Rate (CoV FR) for experimental data over time windows.

    Parameters:
    - recording_object: A recording extractor object containing the experiment's total duration.
    - sorting_object: A sorting extractor object with spike times for each unit.
    - sampling_rate: Sampling rate of the recording (in Hz, default: 10000).
    - window_size: Size of the time window for calculating firing rates (default: 1.0 second).

    Returns:
    - CoV_fr: Coefficient of Variation of Firing Rate (float).
    """
    # Get total duration of the recording (in seconds)
    total_duration = recording_object.get_total_duration()

    # Initialize list to hold firing rates across all windows for all units
    firing_rates = []

    # Get unit IDs from the sorting object
    units = sorting_object.get_unit_ids()

    # Process each unit's spike train
    for unit in units:
        spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate  # Convert spike times to seconds
        assert all(spike_time >= 0 for spike_time in spike_train), f"Spike times for unit {unit} contain negative values"

        if len(spike_train) > 1:
            # Divide the total duration into non-overlapping time windows
            num_windows = int(np.ceil(total_duration / window_size))
            window_edges = np.linspace(0, total_duration, num_windows + 1)

            # Count spikes in each window for the current unit
            spike_counts = np.histogram(spike_train, bins=window_edges)[0]

            # Convert spike counts to firing rates (spikes per second) for this unit
            unit_firing_rates = spike_counts / window_size

            # Append non-NaN and non-inf firing rates to the global list
            firing_rates.extend([fr for fr in unit_firing_rates if not np.isnan(fr) and not np.isinf(fr)])

    # Calculate CoV FR if there are valid firing rates
    if len(firing_rates) > 1:
        mean_fr = np.mean(firing_rates)
        std_fr = np.std(firing_rates)
        CoV_fr = std_fr / mean_fr
        return CoV_fr
    else:
        # Return None if there are no valid firing rates
        return None

def get_mean_isi(recording_object, sorting_object, sampling_rate=10000):
    all_means = []
    units = sorting_object.get_unit_ids()
    for unit in units:
        spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
        assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
        if len(spike_train) > 1:
            isi = np.diff(spike_train)
            mean_isi = np.mean(isi)
            if not np.isnan(mean_isi) and not np.isinf(mean_isi):
                all_means.append(mean_isi)
    
    overall_mean = np.mean(all_means) if all_means else None
    return overall_mean if overall_mean is not None else None  # already in seconds

def get_CoV_isi(recording_object, sorting_object, sampling_rate=10000):
    all_CoVs = []
    units = sorting_object.get_unit_ids()
    for unit in units:
        spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
        assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
        if len(spike_train) > 1:
            isi = np.diff(spike_train)  # already in seconds
            if len(isi) > 1:
                CoV = np.std(isi) / np.mean(isi)
                if not np.isnan(CoV) and not np.isinf(CoV):
                    all_CoVs.append(CoV)
    
    overall_mean_CoV = np.mean(all_CoVs) if all_CoVs else None
    return overall_mean_CoV

def get_unit_fr(recording_object, sorting_object, unit, sampling_rate=10000):
    spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
    assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
    if len(spike_train) > 1:
        duration = recording_object.get_total_duration()  # seconds
        fr = len(spike_train) / duration  # spikes per second
        return fr
    else:
        return 0

def get_unit_fr_CoV(recording_object, sorting_object, unit, sampling_rate=10000):
    spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
    assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
    if len(spike_train) > 1:
        isi = np.diff(spike_train)  # inter-spike intervals
        if len(isi) > 1:
            CoV = np.std(isi) / np.mean(isi)
            return CoV
        else:
            return None
    else:
        return None
    
def get_unit_mean_isi(recording_object, sorting_object, unit, sampling_rate=10000):
    spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
    assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
    if len(spike_train) > 1:
        isi = np.diff(spike_train)
        mean_isi = np.mean(isi)
        return mean_isi  # already in seconds
    else:
        return None
    
def get_unit_isi_CoV(recording_object, sorting_object, unit, sampling_rate=10000):
    spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
    assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
    if len(spike_train) > 2:
        isi = np.diff(spike_train)  # already in seconds
        if len(isi) > 1:
            CoV = np.std(isi) / np.mean(isi)
            return CoV
        else:
            return None
    else:
        return None

import numpy as np

def get_unit_fano_factor(recording_object, sorting_object, unit, sampling_rate=10000, bin_size=6.0):
    """
    Computes the Fano Factor (FF) for a given unit, which is the variance-to-mean ratio
    of spike counts in a 6-second time window.

    Parameters:
    - recording_object: SpikeInterface recording object.
    - sorting_object: SpikeInterface sorting object.
    - unit: The unit ID to analyze.
    - sampling_rate (int): Sampling rate in Hz (default: 10000 Hz).
    - bin_size (float): Window size in seconds for spike count bins (default: 6s).

    Returns:
    - Fano Factor (float) or None if not enough spikes.
    """

    # Get spike train and convert to seconds
    spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate  
    assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'

    # Ensure sufficient spikes
    if len(spike_train) < 2:
        return None  # Not enough data

    # Get total duration of recording (in seconds)
    duration = recording_object.get_total_duration()
    
    # Define bin edges for 6-second bins
    num_bins = max(1, int(duration / bin_size))  # Ensure at least 1 bin
    bin_edges = np.linspace(0, duration, num_bins + 1)
    
    # Compute spike counts in bins
    spike_counts, _ = np.histogram(spike_train, bins=bin_edges)
    
    # Compute mean and variance of spike counts
    mean_count = np.mean(spike_counts)
    var_count = np.var(spike_counts)

    # Compute Fano Factor (variance-to-mean ratio)
    fano_factor = var_count / mean_count if mean_count > 0 else None

    return fano_factor

def extract_metrics_from_experimental_data(spike_times, timeVector, spike_times_by_unit, **kwargs):
    
    #extract spiking metrics from experimental data
    recording_object = kwargs['recording_object']
    sorting_object = kwargs['sorting_object']
    sampling_rate = recording_object.get_sampling_frequency()    
    
    #add immediately available data to network_data
    network_data['source'] = 'experimental'
    network_data['timeVector'] = timeVector
    network_data['spiking_data']['spike_times'] = spike_times
    network_data['spiking_data']['spiking_times_by_unit'] = spike_times_by_unit
    
    #overall mean firing rate
    network_data['spiking_data']['spiking_summary_data']['MeanFireRate'] = get_mean_fr(
        recording_object, sorting_object, sampling_rate=sampling_rate
        )
    
    #overall CoV firing rate
    network_data['spiking_data']['spiking_summary_data']['CoVFireRate'] = get_CoV_fr_experimental(
        recording_object, sorting_object, sampling_rate=sampling_rate
        )
    
    #overall mean ISI
    network_data['spiking_data']['spiking_summary_data']['MeanISI'] = get_mean_isi(
        recording_object, sorting_object, sampling_rate=sampling_rate
    )
    
    #overall CoV ISI
    network_data['spiking_data']['spiking_summary_data']['CoV_ISI'] = get_CoV_isi(
        recording_object, sorting_object, sampling_rate=sampling_rate
    )
    
    #spiking data by unit - just go through the simulated version of the data but de-identify pop
    #simulated_spiking_data_by_unit = network_data['simulated_data']['spiking_data_by_unit']
    spiking_data_by_unit = {}
    units = sorting_object.get_unit_ids()
    for unit in units:
        spiking_data_by_unit[unit] = {
            #'unitProperty': sorting_object.get_unit_property(unit),
            #'SpikeTrain': sorting_object.get_unit_spike_train(unit),
            'FireRate': get_unit_fr(recording_object, sorting_object, unit),
            'fr_CoV': get_unit_fr_CoV(recording_object, sorting_object, unit),
            'meanISI': get_unit_mean_isi(recording_object, sorting_object, unit),
            'isi_CoV': get_unit_isi_CoV(recording_object, sorting_object, unit),
            'fano_factor': get_unit_fano_factor(recording_object, sorting_object, unit), # HACK: Hardcoded 6 second bin size
                                                                                        # consistent with the 6 second bin size
                                                                                        # user here: https://pmc.ncbi.nlm.nih.gov/articles/PMC3434456/
            'spike_times': spike_times_by_unit[unit],
        }
    network_data['spiking_data']['spiking_data_by_unit'] = spiking_data_by_unit

def get_spike_times(recording_object, sorting_object, sampling_rate=10000):
    spike_times = []
    units = sorting_object.get_unit_ids()
    total_duration = recording_object.get_total_duration() #seconds
    for unit in units:
        spike_train = sorting_object.get_unit_spike_train(unit) #in samples
        spike_times.extend(spike_train)
    spike_times.sort()
    spike_times = np.array(spike_times) / sampling_rate #convert to seconds
    
    # quality control
    # spike_times_max need to be rounded to the nearest thousandths place to avoid floating point errors
    spike_times_max = np.round(max(spike_times), 3)
    
    print(f'Max Spike Time: {max(spike_times)}')
    #assert max(spike_times) <= total_duration, 'Spike times are not in seconds'
    assert spike_times_max <= total_duration, 'Spike times are not in seconds'
    assert all(spike_time >= 0 for spike_time in spike_times), 'Spike times contain negative values'
    
    #return
    return spike_times

def get_time_vector(recording_object, sampling_rate=10000):
    duration = recording_object.get_total_duration() #seconds
    time_vector = np.linspace(0, duration, int(duration * sampling_rate))
    assert len(time_vector) == int(duration * sampling_rate), 'Time vector length mismatch'
    return time_vector

def get_spike_times_by_unit(sorting_object, sampling_rate=10000):
    spike_times_by_unit = {}
    units = sorting_object.get_unit_ids()
    for unit in units:
        spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
        assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
        spike_times_by_unit[unit] = spike_train
    return spike_times_by_unit

def init_network_data_dict():
    # global network_data
    # from RBS_network_models.network_data_struct import network_data_dict
    # empty_network_data = network_data_dict
    # network_data = empty_network_data
    # return empty_network_data
    
    empty_network_data = {
        #General Data
        'source': None, # 'simulated' or 'experimental'
        'timeVector': None,
        
        #Simulated Data
        'simulated_data': {
            'soma_voltage': None,
            'E_Gids': None,
            'I_Gids': None,
            'MeanFireRate_E': None,
            'CoVFireRate_E': None,
            'MeanFireRate_I': None,
            'CoVFireRate_I': None,
            'MeanISI_E': None,
            'MeanISI_I': None,
            'CoV_ISI_E': None,
            'CoV_ISI_I': None,
            'spiking_data_by_unit': None, 
        },
        
        #Spiking Data
        'spiking_data': {
            'spike_times': None,
            'spiking_summary_data': {
                #'spike_times': None,
                'MeanFireRate': None,
                'CoVFireRate': None,
                'MeanISI': None,
                'CoV_ISI': None,         
            },
            'spiking_times_by_unit': None,
            'spiking_data_by_unit': None,
        },
        
        #Bursting Data
        'bursting_data': {
            'bursting_summary_data': {
                'baseline': None,
                'MeanWithinBurstISI': None,
                'CovWithinBurstISI': None,
                'MeanOutsideBurstISI': None,
                'CoVOutsideBurstISI': None,
                'MeanNetworkISI': None,
                'CoVNetworkISI': None,
                'NumUnits': None,
                'Number_Bursts': None,
                'mean_IBI': None,
                'cov_IBI': None,
                'mean_Burst_Peak': None,
                'cov_Burst_Peak': None,
                'fano_factor': None,
            },
            'bursting_data_by_unit': None,
        },
        
        # Neuron Classification Data
        'classification_data': {
            'waveform_data': None,
            'classification_data': None,
            'excitatory_neurons': None,
            'inhibitory_neurons': None,            
        },
    }
    global network_data
    network_data = empty_network_data
    return empty_network_data
    
    ## Keep for reference
    
            # return {
        #     # 'burstPeakValues': None,
        #     # 'IBIs': None,
        #     # 'baseline': None,
        #     # 'peakFreq': None,
        #     # 'firingRate': None,
        #     # 'burstPeakTimes': None,
        #     # 'timeVector': None,
        #     # 'threshold': None,
            
        #     'Number_Bursts': None,
        #     'mean_IBI': None,
        #     'cov_IBI': None,
        #     'mean_Burst_Peak': None,
        #     'cov_Burst_Peak': None,
        #     'fano_factor': None,
        #     'MeanWithinBurstISI': None,
        #     'CoVWithinBurstISI': None,
        #     'MeanOutsideBurstISI': None,
        #     'CoVOutsideBurstISI': None,
        #     'MeanNetworkISI': None,
        #     'CoVNetworkISI': None,
        #     'NumUnits': None,
        #     #'fileName': None
        # }

def get_experimental_network_activity_metrics(sorting_object, recording_segment_object, conv_params, mega_params, **kwargs):
    
    assert sorting_object is not None, 'No sorting object provided'
    
    #initialize network_data
    network_data = init_network_data_dict()
    
    #network data pathing info
    network_data['recording_path'] = recording_segment_object._kwargs.get('file_path', None)
    network_data['sorting_output'] = sorting_object._kwargs.get('folder_path', None)
    assert network_data['sorting_output'] is not None, 'No sorting object source found'
    
    # waveform info
    network_data['waveform_output'] = network_data['sorting_output'].replace('sorted', 'waveforms').replace('sorter_output', '')
    
    #get data
    sorting_object = sorting_object
    recording_object = recording_segment_object
    sampling_rate = recording_object.get_sampling_frequency()
    
    # analyze each neuron. classify if possible. Excitatory (RS) or Inhibitory (FS)
    #from RBS_network_models.experimental import classify_neurons
    #network_data['waveform_output'] = network_data['sorting_output'].replace('sorted', 'neurons').replace('sorter_output', '')
    #network_data['waveform_output'] = network_data['waveform_output'].replace('sorter', 'waveform')
    #classifications = classify_neurons(sorting_object, recording_object, sampling_rate=sampling_rate, output_path=network_data['waveform_output']) 
    
    #convert time to seconds - get initially available data
    spike_times = get_spike_times(recording_object, sorting_object, sampling_rate=sampling_rate) #seconds
    timeVector = get_time_vector(recording_object, sampling_rate=sampling_rate) #seconds
    spike_times_by_unit = get_spike_times_by_unit(sorting_object, sampling_rate=sampling_rate) 
    
    #add conv_params to kwargs
    kwargs['conv_params'] = conv_params
    kwargs['mega_params'] = mega_params
    
    #add sorting_object and recording_object to kwargs if not already present
    kwargs['sorting_object'] = sorting_object
    kwargs['recording_object'] = recording_object
    
    #extract spiking metrics from simulated data
    try: 
        extract_metrics_from_experimental_data(spike_times, timeVector, spike_times_by_unit, **kwargs)
    except Exception as e:
        print(f'Error extracting metrics from simulated data: {e}')
        pass
    
    #extract bursting metrics from simulated data (but this one works for both simulated and experimental data)
    try: 
        extract_bursting_activity_data(spike_times, spike_times_by_unit, **kwargs)
    except Exception as e:
        print(f'Error calculating bursting activity: {e}')
        pass
    
    #classify neurons
    from RBS_network_models.experimental import classify_neurons
    try:
        network_data = classify_neurons(network_data, **kwargs)
    except Exception as e:
        print(f'Error classifying neurons: {e}')
        pass
    
    return network_data

def get_simulated_network_activity_metrics_nvm(conv_params, mega_params, **kwargs):
    
    #initialize network_data
    network_data = init_network_data_dict()
    
    #add conv_params to kwargs
    kwargs['conv_params'] = conv_params
    kwargs['mega_params'] = mega_params
    
    #extract spiking metrics from simulated data
    try: 
        extract_metrics_from_experimental_data(spike_times, timeVector, spike_times_by_unit, **kwargs)
    except Exception as e:
        print(f'Error extracting metrics from simulated data: {e}')
        pass
    
    #extract bursting metrics from simulated data (but this one works for both simulated and experimental data)
    try: 
        extract_bursting_activity_data(spike_times, spike_times_by_unit, **kwargs)
    except Exception as e:
        print(f'Error calculating bursting activity: {e}')
        pass
    
    #classify neurons
    from RBS_network_models.experimental import classify_neurons
    try:
        network_data = classify_neurons(network_data, **kwargs)
    except Exception as e:
        print(f'Error classifying neurons: {e}')
        pass
    
    return network_data
