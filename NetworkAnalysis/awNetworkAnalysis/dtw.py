import time
from fastdtw import fastdtw
from scipy.spatial.distance import euclidean
import numpy as np
import pandas as pd
import time
import random
import multiprocessing
import os
import scipy.stats as stats
import traceback
# =================== #
# Function Definitions

def build_sequence_stacks(network_metrics, burst_metrics, **kwargs):
    #### 
    #burst_metrics = network_metrics['bursting_data']['burst_metrics']
    #HACK
    source = kwargs['source']
    if source == 'experimental':
        #classified_units = network_metrics['classification_output']['classified_units']
        unit_types = network_metrics['classification_output']['classified_units']
    elif source == 'simulated':
        #classified_units = network_metrics['unit_types']
        unit_types = network_metrics['unit_types']
    
    #### Burst Participation Metrics ####
    burst_parts = burst_metrics['burst_parts']
    sequence_stacks = {}
    # time_sequence_mat = []
    # cat_sequence_mat = []
    for burst_id in burst_parts:
        
        #
        burst = burst_parts[burst_id]
        try: unit_sequence = burst['unit_sequence']
        except: unit_sequence = burst['unit_seqeunce'] #TODO #stupid typo in analysis code - fix and rerun
        time_sequence = burst['time_sequence']
        relative_time_sequence = burst['relative_time_sequence']
        
        #classification output
        #cross reference participating units with classification data to get excitatory and inhibitory units
        #participating_units = burst_parts[burst_id]['participating_units']
        #classified_units = network_metric_targets['classification_output']['classified_units']
        #classified_units = network_metrics['classification_output']['classified_units']
        #HACK
        source = kwargs['source']
        if source == 'experimental':
            I_units = [unit for unit in unit_types if unit_types[unit]['desc'] == 'inhib']
            E_units = [unit for unit in unit_types if unit_types[unit]['desc'] == 'excit']
        elif source == 'simulated':
            I_units = [unit for unit in unit_types if unit_types[unit] == 'I']
            E_units = [unit for unit in unit_types if unit_types[unit] == 'E']
        classified_sequence = []
        for unit in unit_sequence:
            if unit in I_units:
                classified_sequence.append('I')
            elif unit in E_units:
                classified_sequence.append('E')
            else:
                classified_sequence.append('U')
        burst_parts[burst_id]['classified_sequence'] = classified_sequence
        #print('classified_sequence:', classified_sequence)
        
        #stack sequences for each unit
        cat_sequence = [0 if x == 'I' else 1 if x == 'E' else 2 for x in classified_sequence]
        burst_parts[burst_id]['cat_sequence'] = cat_sequence
        sequence_stack = np.column_stack((
            #unit_sequence, 
            #time_sequence, 
            cat_sequence,
            relative_time_sequence,)).astype(float)
            #cat_sequence)).astype(float)
            
        # cat_sequence_mat.append(cat_sequence)
        # time_sequence_mat.append(relative_time_sequence)            
        
        #store for DTW analysis
        sequence_stacks[burst_id] = np.array(sequence_stack)
        
        #store for DTW analysis
        burst_parts[burst_id]['sequence_stack'] = sequence_stack
        
    #update network metrics
    network_metrics['bursting_data']['burst_metrics']['burst_parts'] = burst_parts
    
    #
    return sequence_stacks, network_metrics

def dtw_burst_analysis(network_metrics, kwargs):
    print('Running DTW Burst Analysis...')
    
    # Subfuncs ========================================
    def prep_outputs(data_path, results, mega_results):
        # get the most representative sequences - get the column with the lowest sum indicating lowest difference to all other sequences
        # regular bursts
        try:
            min_dtw_sum = np.sum(dtw_distance_matrix, axis=0)
        except:
            min_dtw_sum = np.nan
        try:
            min_dtw_idx = np.argmin(min_dtw_sum)
        except:
            min_dtw_idx = np.nan
        try:
            min_dtw_seq = sequence_stacks[min_dtw_idx]
        except:
            min_dtw_seq = np.nan
        try:
            min_dtw_diff_vector = dtw_distance_matrix[min_dtw_idx]
        except:
            min_dtw_diff_vector = np.nan
        try:
            min_dtw_diff_vector_sorted = np.sort(min_dtw_diff_vector) # get second min diff value - the first will always be 0
        except:
            min_dtw_diff_vector_sorted = np.nan
        try:
            diff_min = min_dtw_diff_vector_sorted[1]
        except:
            diff_min = np.nan
        try:
            diff_max = np.max(min_dtw_diff_vector)
        except:
            diff_max = np.nan
        try:
            diff_mean = np.mean(min_dtw_diff_vector)
        except:
            diff_mean = np.nan
        try:
            diff_median = np.median(min_dtw_diff_vector)
        except:
            diff_median = np.nan
        try:
            diff_std = np.std(min_dtw_diff_vector)
        except:
            diff_std = np.nan
        try:
            diff_cov = diff_std / diff_mean if diff_mean > 0 and diff_std > 0 else 0
        except:
            diff_cov = np.nan
        
        # mega bursts
        try:
            min_mega_dtw_sum = np.sum(mega_dtw_distance_matrix, axis=0)
        except:
            min_mega_dtw_sum = np.nan
        try:
            min_mega_dtw_idx = np.argmin(min_mega_dtw_sum)
        except:
            min_mega_dtw_idx = np.nan
        try:
            min_mega_dtw_seq = mega_sequence_stacks[min_mega_dtw_idx]
        except:
            min_mega_dtw_seq = np.nan
        try:
            min_mega_dtw_diff_vector = mega_dtw_distance_matrix[min_mega_dtw_idx]
        except:
            min_mega_dtw_diff_vector = np.nan
        try:
            min_mega_dtw_diff_vector_sorted = np.sort(min_mega_dtw_diff_vector) # get second min diff value - the first will always be 0
        except:
            min_mega_dtw_diff_vector_sorted = np.nan
        try:
            mega_diff_min = min_mega_dtw_diff_vector_sorted[1]
        except:
            mega_diff_min = np.nan
        try:
            mega_diff_max = np.max(min_mega_dtw_diff_vector)
        except:
            mega_diff_max = np.nan
        try:
            mega_diff_mean = np.mean(min_mega_dtw_diff_vector)
        except:
            mega_diff_mean = np.nan
        try:
            mega_diff_median = np.median(min_mega_dtw_diff_vector)
        except:
            mega_diff_median = np.nan
        try:
            mega_diff_std = np.std(min_mega_dtw_diff_vector)
        except:
            mega_diff_std = np.nan
        try:
            mega_diff_cov = mega_diff_std / mega_diff_mean if mega_diff_mean > 0 and mega_diff_std > 0 else 0
        except:
            mega_diff_cov = np.nan
        
        # update result dicts
        results['mindiff_sequence_analysis'] = {
            'min_dtw_idx': min_dtw_idx,
            'min_dtw_seq': min_dtw_seq,
            'min_dtw_diff_vector': min_dtw_diff_vector,
            'diff_min': diff_min,
            'diff_max': diff_max,
            'diff_mean': diff_mean,
            'diff_median': diff_median,
            'diff_std': diff_std,
            'diff_cov': diff_cov
        }
        
        mega_results['mindiff_sequence_analysis'] = {
            'min_mega_dtw_idx': min_mega_dtw_idx,
            'min_mega_dtw_seq': min_mega_dtw_seq,
            'min_mega_dtw_diff_vector': min_mega_dtw_diff_vector,
            'mega_diff_min': mega_diff_min,
            'mega_diff_max': mega_diff_max,
            'mega_diff_mean': mega_diff_mean,
            'mega_diff_median': mega_diff_median,
            'mega_diff_std': mega_diff_std,
            'mega_diff_cov': mega_diff_cov
        }
        
        return results, mega_results
    
    # Main ========================================
    # check if doing dtw computations is possible.
    # check if burst parts are available to compute dtw
    bursting_data = network_metrics.get('bursting_data', None)
    if bursting_data is not None:
        burst_metrics = bursting_data.get('burst_metrics', None)
        if burst_metrics is not None:
            burst_parts = burst_metrics.get('burst_parts', None)
            if 'Burst sequencing not enabled' in burst_parts:
                network_metrics['dtw_output'] = 'Burst sequencing not enabled'
                network_metrics['mega_dtw_output'] = 'Burst sequencing not enabled'
                network_metrics['dtw_burst_analysis'] = 'Burst sequencing not enabled'
                network_metrics['dtw_mega_burst_analysis'] = 'Burst sequencing not enabled'
                return network_metrics
    
    #init
    max_workers = kwargs.get('max_workers', 4)  

    # prep sequence stacks
    burst_metrics = network_metrics['bursting_data']['burst_metrics']
    sequence_stacks, network_metrics = build_sequence_stacks(network_metrics, burst_metrics, **kwargs)
    
    mega_burst_metrics = network_metrics['mega_bursting_data']['burst_metrics']
    mega_sequence_stacks, network_metrics = build_sequence_stacks(network_metrics, mega_burst_metrics, **kwargs)
    
    # dtw_path = network_metrics['output_paths']['dtw_output']
    # mega_dtw_path = network_metrics['output_path']['mega_dtw_output']
    dtw_temp_dir = kwargs['dtw_temp']
    dtw_path = os.path.join(dtw_temp_dir, 'dtw_output')
    mega_dtw_path = os.path.join(dtw_temp_dir, 'mega_dtw_output')
    
    # save network metrics dirs
    network_metrics['dtw_output'] = dtw_path
    network_metrics['mega_dtw_output'] = mega_dtw_path
    
    # set data path - base on sorter_output
    # sorting_output_dir = network_metrics['sorting_output']
    # # replace 'sorted' with 'dtw' in path
    # data_path = sorting_output_dir.replace('sorted', 'dtw')
    # # replace 'sorter_output' with 'dtw_output' in path
    # data_path = data_path.replace('sorter_output', 'dtw_output')
    # mega_data_path = data_path.replace('dtw_output', 'mega_dtw')
    
    # create data paths if they dont exist
    #import os
    if not os.path.exists(dtw_path):
        os.makedirs(dtw_path)
    if not os.path.exists(mega_dtw_path):
        os.makedirs(mega_dtw_path)
    
    # run dtw analysis - regular and mega bursts
    #num_workers = 50    # NOTE: aw 2025-02-24 01:50:50 - good number for login node where I'll usually run this
    #num_workers = kwargs.get('max_workers', 50)
    batch_size = 1000 # this will usually work for regualr bursting, but this might be too large for mega bursting.
    mega_batch_size = 5 # based on the experimental data set, 23 bursts. 23^2 computations = 529 - so to use 50 workers, we need to use a batch size of, at most, 529 // 50 = 10
    num_workers = max_workers
    
    # HACK
    if len(sequence_stacks) < 1000:
        batch_size = len(sequence_stacks)//4
        
    if len(mega_sequence_stacks) < 5:
        #mega_batch_size = len(mega_sequence_stacks)
        mega_batch_size = 1
    
    # Do Mega first - it's probably quicker to debug this way.
    #if len(mega_sequence_stacks)>1: # cant cross compare sequences if there's only one...
    mega_results = dtw_analysis_dynamic_v2(
        mega_sequence_stacks, 
        data_path=mega_dtw_path, 
        num_workers=num_workers,
        batch_size=mega_batch_size
        )    
    
    # Do regular bursts - this should take quite a bit longer...
    results = dtw_analysis_dynamic_v2(
        sequence_stacks, 
        data_path=dtw_path, 
        num_workers=num_workers,
        batch_size=batch_size        
        )
    
    # get sequences with the lowest DTW distances - i.e. the most representative sequences regular and mega bursts
    dtw_distance_matrix = results['dtw_distance_matrix']
    mega_dtw_distance_matrix = mega_results['dtw_distance_matrix']
    
    # prep outputs
    results, mega_results = prep_outputs(dtw_path, results, mega_results)
        
    # update network metrics
    network_metrics['dtw_burst_analysis'] = results
    network_metrics['dtw_mega_burst_analysis'] = mega_results
    
    return network_metrics   
    
def init_worker(l, c, sv, ssv, nv, sf, pm, er, st, gmvl):
    """Each worker initializes shared variables (avoiding pickling issues)."""
    global lock, counter, sum_value, sum_sq_value, n_value, stop_flag, prev_means, effective_rate, start_time, global_matrix_value_list
    lock = l
    counter = c
    sum_value = sv
    sum_sq_value = ssv
    n_value = nv
    stop_flag = sf
    prev_means = pm
    effective_rate = er
    start_time = st
    global_matrix_value_list = gmvl

def compute_dtw_dynamic_v2(pairs_sample, sequence_stacks, data_path=None, cv_threshold=0.05, moving_avg_window=500, moving_avg_threshold=0.01, confidence_threshold=0.99):
    """Computes DTW distances and updates running statistics in parallel."""
    
    global lock, counter, sum_value, sum_sq_value, n_value, stop_flag, prev_means, effective_rate, start_time, global_matrix_value_list

    proc_start_time = time.time()
    print(f"Worker started processing at {time.strftime('%H:%M:%S', time.gmtime(proc_start_time))}")

    # get batch_size for progress updates
    batch_size = len(pairs_sample)
    
    for i, j in pairs_sample:
        if not np.isnan(global_matrix_value_list[i * len(sequence_stacks) + j]) and \
           not np.isnan(global_matrix_value_list[j * len(sequence_stacks) + i]):
            continue  # Skip if already computed

        seq_i = sequence_stacks[i]
        seq_j = sequence_stacks[j]
        distance, _ = fastdtw(seq_i, seq_j, dist=euclidean)

        with lock:  # 🔒 Ensure only one process updates shared resources at a time
            global_matrix_value_list[i * len(sequence_stacks) + j] = distance
            global_matrix_value_list[j * len(sequence_stacks) + i] = distance
            sum_value.value += distance
            sum_sq_value.value += distance ** 2
            n_value.value += 1
            counter.value += 1

            mean = sum_value.value / n_value.value
            prev_means.append(mean)
            if len(prev_means) > 500:
                prev_means.pop(0)

            prog_count = batch_size // 5
            if prog_count == 0:
                prog_count = 1
            #if counter.value % 100 == 0:
            if counter.value % prog_count == 0:
                variance = (sum_sq_value.value / n_value.value) - mean ** 2
                std_dev = np.sqrt(variance) if variance > 0 else 0
                effective_rate.value = counter.value / (time.time() - start_time)

                #print(f"Progress: {n_value.value} | Mean: {mean:.3f} | Std: {std_dev:.3f} | Rate: {effective_rate.value:.3f} comp/s")
                
                # get z-score for confidence threshold

                z = stats.norm.ppf(confidence_threshold)
                margin_error = z * (std_dev / np.sqrt(n_value.value))               
                cov = std_dev / mean if mean > 0 and std_dev > 0 else 0
                moving_average_diff = abs(prev_means[-1] - prev_means[0])
                moving_threshold = moving_avg_threshold * prev_means[-1]
                
                # Print Progress
                total_possible_computations = len(sequence_stacks) * (len(sequence_stacks) - 1) / 2
                estimated_time_remaining = (total_possible_computations - n_value.value) / effective_rate.value
                estimated_time_remaining = time.strftime('%H:%M:%S', time.gmtime(estimated_time_remaining))
                print(f"Progress: {n_value.value}/{total_possible_computations}"
                      F"| Mean: {mean:.3f} | Std: {std_dev:.3f} | "
                      #f"CoV: {cov:.3f} | "
                      f"Margin Error: {margin_error:.3f} --> Target: {confidence_threshold:.3f} | "
                      f"CoV: {cov:.3f} --> Target: {cv_threshold:.3f} | "
                      #f"Moving Avg Diff: {moving_average_diff:.3f} --> Target: {moving_threshold:.3f} | "
                      #f"Moving Avg Diff: {moving_average_diff:.3f} | "
                      #f"Threshold: {moving_threshold:.3f} | "
                      f"Effective Rate: {effective_rate.value:.3f} computations/s | "
                      f"ETA: {estimated_time_remaining}"
                      )
                
                if margin_error < confidence_threshold:
                    switch = True
                    if switch:
                        print(f"\tConfidence Stabilized: {margin_error:.3f} < {confidence_threshold:.3f}. "
                              "Confidence has stabilized. "
                              "This stopping condition is enabled. "
                              "Stopping computation.")
                        stop_flag.value = 1
                        break
                    else:
                        print(f"\tConfidence Stabilized: {margin_error:.3f} < {confidence_threshold:.3f}. "
                              "Confidence has stabilized. "
                              "This stopping condition is disabled. "
                              "Continuing computation.")
                
                if cov < cv_threshold:
                    switch = False
                    if switch:
                        print(f"\tCV Stabilized: {cov:.3f} < {cv_threshold:.3f}. "
                              "CV has stabilized. "
                              "This stopping condition is enabled. "
                              "Stopping computation.")
                        stop_flag.value = 1
                        break
                    else:
                        print(f"\tCV Stabilized: {cov:.3f} < {cv_threshold:.3f}. "
                              "CV has stabilized. "
                              "This stopping condition is disabled. "
                              "Continuing computation.")
                        
                if moving_average_diff < moving_threshold:
                    switch = False
                    if switch:
                        print(f"\tMoving Average Stabilized: {moving_average_diff:.3f} < {moving_threshold:.3f}. "
                              "Moving average has stabilized. "
                              "This stopping condition is enabled. "
                              "Stopping computation.")
                        stop_flag.value = 1
                        break
                    else:
                        verbose = False
                        if verbose:
                            print(f"\tMoving Average Stabilized: {moving_average_diff:.3f} < {moving_threshold:.3f}. "
                                "Moving average has stabilized. "
                                "This stopping condition is disabled. "
                                "Continuing computation.")
                            
            # save global matrix to file every 5000 computations
            if data_path is not None:
                save_count = batch_size * 10
                #if counter.value % 100000 == 0:
                if counter.value % save_count == 0:
                    # Save global matrix to file
                    print("\tSaving global matrix to file...")
                    import pandas as pd
                    global_matrix = np.array(global_matrix_value_list).reshape(len(sequence_stacks), len(sequence_stacks))
                    csv_path = os.path.join(data_path, 'dtw_analysis_dynamic_results.csv')
                    pd.DataFrame(global_matrix).to_csv(csv_path)
                    print(f"\tGlobal matrix saved to {csv_path}")
                    
                    #pd.DataFrame(global_matrix).to_csv(data_path + 'dtw_analysis_dynamic_results.csv')
        
    # save global matrix to file workers complete work
    # with lock:
    #     print("\tWorking finishing processing. Saving global matrix to file...")
        
    #     import pandas as pd
    #     global_matrix = np.array(global_matrix_value_list).reshape(len(sequence_stacks), len(sequence_stacks))       
    #     pd.DataFrame(global_matrix).to_csv(data_path + 'dtw_analysis_dynamic_results.csv')
        
    # print progress
    print(f"\tWorker finished processing in {time.strftime('%H:%M:%S', time.gmtime(time.time() - proc_start_time))}")
    
    # import sys
    # sys.exit()
        
def dtw_analysis_dynamic_v2(sequence_stacks, data_path=None, 
                            cv_threshold=0.05, moving_avg_window=500, 
                            moving_avg_threshold=0.01, confidence_threshold=0.99,
                            batch_size=None, num_workers=None,
                            ):
    """Runs DTW analysis with real-time confidence-based stopping."""

    # Subfuncs ========================================
    
    def try_load_data(data_path):
        global global_matrix_value_list, sum_value, sum_sq_value, n_value
        try:
            if data_path is not None:
                csv_path = os.path.join(data_path, 'dtw_analysis_dynamic_results.csv')
                if os.path.exists(csv_path):
                    df = pd.read_csv(csv_path)
                    global_matrix = df.to_numpy()
                    global_matrix = df.to_numpy()[:, 1:]  # Remove index column

                    for idx, val in np.ndenumerate(global_matrix):
                        if not np.isnan(val):
                            global_matrix_value_list[idx[0] * num_sequences + idx[1]] = val
                            global_matrix_value_list[idx[1] * num_sequences + idx[0]] = val

                    # Update Shared Stats Based on Loaded Data
                    prev_computed = [global_matrix_value_list[i * num_sequences + j] 
                                    for i in range(num_sequences) for j in range(i + 1, num_sequences) 
                                    if not np.isnan(global_matrix_value_list[i * num_sequences + j])]
                    
                    sum_value.value = sum(prev_computed)
                    sum_sq_value.value = sum([x**2 for x in prev_computed])
                    n_value.value = len(prev_computed)

                    print(f"Loaded {n_value.value} precomputed distances.")
                else:
                    pass # data_path will be used to save results later
        except Exception as e:
            traceback.print_exc()
            print(f"Error loading data: {e}")
            
            # return global_matrix_value_list, sum_value, sum_sq_value, n_value
        
    # Main ========================================
    global lock, counter, sum_value, sum_sq_value, n_value, stop_flag, prev_means, effective_rate, start_time, global_matrix_value_list

    # 🔹 Initialize Shared Memory Resources
    lock = multiprocessing.Lock()
    counter = multiprocessing.Value('i', 0)
    sum_value = multiprocessing.Value('d', 0.0)
    sum_sq_value = multiprocessing.Value('d', 0.0)
    n_value = multiprocessing.Value('i', 0)
    stop_flag = multiprocessing.Value('i', 0)
    prev_means = multiprocessing.Manager().list()
    effective_rate = multiprocessing.Value('f', 0.0)
    start_time = time.time()
    
    num_sequences = len(sequence_stacks)
    global_matrix_value_list = multiprocessing.Array('d', num_sequences * num_sequences)
    global_matrix_value_list[:] = [np.nan] * (num_sequences * num_sequences)

    # 🔹 Load Previously Computed Data from temp files (if available) - with larger datasets this function can take a long time to run - saving partial computations helps
    if data_path is not None: try_load_data(data_path)
    
    # 🔹 Generate Pair Combinations & Shuffle for Random Sampling
    total_pairs = [(i, j) for i in range(num_sequences) for j in range(i + 1, num_sequences)]
    random.shuffle(total_pairs)
    
    # HACK - I should set this batch size in run script - fix later
    if batch_size is None:
        batch_size = 1000
    if num_workers is None:
        num_workers = min(multiprocessing.cpu_count(), 50, batch_size)  # Limit to 50 workers
        
    pair_batches = [total_pairs[i:i + batch_size] for i in range(0, len(total_pairs), batch_size)]

    print(f"Using {num_workers} workers with batch size {batch_size}.")

    # 🔹 Run Multiprocessing Pool
    
    if len(pair_batches) > 1:
        print(f"Running {len(pair_batches)} batches across {num_workers} workers.")
        with multiprocessing.Pool(
            processes=num_workers, initializer=init_worker, 
            initargs=(lock, counter, sum_value, sum_sq_value, n_value, stop_flag, prev_means, effective_rate, start_time, global_matrix_value_list)
        ) as pool:
            pool.starmap(compute_dtw_dynamic_v2, [(batch, sequence_stacks, data_path,
                                                cv_threshold, moving_avg_window, moving_avg_threshold, confidence_threshold,                                               
                                                ) for batch in pair_batches])

        # 🔹 Compute Final Statistics
        mean_dtw = sum_value.value / n_value.value
        variance_dtw = (sum_sq_value.value / n_value.value) - mean_dtw ** 2
        std_dtw = np.sqrt(variance_dtw)
        cov_dtw = std_dtw / mean_dtw if mean_dtw > 0 and std_dtw > 0 else 0 # Coefficient of Variation
        median_dtw = np.median(global_matrix_value_list)
        min_dtw = np.min(global_matrix_value_list)
        max_dtw = np.max(global_matrix_value_list)
        
        # 🔹 Convert Computed Distances into a Matrix
        global_matrix = np.array(global_matrix_value_list).reshape(num_sequences, num_sequences)
        
        # if value is np.nan and in diagonal position, replace with 0
        for i in range(num_sequences):
            global_matrix[i, i] = 0

        print(f"Final DTW Analysis:")
        print(f"Mean: {mean_dtw:.3f} | Std: {std_dtw:.3f} | Variance: {variance_dtw:.3f} | Computations: {n_value.value}")
    else:
        print(f"Only one batch to process. Cannot compute dtw_distance_matrix with only one batch.")
        
        # 🔹 Compute Final Statistics
        # mean_dtw = sum_value.value / n_value.value
        # variance_dtw = (sum_sq_value.value / n_value.value) - mean_dtw ** 2
        # std_dtw = np.sqrt(variance_dtw)
        # cov_dtw = std_dtw / mean_dtw if mean_dtw > 0 and std_dtw > 0 else 0
        
        # Final Statistics
        mean_dtw = np.nan
        std_dtw = np.nan
        variance_dtw = np.nan
        cov_dtw = np.nan
        n_value.value = 0       
        
        # 🔹 Convert Computed Distances into a Matrix
        global_matrix = np.zeros((num_sequences, num_sequences))
        
        # print
        print(f"Final DTW Analysis:")
        print(f"Mean: {mean_dtw:.3f} | Std: {std_dtw:.3f} | Variance: {variance_dtw:.3f} | Computations: {n_value.value}")
    
    # 🔹 Save Results
    results_dict = {
        'dtw_distance_matrix': global_matrix,
        'mean_dtw': mean_dtw,
        'std_dtw': std_dtw,
        'variance_dtw': variance_dtw,
        'cov_dtw': cov_dtw,
        'median_dtw': np.median(global_matrix),
        'min_dtw': np.min(global_matrix),
        'max_dtw': np.max(global_matrix),

        'num_computations': n_value.value,
        #'effective_rate': effective_rate.value
    }

    #np.save(data_path + 'dtw_analysis_dynamic_results.npy', results_dict)
    npy_path = os.path.join(data_path, 'dtw_analysis_dynamic_results.npy')
    np.save(npy_path, results_dict)
    print(f"Results saved to {npy_path}")
    
    # Save global matrix to file
    csv_path = os.path.join(data_path, 'dtw_analysis_dynamic_results.csv')
    pd.DataFrame(global_matrix).to_csv(csv_path)
    #pd.DataFrame(global_matrix).to_csv(data_path + 'dtw_analysis_dynamic_results.csv')

    return results_dict
    
# =================== #
# old funcs:
    # def compute_dtw_dynamic(
    #     pairs_sample, sequence_stacks, 
    #     #lock, 
    #     counter, total, 
    #     sum_value, sum_sq_value, n_value, stop_flag, prev_means,
    #     effective_rate, start_time, 
    #     global_matrix_value_list,
    #     cv_threshold = 0.05, 
    #     moving_avg_window = 500,                        
    #     confidence_threshold=0.99, 
    #     cov_threshold=0.05
    #     ):
        
    #     global lock  # Use the global lock
        
    #     """Computes DTW distances and updates running statistics in parallel."""
            
    #     proc_start_time = time.time()
    #     #with lock:
    #     print(f"Worker started processing at {time.strftime('%H:%M:%S', time.gmtime(proc_start_time))}")
        
    #     # local_sum = 0
    #     # local_sum_sq = 0
    #     # local_count = 0
    #     #local_matrix = np.zeros((len(sequence_stacks), len(sequence_stacks)))
    #     #local_matrix = np.nans((len(sequence_stacks), len(sequence_stacks)))
        
    #     for i, j in pairs_sample:
            
    #         # if global_matrix_list_value already has a value, skip computation
    #         if not np.isnan(global_matrix_value_list[i*len(sequence_stacks) + j]) and not np.isnan(global_matrix_value_list[j*len(sequence_stacks) + i]):
    #             continue
            
    #         seq_i = sequence_stacks[i]
    #         seq_j = sequence_stacks[j]
    #         distance, _ = fastdtw(seq_i, seq_j, dist=euclidean)
    #         # local_matrix[i, j] = distance
    #         # local_matrix[j, i] = distance  # Symmetric property
            
    #         # Update local statistics
    #         # local_sum += distance
    #         # local_sum_sq += distance**2
    #         # local_count += 1

    #         with lock:
    #             # update global matrix list with distance
    #             global_matrix_value_list[i*len(sequence_stacks) + j] = distance
    #             global_matrix_value_list[j*len(sequence_stacks) + i] = distance            
                
    #             # Update shared statistics
    #             #sum_value.value += local_sum
    #             #sum_sq_value.value += local_sum_sq
    #             #n_value.value += local_count
    #             sum_value.value += distance
    #             sum_sq_value.value += distance**2
    #             n_value.value += 1
    #             counter.value += 1
    #             mean = sum_value.value / n_value.value
                
    #             # NOTE: Moving Average Stabilization doesnt have memory of previous means on a new run - idk if this is a problem
    #             prev_means.append(mean)
    #             if len(prev_means) > moving_avg_window:
    #                 prev_means.pop(0)
                
    #             # Compute running statistics
    #             if n_value.value > 1:
    #                 # Print progress
    #                 if counter.value % 100 == 0 or n_value.value == total:
    #                     #mean = sum_value.value / n_value.value
    #                     variance = (sum_sq_value.value / n_value.value) - mean**2
    #                     std_dev = np.sqrt(variance) if variance > 0 else 0

    #                     # Confidence Interval Calculation (99%)
    #                     # Z = 2.576  # Critical Z-score for 99% confidence
    #                     # Critical Z-score for 95% confidence: 1.96
    #                     #Z = 1.96
    #                     #margin_error = Z * (std_dev / np.sqrt(n_value.value))
                        
    #                     if mean > 0 and std_dev > 0:
    #                         cov = std_dev / mean # Coefficient of Variation
                            
    #                     # Moving Average Stabilization
    #                     if len(prev_means) > 1:
    #                         moving_average_diff = abs(prev_means[-1] - prev_means[0])
    #                         moving_threshold = 0.01 * prev_means[-1]  # 1% stability
    #                         #50% stability for quick testing
    #                         #moving_threshold = 0.5 * prev_means[-1]
    #                     else:
    #                         moving_average_diff = 0
    #                         moving_threshold = 0
                            
    #                     effective_rate.value = counter.value / (time.time() - start_time)
    #                     estimated_time_remaining = (total - counter.value) / effective_rate.value
    #                     estimated_time_remaining = time.strftime('%H:%M:%S', time.gmtime(estimated_time_remaining))
                        
    #                     #print(f"Progress: {counter.value}/{total}"
    #                     print(f"Progress: {n_value.value}/{total}"
    #                             f" | Mean: {mean:.3f} | Std: {std_dev:.3f}"
    #                             f" | CoV: {cov:.3f}"
    #                             f" | Moving Avg Diff: {moving_average_diff:.3f}"
    #                             f" | Threshold: {moving_threshold:.3f}"
    #                             #f" | Margin Error: {margin_error:.3f}"
    #                             #f" | Confidence: {1 - margin_error:.3f}"
    #                             f" | Effective Rate: {effective_rate.value:.3f} computations/s"
    #                             f" | ETA: {estimated_time_remaining}")

    #                     # # **Break Condition**: Stop if confidence is reached
    #                     # if margin_error < confidence_threshold:
    #                     #     stop_flag.value = 1
    #                     #     break  # Stop further processing in this worker
                        
    #                     # # **Break Condition**: Stop if coefficient of variation is reached
    #                     # if cov < cov_threshold:
    #                     #     stop_flag.value = 1
    #                     #     break  # Stop further processing in this worker
                        
    #                     # **Check Coefficient of Variation (CV)**
    #                     if mean > 0 and std_dev > 0:
    #                         #cv = std_dev / mean  # CV = σ / μ
    #                         if cov < cv_threshold:
    #                             # print(f"CV Stabilized: {cov:.3f} < {cv_threshold:.3f}")
    #                             switch = True
    #                             if switch:
    #                                 print(f"\tCV Stabilized: {cov:.3f} < {cv_threshold:.3f}. "
    #                                       "CV has stabilized. "
    #                                       "This stopping condition is enabled. "
    #                                       "Stopping computation.")
    #                                 stop_flag.value = 1
    #                                 break
    #                             else:
    #                                 print(f"\tCV Stabilized: {cov:.3f} < {cv_threshold:.3f}. "
    #                                         "CV has stabilized. "
    #                                         "This stopping condition is disabled. "
    #                                         "Continuing computation.")
                                    
    #                     # **Check Moving Average Stabilization**
    #                     #prev_means.append(mean)
    #                     #if len(prev_means) > moving_avg_window:
    #                         #prev_means.pop(0)  # Keep window size fixed
    #                         #if abs(prev_means[-1] - prev_means[0]) < 0.01 * prev_means[-1]:  # 1% stability
    #                     if moving_average_diff < moving_threshold:
    #                         # print(f"Moving Average Stabilized: {moving_average_diff:.3f} < {moving_threshold:.3f}")
    #                         switch = False
    #                         if switch:
    #                             print(f"\tMoving Average Stabilized: {moving_average_diff:.3f} < {moving_threshold:.3f}. " 
    #                                   "Moving average has stabilized. "
    #                                   "This stopping condition is enabled. "
    #                                   "Stopping computation.")
    #                             stop_flag.value = 1
    #                             break
    #                         else:
    #                             print(f"\tMoving Average Stabilized: {moving_average_diff:.3f} < {moving_threshold:.3f}. "
    #                                   "Moving average has stabilized. "
    #                                   "This stopping condition is disabled. "
    #                                   "Continuing computation.")
        
    #     proc_time_taken = time.time() - proc_start_time
    #     #format time taken
    #     proc_time_taken = time.strftime('%H:%M:%S', time.gmtime(proc_time_taken))
    #     print(f"Worker finished processing in {proc_time_taken}")
    #     #return None  # No need to return results since statistics are shared globally
        
    #     #for debug, print moving average window just to see what it looks like
    #     #print(f"Moving Average Window: {prev_means}")
    #     #print(f"Length of Moving Average Window: {len(prev_means)}")
    #     #return local_matrix

    # def dtw_analysis_dynamic(sequence_stacks, data_path=None,
    #                          #confidence_threshold=0.99,
    #                          cv_threshold=0.05, 
    #                          min_samples=500,
    #                          moving_avg_window=500,
    #                          ):
    #     """Runs DTW analysis with real-time confidence-based stopping."""

    #     # init lock object and shared variables
    #     global lock # Use the global lock
    #     lock = multiprocessing.Lock()  # Lock for updating counter safely
    #     effective_rate = multiprocessing.Value('f', 0.0)
    #     start_time = time.time()
    #     counter = multiprocessing.Value('i', 0)
    #     sum_value = multiprocessing.Value('d', 0.0)
    #     sum_sq_value = multiprocessing.Value('d', 0.0)
    #     n_value = multiprocessing.Value('i', 0)
    #     stop_flag = multiprocessing.Value('i', 0)  # Shared flag to indicate stopping
    #     prev_means = multiprocessing.Manager().list()  # Shared list for moving average
    #     #initializer_matrix = np.zeros((num_sequences, num_sequences))
    #     #global_matrix = multiprocessing.Array('d', initializer_matrix)
        
    #     # init global distance matrix - populate with nan values
    #     num_sequences = len(sequence_stacks)
    #     global_matrix_value_list = multiprocessing.Array('d', num_sequences*num_sequences)
    #     global_matrix_value_list[:] = [np.nan for i in range(num_sequences*num_sequences)]
        
    #     #add zeros to diagonal values, this is the distance between the same sequence - so we can get it out of the way.
    #     # also this provides structure to the initial matrix
    #     for i in range(num_sequences):
    #         global_matrix_value_list[i*num_sequences + i] = 0    
        
    #     # if data_path is not None, load pandas dataframe
    #     if data_path is not None:
    #         # load data
    #         import pandas as pd
    #         df = pd.read_csv(data_path + 'dtw_analysis_dynamic_results.csv')
    #         global_matrix = df.to_numpy()
            
    #         # remove first column from global_matrix, these are just the indices
    #         global_matrix = global_matrix[:, 1:]
            
    #         # load data onto global_matrix_value_list
    #         print(f"Loading computation history matrix from file: {global_matrix.shape}")
    #         for idx, val in np.ndenumerate(global_matrix):
    #             if np.isnan(val):
    #                 continue # dont waste time on nans
    #             #print(f"Index: {idx} | Value: {val}")
    #             global_matrix_value_list[idx[0]*num_sequences + idx[1]] = val
    #             global_matrix_value_list[idx[1]*num_sequences + idx[0]] = val
                
    #         # reformat global_matrix_value_list for calcs
    #         # remove recieprocal values and remove diagonal values
    #         print(f"Updating global matrix statistics from loaded data...")
    #         true_sample_values = []
    #         for i in range(num_sequences):
    #             for j in range(i+1, num_sequences):
    #                 if i == j:
    #                     continue # skip diagonal values
    #                 if i > j:
    #                     continue # skip reciprocal values
    #                 if np.isnan(global_matrix_value_list[i*num_sequences + j]):
    #                     continue                
    #                 true_sample_values.append(global_matrix_value_list[i*num_sequences + j])
            
    #         # Update shared statistics based on global matrix list
    #         # sum_value.value = np.nansum(global_matrix_value_list)
    #         # sum_sq_value.value = np.nansum(np.square(global_matrix_value_list))
    #         # n_value.value = np.count_nonzero(~np.isnan(global_matrix_value_list))
    #         sum_value.value = sum(true_sample_values)
    #         sum_sq_value.value = sum([x**2 for x in true_sample_values])
    #         n_value.value = len(true_sample_values)
    #         #counter.value = n_value.value
    #         mean = sum_value.value / n_value.value
            
    #         # current confidence stats
    #         variance = (sum_sq_value.value / n_value.value) - mean**2
    #         std_dev = np.sqrt(variance) if variance > 0 else 0
    #         cov = std_dev / mean if mean > 0 and std_dev > 0 else 0
    #         print(f"Current Confidence Stats:")
    #         print(f"Mean: {mean:.3f} | Std: {std_dev:.3f} | CoV: {cov:.3f}")        
    #         print(f"Total Computations Loaded: {n_value.value}")       
            
    #         print(f"Loaded computation history matrix from file: {global_matrix.shape}")
            
    #         # some debuggery to check if values are being loaded correctly
    #         # import sys
    #         # #print(f'true_sample_values: {true_sample_values}')
    #         # #print(f'len(true_sample_values): {len(true_sample_values)}')
    #         # #print(f'sum_value: {sum_value.value}')
    #         # #print(f'sum_sq_value: {sum_sq_value.value}')
    #         # #print(f'n_value: {n_value.value}')
    #         # #print(f'counter: {counter.value}')
    #         # #print(f'mean: {mean}')
    #         # sys.exit()
        
    #     # Run DTW Analysis
    #     num_sequences = len(sequence_stacks)
    #     total_pairs = [(i, j) for i in range(num_sequences) for j in range(i + 1, num_sequences)]
    #     print(f"Total possible pairs: {len(total_pairs)}")

    #     # Shuffle pairs for better randomness in early sampling
    #     random.shuffle(total_pairs)

    #     #num_workers = min(multiprocessing.cpu_count(), 50)  # Limit to 50 CPUs
    #     # use all possible cpus. 1 worker per cpu
    #     num_workers = multiprocessing.cpu_count()
    #     #batch_size = len(total_pairs) // (num_workers*2)
    #     #num_workers = 50
    #     batch_size = 1000
    #     pair_batches = [total_pairs[i:i + batch_size] for i in range(0, len(total_pairs), batch_size)]

    #     print(f"Using {num_workers} workers with batch size {batch_size}.")
    #     #
    #     # import sys
    #     # sys.exit()

    #     # Run parallel computation
    #     results = []
    #     with multiprocessing.Pool(processes=num_workers) as pool:
            
    #         # # aw 2025-02-23 16:56:44 - apprently this is running in serial, not parallel
    #         # for pair_batch in pair_batches:
    #         #     distance_matrix = compute_dtw_dynamic(
    #         #         pair_batch, sequence_stacks, lock, counter, len(total_pairs),
    #         #         sum_value, sum_sq_value, n_value, stop_flag, prev_means,
    #         #         effective_rate, start_time, global_matrix_value_list,
    #         #         #cv_threshold, 
    #         #         #moving_avg_window, 
    #         #         #confidence_threshold
    #         #         )
    #         #     results.append(distance_matrix)
    #         #     if stop_flag.value:
    #         #         break
            
    #         pool.starmap(
    #             compute_dtw_dynamic,
    #             [(batch, sequence_stacks, 
    #               #lock, 
    #               counter, len(total_pairs),
    #             sum_value, sum_sq_value, n_value, stop_flag, prev_means,
    #             effective_rate, start_time, global_matrix_value_list)
    #             for batch in pair_batches]
    #         )


    #     # Compute final statistics
    #     mean_dtw = sum_value.value / n_value.value
    #     variance_dtw = (sum_sq_value.value / n_value.value) - mean_dtw**2
    #     std_dtw = np.sqrt(variance_dtw)

    #     print(f"\nFinal DTW Analysis (Dynamic Stopping):")
    #     print(f"Mean DTW Distance: {mean_dtw:.3f}")
    #     print(f"Standard Deviation: {std_dtw:.3f}")
    #     print(f"Variance: {variance_dtw:.3f}")
    #     print(f"Total Computations: {n_value.value}")
        
    #     # replace zeros in distance matrix with nan, then flatten matrices - they all represent differnt coordinates of the same matrix
    #     # results = [np.where(matrix == 0, np.nan, matrix).flatten() for matrix in results]
    #     # results = np.concatenate(results)
    #     #global_matrix = np.where(results == 0, np.nan, results)
        
    #     #convert global_matrix_list to matrix with 2 dimensions
    #     global_matrix = np.array(global_matrix_value_list).reshape(num_sequences, num_sequences)
    #     print(f"Global Matrix Shape: {global_matrix.shape}")
        
    #     results_dict = {
    #         'mean_dtw': mean_dtw,
    #         'std_dtw': std_dtw,
    #         'variance_dtw': variance_dtw,
    #         'global_matrix': global_matrix,
    #         'num_computations': n_value.value,
    #         'effective_rate': effective_rate.value,
    #         #'confidence_threshold': confidence_threshold        
    #     }
        
    #     # save results to file
    #     path = '/pscratch/sd/a/adammwea/workspace/RBS_network_models/tests/dtw_test_data/'
    #     np.save(path + 'dtw_analysis_dynamic_results.npy', results_dict)   
        
    #     #save matrix to pandas csv so it can be read in as a dataframe    
    #     import pandas as pd
    #     df = pd.DataFrame(global_matrix)
    #     df.to_csv(path + 'dtw_analysis_dynamic_results.csv')
        
    #     return mean_dtw, std_dtw, variance_dtw, global_matrix

    # def compute_dtw_batch_3(pair_batch, sequence_stacks, lock, counter, total, effective_rate, start_time):
    #     """Computes DTW distances for a batch of pairs (vectorized approach)."""
    #     results = []
        
    #     # loop through pairs
    #     for i, j in pair_batch:
    #         # seq_i = np.array(time_sequence_dict[i])
    #         # seq_j = np.array(time_sequence_dict[j])
            
    #         # # Downsample for speed (adjust this if accuracy is a concern)
    #         # aw 2025-02-22 22:39:34
    #         # apparently fastdtw already downsamples by half internally. That's why it's so fast lol.
    #         # accuracy is... vaugely concerning, but not really since we're just trying to decipher a 
    #         # general trend to tell simulations to target.
    #         # Just going to try downsampling by factor of 4 to see if it speeds things up.
            
    #         #distance_test, _ = fastdtw(seq_i, seq_j, dist=euclidean)  # Compute DTW

            
    #         # seq_i = downsample_sequence(seq_i, downsample_factor=4)
    #         # seq_j = downsample_sequence(seq_j, downsample_factor=4)
            
    #         seq_i = sequence_stacks[i]
    #         seq_j = sequence_stacks[j]

    #         distance, _ = fastdtw(seq_i, seq_j, dist=euclidean)  # Compute DTW
            
    #         results.append((i, j, distance))

    #         with lock:  # Safely update progress counter
    #             counter.value += 1
    #             effective_rate.value = counter.value / (time.time() - start_time)
    #             # if counter.value % 25 == 0 or counter.value == total:
    #             #     #print(f"Progress: {counter.value}/{total} distances computed")
    #             #     print(f"EFFECTIVE RATE: {effective_rate.value} distances computed")
    #             #     estimated_time_remaining = (total - counter.value) / effective_rate.value
    #             #     # convert to minutes
    #             #     estimated_time_remaining = estimated_time_remaining / 60
    #             #     print(f"Estimated time remaining: {estimated_time_remaining} minutes")
    #             if counter.value % 100 == 0 or counter.value == total:  
    #                 #print(f"Progress: {counter.value}/{total} distances computed")
    #                 #put all messages together on one line
    #                 estimated_time_remaining = (total - counter.value) / effective_rate.value
    #                 # convert to HH:MM format
    #                 estimated_time_remaining = time.strftime('%H:%M', time.gmtime(estimated_time_remaining))
    #                 print(f"Progress: {counter.value}/{total} distances computed | "
    #                       f"Effective Rate: {effective_rate.value} computations/s | "
    #                       f"Estimated time remaining: {estimated_time_remaining}")
        
    #     return results

    # def dtw_analysis_v3(time_sequence_dict, cat_sequence_dict, sequence_stacks, bursting_data):
        
    #     ## =================== ##
    #     ## save test data for easier debugging
    #     # path = '/pscratch/sd/a/adammwea/workspace/RBS_network_models/tests/dtw_test_data/'
    #     # time_sequence_dict = {i: time_sequence_mat[i] for i in range(len(time_sequence_mat))}
    #     # cat_sequence_dict = {i: cat_sequence_mat[i] for i in range(len(cat_sequence_mat))}
    #     # np.save(path + 'time_sequence_dict.npy', time_sequence_dict)
    #     # np.save(path + 'cat_sequence_dict.npy', cat_sequence_dict)
    #     # np.save(path + 'sequence_stacks.npy', sequence_stacks)
    #     # np.save(path + 'bursting_data.npy', bursting_data)
        
    #     ## =================== ##
    #     ## DTW Analysis    
        
    #     # print('Length of cat_sequence_mat:', len(cat_sequence_mat))
    #     # print('Length of time_sequence_mat:', len(time_sequence_mat))
    #     print('Length of cat_sequence_dict:', len(cat_sequence_dict))
    #     print('Length of time_sequence_dict:', len(time_sequence_dict))

    #     #num_sequences = len(time_sequence_mat)
    #     num_sequences = len(time_sequence_dict)
    #     distance_matrix = np.zeros((num_sequences, num_sequences))  # Initialize matrix

    #     # Create a list of unique index pairs (i, j)
    #     pairs = [(i, j) for i in range(num_sequences) for j in range(i + 1, num_sequences)]
    #     total_distances = len(pairs)
    #     print('total_distances to compute:', total_distances)
        
    #     # Setup multiprocessing resources
    #     counter = multiprocessing.Value('i', 0)  # Shared counter for tracking progress
    #     lock = multiprocessing.Lock()  # Lock for updating counter safely
    #     start_time = time.time()
    #     effective_rate = multiprocessing.Value('f', 0.0)
        
    #     # # one dtw at a time for now
    #     # for i in range(num_sequences):
    #     #     for j in range(i + 1, num_sequences):
    #     #         # seq_i = np.array(time_sequence_dict[i])
    #     #         # seq_j = np.array(time_sequence_dict[j])
                
    #     #         # correct way to do it - sensitive to spike times and i/e categorization
    #     #         seq_i = sequence_stacks[i]
    #     #         seq_j = sequence_stacks[j]
    #     #         distance, _ = fastdtw(seq_i, seq_j, dist=euclidean)
    #     #         print(f'Shape of seq_i: {seq_i.shape}')
    #     #         print(f'Shape of seq_j: {seq_j.shape}')
    #     #         print(f'DTW distance between {i} and {j}: {distance}')
                
    #     #         # test - time data alone
    #     #         # # seq_i_test = time_sequence_dict[i]
    #     #         # # seq_j_test = time_sequence_dict[j]
    #     #         # seq_i_test = np.array(time_sequence_dict[i])
    #     #         # seq_j_test = np.array(time_sequence_dict[j])
    #     #         # distance_test, _ = fastdtw(seq_i_test, seq_j_test, dist=euclidean)
    #     #         # print(f'Shape of seq_i_test: {seq_i_test.shape}')
    #     #         # print(f'Shape of seq_j_test: {seq_j_test.shape}')
    #     #         # print(f'DTW distance between {i} and {j}: {distance_test}')
                
    #     #         distance_matrix[i, j] = distance
    #     #         distance_matrix[j, i] = distance  # Symmetric property
    #     #         with lock:  # Safely update progress counter
    #     #             counter.value += 1
    #     #             if counter.value % 100 == 0 or counter.value == total_distances:
    #     #                 print(f"Progress: {counter.value}/{total_distances} distances computed")
                        
    #     # parallelized dtw and batched
    #     # Batch the pairs to reduce overhead
    #     #batch_size = 500  # Adjust this for better performance (100-500 is a good range)
    #     #batch_size = 1000
        
    #     #num cpus byos count
    #     import os
    #     num_cpus = os.cpu_count()
    #     print(f'Number of cpus by os: {num_cpus}')
        
    #     #num_workers = 50
    #     # use all possible cpus. 1 worker per cpu
    #     num_workers = multiprocessing.cpu_count()
    #     print(f'Number of cpus: {multiprocessing.cpu_count()}')
    #     print(f'Number of workers: {num_workers}')
        
    #     import sys
    #     sys.exit()
    #     # override test
    #     # num_workers = 256
    #     # print(f'Number of workers (overridden): {num_workers}')
        
    #     # split pairs into 1 batch per worker
    #     # split batches squarely in to num_workers-1 batches
    #     # then use the last worker-batch to fill in the rest
    #     batch_size = len(pairs)//(num_workers)
    #     #pair_batches = [pairs[i:i + batch_size] for i in range(0, len(pairs), batch_size)]
        
    #     # randomize pairs for better approximation of time
    #     np.random.shuffle(pairs)
    #     pair_batches = [pairs[i:i + batch_size] for i in range(0, len(pairs), batch_size)]    

    #     #num_workers = min(multiprocessing.cpu_count(), len(pair_batches))  # Limit workers to CPU count
    #     #num_workers = 50
    #     with multiprocessing.Pool(processes=num_workers) as pool:
    #         # Use partial function to pass additional arguments
    #         results_batches = [compute_dtw_batch_3(pair_batch, sequence_stacks, 
    #                                                lock, counter, total_distances, 
    #                                                effective_rate, start_time) for pair_batch in pair_batches]
    #         #results_batches = [func, pair_batch for pair_batch in pair_batches]
            
    #     # try multi threading instead
        
        
        
            
    #     # Flatten results and store them in the distance matrix
    #     for results in results_batches:
    #         for i, j, distance in results:
    #             distance_matrix[i, j] = distance
    #             distance_matrix[j, i] = distance  # Symmetric property

    #     print("DTW distance matrix computation complete!")
    #     return distance_matrix
