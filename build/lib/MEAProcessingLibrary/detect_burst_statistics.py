import numpy as np
def detect_bursts_statistics(spike_times, isi_threshold):
    """
    Detect bursts in a spike train and calculate burst statistics.

    Parameters  
    ----------
    spike_times : dict  
        A dictionary of spike times for each unit.
    isi_threshold : float
        The threshold for detecting bursts (in seconds).

    Returns

    -------
    bursts : dict   
        A dictionary of burst times for each unit.
    
    """
# Dictionary to store results
    results = {}

    for unit, times in spike_times.items():
        
        # Step 1: Calculate the ISIs
        isis = np.diff(times)
        
        # Step 2: Identify where ISIs are below the threshold
        burst_mask = isis < isi_threshold
        
        # Step 3: Find the start and end indices of bursts using diff and boolean indexing
        burst_starts = np.where(np.diff(np.insert(burst_mask, 0, False).astype(int)) == 1)[0]
        burst_ends = np.where(np.diff(np.append(burst_mask, False).astype(int)) == -1)[0]
        
        # Step 4: Group spikes into bursts
        bursts = [times[start:end + 1] for start, end in zip(burst_starts, burst_ends)]
        
        # Extract ISIs within bursts using list comprehension and vectorized operations
        isis_within_bursts = np.concatenate([isis[start:end] for start, end in zip(burst_starts, burst_ends)]) if len(burst_starts) > 0 else np.array([])
        isis_outside_bursts = isis[~burst_mask]

        # Calculate statistics using vectorized numpy operations
        mean_isi_within = np.mean(isis_within_bursts) if isis_within_bursts.size > 0 else np.nan
        cov_isi_within = np.cov(isis_within_bursts) if isis_within_bursts.size > 0 else np.nan
        
        mean_isi_outside = np.mean(isis_outside_bursts) if isis_outside_bursts.size > 0 else np.nan
        cov_isi_outside = np.cov(isis_outside_bursts) if isis_outside_bursts.size > 0 else np.nan
        
        mean_isi_all = np.mean(isis) if isis.size > 0 else np.nan
        cov_isi_all = np.cov(isis) if isis.size > 0 else np.nan

        results[unit] = {
            "bursts": bursts,
            "mean_isi_within": mean_isi_within,
            "cov_isi_within": cov_isi_within,
            "mean_isi_outside": mean_isi_outside,
            "cov_isi_outside": cov_isi_outside,
            "mean_isi_all": mean_isi_all,
            "cov_isi_all": cov_isi_all,
            "isis_within_bursts": isis_within_bursts,
            "isis_outside_bursts": isis_outside_bursts,
            "isis_all": isis
        }

    return results