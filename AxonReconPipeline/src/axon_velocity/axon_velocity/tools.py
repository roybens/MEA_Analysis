import numpy as np
from scipy.stats import linregress


def compute_peak_time_stds(template, locations, fs, neighbor_distance=30):
    peak_times = np.argmin(template, 1) / fs
    peak_time_stds = np.zeros(len(locations))
    x, y = locations.T
    for i in np.arange(len(locations)):
        neighboring_chs_idxs = np.where((x > x[i] - neighbor_distance) & (x < x[i] + neighbor_distance)
                                        & (y > y[i] - neighbor_distance) & (y < y[i] + neighbor_distance))[0]
        if len(neighboring_chs_idxs) > 1:
            peak_time_stds[i] = np.std(peak_times[neighboring_chs_idxs])

    return peak_time_stds


def compute_axon_velocity_on_branches(template, gtr_original):
    """
    Computes axon velocities on the same branches as 'gtr_original'

    Parameters
    ----------
    template
    gtr_original

    Returns
    -------

    """
    from .tracking_classes import AxonTracking
    fs = gtr_original.fs / gtr_original._upsample
    gtr_new = AxonTracking(template, gtr_original.locations, fs, upsample=gtr_original._upsample)

    # Create branch dictionary
    branches = []
    gtr_new._paths_clean = []
    for p_i, branch in enumerate(gtr_original.branches):
        path_rev = branch["channels"]
        peaks = (gtr_new.peak_times[path_rev] - gtr_new.peak_times[path_rev[0]]) / gtr_new.fs * 1000
        dists = []
        cum_dist = 0
        for i, p in enumerate(path_rev):
            if i == 0:
                cum_dist += 0
            else:
                cum_dist += np.linalg.norm(gtr_new.locations[p] - gtr_new.locations[path_rev[i - 1]])
            dists.append(cum_dist)
        peaks = np.array(peaks)
        dists = np.array(dists)

        velocity, offset, r_value, p_value, std_err = linregress(peaks, dists)
        r2 = r_value ** 2

        gtr_new._paths_clean.append(path_rev)
        branch_dict = {}
        branch_dict['channels'] = path_rev.astype(int)
        branch_dict['velocity'] = velocity
        branch_dict['offset'] = offset
        branch_dict['r2'] = r2
        branch_dict['pval'] = p_value
        branch_dict['distances'] = np.array(dists)
        branch_dict['peak_times'] = np.array(peaks)
        branch_dict['raw_path_idx'] = p_i
        # branch_dict['avg_heur'] = gtr._compute_path_avg_heur(path)
        branches.append(branch_dict)

    gtr_new.branches = branches

    return gtr_new


def distance_numpy(A, B, P):
    """ segment line AB, point P, where each one is an array([x, y]) """
    if np.linalg.norm(P - A) < 1e-5 or np.linalg.norm(P - B) < 1e-5:
        return 0
    if np.arccos(np.round(np.dot((P - A) / np.linalg.norm(P - A), (B - A) / np.linalg.norm(B - A)), 5)) > np.pi / 2:
        return np.linalg.norm(P - A)
    if np.arccos(np.round(np.dot((P - B) / np.linalg.norm(P - B), (A - B) / np.linalg.norm(A - B)), 5)) > np.pi / 2:
        return np.linalg.norm(P - B)
    # return np.linalg.norm(np.cross(A - B, A - P)) / np.linalg.norm(B - A)
    theta = np.arccos(np.round(np.dot((P - A) / np.linalg.norm(P - A), (B - A) / np.linalg.norm(B - A)), 5))
    return np.abs(np.linalg.norm(P - A) * np.sin(theta))

