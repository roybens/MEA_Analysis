import numpy as np
from scipy.signal import resample_poly
from scipy.stats import linregress
from .tools import distance_numpy


def extract_ground_truth_velocity(cell, segment_sections=None, min_segs=3, min_length=50,
                                  upsample=10, min_vel=100):
    """
    Extract information about ground truth branches.

    Parameters
    ----------
    cell: LFPy.Cell
        The cell model
    segment_sections: array-like
        String array specifying the section name of each segment (if not present in the cell model)
    min_segs: int
        Minimum number of segments for an axonal branch
    min_length: float
        Minimum length in um for an axonal branch
    upsample: int
        Upsampling factor for membrane potential (default 10)

    Returns
    -------
    axon_branches: list
        List of dictionaries containing:
            - idxs: indexes of the cell model belonging to the axon
            - section: section name
            - velocity: ground-truth velocity
            - peaks: peak latency for in ms for each axonal segment
            - distances: distances for each axonal segment
    """
    fs = 1 / cell.dt * 1e3  # hz
    valid_axon_sections = []
    valid_axon_idxs = []

    if upsample > 1:
        vmem = resample_poly(cell.vmem, up=upsample, down=1, axis=1)
        fs_up = fs * upsample
    else:
        vmem = cell.vmem
        fs_up = fs

    axon_branches = []
    if segment_sections is None:
        assert cell.allseclist is not None, "Could not find section information!"
        for sec in cell.allseclist:
            if 'axon' in sec.name() and sec.nseg > min_segs:
                valid_axon_sections.append(sec.name())
                valid_axon_idxs.append(cell.get_idx(section=sec.name()))
    else:
        assert len(segment_sections) == len(cell.x), "'segment_sections' has doesn't have the same length as " \
                                                          "the number of segments!"
        axon_segs = []
        for seg_sec in segment_sections:
            if 'axon' in seg_sec:
                axon_segs.append(seg_sec)
        axon_segs = np.array(axon_segs)
        axon_sections, nsegs = np.unique(axon_segs, return_counts=True)
        for sec, count in zip(axon_sections, nsegs):
            if count > min_segs:
                valid_axon_sections.append(sec)
                valid_axon_idxs.append(np.where(segment_sections == sec)[0])

    for axon, axon_idxs in zip(valid_axon_sections, valid_axon_idxs):
        cum_distance = np.zeros(len(axon_idxs) - 1)
        peaks = np.zeros(len(axon_idxs) - 1)
        for i, (idx, idx_prev) in enumerate(zip(axon_idxs[1:], axon_idxs[:-1])):
            peak_time = np.argmax(vmem[idx]) / fs_up * 1000
            dist = np.linalg.norm(np.array([np.mean(cell.x[idx]), np.mean(cell.y[idx]), np.mean(cell.z[idx])]) -
                                  np.array([np.mean(cell.x[idx_prev]), np.mean(cell.y[idx_prev]),
                                            np.mean(cell.z[idx_prev])]))
            if i > 0:
                cum_distance[i] = cum_distance[i - 1] + dist
            else:
                cum_distance[i] = dist
            peaks[i] = peak_time

        if cum_distance[-1] > min_length:
            velocity, offset, r_value, p_value, std_err = linregress(peaks, cum_distance)
            if velocity > min_vel:
                axon_dict = {}
                axon_dict['idxs'] = axon_idxs
                axon_dict['section'] = axon
                axon_dict['velocity'] = velocity
                axon_dict['peaks'] = peaks
                axon_dict['distances'] = cum_distance
                axon_branches.append(axon_dict)
        else:
            print(f"axon {axon} is too short: {cum_distance[-1]} um")

    return axon_branches


def evaluate_tracking_accuracy(branches, gt_branches, cell, locations, max_median_dist_for_match=30,
                               min_dist_for_overlap=50, min_est_seg_for_match=3, max_overlap=0.2):
    """
    Match extracted branch with ground-truth axon and computes errors

    Parameters
    ----------
    branches: dict
        Extracted branches from tracking algorithm
    gt_branches: dict
        Ground-truth axons from 'extract_ground_truth_velocity'
    cell: LFPy.Cell
        The cell model
    locations: np.array
        x/y locations
    max_median_dist_for_match: float
        Maximum median distance for extracted branch to gt branch to consider a match
    min_dist_for_overlap: float
        Minimum distance to consider an overlapping branch (EXPLAIN)
    min_est_seg_for_match: int
        Minimum number of segments in an estimated path mathed to a gt branch for the gt branch to be considered
    max_overlap: float
        (EXPLAIN)


    Returns
    -------
    evaluation: list
        List of dictionaries with evaluation results for each matched axon containing:
            - branch_idx: index of matched estimated branch
            - gt_matches: list of ground truth branches matched to estimated branch
            - axon_idxs: list of indexes of matched ground-truth axons
            - velocity_ground_truth
            - velocity_estimated
            - errors: array with error with respect to closest ground-truth index
            - mean_error
            - std_error
    """
    evaluation = []
    matches = np.zeros((len(gt_branches), len(branches)))
    num_est_segments = np.zeros((len(gt_branches), len(branches)), dtype=int)

    cell_locations = np.array([[np.mean(x), np.mean(y)] for (x, y) in zip(cell.x, cell.y)])

    # find median distance between GT and est branches
    for i, (gt_branch) in enumerate(gt_branches):
        median_distances_to_estimated_branch = []
        est_segments_distance_ids = []
        for br, branch_dict in enumerate(branches):
            branch_channels = branch_dict['channels']
            dists = []
            ids_fo_gt_branch = []
            for idx in gt_branch['idxs']:
                seg_loc = cell_locations[idx]
                dists_to_estimated_branch = []
                for i_ch, ch in enumerate(branch_channels[:-1]):
                    loc_0 = locations[ch]
                    loc_1 = locations[branch_channels[i_ch + 1]]
                    if ch == branch_channels[i_ch + 1]:
                        raise Exception
                    d = distance_numpy(loc_0, loc_1, seg_loc)
                    dists_to_estimated_branch.append(d)
                ids_fo_gt_branch.append(np.argmin(dists_to_estimated_branch))
                dists.append(np.min(dists_to_estimated_branch))
            est_segments_distance_ids.append(len(np.unique(ids_fo_gt_branch)))
            median_distances_to_estimated_branch.append(np.median(dists))
        matches[i] = median_distances_to_estimated_branch
        num_est_segments[i] = est_segments_distance_ids

    # assign gt_branch to detected branch
    gt_matches = np.argmin(matches, axis=1)

    for i, gt_match in enumerate(gt_matches):
        if matches[i, gt_match] > max_median_dist_for_match or num_est_segments[i, gt_match] < min_est_seg_for_match:
            gt_matches[i] = -1

    # check that gt_matches are not overlapping branches (in case choose best one)
    matched_idxs = np.unique(gt_matches[gt_matches != -1])

    match_dict = {}
    for matched_idx in matched_idxs:
        matched_branches = np.nonzero(gt_matches == matched_idx)[0]

        print(f"Matched gt branches for {matched_idx}: {matched_branches}")

        if len(matched_branches) > 1:
            axon_locations = [cell_locations[gt_branch['idxs']]
                              for gt_branch in np.array(gt_branches)[matched_branches]]

            branches_to_delete = []
            for i1, loc_1 in enumerate(axon_locations):
                for i2, loc_2 in enumerate(axon_locations):
                    if i1 > i2:
                        shorter_axon_idx = np.argmin([len(loc_1), len(loc_2)])
                        if shorter_axon_idx == 0:
                            shorter_axon = loc_1
                            longer_axon = loc_2
                        else:
                            shorter_axon = loc_2
                            longer_axon = loc_1

                        n_overlaps = 0
                        for loc1 in shorter_axon:
                            dists = np.array([np.linalg.norm(loc1 - loc2) for loc2 in longer_axon])
                            min_dist_longer = np.min(dists)
                            if min_dist_longer < min_dist_for_overlap:
                                n_overlaps += 1
                        overlap = n_overlaps / len(shorter_axon)

                        if overlap > max_overlap:
                            # Remove shortest axon match
                            branch_dict = branches[matched_idx]
                            branch_channels = branch_dict['channels']

                            cum_dist1 = np.sum([np.min([np.linalg.norm(seg_loc - locations[ch])])
                                                        for ch in branch_channels
                                                for seg_loc in loc_1]) / len(loc_1)
                            cum_dist2 = np.sum([np.min([np.linalg.norm(seg_loc - locations[ch])])
                                                        for ch in branch_channels
                                                for seg_loc in loc_1]) / len(loc_2)
                            if cum_dist1 > cum_dist2:
                                branches_to_delete.append(matched_branches[i1])
                            else:
                                branches_to_delete.append(matched_branches[i2])

            if len(branches_to_delete) > 0:
                branches_to_delete = np.unique(branches_to_delete)
                for br in branches_to_delete:
                    matched_branches = np.delete(matched_branches, np.where(matched_branches == br))

            print(f"Matched branches for {matched_idx}: {matched_branches} - "
                  f"velocities {[gt_branches[m]['velocity'] for m in matched_branches]}")

        match_dict[matched_idx] = matched_branches

    for br, gt_match_idxs in match_dict.items():
        eval_dict = {}
        branch_dict = branches[br]
        channels = branch_dict["channels"]

        error_dist = []
        velocities = []
        lens_gt_branch = []
        axon_idxs = []

        for gt_match_idx in gt_match_idxs:
            gt_branch = gt_branches[gt_match_idx]
            gt_idxs = gt_branch["idxs"]
            axon_idxs.append(gt_idxs)

            lens_gt_branch.append(len(gt_branch['idxs']))
            velocities.append(gt_branch['velocity'])

        for ch in channels:
            channel_loc = locations[ch]
            dists = []
            for axon_idx in axon_idxs:
                dists += [np.linalg.norm(cell_locations[idx] - channel_loc) for idx in axon_idx]
            dist = np.min(dists)
            error_dist.append(dist)

        eval_dict['branch_idx'] = br
        eval_dict['gt_matches'] = gt_match_idxs
        eval_dict['axon_idxs'] = axon_idxs
        eval_dict['channels'] = branch_dict["channels"]
        eval_dict['velocity_estimated'] = branch_dict['velocity']
        eval_dict['velocity_ground_truth'] = np.sum([l * v for (l, v) in zip(lens_gt_branch, velocities)]) / \
                                             np.sum(lens_gt_branch)
        eval_dict['errors'] = error_dist
        eval_dict['mean_error'] = np.mean(error_dist)
        eval_dict['std_error'] = np.std(error_dist)
        evaluation.append(eval_dict)

    return evaluation
