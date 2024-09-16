from copy import deepcopy

from .tracking_classes import GraphAxonTracking

def compute_graph_propagation_velocity(template, locations, fs, upsample=1, init_delay=0,
                                       detect_threshold=0.1, kurt_threshold=0.3,
                                       peak_std_threshold=None, peak_std_distance=30, remove_isolated=True,
                                       detection_type="relative",  min_selected_points=30, min_path_points=5,
                                       mad_threshold=5, r2_threshold=0.9, max_distance_for_edge=300,
                                       max_distance_to_init=200, n_neighbors=3, min_path_length=100,
                                       init_amp_peak_ratio=0.2, edge_dist_amp_ratio=0.3,
                                       min_points_after_branching=3, max_peak_latency_for_splitting=1,
                                       theilsen_maxiter=2000, split_paths=True,
                                       distance_exp=2, neighbor_radius=100,
                                       r2_threshold_for_outliers=0.98, min_outlier_tracking_error=50,
                                       verbose=False, debug_mode=False):
    """
    Compute velocity of branches using a graph-based approach.
    After an initial channel selection based on amplitude and kurtosis, each channel make up a node. Starting from
    nodes with a late peak, a graph is build so that each node has edges to the closest channels whose peak comes
    earlier. Branches are found as the shortest path between nodes and the initial channel (the one with the absolute
    minimum). Finally, branches are cleand for duplicates.

    Parameters
    ----------
    template: np.array (num_channels x num_timepoints)
        The template to compute spread on
    locations: np.array (num_channels x 2)
        The channel locations
    fs: float
        Sampling frequency in Hz
    upsample: int
        Upsanpling factor (default 10)
    init_delay: int
        Delay in ms from starting frame to discard channels (default 0)
    min_selected_points: int
        Minimum number of channels to run tracking algorithm
    detect_threshold: float
        Detection threshold relative to max amplitude (default 0.1)
    kurt_threshold: float
        Threshold on kurtosis for single channel templates (default 0.3). The channels containing a spike are expected
        to be super-gaussian (kurtosis > 0)
    peak_std_threshold: float
        Threshold to remove channels based on peak time standard deviation
    peak_std_distance: float
        Distance to compute peak time std (default 30)
    remove_isolated: bool
        If True (default), isolated channels are removed
    detection_type: str
        "relative": detection threshold is relative to max amplitude
        "absolute": detection threshold is absolute
    r2_threshold: float
        Threshold on the r2 values of velocity estimation to keep a branch (default 0.9)
    mad_threshold: float
        Threshold on median absolute deviation on the fit error to consider points as outliers in the velocity
        estimation (default 5)
    max_distance_for_edge: float
        Maximum distance to create an edge between channels (same units as locations - default 300)
    max_distance_to_init: float
        Maximum distance to create an edge between channels and the init channel (same units as locations - default 200)
    n_neighbors: int
        Maximum number of neighbors (edges) for each node (default 3)
    min_path_length: float
        Minimum length of a path to keep the reconstructed branch (default 100)
    min_path_points: int
        Minimum number of channels to in an branch (default 5)
    min_points_after_branching: int
        Minimum number of points to keep a sub-path after a branching point
    max_peak_latency_for_splitting: float
        If there is a jump in the peak latency above this value, the branch is split and sub-branches
        are fitted separately
    init_amp_peak_ratio: float
        Ratio between amplitude and peak latency to sort channels for path search (default 0.2)
    edge_dist_amp_ratio: float
        Ratio between distance and amplitude to select nearest neighbors for edges
        (default 0.3 - i.e. proximity weighted 30%, amplitude 70%)
    distance_exp: int
        Exponent of distance to penalize long jumps (default 2)
    neighbor_radius: float
        Radius to exclude neighboring channels around an identified path (default 50)
    r2_threshold_for_outliers: float
        R2 value below which the algorithm looks for outliers in the fitting (default 0.98)
    min_outlier_tracking_error: float
        Minimum tracking error for a channel to be considered an outlier
    theilsen_maxiter: int
        Maximum number of iteration for TheilSen regressor
    split_paths: bool
        If True (defualt), branches can be split in the clean_paths phase to achieve a better velocoty fit
    verbose: bool
        If True the output is verbose
    debug_mode: bool
        If True, several plots are created to debug

    Returns
    -------
    gtr: GraphAxonTracking
        Output object containing the following fields:
        - branches
            List of dictionaries containing the following fields:
                - 'channels': selected channels in the path
                - 'velocity': velocity estimate in mm/s (if locations in um and fs in Hz)
                - 'offset': offset of velocity estimate (same units as locations)
                - 'r2': r-squared of the velocity fit
                - 'error': standard error of the linear fit
                - 'pval': p_value of the fit
                - 'distances': array with distances computed along the branch
                - 'peak_times': array with peak differences with initial channel
        - selected_channels
            List of selected channels
        - graph
            NetworkX directed graph
    """
    if debug_mode:
        print("\n\n##############################")
        print("######### DEBUG MODE #########")
        print("##############################\n\n")
        verbose = True

    gtr = GraphAxonTracking(template, locations, fs, min_selected_points=min_selected_points,
                            init_delay=init_delay,
                            detect_threshold=detect_threshold,
                            kurt_threshold=kurt_threshold,
                            peak_std_threshold=peak_std_threshold,
                            r2_threshold=r2_threshold,
                            remove_isolated=remove_isolated,
                            detection_type=detection_type,
                            max_distance_for_edge=max_distance_for_edge,
                            max_distance_to_init=max_distance_to_init,
                            mad_threshold=mad_threshold,
                            n_neighbors=n_neighbors,
                            upsample=upsample,
                            peak_std_distance=peak_std_distance,
                            min_path_length=min_path_length,
                            min_path_points=min_path_points,
                            min_points_after_branching=min_points_after_branching,
                            max_peak_latency_for_splitting=max_peak_latency_for_splitting,
                            distance_exp=distance_exp,
                            neighbor_radius=neighbor_radius,
                            r2_threshold_for_outliers=r2_threshold_for_outliers,
                            min_outlier_tracking_error=min_outlier_tracking_error,
                            init_amp_peak_ratio=init_amp_peak_ratio,
                            edge_dist_amp_ratio=edge_dist_amp_ratio,
                            theilsen_maxiter=theilsen_maxiter,
                            split_paths=split_paths,
                            verbose=verbose)
    gtr.track_axons()

    if debug_mode:
        _ = gtr.plot_channel_selection()
        _ = gtr.plot_graph()
        _ = gtr.plot_branches()
        _ = gtr.plot_velocities()

    return gtr


def get_default_graph_velocity_params():
    return deepcopy(GraphAxonTracking.default_params)
