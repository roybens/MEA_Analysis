import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import numpy as np
import MEAutility as mu
import probeinterface as pi
from probeinterface import plotting
from tqdm import tqdm

from .tools import compute_peak_time_stds


def plot_branch_neurites(branches, x, y, ax=None, cmap='rainbow', **kwargs):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    cm = plt.get_cmap(cmap)
    n_branches = len(branches)

    for i, branch in enumerate(branches):
        sel_idxs = branch['selected_channels']
        ax.plot(x[sel_idxs], y[sel_idxs], marker='o', ls='', color=cm(i / n_branches), **kwargs)

    return ax


def plot_branch_velocities(branches, ax=None, cmap='rainbow', alpha_marker=0.7, alpha_ln=1, fontsize=12,
                           legend=True):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    cm = plt.get_cmap(cmap)
    n_branches = len(branches)

    for i, branch in enumerate(branches):
        velocity = branch['velocity']
        offset = branch['offset']
        r2 = branch['r2']
        distances = branch['distances']
        peak_times = branch['peak_times']

        ax.plot(peak_times, distances, marker='o', ls='', color=cm(i / n_branches), alpha=alpha_marker)
        xs = np.array([peak_times[0], peak_times[-1]])
        ys = velocity * xs + offset
        ax.plot(xs, ys, ls='--', color=cm(i / n_branches), alpha=alpha_ln,
                label=f"vel.: {np.round(velocity)} mm/s\n"
                      f"$r^2$: {np.round(r2, 2)}")
        ax.set_xlabel("Peak time (ms)", fontsize=fontsize)
        ax.set_ylabel("Distance ($\mu$m)", fontsize=fontsize)
        if legend:
            ax.legend()
    return ax


def plot_velocity(peak_times, distances, velocity, offset, color=None, r2=None, ax=None,
                  extend_line=0.2, alpha_markers=0.3, alpha_line=0.8, lw=1, fs=15, plot_markers=True, **kwargs):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    if color is None:
        color = 'C0'

    if plot_markers:
        ax.plot(peak_times, distances, marker='o', ls='', color=color, alpha=alpha_markers, **kwargs)
    pt_ptp = np.ptp(peak_times)
    xs = np.array([np.min(peak_times) - extend_line * pt_ptp, np.max(peak_times) + extend_line * pt_ptp])
    ys = velocity * xs + offset
    if r2 is not None:
        label = f"velocity: {np.round(velocity, 1)} mm/s\nr2: {np.round(r2, 2)}"
    else:
        label = f"velocity: {np.round(velocity, 1)} mm/s"
    ax.plot(xs, ys, ls='--', color=color, alpha=alpha_line, label=label, lw=lw)
    ax.set_xlabel("Peak time (ms)", fontsize=fs)
    ax.set_ylabel("Distance ($\mu$m)", fontsize=fs)
    ax.legend(fontsize=fs)
    return ax


def plot_template_propagation(template, locations, selected_channels, sort_templates=False,
                              color=None, color_marker=None, ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    if color is None:
        color = 'C0'
    if color_marker is None:
        color_marker = 'r'

    template_selected = template[selected_channels]
    locs = locations[selected_channels]
    if sort_templates:
        peaks = np.argmin(template_selected, 1)
        sorted_idxs = np.argsort(peaks)
        template_sorted = template_selected[sorted_idxs]
        locs_sorted = locs[sorted_idxs]
    else:
        template_sorted = template_selected
        locs_sorted = locs
    dist_peaks_sorted = np.array([np.linalg.norm(loc - locs_sorted[0]) for loc in locs_sorted])
    dist_peaks_sorted /= np.max(dist_peaks_sorted)
    dist_peaks_sorted *= len(template_sorted)

    ptp_glob = np.max(np.ptp(template_sorted, 1))
    for i, temp in enumerate(template_sorted):
        temp_shifted = temp + i * 1.5 * ptp_glob  # dist_peaks_sorted[i]
        min_t = np.min(temp_shifted)
        min_idx = np.argmin(temp_shifted)

        ax.plot(temp_shifted, color=color)
        ax.plot(min_idx, min_t, marker='o', color=color_marker)

    ax.axis('off')
    return ax


def plot_peak_latency_map(template, locations, fs, cmap='viridis', log=False,
                          elec_size=8, alpha=0.9, ax=None, colorbar=False,
                          colorbar_orientation="vertical", colorbar_shrink=0.5,
                          plot_image=True, **to_image_kwargs):
    """
    Plots peak latency map for extracellular template
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()

    temp_map = np.abs(np.argmin(template, axis=1)).astype('float')
    temp_map = temp_map / fs * 1000
    probe = get_probe(locations, elec_size)
    ax, im = _plot_image(temp_map, probe, log=log, alpha=alpha, cmap=cmap, ax=ax, plot_image=plot_image,
                         **to_image_kwargs)

    if colorbar:
        if colorbar:
            cbar = fig.colorbar(im,  ax=ax, use_gridspec=True, orientation=colorbar_orientation,
                                shrink=colorbar_shrink)
            max_value = np.max(temp_map)
            ticks = np.arange(int(max_value + 1)).astype("int")
            cbar.set_ticks(ticks)
        if log:
            min_value = np.min(temp_map)
            max_value = np.max(temp_map)
            ticks = np.arange(int(max_value + 1)).astype("int")
            log_ticks = np.log(ticks - min_value + 1)
            cbar.set_ticks(log_ticks)
        cbar.set_ticklabels(ticks)
        cbar.set_label('Latency (ms)')

    return ax


def plot_amplitude_map(template, locations, cmap='viridis', log=False,
                       elec_size=8, alpha=0.9, ax=None, colorbar=False,
                       colorbar_orientation="vertical", colorbar_shrink=0.5,
                       plot_image=True, **to_image_kwargs):
    """
    Plots amplitude map of extracellular template
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()

    temp_map = np.max(np.abs(template), axis=1)
    probe = get_probe(locations, elec_size)
    ax, im = _plot_image(temp_map, probe, log=log, alpha=alpha, cmap=cmap, ax=ax, plot_image=plot_image,
                         **to_image_kwargs)

    if colorbar:
        if colorbar:
            cbar = fig.colorbar(im,  ax=ax, use_gridspec=True, orientation=colorbar_orientation,
                                shrink=colorbar_shrink)
        if log:
            min_value = np.min(temp_map)
            max_value = np.max(temp_map)

            max_exp = int(np.floor(np.log10(max_value)))

            ticks = [10**i for i in range(max_exp + 1)]

            if max_value - ticks[-1] > 30:
                ticks = ticks + [np.round(max_value - 1)]
            ticks = np.array(ticks).astype(int)

            log_ticks = np.log(ticks - min_value + 1)
            cbar.set_ticks(log_ticks)
            cbar.set_ticklabels(ticks)
        cbar.set_label('Amplitude ($\mu$V)')

    return ax


def plot_peak_std_map(template, locations, fs, cmap='viridis', neighbor_distance=30, log=False,
                      elec_size=8, alpha=0.9, ax=None, plot_image=True, **to_image_kwargs):
    """
    Plots peak time standard deviation map
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    peak_time_stds = compute_peak_time_stds(template, locations, fs, neighbor_distance)
    temp_map = peak_time_stds

    probe = get_probe(locations, elec_size)
    ax, im = _plot_image(temp_map, probe, log=log, alpha=alpha, cmap=cmap, ax=ax, plot_image=plot_image,
                         **to_image_kwargs)

    return ax


def play_template_map(template, locations, gtr=None, elec_size=8, cmap='viridis', log=False, ax=None, skip_frames=1,
                      interval=10, **to_image_kwargs):
    """
    Plays animation of extracellular template
    """
    from matplotlib import animation

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()

    if gtr is not None:
        for br in gtr.branches:
            ax.plot(gtr.locations[br["channels"], 0], gtr.locations[br["channels"], 1], color="k", alpha=0.2)

    probe = get_probe(locations, elec_size)

    if log:
        min_value = np.min(template).copy()
        template = np.log(template - np.min(template) + 1)
    else:
        min_value = None
    vmin, vmax = [np.min(template), np.max(template)]
    ims = []
    template_skip = template[:, ::skip_frames]
    for i in tqdm(np.arange(template_skip.shape[1]), desc="Generating frames"):
        t = template_skip[:, i]
        img, xlims, ylims = _get_image(t, probe, log=False, compress=False, min_value=min_value, **to_image_kwargs)
        im = ax.imshow(img, extent=xlims + ylims, origin="lower", vmin=vmin, vmax=vmax, cmap=cmap)

        ims.append([im])

    ax.axis("off")

    ani = animation.ArtistAnimation(fig, ims, interval=interval, blit=True,
                                    repeat_delay=1000)
    return ani


# def play_contour(template, locations, elec_size=8, cmap='viridis', log=False, ax=None, skip_frames=1,
#                  interval=2, **to_image_kwargs):
#     from matplotlib import animation
#
#     if ax is None:
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
#     else:
#         fig = ax.get_figure()
#
#     probe = _get_probe(locations, elec_size)
#
#     if log:
#         min_value = np.min(template).copy()
#         template = np.log(template - np.min(template) + 1)
#     else:
#         min_value = None
#     vmin, vmax = [np.min(template), np.max(template)]
#     ims = []
#     template_skip = template[:, ::skip_frames]
#     for i in tqdm(np.arange(template_skip.shape[1]), desc="Generating frames"):
#         t = template_skip[:, i]
#         img, xlims, ylims = _get_image(t, probe, log=False, compress=False, min_value=min_value, **to_image_kwargs)
#         im = ax.contourf(img, extent=xlims + ylims, origin="lower", vmin=vmin, vmax=vmax, cmap=cmap)
#         ims.append([im])
#
#     ani = animation.ArtistAnimation(fig, ims, interval=interval, blit=True,
#                                     repeat_delay=1000)
#     return ani


def plot_template(template, locations, channels=None, ax=None, pitch=None, **kwargs):
    """
    Plots extracellular template
    """
    probe = mu.return_mea(info={'pos': locations, 'center': False, 'pitch': pitch})
    ax = mu.plot_mea_recording(template, probe, channels=channels, ax=ax, **kwargs)

    return ax


def plot_axon_summary(gtr, ax=None, fig=None, figsize=(10, 7),
                      cmap='rainbow', alpha=0.8, marker='.', markersize=5):
    """
    Plots as summary image with:
        - template amplitude map with axonal branches
        - propagation of each branch
        - velocity of each branch

    Parameters
    ----------
    gtr: GraphAxonTracking

    Returns
    -------
    ax: matplotlib axis

    """
    if ax is None:
        if fig is None:
            fig = plt.figure(figsize=figsize)
            ax1 = fig.add_subplot(2, 2, 1)
            ax2 = fig.add_subplot(2, 2, 2)
            ax3 = fig.add_subplot(2, 2, 3)
            ax4 = fig.add_subplot(2, 2, 4)
        else:
            ax1 = fig.add_subplot(2, 2, 1)
            ax2 = fig.add_subplot(2, 2, 2)
            ax3 = fig.add_subplot(2, 2, 3)
            ax4 = fig.add_subplot(2, 2, 4)
    else:
        gs_global = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=ax)
        fig = ax.get_figure()
        ax1 = fig.add_subplot(gs_global[0, 0])
        ax2 = fig.add_subplot(gs_global[0, 1])
        ax3 = fig.add_subplot(gs_global[1, 0])
        ax4 = fig.add_subplot(gs_global[1, 1])
        ax.axis('off')

    axes = [ax1, ax2, ax3, ax4]

    cm = plt.get_cmap(cmap)

    ax3.axis("off")
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    template = gtr.template
    locations = gtr.locations
    init = locations[gtr.init_channel]

    plot_amplitude_map(template, locations, log=True, ax=ax1)
    ax1.plot(init[0], init[1], marker='o', color='r', markersize=5)
    plot_peak_latency_map(template, locations, gtr.fs, ax=ax2)
    ax2.plot(init[0], init[1], marker='o', color='r', markersize=5)
    gs = gridspec.GridSpecFromSubplotSpec(1, len(gtr.branches), subplot_spec=ax3)
    ax1.set_title("amplitude", fontsize=15)
    ax2.set_title("peak latency", fontsize=15)

    for i, branch in enumerate(gtr.branches):
        color = cm(i / len(gtr.branches))
        sel_idxs = branch['channels']

        axpr = fig.add_subplot(gs[0, i])
        _ = plot_template_propagation(template, locations, sel_idxs,
                                      color='k', sort_templates=False,
                                      color_marker=color, ax=axpr)
        axpr.text(0.45, -0.05, f'br.{i}', horizontalalignment='center',
                  verticalalignment='center', transform=axpr.transAxes, fontsize=10)
        for i, sel in enumerate(sel_idxs[:-1]):
            ax1.plot([locations[sel, 0], locations[sel_idxs[i + 1], 0]],
                     [locations[sel, 1], locations[sel_idxs[i + 1], 1]],
                     color=color, lw=1, marker=marker, markersize=markersize, alpha=alpha)
            ax2.plot([locations[sel, 0], locations[sel_idxs[i + 1], 0]],
                     [locations[sel, 1], locations[sel_idxs[i + 1], 1]],
                     color=color, lw=1, marker=marker, markersize=markersize, alpha=alpha)
    plot_branch_velocities(gtr.branches, ax=ax4, cmap=cmap, fontsize=10)
    ax3.set_title("propagation", fontsize=15)
    ax4.set_title("velocity", fontsize=15)

    fig.subplots_adjust(wspace=0.5, hspace=0.2)

    return fig, axes


def get_probe(locations, width=10):
    """
    Returns a probeinterface probe with swuare electrodes.

    Parameters
    ----------
    locations: np.array
        Locations of electrodes
    width: float
        Width of electrodes (default=10)

    Returns
    -------
    print: probeinterface.Probe
        The probe object
    """
    shapes = "square"
    shape_params = {'width': width}

    probe = pi.Probe(ndim=2, si_units='um')
    probe.set_contacts(positions=locations,
                       shapes=shapes, shape_params=shape_params)
    probe.create_auto_shape(probe_type="rect")
    return probe


def _get_image(values, probe, **to_image_kwargs):
    img, xlims, ylims = probe.to_image(values, **to_image_kwargs)

    return img, xlims, ylims


def _plot_image(values, probe, log=False, compress=False, cmap="viridis", alpha=1., ax=None,  plot_image=False,
                **to_image_kwargs):
    if ax is None:
        fig, ax = plt.subplots()

    if log:
        min_value = np.min(values)
        values = np.log(values - min_value + 1)
    if compress:
        values = 1 / (1 + np.exp(-values))

    if plot_image:
        img, xlims, ylims = _get_image(values, probe, **to_image_kwargs)
        im = ax.imshow(img, extent=xlims + ylims, origin="lower", cmap=cmap, alpha=alpha)
    else:
        probe_shape_kwargs = dict(facecolor='green', edgecolor="k", lw=0, alpha=0.0)
        contacts_kargs = dict(edgecolor=None, lw=0)
        _ = plotting.plot_probe(probe, contacts_values=values, ax=ax, probe_shape_kwargs=probe_shape_kwargs,
                                contacts_kargs=contacts_kargs, cmap=cmap)
        im = None
    ax.axis("off")

    return ax, im
