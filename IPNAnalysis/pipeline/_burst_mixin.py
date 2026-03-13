# ==========================================================
# pipeline/_burst_mixin.py
# BurstMixin — network burst analysis and raster-sort helpers.
# ==========================================================

import json
import os
import traceback

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for headless environments
import matplotlib.pyplot as plt

try:
    from ..parameter_free_burst_detector import compute_network_bursts
    from .. import helper_functions as helper
except ImportError:
    try:
        from MEA_Analysis.IPNAnalysis.parameter_free_burst_detector import compute_network_bursts
        from MEA_Analysis.IPNAnalysis import helper_functions as helper
    except ImportError:
        from parameter_free_burst_detector import compute_network_bursts
        import helper_functions as helper


class BurstMixin:
    """Handles network burst analysis and raster-plot generation."""

    # ------------------------------------------------------------------
    # Network Burst Analysis
    # ------------------------------------------------------------------

    def _run_burst_analysis(
        self,
        ids_list=None,
        plot_mode='separate',
        plot_debug=False,
        raster_sort='none',
        fixed_y=True,
    ):
        self.logger.info("Running Network Burst Analysis...")

        spike_times = {}

        # Load spike times — from sorter or from saved .npy file.
        if self.sorting:
            fs = self.recording.get_sampling_frequency()
            if ids_list is None:
                ids_list = self.analyzer.unit_ids

            missing_unit_ids = []
            for uid in ids_list:
                try:
                    spike_times[uid] = self.sorting.get_unit_spike_train(uid) / fs
                except KeyError:
                    missing_unit_ids.append(uid)

            if missing_unit_ids:
                self.logger.warning(
                    "Skipping %d unit(s) not present in active sorting during burst analysis: %s",
                    len(missing_unit_ids),
                    missing_unit_ids[:20],
                )

            if not spike_times:
                self.logger.error(
                    "No valid units left for burst analysis after filtering missing unit IDs."
                )
                return

            np.save(self.output_dir / "spike_times.npy", spike_times)
        else:
            npy_spike_file = self.output_dir / "spike_times.npy"
            if npy_spike_file.exists():
                spike_times = np.load(npy_spike_file, allow_pickle=True).item()
            else:
                self.logger.error("No spike times found for burst analysis.")
                return

        if not spike_times:
            return

        try:
            # Network burst computation
            network_data = compute_network_bursts(SpikeTimes=spike_times, plot=False)

            network_data_clean = helper.recursive_clean(network_data)
            network_data_clean['n_units'] = len(spike_times)

            # Atomic write (temp → rename) to prevent corruption.
            temp_file = self.output_dir / "network_results.tmp.json"
            final_file = self.output_dir / "network_results.json"
            with open(temp_file, 'w') as f:
                json.dump(network_data_clean, f, indent=2)
            if temp_file.exists():
                os.replace(temp_file, final_file)
                self.logger.info(f"Successfully saved: {final_file}")

            sorted_units = self._sort_units_for_raster(spike_times, raster_sort)

            # Build figure according to requested plot mode.
            if plot_mode == 'separate':
                fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
                ax_raster, ax_network = axs
                helper.plot_clean_raster(
                    ax_raster,
                    spike_times,
                    sorted_units,
                    color='gray',
                    markersize=4,
                    markeredgewidth=0.5,
                    alpha=1.0,
                )
                helper.plot_clean_network(ax_network, **network_data["plot_data"])

            elif plot_mode == 'merged':
                fig, ax = plt.subplots(figsize=(12, 5))
                ax_raster = ax
                ax_network = ax.twinx()
                helper.plot_clean_raster(
                    ax_raster,
                    spike_times,
                    sorted_units,
                    color='gray',
                    markersize=4,
                    markeredgewidth=0.5,
                    alpha=1.0,
                )
                helper.plot_clean_network(ax_network, **network_data["plot_data"])
                ax.spines["right"].set_visible(True)
            else:
                self.logger.warning(f"Unknown plot mode: {plot_mode}")
                return

            plt.tight_layout()
            plt.subplots_adjust(hspace=0.05)

            if plot_debug:
                nb_events = network_data["network_bursts"]["events"]
                burst_intervals = [(ev["start"], ev["end"]) for ev in nb_events]
                sb_events = network_data["superbursts"]["events"]
                sb_intervals = [(ev["start"], ev["end"]) for ev in sb_events]
                for start, end in burst_intervals:
                    ax_network.axvspan(start, end, color='gray', alpha=0.1)
                for start, end in sb_intervals:
                    ax_network.axvspan(start, end, color='gray', alpha=0.2)

            plt.savefig(self.output_dir / "raster_burst_plot.svg")

            # 60 s zoom
            ax_raster.set_xlim(0, 60)
            ax_network.set_xlim(0, 60)
            plt.savefig(self.output_dir / "raster_burst_plot_60s.svg")

            # 30 s zoom
            ax_raster.set_xlim(0, 30)
            ax_network.set_xlim(0, 30)
            ax_network.set_xlabel("Time (s)")
            plt.savefig(self.output_dir / "raster_burst_plot_30s.svg")

            plt.savefig(self.output_dir / "raster_burst_plot.png", dpi=300)
            plt.close(fig)

            # Optional fixed-y re-plot using a pre-computed project-wide y-max.
            if fixed_y:
                self._plot_fixed_y(
                    spike_times,
                    network_data,
                    sb_intervals if plot_debug else [],
                )

        except Exception as e:
            self.logger.error(f"Burst analysis error: {e}")
            traceback.print_exc()
            raise

    def _plot_fixed_y(self, spike_times, network_data, sb_intervals):
        """Re-plot raster with a fixed y-axis derived from the project summary."""
        summary_file = (
            self.output_root / self.project_name / f"{self.project_name}_y_max_summary.json"
        )
        if not summary_file.exists():
            self.logger.error(
                f"No y-max summary found at {summary_file}. Run without --fixed-y first."
            )
            return

        with open(summary_file, 'r') as f:
            summary = json.load(f)

        all_maxima = [
            v
            for date in summary.values()
            for chip in date.values()
            for v in chip.values()
        ]
        global_max = max(all_maxima)
        self.logger.info(f"Applying fixed y-max: {global_max:.4f}")

        fig2, axs2 = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        ax_raster2, ax_network2 = axs2
        helper.plot_clean_raster(
            ax_raster2, spike_times, color='gray', markersize=4, markeredgewidth=0.5, alpha=1.0
        )
        helper.plot_clean_network(ax_network2, **network_data["plot_data"])
        ax_network2.set_ylim(0, global_max)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.05)
        for start, end in sb_intervals:
            ax_network2.axvspan(start, end, color='gray', alpha=0.3)

        plt.savefig(self.output_dir / "fixed_y_raster_burst_plot.svg")
        plt.savefig(self.output_dir / "fixed_y_raster_burst_plot.png", dpi=300)
        ax_raster2.set_xlim(0, 60)
        ax_network2.set_xlim(0, 60)
        plt.savefig(self.output_dir / "fixed_y_raster_burst_plot_60s.svg")
        ax_raster2.set_xlim(0, 30)
        ax_network2.set_xlim(0, 30)
        ax_network2.set_xlabel("Time (s)")
        plt.savefig(self.output_dir / "fixed_y_raster_burst_plot_30s.svg")
        plt.savefig(self.output_dir / "fixed_y_raster_burst_plot_30s.png", dpi=300)
        plt.close(fig2)

    # ------------------------------------------------------------------
    # Raster sort
    # ------------------------------------------------------------------

    def _sort_units_for_raster(self, spike_times, raster_sort):
        """Return an ordered list of unit keys for the raster y-axis."""
        if raster_sort == 'none' or raster_sort is None:
            return None  # plot_clean_raster handles default ordering

        if raster_sort == 'firing_rate':
            return sorted(spike_times.keys(), key=lambda uid: len(spike_times[uid]))

        if raster_sort == 'unit_id':
            return sorted(spike_times.keys())

        self.logger.warning(f"Unknown raster_sort: {raster_sort}. Falling back to none.")
        return None
