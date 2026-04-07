# ==========================================================
# pipeline/_reports_mixin.py
# ReportsMixin — Phase 4: reports, unit curation, probe and waveform plots.
# ==========================================================

import traceback
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for headless environments
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf
import spikeinterface.full as si

from .stages import ProcessingStage

try:
    from ..scalebury import add_scalebar
except ImportError:
    try:
        from MEA_Analysis.IPNAnalysis.scalebury import add_scalebar
    except ImportError:
        from scalebury import add_scalebar


class ReportsMixin:
    """Generates quality-metric reports, curates units, and produces static plots."""

    # ------------------------------------------------------------------
    # Phase 4: Reports & Curation
    # ------------------------------------------------------------------

    def generate_reports(
        self,
        thresholds=None,
        no_curation=False,
        export_phy=False,
        plot_mode="separate",
        plot_debug=False,
        raster_sort=None,
        fixed_y=False,
    ):
        if self.state['stage'] == ProcessingStage.REPORTS_COMPLETE.value:
            return

        self.logger.info("--- [Phase 4] Reports & Curation ---")
        try:
            q_metrics = self.analyzer.get_extension("quality_metrics").get_data()
            t_metrics = self.analyzer.get_extension("template_metrics").get_data()
            locations = self.analyzer.get_extension("unit_locations").get_data()

            q_metrics['loc_x'] = locations[:, 0]
            q_metrics['loc_y'] = locations[:, 1]

            q_metrics.to_excel(self.output_dir / "qm_unfiltered.xlsx")
            t_metrics.to_excel(self.output_dir / "tm_unfiltered.xlsx")
            self._plot_probe_locations(q_metrics.index.values, locations, "locations_unfiltered.pdf")

            if no_curation:
                self.logger.info("Skipping curation.")
                clean_units = q_metrics.index.values
            else:
                self.logger.info("Applying curation.")
                clean_metrics, rejection_log = self._apply_curation_logic(q_metrics, thresholds)
                clean_units = clean_metrics.index.values
                clean_metrics.to_excel(self.output_dir / "metrics_curated.xlsx")
                rejection_log.to_excel(self.output_dir / "rejection_log.xlsx")
                t_metrics.loc[clean_units].to_excel(self.output_dir / "tm_curated.xlsx")

            if len(clean_units) == 0:
                self.logger.warning("No units passed curation.")
                self._save_checkpoint(ProcessingStage.REPORTS_COMPLETE, n_units=0)
                return

            mask = np.isin(self.analyzer.unit_ids, clean_units)
            self._plot_probe_locations(
                clean_units, locations[mask], f"locations_{len(clean_units)}_units.pdf"
            )
            self._plot_waveforms_grid(clean_units)
            self._run_burst_analysis(
                clean_units,
                plot_mode=plot_mode,
                plot_debug=plot_debug,
                raster_sort=raster_sort,
                fixed_y=fixed_y,
            )

            if export_phy:
                si.export_to_phy(
                    self.analyzer.select_units(clean_units),
                    output_folder=self.output_dir / "phy_output",
                    remove_if_exists=True,
                    copy_binary=False,
                )

            self._save_checkpoint(
                ProcessingStage.REPORTS_COMPLETE,
                n_units=len(clean_units),
                failed_stage=None,
                error=None,
            )
        except Exception as e:
            err = {
                "failed_stage": ProcessingStage.REPORTS.name,
                "exception": type(e).__name__,
                "message": str(e),
                "traceback": traceback.format_exc(),
                "time": str(datetime.now()),
            }
            self.logger.error(err["traceback"])
            self._save_checkpoint(ProcessingStage.ANALYZER_COMPLETE, error=err)
            raise

    # ------------------------------------------------------------------
    # Curation
    # ------------------------------------------------------------------

    def _apply_curation_logic(self, metrics, user_thresholds):
        defaults = {
            'presence_ratio': 0.75,
            'rp_contamination': 0.15,
            'firing_rate': 0.05,
            'amplitude_median': -20,
            'amplitude_cv_median': 0.5,
        }
        if user_thresholds:
            defaults.update(user_thresholds)

        keep_mask = np.ones(len(metrics), dtype=bool)
        rejections = []
        for idx, row in metrics.iterrows():
            reasons = []
            if row.get('presence_ratio', 1) < defaults['presence_ratio']:
                reasons.append("Low Presence")
            if row.get('rp_contamination', 0) > defaults['rp_contamination']:
                reasons.append("High Contam")
            if row.get('firing_rate', 0) < defaults['firing_rate']:
                reasons.append("Low FR")
            if row.get('amplitude_median', -100) > defaults['amplitude_median']:
                reasons.append("Low Amp")
            # TODO: add cv_median logic after confirming metric exists in SI

            if reasons:
                keep_mask[metrics.index.get_loc(row.name)] = False
                rejections.append({"unit_id": row.name, "reasons": "; ".join(reasons)})

        return metrics[keep_mask], pd.DataFrame(rejections)

    # ------------------------------------------------------------------
    # Probe and waveform plots
    # ------------------------------------------------------------------

    def _plot_probe_locations(self, unit_ids, locations, filename):
        fig, ax = plt.subplots(figsize=(10.5, 6.5))
        si.plot_probe_map(self.recording, ax=ax, with_channel_ids=False)
        ax.scatter(locations[:, 0], locations[:, 1], s=10, c='blue', alpha=0.6)
        ax.invert_yaxis()
        fig.savefig(self.output_dir / filename)
        plt.close(fig)

    def _plot_waveforms_grid(self, unit_ids):
        pdf_path = self.output_dir / "waveforms_grid.pdf"
        self.logger.info(f"Generating PDF: {pdf_path}")

        wf_ext = self.analyzer.get_extension("waveforms")
        fs = self.recording.get_sampling_frequency()

        with pdf.PdfPages(pdf_path) as pdf_doc:
            units_per_page = 12
            for i in range(0, len(unit_ids), units_per_page):
                batch = unit_ids[i: i + units_per_page]
                fig, axes = plt.subplots(3, 4, figsize=(12, 9))
                axes = axes.flatten()

                for ax, uid in zip(axes, batch):
                    wf = wf_ext.get_waveforms_one_unit(uid)
                    mean_wf = np.mean(wf, axis=0)
                    best_ch = np.argmin(np.min(mean_wf, axis=0))

                    time_ms = np.arange(wf.shape[1]) / fs * 1000

                    n_spikes = wf.shape[0]
                    if n_spikes > 100:
                        indices = np.random.choice(n_spikes, 100, replace=False)
                        spikes_to_plot = wf[indices, :, best_ch]
                    else:
                        spikes_to_plot = wf[:, :, best_ch]

                    ax.plot(time_ms, spikes_to_plot.T, c='gray', lw=0.5, alpha=0.3)
                    ax.plot(time_ms, mean_wf[:, best_ch], c='red', lw=1.5)
                    ax.set_title(f"Unit {uid} | Ch {best_ch}", fontsize=10)

                    try:
                        add_scalebar(
                            ax,
                            matchx=False,
                            matchy=False,
                            sizex=1.0,
                            labelx='1 ms',
                            sizey=50,
                            labely='50 µV',
                            loc='lower right',
                            hidex=True,
                            hidey=True,
                        )
                    except Exception:
                        ax.spines['top'].set_visible(False)
                        ax.spines['right'].set_visible(False)

                for j in range(len(batch), len(axes)):
                    axes[j].axis('off')

                pdf_doc.savefig(fig)
                plt.close(fig)
