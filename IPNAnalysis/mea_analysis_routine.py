import os
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
os.environ["QT_QPA_PLATFORM"] = "offscreen"  # For headless environments
import sys
import shutil
import gc
import json
import re
import traceback
import logging
import argparse
import configparser
from pathlib import Path
from datetime import datetime
from timeit import default_timer as timer
from enum import Enum
import psutil
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf

# Scientific Libraries
import spikeinterface.full as si
import spikeinterface.preprocessing as spre
# Import your custom modules
try:
    from parameter_free_burst_detector import compute_network_bursts
    import helper_functions as helper
except ImportError:
    print("Warning: Custom helper modules not found. Analysis steps may fail.")

# --- Enums for Checkpointing ---
class ProcessingStage(Enum):
    NOT_STARTED = 0
    PREPROCESSING_COMPLETE = 1
    SORTING_COMPLETE = 2
    ANALYZER_COMPLETE = 3
    REPORTS_COMPLETE = 4
    FAILED = -1

# --- The Main Pipeline Class ---
class MEAPipeline:
    """
    
    Encapsulates Preprocessing, Sorting, Analysis, and Curation.
    """

    def __init__(self, file_path, stream_id='well000', recording_num='rec0000', output_root=None, checkpoint_root=None, 
                 sorter='kilosort4', docker_image=None, verbose=True, 
                 cleanup=False, force_restart=False):
        
        self.file_path = Path(file_path).resolve()
        self.stream_id = stream_id
        self.recording_num = recording_num
        self.sorter = sorter
        self.docker_image = docker_image
        self.verbose = verbose
        self.cleanup_flag = cleanup
        self.force_restart = force_restart

        # 1. Parse Metadata & Paths
        self.metadata = self._parse_metadata()
        self.run_id = self.metadata.get('run_id', 'UnknownRun')
        self.project_name = self.metadata.get('project', 'UnknownProject')
        
        # Define Directory Structure: Output / Pattern / Well
        # E.g., ./AnalyzedData/250908/M06844/Network/000029/well000
        self.relative_pattern = self.metadata.get('relative_pattern', 'UnknownPattern')
        self.output_dir = Path(output_root) / self.relative_pattern / self.stream_id
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 2. Logger Setup
        log_file = self.output_dir / f"{self.run_id}_{self.stream_id}_pipeline.log"
        self.logger = self._setup_logger(log_file)
        
        # 3. Checkpointing Setup
        ckpt_root = Path(checkpoint_root) if checkpoint_root else self.output_dir / "checkpoints"
        ckpt_root.mkdir(parents=True, exist_ok=True)
        self.checkpoint_file = ckpt_root / f"{self.project_name}_{self.run_id}_{self.stream_id}_checkpoint.json"
        self.state = self._load_checkpoint()

        # Data placeholders
        self.recording = None
        self.sorting = None
        self.analyzer = None

    # --- Setup Methods ---
    def _setup_logger(self, log_file):
        logger = logging.getLogger(f"mea_{self.stream_id}")
        logger.setLevel(logging.DEBUG if self.verbose else logging.INFO)
        # Prevent duplicate handlers if running in loop
        if not logger.handlers:
            formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s')
            fh = logging.FileHandler(log_file, mode='w') # Overwrite log per run
            fh.setFormatter(formatter)
            logger.addHandler(fh)
            ch = logging.StreamHandler(sys.stdout)
            ch.setFormatter(formatter)
            logger.addHandler(ch)
        return logger

    def _parse_metadata(self):
            """Extracts RunID, ChipID, Date, and specific Well data based on stream_id"""
            meta = {
                'run_id': None, 
                'chip_id': None, 
                'project': None, 
                'relative_pattern': f"{self.file_path.parent.parent.name}/{self.file_path.parent.name}/{self.file_path.name}",
                'date': None,
                # These will be populated if the specific well is found
                'well_label': None,
                'plating_date': None,
                'well_id': None
            }
            
            # Strategy A: Regex on path (Fallback)
            try:
                path_str = str(self.file_path)
                match = re.search(r"/(\d+)/data.raw.h5", path_str)
                if match:
                    meta['run_id'] = match.group(1)
                
                parts = path_str.split(os.sep)
                if len(parts) > 5:
                    meta['relative_pattern'] = os.path.join(*parts[-6:-2])
                    meta['project'] = parts[-5]
                    meta['date'] = parts[-4]
                    meta['chip_id'] = parts[-3]
            except Exception:
                pass

            # Strategy B: .metadata file (Overrides regex & gets specific well info)
            meta_file = self.file_path.parent / ".metadata"
            if meta_file.exists():
                try:
                    cfg = configparser.ConfigParser()
                    cfg.read(meta_file, encoding='utf-8') 
                    
                    # 1. Standard Properties
                    if 'properties' in cfg:
                        meta['run_id'] = cfg['properties'].get('runid', meta.get('run_id'))
                        meta['project'] = cfg['properties'].get('project_title', meta.get('project'))
                    
                    if 'runtime' in cfg:
                        meta['chip_id'] = cfg['runtime'].get('chipid', meta.get('chip_id'))

                    # 2. Extract ONLY the well matching self.stream_id
                    if 'wells' in cfg:
                        wells = cfg['wells']
                        num_wells = int(wells.get(r'info\size', 0))
                        
                        # We need to match this specific stream ID
                        target_id = int(self.stream_id)
                        
                        for i in range(1, num_wells + 1):
                            prefix = f"info\\{i}"
                            
                            # Check if this well's ID matches our stream
                            current_well_id = int(wells.get(f"{prefix}\\id", -99))
                            
                            if current_well_id == target_id:
                                # Match found! Extract only this well's data
                                meta['well_id'] = current_well_id
                                meta['well_label'] = wells.get(f"{prefix}\\groupname", '')
                                meta['plating_date'] = wells.get(f"{prefix}\\annotations\\property\\1\\propertyValue", '')
                                
                                # We found our match, no need to parse the rest
                                break

                except Exception as e:
                    # print(f"Metadata parsing warning: {e}")
                    pass
                    
            return meta

    # --- Checkpointing Methods ---
    def _load_checkpoint(self):
        if self.checkpoint_file.exists() and not self.force_restart:
            with open(self.checkpoint_file, 'r') as f:
                return json.load(f)
        return {'stage': ProcessingStage.NOT_STARTED.value}

    def _save_checkpoint(self, stage, **kwargs):
        self.state['stage'] = stage.value
        self.state['last_updated'] = str(datetime.now())
        self.state.update(kwargs)
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.state, f, indent=2)
        self.logger.info(f"Checkpoint Saved: {stage.name}")

    def should_skip(self):
        # Skip if already done and not forcing restart
        if self.state['stage'] == ProcessingStage.REPORTS_COMPLETE.value and not self.force_restart:
            self.logger.info("Pipeline already completed. Skipping.")
            return True
        return False

    # --- Phase 1: Preprocessing ---
    def run_preprocessing(self):
        if self.state['stage'] >= ProcessingStage.PREPROCESSING_COMPLETE.value:
            self.logger.info("Resuming: Loading preprocessed data from cache.")
            zarr_path = self.output_dir / "preprocessed.zarr"
            if zarr_path.exists():
                self.recording = si.load_extractor(zarr_path)
                return

        self.logger.info("--- [Phase 1] Preprocessing ---")
        
        # 1. Universal Loading
        raw_rec = self._load_recording_file()

        # 2.  Preprocessing Chain
        # A. Phase Shift (Maxwell Hardware correction)
        if raw_rec.has_property('inter_sample_shift'):
            self.logger.info("Applying Phase Shift...")
            rec = spre.phase_shift(raw_rec)
        else:
            rec = raw_rec

        # B. Unsigned to Signed conversion (if needed)
        if rec.get_dtype().kind == 'u':
            rec = spre.unsigned_to_signed(rec)

        # C. Bad Channel Interpolation
        bad_channels, _ = spre.detect_bad_channels(rec, method="coherence+psd", std_dmad_threshold=5)
        if len(bad_channels) > 0:
            self.logger.info(f"Interpolating {len(bad_channels)} bad channels")
            rec = spre.interpolate_bad_channels(rec, bad_channels)

        # D. Highpass Filter (Remove DC drift, keep spikes)
        rec = spre.highpass_filter(rec, freq_min=300)

        # E. Local Common Median Reference (Better for bursting)
        try:
            rec = spre.common_reference(rec, reference='local', operator='median', local_radius=(250, 250))
        except:
            self.logger.warning("Local CMR failed (missing locations?), using Global CMR.")
            rec = spre.common_reference(rec, reference='global', operator='median')

        # F. Safe Scaling to Int16
        if rec.get_dtype() != 'int16':
            rec = spre.scale(rec, dtype='int16')

        # 3. Save to Zarr (Cache)
        zarr_path = self.output_dir / "preprocessed.zarr"
        if zarr_path.exists(): shutil.rmtree(zarr_path) # Clean start
        
        self.logger.info("Caching preprocessed data to Zarr...")
        self.recording = rec.save(folder=zarr_path, format="zarr", compressor="blosc", 
                                  n_jobs=16, chunk_duration="1s", progress_bar=self.verbose)
        
        self._save_checkpoint(ProcessingStage.PREPROCESSING_COMPLETE)

    def _load_recording_file(self):
        fpath = str(self.file_path)
        if fpath.endswith(".h5"):
            return si.read_maxwell(fpath, stream_id=self.stream_id,rec_name=self.recording_num, load_if_exists=True)
        elif fpath.endswith(".nwb"):
            return si.read_nwb(fpath, load_if_exists=True)
        elif self.file_path.is_dir():
            return si.load_extractor(self.file_path)
        raise ValueError(f"Unknown format: {fpath}")

    # --- Phase 2: Sorting ---
    def run_sorting(self):
        if self.state['stage'] >= ProcessingStage.SORTING_COMPLETE.value:
            self.logger.info("Resuming: Loading existing sorting.")
            sorter_folder = self.output_dir / "sorter_output"
            # Try generic loader, fallback to kilosort specific
            try:
                self.sorting = si.read_sorter_folder(sorter_folder)
            except:
                self.sorting = si.read_kilosort(sorter_folder)
            return

        self.logger.info(f"--- [Phase 2] Spike Sorting ({self.sorter}) ---")
        sorter_folder = self.output_dir / "sorter_output"

        #Parameters for 2D Static Cultures (Kilosort4)
        ks_params = {
            'nblocks': 0,           # DISABLE drift correction (Crucial for 2D)
            'Th_universal': 9,      # Higher threshold = cleaner units
            'min_template_size': 10,
            'batch_size': 60000,
            'do_correction': False
        }

        # Unified execution (Docker vs Local handled by argument)
        start = timer()
        self.sorting = si.run_sorter(
            sorter_name=self.sorter,
            recording=self.recording,
            output_folder=sorter_folder,
            remove_existing_folder=True,
            verbose=self.verbose,
            docker_image=self.docker_image, # None = Local, String = Docker
            **ks_params
        )
        self.logger.info(f"Sorting finished in {timer()-start:.2f}s")

        # Post-processing: Remove Redundant Units (Fix over-splitting)
        self.logger.info("Running Redundancy Removal...")
        self.sorting = si.remove_redundant_units(self.sorting, duplicate_threshold=0.9)
        
        self._save_checkpoint(ProcessingStage.SORTING_COMPLETE)

    def _spike_detection_only(self):
            """
            Detects spikes (threshold crossings) on each channel without clustering.
            Used when --skip-spikesorting is active.
            """
            self.logger.info("--- [Phase 2-Alt] Spike Detection (No Sorting) ---")
            
            # 1. Detect Peaks using SpikeInterface (Memory Efficient)
            # method='by_channel' detects peaks independently on each trace
            self.logger.info("Running threshold detection (std=5, sign=neg)...")
            
            # Use simple dictionary for job kwargs inside this method
            job_kwargs = {'n_jobs': 16, 'chunk_duration': '1s', 'progress_bar': self.verbose}
            
            peaks = si.detect_peaks(
                self.recording,
                method='by_channel',
                peak_sign='neg',
                detect_threshold=5,
                exclude_sweep_ms=0.1,
                **job_kwargs
            )
            
            self.logger.info(f"Detected {len(peaks)} total spikes across all channels.")

            # 2. Convert to Dictionary format {channel_id: times_in_seconds}
            fs = self.recording.get_sampling_frequency()
            channel_ids = self.recording.get_channel_ids()
            spike_times = {}

            # peaks is a structured numpy array with 'sample_index' and 'channel_index'
            for ch_index, ch_id in enumerate(channel_ids):
                # Extract spikes for this specific channel index
                mask = peaks['channel_index'] == ch_index
                if np.any(mask):
                    # Convert samples to seconds
                    spike_times[ch_id] = peaks['sample_index'][mask] / fs
            
            # 3. Save Results
            # We save this as 'spike_times.npy' so the Burst Analysis can load it 
            # just like it would for sorted units.
            np.save(self.output_dir / "spike_times.npy", spike_times)
            
            # 4. Generate a Channel-based Raster Plot immediately
            self.logger.info("Generating Channel Raster Plot...")
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Helper to plot raster from dict
            # Sort channels by ID for cleaner visualization
            sorted_chs = sorted(spike_times.keys())
            for i, ch in enumerate(sorted_chs):
                times = spike_times.get(ch, [])
                if len(times) > 0:
                    ax.plot(times, np.ones_like(times) * i, '|', markersize=1, color='black', alpha=0.5)
            
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Channel Index")
            ax.set_title(f"Multi-Unit Activity (Threshold Detection) - {self.stream_id}")
            plt.tight_layout()
            plt.savefig(self.output_dir / "mua_channel_raster.png", dpi=300)
            plt.close(fig)

            # Return the list of "ids" (channels) for the next step
            return list(spike_times.keys())

    # --- Phase 3: Analyzer ---
    def run_analyzer(self):
            analyzer_folder = self.output_dir / "analyzer_output"

            if self.state['stage'] >= ProcessingStage.ANALYZER_COMPLETE.value:
                if analyzer_folder.exists():
                    self.logger.info("Resuming: Loading Sorting Analyzer.")
                    self.analyzer = si.load_sorting_analyzer(analyzer_folder)
                    return
            
            self.logger.info("--- [Phase 3] Computing Sorting Analyzer ---")
            if analyzer_folder.exists(): shutil.rmtree(analyzer_folder)

            # 1. Sparsity (Radius 50um optimized for Maxwell, matches old script)
            sparsity = si.estimate_sparsity(
                self.recording, self.sorting, 
                method="radius", radius_um=50, peak_sign='neg'
            )

            # 2. Create Object
            self.analyzer = si.create_sorting_analyzer(
                self.sorting, self.recording, format="binary_folder", 
                folder=analyzer_folder, sparsity=sparsity, return_in_uV=True
            )

            # 3. Compute Extensions (Robust Loop)
            # Added 'template_metrics' which is required for 'amplitude_cv_median' curation
            ext_list = ["random_spikes", "waveforms", "templates", "noise_levels", 
                        "quality_metrics", "template_metrics", "unit_locations"]
            
            ext_params = {
                "waveforms": {"ms_before": 1.0, "ms_after": 2.0},
                "unit_locations": {"method": "monopolar_triangulation"}
            }
            
            for ext in ext_list:
                try:
                    self.logger.info(f"Computing extension: {ext}...")
                    start_ext = timer()
                    self.analyzer.compute(ext, extension_params=ext_params, verbose=self.verbose)
                    self.logger.info(f"  ✓ {ext} done ({timer()-start_ext:.2f}s)")
                except Exception as e:
                    self.logger.error(f"  ✗ FAILED to compute {ext}: {e}")
                    # Crash only on critical extensions
                    if ext in ["waveforms", "templates"]:
                        raise e
            
            self._save_checkpoint(ProcessingStage.ANALYZER_COMPLETE)

    # --- Phase 4: Reports & Curation ---
    def generate_reports(self, thresholds=None, no_curation=False, export_phy=False):
        if self.state['stage'] == ProcessingStage.REPORTS_COMPLETE.value:
            return

        self.logger.info("--- [Phase 4] Reports & Curation ---")
        
        # 1. Load Metrics
        q_metrics = self.analyzer.get_extension("quality_metrics").get_data()
        t_metrics = self.analyzer.get_extension("template_metrics").get_data()
        locations = self.analyzer.get_extension("unit_locations").get_data()
        
        # Attach locations to metrics for Excel
        q_metrics['loc_x'] = locations[:, 0]
        q_metrics['loc_y'] = locations[:, 1]
        
        # Save Unfiltered
        q_metrics.to_excel(self.output_dir / "qm_unfiltered.xlsx")
        t_metrics.to_excel(self.output_dir / "tm_unfiltered.xlsx")

        #plot the unfiltered Probe Map
        self._plot_probe_locations(q_metrics.index.values, locations, "locations_unfiltered.pdf")
        # 2. Automatic Curation
        if no_curation:
            self.logger.info("Skipping curation (User Flag).")
            clean_units = q_metrics.index.values
        else:
            self.logger.info("Applying Curation Rules...")
            clean_metrics, rejection_log = self._apply_curation_logic(q_metrics, thresholds)
            clean_units = clean_metrics.index.values
            
            # Save Logs
            clean_metrics.to_excel(self.output_dir / "metrics_curated.xlsx")
            rejection_log.to_excel(self.output_dir / "rejection_log.xlsx")
            self.logger.info(f"Units passing curation: {len(clean_units)} / {len(q_metrics)}")
            t_metrics_filtered = t_metrics.loc[clean_units]
            t_metrics_filtered.to_excel(self.output_dir / "tm_curated.xlsx")

        if len(clean_units) == 0:
            self.logger.warning("No units passed curation! Skipping plots/bursts.")
            self._save_checkpoint(ProcessingStage.REPORTS_COMPLETE, n_units=0)
            return
        # Filter locations for plot
        mask = np.isin(self.analyzer.unit_ids, clean_units)
        self._plot_probe_locations(clean_units, locations[mask], f"locations_{len(clean_units)}_units.pdf")

        # 3. Plotting (Waveforms Grid PDF)
        self._plot_waveforms_grid(clean_units)

        # 4. Burst Analysis
        self._run_burst_analysis(clean_units)

        # 5. Export to Phy
        if export_phy:
            self.logger.info("Exporting to Phy...")
            phy_dir = self.output_dir / "phy_output"
            # We slice the analyzer to only export "Good" units
            clean_analyzer = self.analyzer.select_units(clean_units)
            si.export_to_phy(clean_analyzer, output_folder=phy_dir, remove_if_exists=True, copy_binary=False)

        self._save_checkpoint(ProcessingStage.REPORTS_COMPLETE, n_units=len(clean_units))

    def _apply_curation_logic(self, metrics, user_thresholds):
        """curation logic with rejection tracking"""
        # Default Thresholds
        defaults = {
            'presence_ratio': 0.75, 'rp_contamination': 0.15,
            'firing_rate': 0.05, 'amplitude_median': -20 ,# uV (negative peak)
            'amplitude_cv_median': 0.5
        }
        if user_thresholds: defaults.update(user_thresholds)
        
        keep_mask = np.ones(len(metrics), dtype=bool)
        rejections = []

        for idx, row in metrics.iterrows():
            reasons = []
            unit_id = row.name # unit_id is index

            # Apply Rules
            if 'presence_ratio' in row and row['presence_ratio'] < defaults['presence_ratio']:
                reasons.append(f"Presence < {defaults['presence_ratio']}")
            if 'rp_contamination' in row and row['rp_contamination'] > defaults['rp_contamination']:
                reasons.append(f"ISI Contam > {defaults['rp_contamination']}")
            if 'firing_rate' in row and row['firing_rate'] < defaults['firing_rate']:
                reasons.append(f"FR < {defaults['firing_rate']}")
            if 'amplitude_median' in row and row['amplitude_median'] > defaults['amplitude_median']:
                reasons.append(f"Amp > {defaults['amplitude_median']} (too small)")
            if 'amplitude_cv_median' in row and row['amplitude_cv_median'] >= defaults['amplitude_cv_median']:
                reasons.append(f"Amp CV >= {defaults['amplitude_cv_median']}")
            if reasons:
                keep_mask[metrics.index.get_loc(unit_id)] = False
                rejections.append({"unit_id": unit_id, "reasons": "; ".join(reasons)})

        return metrics[keep_mask], pd.DataFrame(rejections)

    def _plot_probe_locations(self, unit_ids, locations, filename):
            """spatial map visualization"""
            fig, ax = plt.subplots(figsize=(10.5, 6.5))
            # Plot the MEA hardware map
            si.plot_probe_map(self.recording, ax=ax, with_channel_ids=False)
            
            # Scatter plot of units
            # locations is (N, 2) or (N, 3)
            ax.scatter(locations[:, 0], locations[:, 1], s=10, c='blue', alpha=0.6)
            
            ax.invert_yaxis() # Maxwell convention
            ax.set_title(f"Unit Locations (n={len(unit_ids)})")
            fig.savefig(self.output_dir / filename)
            plt.close(fig)

    def _plot_waveforms_grid(self, unit_ids):
        """Generates the multi-page PDF of waveforms"""
        pdf_path = self.output_dir / "waveforms_grid.pdf"
        self.logger.info(f"Generating PDF: {pdf_path}")
        
        wf_ext = self.analyzer.get_extension("waveforms")
        
        with pdf.PdfPages(pdf_path) as pdf_doc:
            units_per_page = 12 # 4 cols x 3 rows
            for i in range(0, len(unit_ids), units_per_page):
                batch = unit_ids[i : i + units_per_page]
                fig, axes = plt.subplots(3, 4, figsize=(12, 9))
                axes = axes.flatten()
                
                for ax, uid in zip(axes, batch):
                    # Get waveforms on best channel
                    # Note: Simplified plotting for robustness
                    wf = wf_ext.get_waveforms_one_unit(uid)
                    mean_wf = np.mean(wf, axis=0) # (time, channels)
                    best_ch = np.argmin(np.min(mean_wf, axis=0)) # neg peak
                    
                    # Plot 50 random traces + mean
                    ax.plot(wf[:50, :, best_ch].T, c='gray', alpha=0.3, lw=0.5)
                    ax.plot(mean_wf[:, best_ch], c='red', lw=1.5)
                    ax.set_title(f"Unit {uid} (Ch {best_ch})")
                    ax.axis('off')
                
                # Clear unused axes
                for j in range(len(batch), len(axes)): axes[j].axis('off')
                
                pdf_doc.savefig(fig)
                plt.close(fig)

    def _run_burst_analysis(self, ids_list=None):
            """
            Runs burst analysis. 
            Adapted to work with either Unit IDs (Sorting) or Channel IDs (Detection).
            """
            self.logger.info("Running Network Burst Analysis...")
            
            fs = self.recording.get_sampling_frequency()
            spike_times = {}
            
            # Logic Branch: Sorter vs Detection Only
            if self.sorting is not None:
                # Case A: We have sorted units
                if ids_list is None: ids_list = self.analyzer.unit_ids
                for uid in ids_list:
                    spike_times[uid] = self.sorting.get_unit_spike_train(uid) / fs
                # Save for record
                np.save(self.output_dir / "spike_times.npy", spike_times)
            else:
                # Case B: We skipped sorting, check for saved detection file
                spike_file = self.output_dir / "spike_times.npy"
                if spike_file.exists():
                    self.logger.info("Loading spike times from detection step...")
                    # Allow pickle=True because we saved a dictionary
                    spike_times = np.load(spike_file, allow_pickle=True).item()
                else:
                    self.logger.error("No sorting object and no spike_times.npy found.")
                    return

            # Check if we have data
            if not spike_times:
                self.logger.warning("No spikes found to analyze.")
                return

            # --- The rest is identical to your original logic ---
            try:
                # A. Basic Burst Stats
                burst_stats = helper.detect_bursts_statistics(spike_times, isi_threshold=0.1)
                bursts = [x['bursts'] for x in burst_stats.values()]
                
                # B. Plot Raster (Re-using helper)
                fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
                
                title_suffix = "(Sorted Units)" if self.sorting else "(Channel MUA)"
                helper.plot_raster_with_bursts(axs[0], spike_times, bursts, title_suffix=title_suffix)
                
                # C. Network Burst (Parameter Free)
                # 
                network_data = compute_network_bursts(
                    ax_raster=None, 
                    ax_macro=axs[1], 
                    SpikeTimes=spike_times, 
                    plot=True
                )
                
                # Add metadata
                network_data['file'] = str(self.relative_pattern)
                network_data['well'] = self.stream_id
                network_data['analysis_type'] = "sorted_units" if self.sorting else "channel_mua"
                
                plt.tight_layout()
                plt.savefig(self.output_dir / "raster_burst_plot.svg")
                plt.close(fig)

                # D. Save JSON
                with open(self.output_dir / "network_results.json", 'w') as f:
                    class NpEncoder(json.JSONEncoder):
                        def default(self, obj):
                            if isinstance(obj, np.integer): return int(obj)
                            if isinstance(obj, np.floating): return float(obj)
                            if isinstance(obj, np.ndarray): return obj.tolist()
                            return super(NpEncoder, self).default(obj)
                    json.dump(network_data, f, cls=NpEncoder, indent=2)

            except Exception as e:
                self.logger.error(f"Burst Analysis Failed: {e}")
                self.logger.error(traceback.format_exc())

    def cleanup(self):
        """Frees memory and deletes temp files"""
        if self.cleanup_flag:
            self.logger.info("Cleaning up temp files...")
            shutil.rmtree(self.output_dir / "binary", ignore_errors=True)
            shutil.rmtree(self.output_dir / "sorter_output", ignore_errors=True)
            shutil.rmtree(self.output_dir / "preprocessed.zarr", ignore_errors=True)

        # Release RAM
        self.recording = None
        self.sorting = None
        self.analyzer = None
        gc.collect()
        # GPU  Memeory Clear
        try:
            import torch
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        except ImportError:
            pass
        
        try:
            import cupy
            cupy.get_default_memory_pool().free_all_blocks()
        except ImportError:
            pass


# --- CLI Entry Point ---
def main():
    parser = argparse.ArgumentParser(description="SOTA 2D MEA Pipeline")
    parser.add_argument('file_path', help='Path to .h5 or .nwb file')
    parser.add_argument('--well', required=True, help='Stream ID (e.g. well000)')
    parser.add_argument('--output-dir', required=True, help='Root output directory')
    parser.add_argument('--checkpoint-dir', help='Custom checkpoint location')
    parser.add_argument('--docker', help='Docker image name (e.g. spikeinterface/kilosort4-base)')
    parser.add_argument('--sorter', default='kilosort4', help='Sorter algorithm')
    parser.add_argument('--params', help='JSON string/file for curation thresholds')
    parser.add_argument('--rec', default='rec0000', help='Recording number within file')

    
    # Flags
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--clean-up', action='store_true')
    parser.add_argument('--force-restart', action='store_true')
    parser.add_argument('--export-to-phy', action='store_true')
    parser.add_argument('--no-curation', action='store_true')
    parser.add_argument('--skip-spikesorting', action='store_true')
    

    args = parser.parse_args()

    # Load params
    thresholds = None
    if args.params:
        if os.path.exists(args.params):
            with open(args.params) as f: thresholds = json.load(f)
        else:
            thresholds = json.loads(args.params)

    try:
        # Instantiate Pipeline
        pipeline = MEAPipeline(
            file_path=args.file_path,
            stream_id=args.well,
            recording_num=args.rec,
            output_root=args.output_dir,
            checkpoint_root=args.checkpoint_dir,
            sorter=args.sorter,
            docker_image=args.docker,
            verbose=args.debug,
            cleanup=args.clean_up,
            force_restart=args.force_restart
        )

        if pipeline.should_skip():
            sys.exit(0)

        # Execute Stages
        pipeline.run_preprocessing()

        if not args.skip_spikesorting:
            pipeline.run_sorting()
            pipeline.run_analyzer()
            pipeline.generate_reports(
                thresholds=thresholds,
                no_curation=args.no_curation,
                export_phy=args.export_to_phy
            )
        else :
            # Alternate Path: Spike Detection Only with electrode-level analysis
            # 1. Detect spikes on channels
            channel_ids = pipeline._spike_detection_only()
            
            # 2. Run burst analysis using those channels as "units"
            pipeline._run_burst_analysis(ids_list=channel_ids)



        
        # Cleanup
        pipeline.cleanup()
        print(f"Processing Complete for {args.well}")

    except Exception as e:
        print(f"CRITICAL FAILURE in {args.well}: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()


