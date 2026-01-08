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
from math import floor
import psutil
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf

from spikeinterface.sortingcomponents.peak_detection import detect_peaks

# Scientific Libraries
import spikeinterface.full as si
import spikeinterface.preprocessing as spre

# Import your custom modules
try:
    from parameter_free_burst_detector import compute_network_bursts
    import helper_functions as helper
    from scalebury import add_scalebar
except ImportError:
    print("Warning: Custom helper modules not found. Analysis steps may fail.")

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer): return int(obj)
        if isinstance(obj, np.floating): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        return super(NpEncoder, self).default(obj)
    
# --- Enums for Checkpointing ---
class ProcessingStage(Enum):
    NOT_STARTED = 0
    PREPROCESSING = 1
    PREPROCESSING_COMPLETE = 2
    SORTING = 3
    SORTING_COMPLETE = 4
    ANALYZER = 5
    ANALYZER_COMPLETE = 6
    REPORTS = 7
    REPORTS_COMPLETE = 8



# --- The Main Pipeline Class ---
class MEAPipeline:
    """
    SOTA Pipeline for 2D MEA Analysis (Maxwell Biosystems).
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
        self.well = self.metadata.get('well', 'UnknownWell')
        self.chip_id = self.metadata.get('chip_id', 'UnknownChip')
        self.date = self.metadata.get('date', 'UnknownDate')

        
        # Define Directory Structure: Output / Pattern / Well
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
            fh = logging.FileHandler(log_file, mode='a') # Append to log per run
            #add a big header
            fh.stream.write("\n" + "="*80 + "\n")
            fh.setFormatter(formatter)
            logger.addHandler(fh)
            ch = logging.StreamHandler(sys.stdout)
            ch.setFormatter(formatter)
            logger.addHandler(ch)
        return logger

    def _parse_metadata(self):
    
        meta = {
            'run_id': None, 'chip_id': None, 'project': None, 
            'relative_pattern': f"{self.file_path.parent.parent.name}/{self.file_path.parent.name}/{self.file_path.name}",
            'date': None, 'well': None
        }
        
        # Strategy A: Regex on path (Fallback)
        try:
            path_str = str(self.file_path)
            match = re.search(r"/(\d+)/data.raw.h5", path_str)
            if match: meta['run_id'] = match.group(1)
            parts = path_str.split(os.sep)
            if len(parts) > 5:
                meta['relative_pattern'] = os.path.join(*parts[-6:-1])
                meta['project'] = parts[-6]
                meta['date'] = parts[-5]
                meta['chip_id'] = parts[-4]
                meta['well'] = self.stream_id
        except Exception: pass

        # Strategy B: .metadata file (Overrides regex)
        meta_file = self.file_path.parent / ".metadata"
        if meta_file.exists():
            try:
                cfg = configparser.ConfigParser()
                cfg.read(meta_file, encoding='utf-8') 
                if 'properties' in cfg:
                    meta['run_id'] = cfg['properties'].get('runid', meta.get('run_id'))
                    meta['project'] = cfg['properties'].get('project_title', meta.get('project'))
                if 'runtime' in cfg:
                    meta['chip_id'] = cfg['runtime'].get('chipid', meta.get('chip_id'))
            except: pass
        return meta

    # --- Checkpointing Methods ---
    def _load_checkpoint(self):
        if self.checkpoint_file.exists() and not self.force_restart:
            with open(self.checkpoint_file, 'r') as f: return json.load(f)
        return {'stage': ProcessingStage.NOT_STARTED.value, #last stage completed
                'failed_stage': None,
                'last_updated': None,
                'run_id': self.run_id,
                'chip_id': self.chip_id,
                'well': self.well,
                'project': self.project_name,
                'date': self.date,
                'output_dir': str(self.output_dir),
                'error': None,
                }

    def _save_checkpoint(self, stage, **kwargs):
        self.state['stage'] = stage.value
        self.state['last_updated'] = str(datetime.now())
        self.state.update(kwargs)
        with open(self.checkpoint_file, 'w') as f: json.dump(self.state, f, indent=2)
        self.logger.info(f"Checkpoint Saved: {stage.name}")

    def should_skip(self):
        if self.state['stage'] == ProcessingStage.REPORTS_COMPLETE.value and not self.force_restart:
            self.logger.info("Pipeline already completed. Skipping.")
            return True
        return False

    def _load_recording_file(self):
        fpath = str(self.file_path)
        if fpath.endswith(".h5"):
            # Maxwell loader does NOT accept load_if_exists
            return si.read_maxwell(fpath, stream_id=self.stream_id, rec_name=self.recording_num)
        elif fpath.endswith(".nwb"):
            return si.read_nwb(fpath, load_if_exists=True)
        elif self.file_path.is_dir():
            return si.load_extractor(self.file_path)
        raise ValueError(f"Unknown format: {fpath}")

    # --- Phase 1: Preprocessing ---
    def run_preprocessing(self):
        binary_folder = self.output_dir / "binary"

        # 1. Resume Check
        if self.state['stage'] >= ProcessingStage.PREPROCESSING_COMPLETE.value and binary_folder.exists():
            self.logger.info("Resuming: Loading preprocessed data from binary cache.")
            try: self.recording = si.load(binary_folder)
            except: self.recording = si.load_extractor(binary_folder)
            return 
        self._save_checkpoint(ProcessingStage.PREPROCESSING)
        self.logger.info("--- [Phase 1] Preprocessing ---")
        
        # 2. Universal Loading
        rec = self._load_recording_file()

        # 3. Time Slicing (Remove last 1s to avoid end-of-file artifacts)
        fs = rec.get_sampling_frequency()
        self.metadata['fs'] = fs
        total_frames = rec.get_num_frames()
        # Calculate new end frame flooring to nearest integer second
        end_frame = floor(total_frames)
        #end_frame = total_frames - int(1.0 * fs)
        if end_frame > 0:
            self.logger.info(f"Trimming recording: {total_frames} -> {end_frame} frames (removed last 1s).")
            rec = rec.frame_slice(start_frame=0, end_frame=end_frame)

        if rec.get_dtype().kind == 'u':
            rec = spre.unsigned_to_signed(rec)


        rec = spre.highpass_filter(rec, freq_min=300)

        # Local Common Median Reference (Preserves Network Bursts)
        try:
            rec = spre.common_reference(rec, reference='local', operator='median', local_radius=(250, 250))
        except:
            self.logger.warning("Local CMR failed (missing locations?), using Global CMR.")
            rec = spre.common_reference(rec, reference='global', operator='median')
        
        rec.annotate(is_filtered=True)

        # F. Force Float32 (Prevents Signal Crushing)
        # Int16 scaling caused Kilosort crash on artifacts. Float32 is safer.
        if rec.get_dtype() != 'float32':
            self.logger.info("Converting to float32 to preserve signal fidelity...")
            rec = spre.astype(rec, 'float32')

        # 5. Save to BINARY
        if binary_folder.exists():
            shutil.rmtree(binary_folder)
        
        self.logger.info(f"Saving binary recording to {binary_folder}...")
        
        rec.save(
            folder=binary_folder, 
            format='binary', 
            overwrite=True, 
            n_jobs=16, 
            chunk_duration='1s', 
            progress_bar=self.verbose
        )
        
        # 6. Reload from Binary (Memory Map)
        self.recording = si.load(binary_folder)
        self._save_checkpoint(ProcessingStage.PREPROCESSING_COMPLETE)

    # --- Phase 2: Sorting ---
    def run_sorting(self):
        sorter_folder = self.output_dir / "sorter_output"
        if self.state['stage'] >= ProcessingStage.SORTING_COMPLETE.value:
            self.logger.info("Resuming: Loading existing sorting.")
            try: self.sorting = si.read_sorter_folder(sorter_folder)
            except: self.sorting = si.read_kilosort(sorter_folder)
            return
        self._save_checkpoint(ProcessingStage.SORTING)
        self.logger.info(f"--- [Phase 2] Spike Sorting ({self.sorter}) ---")
        
        import torch
        if torch.cuda.is_available():
            total_vram = torch.cuda.get_device_properties(0).total_memory / (1024**3) # Convert to GB
            self.logger.info(f"GPU Detected: {torch.cuda.get_device_name(0)} with {total_vram:.2f} GB VRAM")
        else:
            self.logger.warning("No GPU detected! Kilosort4 will likely fail or run extremely slowly on CPU.")
            total_vram = 0

        fs = int(self.recording.get_sampling_frequency())

        # Parameters for 2D Static Cultures (Kilosort4)
        # Th_universal lowered to 8 to improve detection on organoids
        ks_params_high_vram ={
                        'batch_size': int(self.recording.get_sampling_frequency()) * 2,
                        'clear_cache': True,
                        'invert_sign': True,
                        'cluster_downsampling': 20,
                        'max_cluster_subset': None,
                        'nblocks':0,
                        'dmin':17,
                        'do_correction': False,
                    }
        ks_params_low_vram = {
                        'batch_size': int(self.recording.get_sampling_frequency() * 0.5),
                        'clear_cache': True,
                        'invert_sign': True,
                        'cluster_downsampling': 30, 
                        'max_cluster_subset': 50000,
                        'nblocks':0,
                        'do_correction': False,
        }

        if total_vram >= 14:
            ks_params = ks_params_high_vram
        else:
            ks_params = ks_params_low_vram

        start = timer()
        try:
            self.logger.info(f"Running Kilosort4 with parameters: {ks_params}")
            torch.cuda.empty_cache()
            self.sorting = si.run_sorter(
                sorter_name=self.sorter,
                recording=self.recording,
                folder=sorter_folder,
                delete_output_folder=False,
                remove_existing_folder=True, # Safety flag
                verbose=self.verbose,
                docker_image=self.docker_image, 
                **ks_params
            )
            self.logger.info(f"Sorting finished in {timer()-start:.2f}s")

            #self.logger.info("Running Redundancy Removal...")
            #self.sorting = si.remove_redundant_units(self.sorting, duplicate_threshold=0.9)
            self._save_checkpoint(ProcessingStage.SORTING_COMPLETE, failed_stage=None, error=None)
        except Exception as e:
            err = {
                "failed_stage": ProcessingStage.SORTING.name,
                "exception": type(e).__name__,
                "message": str(e),
                "traceback": traceback.format_exc(),
                "time": str(datetime.now())
            }
            self.logger.error(err["traceback"])
            self._save_checkpoint(ProcessingStage.PREPROCESSING_COMPLETE, error=err)
            raise

    def _spike_detection_only(self):
        """Detects spikes (threshold crossings) without sorting."""
        self.logger.info("--- [Phase 2-Alt] Spike Detection (No Sorting) ---")
        job_kwargs = {'n_jobs': 16, 'chunk_duration': '1s', 'progress_bar': self.verbose}
        
        peaks = detect_peaks(
            self.recording, method='by_channel', 
            detect_threshold=5, peak_sign='neg', exclude_sweep_ms=0.1, **job_kwargs
        )

        
        self.logger.info(f"Detected {len(peaks)} total spikes.")
        fs = self.recording.get_sampling_frequency()
        self.metadata['fs'] = fs
        channel_ids = self.recording.get_channel_ids()
        spike_times = {}

        for ch_index, ch_id in enumerate(channel_ids):
            mask = peaks['channel_index'] == ch_index
            if np.any(mask):
                spike_times[ch_id] = peaks['sample_index'][mask] / fs
        
        np.save(self.output_dir / "spike_times.npy", spike_times)
        
        # Quick raster plot
        fig, ax = plt.subplots(figsize=(12, 8))
        sorted_chs = sorted(spike_times.keys())
        for i, ch in enumerate(sorted_chs):
            times = spike_times.get(ch, [])
            if len(times) > 0:
                ax.plot(times, np.ones_like(times) * i, '|', markersize=1, color='black', alpha=0.5)
        ax.set_title(f"MUA Raster - {self.stream_id}")
        plt.savefig(self.output_dir / "mua_channel_raster.svg", dpi=300)
        #also save a 30s zoom
        ax.set_xlim(0,30)
        plt.savefig(self.output_dir / "mua_channel_raster_30s.svg", dpi=300)
        plt.close(fig)

        return list(spike_times.keys())

    # --- Phase 3: Analyzer ---
    def run_analyzer(self):
        analyzer_folder = self.output_dir / "analyzer_output"
        if self.state['stage'] >= ProcessingStage.ANALYZER_COMPLETE.value and  analyzer_folder.exists():
            self.logger.info("Resuming: Loading Sorting Analyzer.")
            self.analyzer = si.load_sorting_analyzer(analyzer_folder)
            return
        
        self._save_checkpoint(ProcessingStage.ANALYZER)
        self.logger.info("--- [Phase 3] Computing Sorting Analyzer ---")
        try:
            if analyzer_folder.exists(): shutil.rmtree(analyzer_folder)

            sparsity = si.estimate_sparsity(
                self.sorting, self.recording, 
                method="radius", radius_um=50, peak_sign='neg'
            )

            self.analyzer = si.create_sorting_analyzer(
                self.sorting, self.recording, format="binary_folder", 
                folder=analyzer_folder, sparsity=sparsity, return_in_uV=True
            )

            ext_list = ["random_spikes","spike_amplitudes", "waveforms", "templates", "noise_levels", 
                        "quality_metrics", "template_metrics", "unit_locations"]
            
            ext_params = {
                "waveforms": {"ms_before": 1.0, "ms_after": 2.0},
                "unit_locations": {"method": "monopolar_triangulation"}
            }
            
            self.analyzer.compute(ext_list, extension_params=ext_params, verbose=self.verbose)
            self._save_checkpoint(ProcessingStage.ANALYZER_COMPLETE,failed_stage=None, error=None)
        except Exception as e:
            err = {
                "failed_stage": ProcessingStage.ANALYZER.name,
                "exception": type(e).__name__,
                "message": str(e),
                "traceback": traceback.format_exc(),
                "time": str(datetime.now())
            }
            self.logger.error(err["traceback"])
            self._save_checkpoint(ProcessingStage.SORTING_COMPLETE, error=err)
            raise

    # --- Phase 4: Reports & Curation ---
    def generate_reports(self, thresholds=None, no_curation=False, export_phy=False):
        if self.state['stage'] == ProcessingStage.REPORTS_COMPLETE.value: return

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
            self._plot_probe_locations(clean_units, locations[mask], f"locations_{len(clean_units)}_units.pdf")
            self._plot_waveforms_grid(clean_units)
            self._run_burst_analysis(clean_units)

            if export_phy:
                si.export_to_phy(self.analyzer.select_units(clean_units), 
                            output_folder=self.output_dir/"phy_output", remove_if_exists=True, copy_binary=False)

            self._save_checkpoint(ProcessingStage.REPORTS_COMPLETE, n_units=len(clean_units),failed_stage=None, error=None)
        except Exception as e:
            err = {
                "failed_stage": ProcessingStage.REPORTS.name,
                "exception": type(e).__name__,
                "message": str(e),
                "traceback": traceback.format_exc(),
                "time": str(datetime.now())
            }
            self.logger.error(err["traceback"])
            self._save_checkpoint(ProcessingStage.ANALYZER_COMPLETE, error=err)
            raise

    def _apply_curation_logic(self, metrics, user_thresholds):
        defaults = {'presence_ratio': 0.75, 'rp_contamination': 0.15, 'firing_rate': 0.05, 
                    'amplitude_median': -20, 'amplitude_cv_median': 0.5}
        if user_thresholds: defaults.update(user_thresholds)
        
        keep_mask = np.ones(len(metrics), dtype=bool)
        rejections = []
        for idx, row in metrics.iterrows():
            reasons = []
            if row.get('presence_ratio', 1) < defaults['presence_ratio']: reasons.append("Low Presence")
            if row.get('rp_contamination', 0) > defaults['rp_contamination']: reasons.append("High Contam")
            if row.get('firing_rate', 0) < defaults['firing_rate']: reasons.append("Low FR")
            if row.get('amplitude_median', -100) > defaults['amplitude_median']: reasons.append("Low Amp")
            
            if reasons:
                keep_mask[metrics.index.get_loc(row.name)] = False
                rejections.append({"unit_id": row.name, "reasons": "; ".join(reasons)})
        
        return metrics[keep_mask], pd.DataFrame(rejections)

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
        
        # [1] Get sampling freq to calculate milliseconds
        fs = self.recording.get_sampling_frequency()
        
        with pdf.PdfPages(pdf_path) as pdf_doc:
            units_per_page = 12 
            for i in range(0, len(unit_ids), units_per_page):
                batch = unit_ids[i : i + units_per_page]
                fig, axes = plt.subplots(3, 4, figsize=(12, 9))
                axes = axes.flatten()
                
                for ax, uid in zip(axes, batch):
                    wf = wf_ext.get_waveforms_one_unit(uid)
                    mean_wf = np.mean(wf, axis=0)
                    best_ch = np.argmin(np.min(mean_wf, axis=0))
                    
                    # [2] Create time vector in milliseconds
                    # (NumSamples / fs) * 1000
                    time_ms = np.arange(wf.shape[1]) / fs * 1000

                    # Random subset logic
                    n_spikes = wf.shape[0]
                    if n_spikes > 50:
                        indices = np.random.choice(n_spikes, 50, replace=False)
                        spikes_to_plot = wf[indices, :, best_ch]
                    else:
                        spikes_to_plot = wf[:, :, best_ch]
                    
                    # [3] Plot against time_ms (X-axis)
                    ax.plot(time_ms, spikes_to_plot.T, c='gray', lw=0.5, alpha=0.3)
                    ax.plot(time_ms, mean_wf[:, best_ch], c='red', lw=1.5)
                    
                    ax.set_title(f"Unit {uid} | Ch {best_ch}", fontsize=10)
                    
                    # [4] Add Scale Bar (1 ms x 50 uV)
                    # hidex/hidey=True removes the standard axis ticks
                    try:
                        add_scalebar(ax, 
                                    matchx=False, matchy=False, 
                                    sizex=1.0, labelx='1 ms', 
                                    sizey=50, labely='50 µV', 
                                    loc='lower right', 
                                    hidex=True, hidey=True)
                    except Exception:
                        # Fallback if scalebury is missing or fails
                        ax.spines['top'].set_visible(False)
                        ax.spines['right'].set_visible(False)

                for j in range(len(batch), len(axes)): 
                    axes[j].axis('off')
                
                pdf_doc.savefig(fig)
                plt.close(fig)

    def _run_burst_analysis(self, ids_list=None):
            self.logger.info("Running Network Burst Analysis...")
            
            spike_times = {}
            
            # 1. Load Spike Times
            if self.sorting:
                fs = self.recording.get_sampling_frequency()
                if ids_list is None: ids_list = self.analyzer.unit_ids
                for uid in ids_list:
                    spike_times[uid] = self.sorting.get_unit_spike_train(uid) / fs
                np.save(self.output_dir / "spike_times.npy", spike_times)
            else:
                npy_spike_file = self.output_dir / "spike_times.npy"
                if npy_spike_file.exists():
                    spike_times = np.load(npy_spike_file, allow_pickle=True).item()
                else:
                    self.logger.error("No spike times found for burst analysis.")
                    return

            if not spike_times: return

            # 2. Main Analysis Block
            try:
                # A. Statistics & Raster Plot
                burst_stats = helper.detect_bursts_statistics(spike_times, isi_threshold=0.1)
                bursts = [x['bursts'] for x in burst_stats.values()]
                
                fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
                #title_suffix = "(Sorted Units)" if self.sorting else "(Channel MUA)"
                title_suffix = f" "
                helper.plot_raster_with_bursts(axs[0], spike_times, bursts, title_suffix=title_suffix)
                
                # B. Network Burst Calculation
                network_data = compute_network_bursts(ax_raster=None, ax_macro=axs[1], SpikeTimes=spike_times, plot=True)
                network_data['file'] = str(self.relative_pattern)
                network_data['well'] = self.stream_id
                
                plt.tight_layout()
                plt.subplots_adjust(hspace=0.05)
                plt.savefig(self.output_dir / "raster_burst_plot.svg")
                # Save 30s zoom
                axs[0].set_xlim(0, 30)
                axs[1].set_xlim(0, 30)
                plt.savefig(self.output_dir / "raster_burst_plot_30s.svg")
                plt.close(fig)
                
                # C. Robust JSON Saving (Handles numpy types AND dictionary keys)
                def recursive_clean(obj):
                    """Recursively converts numpy types and keys to Python standard types."""
                    if isinstance(obj, dict):
                        new_dict = {}
                        for k, v in obj.items():
                            # Force keys to string (JSON requirement)
                            clean_k = str(k)
                            new_dict[clean_k] = recursive_clean(v)
                        return new_dict
                    if isinstance(obj, list):
                        return [recursive_clean(v) for v in obj]
                    if isinstance(obj, (np.integer, int)):
                        return int(obj)
                    if isinstance(obj, (np.floating, float)):
                        return float(obj)
                    if isinstance(obj, np.ndarray):
                        return recursive_clean(obj.tolist())
                    return obj

                self.logger.info(f"Sanitizing JSON data for {len(network_data.get('unit_bursts', []))} units...")
                
                # Clean data in memory
                network_data_clean = recursive_clean(network_data)

                # Atomic Write (Temp -> Rename) to prevent corruption
                temp_file = self.output_dir / "network_results.tmp.json"
                final_file = self.output_dir / "network_results.json"

                with open(temp_file, 'w') as f:
                    json.dump(network_data_clean, f, indent=2)
                
                if temp_file.exists():
                    os.replace(temp_file, final_file)
                    self.logger.info(f"Successfully saved: {final_file}")

            except Exception as e:
                self.logger.error(f"Burst analysis error: {e}")
                traceback.print_exc()
                raise e

    # --- Cleanup ---
    def cleanup(self):
        if self.cleanup_flag:
            self.logger.info("Cleaning up temp files...")
            shutil.rmtree(self.output_dir / "binary", ignore_errors=True)
            shutil.rmtree(self.output_dir / "sorter_output", ignore_errors=True)
        self.recording = None
        self.sorting = None
        self.analyzer = None
        gc.collect()

# --- CLI Entry Point ---
def main():
    parser = argparse.ArgumentParser(description="SOTA 2D MEA Pipeline")
    parser.add_argument('file_path')
    parser.add_argument('--well', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--checkpoint-dir')
    parser.add_argument('--docker')
    parser.add_argument('--sorter', default='kilosort4')
    parser.add_argument('--params')
    parser.add_argument('--rec', default='rec0000')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--clean-up', action='store_true')
    parser.add_argument('--force-restart', action='store_true')
    parser.add_argument('--export-to-phy', action='store_true')
    parser.add_argument('--no-curation', action='store_true')
    parser.add_argument('--skip-spikesorting', action='store_true')
    parser.add_argument("--reanalyze-bursts", action="store_true", help="Only re-run burst analysis on existing spike times.")

    args = parser.parse_args()
    
    thresholds = None
    if args.params:
        if os.path.exists(args.params):
            with open(args.params) as f: thresholds = json.load(f)
        else: thresholds = json.loads(args.params)

    try:
        pipeline = MEAPipeline(
            file_path=args.file_path, stream_id=args.well, recording_num=args.rec,
            output_root=args.output_dir, checkpoint_root=args.checkpoint_dir,
            sorter=args.sorter, docker_image=args.docker, verbose=args.debug,
            cleanup=args.clean_up, force_restart=args.force_restart
        )

        if args.reanalyze_bursts:
            print("Re-analyzing bursts only on existing spike times...")
            pipeline._run_burst_analysis()
            current_stage = ProcessingStage(pipeline.state['stage'])
            pipeline._save_checkpoint(current_stage, note="Burst Re-analysis Performed",last_updated=str(datetime.now()))
            print("Burst Re-analysis Complete.")
            sys.exit(0)


        if pipeline.should_skip(): sys.exit(0)

        pipeline.run_preprocessing()

        if not args.skip_spikesorting:
            pipeline.run_sorting()
            pipeline.run_analyzer()
            pipeline.generate_reports(thresholds, args.no_curation, args.export_to_phy)
        else:
            ids = pipeline._spike_detection_only()
            pipeline._run_burst_analysis(ids)

        pipeline.cleanup()
        print(f"Processing Complete for {args.well}")

    except Exception as e:
        print(f"CRITICAL FAILURE in {args.well}: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()