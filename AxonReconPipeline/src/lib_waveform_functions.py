import sys
import os
import h5py
import spikeinterface.full as si
import shutil
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from MEAProcessingLibrary import mea_processing_library as MPL
import AxonReconPipeline.src.lib_sorting_functions as sorter

# Function to extract waveforms for a specific unit
def extract_unit_waveforms(h5_path, stream_id, segment_sorting, save_root=None, logger=None, te_params={}, **wf_kwargs):
    n_jobs = te_params.get('n_jobs', 4)
    overwrite_wf = te_params.get('overwrite_wf', False)
    
    # Helper function to get assay information
    def get_assay_information(rec_path):
        with h5py.File(rec_path, 'r') as h5:
            pre, post, well_id = -1, -1, 0
            while pre <= 0 or post <= 0:
                well_name = list(h5['wells'].keys())[well_id]
                rec_name = list(h5['wells'][well_name].keys())[well_id]
                pre = h5['wells'][well_name][rec_name]['groups']['routed']['trigger_pre'][0]
                post = h5['wells'][well_name][rec_name]['groups']['routed']['trigger_post'][0]
                well_id += 1
        return [pre, post]

    # Helper function to load waveforms
    def load_waveforms(wf_path, seg_sort):
        try:
            seg_we = si.WaveformExtractor.load(wf_path, with_recording=True, sorting=seg_sort)
            if len(seg_sort.get_unit_ids()) == len(seg_we.unit_ids):
                if logger is not None: logger.debug(f'Waveforms loaded from {wf_path} already exist')
                return wf_path, seg_we
            else:
                raise ValueError('Unit IDs do not match')
        except Exception as e:
            if os.path.exists(wf_path):
                shutil.rmtree(wf_path)
            if logger is not None: logger.error(f"Error loading waveforms: {e}")
            return None, None

    # Main function to extract waveforms
    def extract_waveforms(h5_path, segment_sorting, stream_id, save_root, n_jobs, overwrite_wf):
        cutout = [x / (segment_sorting.get_sampling_frequency() / 1000) for x in get_assay_information(h5_path)]
        h5_details = MPL.extract_recording_details(h5_path)
        date, chip_id, scanType, run_id = h5_details[0]['date'], h5_details[0]['chipID'], h5_details[0]['scanType'], h5_details[0]['runID']
        
        segment_waveforms = {}
        with h5py.File(h5_path, 'r') as h5:
            rec_names = list(h5['wells'][stream_id].keys())
            wf_paths = [os.path.join(save_root, f'waveforms/{date}/{chip_id}/{scanType}/{run_id}/{stream_id}/seg{sel_idx}') for sel_idx in range(len(rec_names))]
            seg_sorts = []
            logger.info(f'Extracting waveforms for stream {stream_id}')
            for sel_idx, wf_path in enumerate(wf_paths):
                try: seg_sorts.append((wf_path, si.SelectSegmentSorting(segment_sorting, sel_idx)))
                except: pass #deal with erroneus recording segments

            if not overwrite_wf:
                if logger is not None: logger.info(f'Loading waveforms for stream {stream_id}')
                with ThreadPoolExecutor(max_workers=n_jobs) as executor:
                    futures = {executor.submit(load_waveforms, seg_sort[0], seg_sort[1]): sel_idx for sel_idx, seg_sort in enumerate(seg_sorts)}

                try:
                    for future in futures:
                        sel_idx = futures[future]
                        try:
                            wf_path, seg_we = future.result()
                            if seg_we is not None:
                                segment_waveforms[rec_names[sel_idx]] = {'path': wf_path, 'waveforms': seg_we}
                        except Exception as e:
                            if logger is not None: logger.error(f"Error loading waveforms for segment {sel_idx}: {e}")
                except Exception as e:
                    if 'Waveform folder does not exist' in str(e):
                        if logger is not None: logger.warning(f'Waveform folder does not exist, extracting waveforms')
                        else: print(f'Waveform folder does not exist, extracting waveforms')
                    if logger is not None: logger.error(f"Error loading waveforms: {e}")
                    else: print(f"Error loading waveforms: {e}")

            for sel_idx, rec_name in enumerate(rec_names):
                if rec_name in segment_waveforms: continue
                wf_path = wf_paths[sel_idx]
                try:
                    rec = si.MaxwellRecordingExtractor(h5_path, stream_id=stream_id, rec_name=rec_name)
                except Exception as e:
                    if logger is not None: logger.error(f"Error extracting recording for segment {sel_idx}: {e}")
                    continue

                chunk_size = min(10000, rec.get_num_samples()) - 100
                rec_centered = si.center(rec, chunk_size=chunk_size)
                seg_sort = si.SelectSegmentSorting(segment_sorting, sel_idx)
                seg_sort = si.remove_excess_spikes(seg_sort, rec_centered)
                seg_sort.register_recording(rec_centered)

                if not os.path.exists(wf_path) or overwrite_wf:
                    if logger is not None: logger.info(f'Extracting waveforms to {wf_path}, n_jobs={n_jobs}')
                    else: print(f'Extracting waveforms to {wf_path}, n_jobs={n_jobs}')
                    os.makedirs(wf_path, exist_ok=True)
                    seg_we = si.WaveformExtractor.create(rec_centered, seg_sort, wf_path, allow_unfiltered=True, remove_if_exists=True)
                    seg_we.set_params(ms_before=cutout[0], ms_after=cutout[1], return_scaled=True)
                    seg_we.run_extract_waveforms(n_jobs=n_jobs)
                    segment_waveforms[rec_name] = {'path': wf_path, 'waveforms': seg_we}
        return segment_waveforms

    waveforms = extract_waveforms(h5_path, segment_sorting, stream_id, save_root, n_jobs, overwrite_wf)
    return waveforms

# Function to select units
def select_units(sorting, min_n_spikes=500, exclude_mua=True, logger=None):
    if exclude_mua:
        ks_label = sorting.get_property('KSLabel')
        mua_idx = ks_label == 'mua'
    else:
        mua_idx = np.full((sorting.get_num_units(),), False, dtype='bool')

    n_spikes = [len(sorting.get_unit_spike_train(x)) for x in sorting.get_unit_ids()]

    bad_n_spikes_idx = np.array(n_spikes) < min_n_spikes
    bad_idx = mua_idx | bad_n_spikes_idx
    bad_id = [i for i, x in enumerate(bad_idx) if x]

    cleaned_sorting = sorting.remove_units(bad_id)
    if logger is not None: logger.debug(f'Removed units: {bad_id}')
    return cleaned_sorting