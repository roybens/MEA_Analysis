import os
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
import csv
import numpy as np
import matplotlib.pyplot as plt
from tsmoothie.smoother import GaussianSmoother
import spikeinterface.full as si
from parameter_free_burst_detector import compute_network_bursts
import helper_functions as helper
from pathlib import Path
from timeit import default_timer as timer
import multiprocessing
import sys
from datetime import datetime
import logging
import re
import pickle
import scipy.io as sio
import argparse
import pandas as pd
from spikeinterface.curation import CurationSorting
import pdb, traceback
from collections import defaultdict
import yaml
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve, find_peaks
from scipy.stats import norm
import json
import h5py
import psutil
from helper_functions import detect_peaks

from enum import Enum



#os.environ['HDF5_PLUGIN_PATH']='/home/mmp/Documents/Maxlab/so/'
# Configure the logger
# Manually clear the log file

def setup_logger(log_file=None):
    if log_file is None:
        log_file =f'./applicaiton_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
    #make dirs if not exist
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    open(log_file, 'w').close()  # Clear log file on each run

    logger = logging.getLogger("mea_pipeline")
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s - %(module)s - %(lineno)d')

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info(f"Logging to {log_file}")
    return logger



args =None

logger =None  #placeholder for global logger

def log_resource_usage():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    logger.debug(f"Memory usage: {mem_info.rss / 1024**3:.2f} GB")
    logger.debug(f"CPU usage: {process.cpu_percent()}%")

##DEFINING GLOBAL VARIABLES

BASE_FILE_PATH =  os.path.dirname(os.path.abspath(__file__))


job_kwargs = dict(n_jobs=32, chunk_duration="1s", progress_bar=False)

## CHECKPOINTING CODE (inspirred from Claude Sonnet)


class ProcessingStage(Enum):
    """Pipeline stages for checkpointing"""
    NOT_STARTED = 0
    SORTING_COMPLETE = 1
    ANALYZER_COMPLETE = 2
    REPORTS_COMPLETE = 3
    FAILED = -1

class PipelineCheckpoint:
    """Manage checkpoints for MEA processing pipeline"""
    
    def __init__(self, base_path, project_name,date,chip_id,run_id, well_id, rec_name):
        self.run_id = run_id
        self.well_id = well_id
        self.rec_name = rec_name
        self.project_name = project_name
        self.date = date
        self.chip_id = chip_id
        self.checkpoint_dir = Path(base_path) / "checkpoints"
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self.checkpoint_file = self.checkpoint_dir / f"{project_name}_{date}_{chip_id}_{run_id}_{well_id}_{rec_name}_checkpoint.json"
        self.state = self.load_checkpoint()
    
    def load_checkpoint(self):
        """Load existing checkpoint or create new one"""
        if self.checkpoint_file.exists():
            with open(self.checkpoint_file, 'r') as f:
                return json.load(f)
        else:
            return {
                'project_name': self.project_name,
                'date': self.date,
                'chip_id': self.chip_id,
                'run_id': self.run_id,
                'well_id': self.well_id,
                'rec_name': self.rec_name,
                'stage': ProcessingStage.NOT_STARTED.value,
                'stage_name': ProcessingStage.NOT_STARTED.name,
                'analyzer_folder': None,
                'last_updated': None,
                'error': None
            }
    
    def save_checkpoint(self, stage, **kwargs):
        """Save checkpoint with current stage and metadata"""
        self.state['stage'] = stage.value
        self.state['stage_name'] = stage.name
        self.state['last_updated'] = str(datetime.now())
        
        # Update any additional metadata
        for key, value in kwargs.items():
            self.state[key] = value
        
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.state, f, indent=2)

        logger.info(f"Checkpoint saved:{self.project_name},{self.chip_id},{self.date},{self.run_id},{self.chip_id} {self.well_id}/{self.rec_name} - {stage.name}")

    def get_stage(self):
        """Get current processing stage"""
        return ProcessingStage(self.state['stage'])
    
    def is_completed(self):
        """Check if processing is fully completed"""
        return self.get_stage() == ProcessingStage.REPORTS_COMPLETE
    
    def should_skip(self, args):
        """Determine if processing should be skipped based on checkpoint"""
        if not args.resume:
            return False
        
        stage = self.get_stage()
        if stage == ProcessingStage.REPORTS_COMPLETE: 
            
            if args.force_restart:
                logger.info(f"Restarting {self.well_id}/{self.rec_name} from scratch")
                return False
            else:
                logger.info(f"Skipping {self.well_id}/{self.rec_name} - already completed")
                return True
        
        return False
        
    def mark_failed(self, error_msg, failed_at_stage=None):
        """
        Mark processing as failed while preserving last successful stage
        
        Args:
            error_msg: Description of the error
            failed_at_stage: Which stage failed (for logging)
        """
        # Keep current stage (last successful one)
        current_stage = self.get_stage()
        
        # Save with current stage preserved, just add error info
        self.save_checkpoint(
            current_stage,  # ✅ Preserve last successful stage
            failed=True,
            failed_at_stage=failed_at_stage,
            error=str(error_msg),
            error_traceback=traceback.format_exc()
        )
        
        logging.error(f"Processing failed at {failed_at_stage }")
        logging.error(f"Last successful stage: {current_stage}")








def get_channel_recording_stats(recording):
    
    channel_ids = recording.get_channel_ids()
    fs = recording.get_sampling_frequency()
    num_chan = recording.get_num_channels()
    num_seg = recording.get_num_segments()
    total_recording = recording.get_total_duration()

    #print('Channel ids:', channel_ids)
    print('Sampling frequency:', fs)
    print('Number of channels:', num_chan)
    print('Number of segments:', num_seg)
    print(f"total_recording: {total_recording} s")
    return fs,num_chan,channel_ids, total_recording

def preprocess(recording):
    global logger
    # 1. Convert to Signed (Essential for Maxwell)
    recording = si.unsigned_to_signed(recording)
    
    # 2. REMOVE DC OFFSET (Highpass Only)
    # Don't lowpass (freq_max). Let KS4 see the sharp spike details.
    # KS4 expects data centered at 0.
    recording = si.highpass_filter(recording, freq_min=300)

    # 3. LOCAL Common Reference (Critical for Superbursts)
    # Instead of 'global', use 'local'.
    # This removes local electrical noise without deleting network-wide bursts.
    # radius=250um is a good balance for Maxwell chips.
    try:
        recording = si.common_reference(
            recording, 
            reference='local', 
            operator='median', 
            local_radius=(250, 250) # radius in um
        )
    except Exception as e:
        # Fallback if channel locations aren't loaded correctly
        logger.warning(f"Local CMR failed ({e}), falling back to Global")
        recording = si.common_reference(recording, reference='global', operator='median')

    recording.annotate(is_filtered=True)

    # 4. Enforce float32 -> int16 conversion 
    # KS4 runs faster on int16. Maxwell Raw is float/int mix sometimes.
    # This scale step ensures the dynamic range fits into int16 without clipping.
    # (Optional but recommended if you aren't already handling dtype elsewhere)
    if recording.get_dtype() != 'int16':
        recording = si.scale(recording, dtype='int16')

    return recording


def resource_cleanup(clear_temp_files = False,file_path_list=None):
            
    import gc
    import shutil
    # Cleanup folders
    if clear_temp_files and file_path_list is not None:
        for file_path in file_path_list:
            if os.path.exists(file_path):
                shutil.rmtree(file_path)
                logger.info(f"Removed folder: {file_path}")
            else:
                logger.info(f"Folder does not exist: {file_path}")
    gc.collect()
    # Optional: free GPU memory
    try:
        import torch
        torch.cuda.empty_cache()
    except:
        pass
    try:
        import cupy
        cupy.get_default_memory_pool().free_all_blocks()
    except:
        pass

    return None,None 
  





def get_channel_locations_mapping(recording):

    channel_locations = recording.get_channel_locations()
    channel_ids = recording.get_channel_ids()
    channel_locations_mappings= {channel_id: location for location, channel_id in zip(channel_locations, channel_ids)}
    return channel_locations_mappings

def get_data_maxwell(file_path,stream_id,rec_num):

    rec_num =  str(rec_num).zfill(4)
    rec_name = 'rec' + rec_num
    recording = si.read_maxwell(file_path,rec_name=rec_name,stream_id=stream_id,install_maxwell_plugin=False)
    return recording,rec_name

def parse_h5_path(file_path):
    pattern = r"/(\d+)/data.raw.h5"
    run_id = int(re.search(pattern, file_path).group(1))
    file_pattern = os.path.dirname(file_path)
    parts = file_pattern.split('/')
    project_name =parts[-5] if len(parts)>=5 else 'UnknownProject'
    date = parts[-4] if len(parts)>=4 else 'UnknownDate'
    chip_id = parts[-3] if len(parts)>=3 else 'UnknownChip'

    desired_pattern = '/'.join(parts[-6:])
    return run_id,project_name,date,chip_id,desired_pattern


def process_block(file_path, time_in_s=None, stream_id='well000', recnumber=0, 
                  clear_temp_files=False, thresholds=None):
    

    global args,logger,BASE_FILE_PATH,job_kwargs

    # Setup paths and checkpoint manager
    if file_path.endswith('.raw.h5'):
        recording, rec_name = get_data_maxwell(file_path, stream_id, recnumber)
        logger.info(f"Processing recording: {file_path}")
        #file pattern extraction
        run_id,project_name,date,chip_id,desired_pattern = parse_h5_path(file_path)

    else:
        rec_name = "binary_recording"   #implement fr binary
        logger.error("Binary file processing not implemented yet")


    # Initialize checkpoint manager
    checkpoint = PipelineCheckpoint(
        base_path=f"{BASE_FILE_PATH}/../AnalyzedData",
        project_name=project_name,
        date=date,
        chip_id=chip_id,
        run_id=run_id,
        well_id=stream_id,
        rec_name=rec_name
    )
    
    # Check if should skip
    if checkpoint.should_skip(args):
        #resource cleanup Cleanup folders
        if clear_temp_files:
            logger.info("Clearing stored files if exists.")
            binary_folder = f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}/binary"
            sorting_folder = f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}/analyzer_output"
            resource_cleanup(clear_temp_files=True,file_path_list=[binary_folder,sorting_folder])
            
                  
    
    if args.force_restart:
        logger.info(f"Restarting {stream_id}/{rec_name} from scratch")
        #manually delete checkpoint file to restart
        if checkpoint.checkpoint_file.exists():
            checkpoint.checkpoint_file.unlink()
            checkpoint.state = checkpoint.load_checkpoint()  # Reset state

    try:
        # Get current stage to resume from
        current_stage = checkpoint.get_stage()
        logger.info(f"Starting from stage: {current_stage.name}")
        
        # ============== STAGE 0: Setup and Preprocessing ==============
        if current_stage.value < ProcessingStage.SORTING_COMPLETE.value:
            logger.info(f"[STAGE 0] Loading and preprocessing recording...")

            # Your existing preprocessing code
            if file_path.endswith('.raw.h5'):
                # Setup output paths
                run_id, project_name, date, chip_id, desired_pattern = parse_h5_path(file_path)
                   
                os.makedirs(
                    f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}",
                    mode=0o777,
                    exist_ok=True
                )

                recording, rec_name = get_data_maxwell(file_path, stream_id, recnumber)
                logger.info(f"Processing recording: {rec_name}")
                fs, num_chan, channel_ids, total_rec_time = get_channel_recording_stats(recording)
                
                if time_in_s is None:
                    time_in_s = total_rec_time - 1
                
                time_start = 0
                time_end = time_start + time_in_s
                recording_chunk = recording.frame_slice(
                    start_frame=int(time_start*fs),
                    end_frame=int(time_end*fs)
                )
                recording_chunk = preprocess(recording_chunk)
                
                recording_chunk.save(
                    folder=f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}/binary",
                    format='binary',
                    overwrite=True,
                    **job_kwargs
                )
                recording_chunk = si.load(
                    f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}/binary"
                )
            else:
                recording_chunk = si.load(file_path)
            
            # ============== STAGE 1:if not spikesorting ==============

            if args.skip_spikesorting:
                logger.info("Skipping spike sorting as per user request.")

                #on all the channels detect peaks and do a raster plot
                channel_list = recording_chunk.get_channel_ids()
                traces = recording_chunk.get_traces(start_frame=0, end_frame=recording_chunk.get_num_frames(), segment_index=0,return_scaled=False)
                spike_times = {}
                for ch in channel_list:
                    x,y = detect_peaks(traces[:,np.where(channel_list == ch)[0][0]], peak_sign="neg", abs_threshold=5.5)
                    spike_times[ch] = x/fs  #convert to seconds

                

                # Save checkpoint after skipping sorting
                checkpoint.save_checkpoint(
                    ProcessingStage.SORTING_COMPLETE,
                    analyzer_folder=analyzer_folder,
                )
         
            analyzer_folder = f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}/analyzer_output"
            output_dir = f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}"
            np.save(f"{output_dir}/spikesorted_spike_times_dict.npy", spike_times)

            
            
    try:
       
        checkpoint.save_checkpoint(
            ProcessingStage.ANALYZER_COMPLETE,
            analyzer_folder=analyzer_folder,
            extensions_computed=[None]
        )



    except Exception as e:
        error_msg = f"ERROR in {project_name},{chip_id} {date} {run_id} {stream_id}/{rec_name}: {str(e)}"
        logger.error(error_msg)
        logger.error(f"TRACEBACK: {traceback.format_exc()}")

        # Mark as failed in checkpoint
        checkpoint.mark_failed(error_msg,"Sorter Analyzing Stage")

        return "error",error_msg

    try:
        # ============== STAGE 3: Generate Reports ==============
        if current_stage.value < ProcessingStage.REPORTS_COMPLETE.value:
            logger.info(f"[STAGE 3] Generating reports and metrics...")
            
            # Load spike times
            spike_times = np.load(f"{output_dir}/spikesorted_spike_times_dict.npy",allow_pickle=True).item()

            ##BURST ANALYSIS CODE
       
            fig, axs = plt.subplots(2, 1, figsize=(8, 8),sharex=True)
            # Define the ISI threshold for burst detection (e.g., 0.1 seconds)
            isi_threshold = 0.1
            # Detect bursts for each unit
            burst_statistics = helper.detect_bursts_statistics(spike_times, isi_threshold)
            bursts = [unit_stats['bursts'] for unit_stats in burst_statistics.values()]
            # Extracting ISIs as combined arrays
            all_isis_within_bursts = np.concatenate([stats['isis_within_bursts'] for stats in burst_statistics.values() if stats['isis_within_bursts'].size > 0])
            all_isis_outside_bursts = np.concatenate([stats['isis_outside_bursts'] for stats in burst_statistics.values() if stats['isis_outside_bursts'].size > 0])
            all_isis = np.concatenate([stats['isis_all'] for stats in burst_statistics.values() if stats['isis_all'].size > 0])

            # Calculate combined statistics
            mean_isi_within_combined = np.mean(all_isis_within_bursts) if all_isis_within_bursts.size > 0 else np.nan
            cov_isi_within_combined = np.cov(all_isis_within_bursts) if all_isis_within_bursts.size > 0 else np.nan

            mean_isi_outside_combined = np.mean(all_isis_outside_bursts) if all_isis_outside_bursts.size > 0 else np.nan
            cov_isi_outside_combined = np.cov(all_isis_outside_bursts) if all_isis_outside_bursts.size > 0 else np.nan

            mean_isi_all_combined = np.mean(all_isis) if all_isis.size > 0 else np.nan
            cov_isi_all_combined = np.cov(all_isis) if all_isis.size > 0 else np.nan

            # Calculate spike counts for each unit
            spike_counts = {unit: len(times) for unit, times in spike_times.items()}

            # Sort units by ascending spike counts
            sorted_units = sorted(spike_counts, key=spike_counts.get)

            axs[0]= helper.plot_raster_with_bursts(axs[0],spike_times, bursts,sorted_units=sorted_units, title_suffix="(Sorted Raster Order)")
            network_data = compute_network_bursts(ax_raster=None,ax_macro=axs[1],SpikeTimes=spike_times,plot=True)

            network_data['NumUnits'] = len(spike_times)
            network_data["fileName"]=f"{desired_pattern}/{stream_id}"

            plt.tight_layout()
            plt.xlim(0, 60)
            plt.savefig(f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}/spike_sorted_raster_plot.svg", format="svg")
        

            # Save the network data to a JSON file
            helper.save_json(f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}/network_data.json", network_data)
            

            # Save final checkpoint
            checkpoint.save_checkpoint(
                ProcessingStage.REPORTS_COMPLETE,
                num_units_filtered=len(spike_times)
            )

            logger.info(f"[STAGE 3] Reports complete - {len(spike_times)} units detected")

            electrodes = None


        # Cleanup folders
        if clear_temp_files and checkpoint.is_completed():
            logger.info("Clearing stored files...")
            binary_folder = f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}/binary"
            if os.path.exists(binary_folder):
                shutil.rmtree(binary_folder)
                logger.info(f"Removed binary folder: {binary_folder}")
            sorting_folder = f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}/analyzer_output"
            if os.path.exists(sorting_folder):
                shutil.rmtree(sorting_folder)
                logger.info(f"Removed sorting folder: {sorting_folder}")
            
        # Final resource cleanup
        
        import gc
        import shutil

        # === Resource cleanup ===
        del recording_chunk
        del sorting_obj
        del sorting_analyzer
        gc.collect()
        # Optional: free GPU memory
        try:
            import torch
            torch.cuda.empty_cache()
        except:
            pass
        try:
            import cupy
            cupy.get_default_memory_pool().free_all_blocks()
        except:
            pass
        return None, checkpoint.state.get('num_units_filtered', 0)
    
    except Exception as e:
        error_msg = f"ERROR in {stream_id}/{rec_name}: {str(e)}"
        logger.error(error_msg)
        logger.error(f"TRACEBACK: {traceback.format_exc()}")
        
        # Mark as failed in checkpoint
        checkpoint.mark_failed(error_msg,"Report Generation Stage")
        
        # if clear_temp_files:
        #     helper.empty_directory(sorting_folder)
        
        return "The sorting failed.", "The sorting failed."
    #helper.dumpdicttofile(failed_sorting_rec,'./failed_recording_sorting')



def main():
    global args,logger

    parser = argparse.ArgumentParser(description="Process single MEA file and well")

    parser.add_argument('file_dir_path', type=str, help='Path to the .h5 file')
    parser.add_argument('--well', type=str, required=True, help='Well IDs to process')
    parser.add_argument('--params', type=str, help='JSON string or path with curation thresholds')
    parser.add_argument('--docker', type=str, help='Docker image (if using containerized sorter)')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    parser.add_argument('--sorter', type=str, default='kilosort4', help='Sorter to use')
    parser.add_argument('--resume', action='store_true', help='Resume from checkpoint')
    parser.add_argument('--force-restart', action='store_true', help='Restart even if completed')
    parser.add_argument('--clean-up', action='store_true', help='Clear temporary files')
    parser.add_argument('--export-to-phy', action='store_true', help='Export results to Phy format')
    parser.add_argument('--checkpoint-dir', type=str, default=f'{BASE_FILE_PATH}/../AnalyzedData/checkpoints', help='Checkpoint folder')

    args = parser.parse_args()
    path = args.file_dir_path
    well = args.well
    #parse file path
    run_id,project_name,date,chip_id,desired_pattern = parse_h5_path(args.file_dir_path)
    logger = setup_logger(log_file=f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/logs/{run_id}_{chip_id}_{date}_{well}.log")
    if args.debug:
        logger.setLevel(logging.DEBUG)


    if not os.path.exists(path):
        logger.error(f"Path does not exist: {path}")
        sys.exit(1)

    # Load thresholds
    thresholds = None
    if args.params:
        try:
            if os.path.isfile(args.params):
                with open(args.params, 'r') as f:
                    thresholds = json.load(f)
                logger.info(f"Loaded parameters from file: {args.params}")
            else:
                thresholds = json.loads(args.params)
                logger.info("Loaded parameters from command line string")
        except Exception as e:
            logger.error(f"Failed to parse params: {e}")
            sys.exit(1)
    try:
        if args.clean_up:
            logger.info("Starting processing with cleanup enabled")
            clear_temp_files = True
        else:
            logger.info("Starting processing without cleanup")
            clear_temp_files = False
        process_block(path, stream_id=well, thresholds=thresholds, clear_temp_files=clear_temp_files)
    except Exception as e:
        logger.error(f"Processing failed for {well}: {e}")
        logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
    


