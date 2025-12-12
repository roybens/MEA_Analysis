import os
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
os.environ["QT_QPA_PLATFORM"] = "offscreen"  # For headless environments
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend
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
import numpy as np
import json
import h5py
import psutil
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






#breakpoint()    
def save_zarr_wrapper(filepath,op_folder):
    """
    saves and compresses the file into the zarr format. 
    IP: file and the corresponding output folder.
    """
    import spikeinterface.full as si
    from flac_numcodecs import Flac
    compressor = Flac(level =8)
    rec_original = si.read_maxwell(filepath)
    rec_int32 = si.scale(rec_original, dtype="int32")
    # remove the 2^15 offset
    rec_rm_offset = si.scale(rec_int32, offset=-2 ** 15)
    # now we can safely cast to int16
    rec_int16 = si.scale(rec_rm_offset, dtype="int16")
    recording_zarr = rec_int16.save(format = "zarr",folder=op_folder,compressor=compressor,
                                    channel_chunk_size =2,n_jobs=64,chunk_duration="1s")
    return recording_zarr





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


def automatic_curation(metrics, thresholds=None):
    """
    Apply thresholds only when values exist (NaN -> skip filtering).
    Also produce a full rejection log showing which rule failed.
    """

    # Default thresholds
    default_thresholds = {
        'presence_ratio': 0.75,
        'rp_contamination': 0.15,
        'firing_rate': 0.05,
        'amplitude_median': -50,          # amplitude must be <= -50 µV
        'amplitude_cv_median': 0.3,       # stable waveform
    }

    if thresholds:
        default_thresholds.update(thresholds)

    # Rejection log list
    rejection_records = []

    # Boolean mask for accepted units
    keep_mask = np.ones(len(metrics), dtype=bool)

    # Iterate rows
    for idx, row in metrics.iterrows():
        reasons = []

        for key, thr in default_thresholds.items():

            if key not in row:
                continue

            val = row[key]

            # Skip filtering if NaN
            if pd.isna(val):
                continue

            #### APPLY RULES ####
            if key == "presence_ratio" and not (val > thr):
                reasons.append(f"{key} <= {thr} ({val:.3f})")

            elif key == "rp_contamination" and not (val < thr):
                reasons.append(f"{key} >= {thr} ({val:.3f})")

            elif key == "firing_rate" and not (val > thr):
                reasons.append(f"{key} <= {thr} ({val:.3f})")

            elif key == "amplitude_median" and not (val <= thr):
                reasons.append(f"{key} > {thr} ({val:.1f} µV)")

            elif key == "amplitude_cv_median" and not (val < thr):
                reasons.append(f"{key} >= {thr} ({val:.3f})")

        # If any reasons collected → mark as rejected
        if len(reasons) > 0:
            keep_mask[idx] = False
            rejection_records.append({
                "unit_id": row.get("unit_id", idx),
                "reasons": "; ".join(reasons)
            })

    # Final outputs
    kept_units = metrics[keep_mask].copy()
    rejection_log = pd.DataFrame(rejection_records)

    return kept_units, rejection_log



""" def automatic_curation(metrics, thresholds=None):

    # Default thresholds
    default_thresholds = {
        'presence_ratio': 0.9,
        'rp_contamination': 0.2,   #used instead of isi_violations_ratio
        'firing_rate': 0.05, # in Hz about 3 spikes in 1 min at least
        'amplitude_median':-20,
        'amplitude_cv_median':0.5
    }

    # Update default thresholds with provided values
    if thresholds:
        default_thresholds.update(thresholds)

    # Function to get query based on key and value
    def get_query_condition(key, value):
        conditions = {
            'num_spikes': f"({key} > {value})",
            'presence_ratio': f"({key} > {value})",
            'rp_contamination': f"({key} < {value})",
            'firing_rate': f"({key} > {value})",
            'amplitude_median': f"({key} <= {value})",
            'amplitude_cv_median': f"({key} < {value})"

        }
        # Return the appropriate condition or a default condition if key is not found
        return conditions.get(key, f"({key} > {value})")

    # Build the query string dynamically using the get_query_condition function
    query_conditions = [get_query_condition(key, value) for key, value in default_thresholds.items()]

    # Join all conditions with ' & '
    our_query = ' & '.join(query_conditions)

    metrics = metrics.query(our_query)
    return metrics
"""


def find_connected_components(pairs):
    # Create a graph using a dictionary
    graph = defaultdict(list)
    
    # Populate the graph with edges
    for x, y in pairs:
        graph[x].append(y)
        graph[y].append(x)  # Assuming undirected graph for transitive relationships

    # Function to perform DFS and find connected components
    def dfs(node, visited, component):
        visited.add(node)
        component.append(node)
        for neighbour in graph[node]:
            if neighbour not in visited:
                dfs(neighbour, visited, component)

    # List to store all components
    components = []
    visited = set()

    # Iterate through all nodes in the graph
    for node in graph:
        if node not in visited:
            component = []
            dfs(node, visited, component)
            components.append(component)

    return components

def merge_similar_templates(sorting,waveforms,sim_score =0.7):

    matrix = si.compute_template_similarity(waveforms,load_if_exists=True)
    metrics = si.compute_quality_metrics(waveforms,load_if_exists=True)
    temp_metrics = si.compute_template_metrics(waveforms,load_if_exists=True)
    n = matrix.shape[0]

    # Find indices where values are greater than 0.5 and not on the diagonal
    indices = [(i,j) for i in range(1,n) for j in range(i-1) if matrix[i,j]>sim_score]

    pairs = []
    # Print the indices
    for ind in indices:
        logger.debug(f"{temp_metrics.index[ind[0]],temp_metrics.index[ind[1]]}")
        pairs.append((temp_metrics.index[ind[0]],temp_metrics.index[ind[1]]))
    connected_components = find_connected_components(pairs)
    cs = CurationSorting(sorting)
    for elements in connected_components:
        cs.merge(elements)
    return cs.sorting

  


def remove_similar_templates(waveforms,sim_score =0.7):

    matrix = si.compute_template_similarity(waveforms,load_if_exists=True)
    metrics = si.compute_quality_metrics(waveforms,load_if_exists=True)
    temp_metrics = si.compute_template_metrics(waveforms,load_if_exists=True)
    n = matrix.shape[0]

    # Find indices where values are greater than 0.5 and not on the diagonal
    indices = [(i,j) for i in range(1,n) for j in range(i-1) if matrix[i,j]>sim_score]

    removables =[]
    # Print the indices
    for ind in indices:
        print(temp_metrics.index[ind[0]],temp_metrics.index[ind[1]])
        if metrics['amplitude_median'].loc[temp_metrics.index[ind[0]]] < metrics['amplitude_median'].loc[temp_metrics.index[ind[1]]]:
            smaller_index = temp_metrics.index[ind[0]]
        else:
            smaller_index = temp_metrics.index[ind[1]]

        removables.append(smaller_index)
    return removables

def analyse_waveforms_sigui(waveforms_folder) :
    import spikeinterface_gui
    waveforms = si.load_waveforms(waveforms_folder)
    # This creates a Qt app
   # waveforms.run_extract_waveforms(**job_kwargs)
    app = spikeinterface_gui.mkQApp() 

    # create the mainwindow and show
    win = spikeinterface_gui.MainWindow(waveforms)
    win.show()
    # run the main Qt6 loop
    app.exec_()
    return


""" 

def get_unique_templates_channels(good_units, waveform):

    
    #get_extremum_channels.
    unit_extremum_channel =spikeinterface.full.get_template_extremum_channel(waveform, peak_sign='neg')
    #Step 1: keep only units that are in good_units 
    unit_extremum_channel = {key:value for key,value in unit_extremum_channel.items() if key in good_units}
    print(f"extremum channel : {unit_extremum_channel}")
    #Step3: get units that correspond to same electrodes.
    output_units = [[key for key, value in unit_extremum_channel.items() if value == v] for v in set(unit_extremum_channel.values()) if list(unit_extremum_channel.values()).count(v) > 1]
    print(f"Units that correspond to same electrode: {output_units}")
    #Step 3: get the metrics
    #missing a step here? - Rohan 

    #Step4: select best units with same electrodes ( based on amp values)
    output=[]
    if output_units :
        for sublist in output_units :
            amp_max = 0 
            for unit in sublist:
                if metrics['amplitude_median'][int(unit)] > amp_max :
                    amp_max = metrics['amplitude_median'][int(unit)]
                    reqd_unit = unit
            output.append(reqd_unit)
    print(f"Best unit among the same electrodes {output}")
    #Step 5 --> unit_extremum_channel - output_units + output
    output_units = [element for sublist in output_units for element in sublist]
    new_list = [ item for item in output_units if item not in output]

    required_templates = {key:value for key,value in unit_extremum_channel.items() if key not in new_list}

    return required_templates
     
 """



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
                  clear_temp_files=False, thresholds=None,checkpoint_root=None,output_root=None):
    

    global args,logger,BASE_FILE_PATH,job_kwargs

    # Setup paths and checkpoint manager
    if file_path.endswith('.raw.h5'):
        recording, rec_name = get_data_maxwell(file_path, stream_id, recnumber)
        logger.info(f"Processing recording: {file_path}")
        #file pattern extraction
        run_id,project_name,date,chip_id,desired_pattern = parse_h5_path(file_path)
        desired_output = os.path.join(output_root,desired_pattern,stream_id)
        os.makedirs(
            desired_output,
            mode=0o777,
            exist_ok=True
        )

    else:
        rec_name = "binary_recording"   #implement fr binary
        logger.error("Binary file processing not implemented yet")
        sys.exit(1)


    # Initialize checkpoint manager
    checkpoint = PipelineCheckpoint(
        base_path=checkpoint_root,
        project_name=project_name,
        date=date,
        chip_id=chip_id,
        run_id=run_id,
        well_id=stream_id,
        rec_name=rec_name
    )
    
    # Check if should skip
    if checkpoint.should_skip(args):
        #resource cleanup
                
        import gc
        import shutil
        # Cleanup folders
        if clear_temp_files:
            logger.info("Clearing stored files if exists.")
            binary_folder =os.path.join(desired_output,"binary")
            if os.path.exists(binary_folder):
                shutil.rmtree(binary_folder)
                logger.info(f"Removed binary folder: {binary_folder}")
            else:
                logger.info(f"Binary folder does not exist: {binary_folder}")
            sorting_folder = os.path.join(desired_output, "analyzer_output")
            if os.path.exists(sorting_folder):
                shutil.rmtree(sorting_folder)
                logger.info(f"Removed sorting folder: {sorting_folder}")    
            else:
                logger.info(f"Sorting folder does not exist: {sorting_folder}")
            
        # Final resource cleanup


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
                desired_output = os.path.join(output_root,desired_pattern,stream_id)
                   
                os.makedirs(
                    desired_output,
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
                    folder=os.path.join(desired_output,"binary"),
                    format='binary',
                    overwrite=True,
                    **job_kwargs
                )
                recording_chunk = si.load(
                    os.path.join(desired_output,"binary")
                )
            else:
                recording_chunk = si.load(file_path)
            

         
            analyzer_folder = os.path.join(desired_output,"analyzer_output")
            
            # ============== STAGE 1: Spike Sorting ==============
            logger.info(f"[STAGE 1] Running spike sorting...")
            if args.debug:
                log_resource_usage()
            start = timer()
            
            if args.docker:
                if args.sorter in ['kilosort2','kilosort2_5', 'kilosort3']:
                    sorting_obj = si.run_sorter(
                        sorter_name=args.sorter,
                        recording=recording_chunk,
                        folder=analyzer_folder,
                        remove_existing_folder=True,
                        delete_output_folder=False,
                        verbose=True,
                        docker_image=args.docker,
                        with_output=True
                    )
                elif args.sorter == 'kilosort4':
                    ks_params ={
                        'batch_size': int(fs) * 2,
                        'clear_cache': True,
                        'invert_sign': True,
                        'cluster_downsampling': 20,
                        'max_cluster_subset': None,
                        'nblocks':0,
                        'dmin':17,
                        'do_correction': False,
                    }
                    sorting_obj = si.run_sorter(
                        sorter_name="kilosort4",
                        recording=recording_chunk,
                        folder=analyzer_folder,
                        delete_output_folder=False,
                        verbose=True,
                        with_output=True,
                        **ks_params
                    )
                else:

                    logger.error(f"Unsupported Mode with sorter with docker: {args.sorter}")
                    raise ValueError(f"Unsupported Mode with sorter with docker: {args.sorter}")
            else:
                if args.sorter == 'kilosort2' or args.sorter == 'kilosort3':
                    sorting_obj = si.run_sorter_local(
                        sorter_name=args.sorter,
                        recording=recording_chunk,
                        folder=analyzer_folder,
                        remove_existing_folder=True,
                        delete_output_folder=False,
                        verbose=True,
                        with_output=True
                    )
                elif args.sorter == 'kilosort4':
                    
                    ks_params ={
                        'batch_size': int(fs) *1,
                        'clear_cache': True,
                        'invert_sign': True,
                        'cluster_downsampling': 20,
                        'max_cluster_subset': None,
                        'nblocks':0,
                        'dmin':35,
                        'dminx':35,
                        'do_correction': False,
                        'nearest_templates': 20,
                        'max_channel_distance': 30,
                        'template_sizes':3,
                        'n_templates': 4,
                        'max_peels':50



                    }
                    
                    sorting_obj = si.run_sorter_local(
                        sorter_name=args.sorter,
                        recording=recording_chunk,
                        folder=analyzer_folder,
                        delete_output_folder=False,
                        verbose=True,
                        with_output=True,
                        **ks_params
                    )
                else:
                    logger.error(f"Unsupported Mode with sorter without docker: {args.sorter}")
                    raise ValueError(f"Unsupported Mode with sorter without docker: {args.sorter}")

            logger.debug(f"Sorting complete in: {timer()-start:.2f}s")

            # Clean up sorting
            sorting_obj = sorting_obj.remove_empty_units()
            sorting_obj = si.remove_excess_spikes(
                sorting_obj,
                recording_chunk
            )
            sorting_obj = si.remove_duplicated_spikes( sorting_obj,censored_period_ms=0.1 )
            #remove redundant units based on template similarity 
            #sorting_obj = si.remove_redundant_units(sorting_obj,duplicate_threshold=0.9,remove_strategy="minimum_shift") #TODO
            #sorting_obj.save(folder=kilosort_output_folder,overwrite=True)  # Overwrite with cleaned sorting
            # Save checkpoint after sorting
            checkpoint.save_checkpoint(
                ProcessingStage.SORTING_COMPLETE,
                analyzer_folder=analyzer_folder,


            )
        
            if args.debug:
                log_resource_usage()
            
        else:
            # Resume from checkpoint - load existing sorting
            logger.info(f"[STAGE 1] Resuming - loading existing sorting...")
            analyzer_folder = checkpoint.state['analyzer_folder']

            logger.info(f"desired_pattern: {desired_pattern}")
            desired_output = os.path.join(output_root,desired_pattern,stream_id)
            # Reload recording
            recording_chunk = si.load(
                os.path.join(desired_output,"binary")
            )
            
            # Try multiple loading strategies
            sorting_obj = None
            
            # Attempt 1: Try read_sorter_folder (if SpikeInterface metadata exists)
            try:
                logger.info("Attempt 1: Using read_sorter_folder")
                sorting_obj = si.read_sorter_folder(folder=analyzer_folder)
                logger.info(f"Successfully loaded with read_sorter_folder: {len(sorting_obj.get_unit_ids())} units")
            except Exception as e:
                logger.debug(f"read_sorter_folder failed: {e}")

            # Attempt 2: Try NumpySorting (if analyzer was created)
            if sorting_obj is None:
                try:
                    logger.info("Attempt 2: Using read_numpy_sorting")
                    sorting_path = f"{analyzer_folder}/sorting"
                    if os.path.exists(sorting_path):
                        sorting_obj = si.read_numpy_sorting_folder(sorting_path)
                        logger.info(f"Successfully loaded with read_numpy_sorting: {len(sorting_obj.get_unit_ids())} units")
                except Exception as e:
                    logger.debug(f"read_numpy_sorting failed: {e}")

            # Attempt 3: Try reading from Kilosort native format
            if sorting_obj is None:
                try:
                    logger.info("Attempt 3: Using read_kilosort from sorter_output")
                    sorter_output = f"{analyzer_folder}/sorter_output"
                    if os.path.exists(sorter_output):
                        sorting_obj = si.read_kilosort(folder=sorter_output)
                        logger.info(f"Successfully loaded with read_kilosort: {len(sorting_obj.get_unit_ids())} units")
                except Exception as e:
                    logger.debug(f"read_kilosort from sorter_output failed: {e}")

            # Attempt 4: Try reading Kilosort from analyzer_folder directly
            if sorting_obj is None:
                try:
                    logger.info("Attempt 4: Using read_kilosort from analyzer_folder")
                    sorting_obj = si.read_kilosort(folder=analyzer_folder)
                    logger.info(f"Successfully loaded with read_kilosort: {len(sorting_obj.get_unit_ids())} units")
                except Exception as e:
                    logger.debug(f"read_kilosort from analyzer_folder failed: {e}")

            # If all attempts failed
            if sorting_obj is None:
                logger.error(f"Failed to load sorting from: {analyzer_folder}")
                logger.error("Tried: read_sorter_folder, read_numpy_sorting_folder, read_kilosort (multiple paths)")
                raise FileNotFoundError(f"Could not load sorting from any format in {analyzer_folder}")
            
            # Re-apply cleaning
            sorting_obj = sorting_obj.remove_empty_units()
            sorting_obj = si.remove_excess_spikes(sorting_obj, recording_chunk)
            sorting_obj = si.remove_duplicated_spikes( sorting_obj,censored_period_ms=0.1 )

    except Exception as e:
        error_msg = f"ERROR in {project_name},{chip_id} {date} {run_id} {stream_id}/{rec_name}: {str(e)}"
        logger.error(error_msg)
        logger.error(f"TRACEBACK: {traceback.format_exc()}")

        # Mark as failed in checkpoint
        checkpoint.mark_failed(error_msg,"Sorting Stage")
        

        
        return "error",error_msg

    try:
        # ============== STAGE 2: Sorting Analyzer ==============
        if current_stage.value < ProcessingStage.ANALYZER_COMPLETE.value:
            logger.info(f"[STAGE 2] Computing sorting analyzer...")
            start = timer()
            
            analyzer_folder = os.path.join(desired_output,"analyzer_output")
            
            # Check if analyzer already exists (partial completion scenario)
            if Path(analyzer_folder).exists():
                logger.info("Found existing analyzer folder, attempting to load...")
                try:
                    sorting_analyzer = si.load_sorting_analyzer(
                        folder=analyzer_folder,
                        recording=recording_chunk
                    )
                    logger.info("Successfully loaded existing analyzer")
                except Exception as e:
                    logger.warning(f"Could not load existing analyzer: {e}. Creating new one...")
                    # Remove corrupted folder and start fresh
                    import shutil
                    shutil.rmtree(analyzer_folder)
                    sorting_analyzer = None
            else:
                sorting_analyzer = None
            
            # Create new analyzer if needed
            if sorting_analyzer is None:
                logger.info("Creating new sorting analyzer...")
                
                # Optimal sparsity for HD-MEA (Maxwell Biosystems)
                sparsity = si.estimate_sparsity(
                    sorting=sorting_obj,
                    recording=recording_chunk,
                    method='radius',      # Best for spatial arrays
                    radius_um=50,         # 40-70 um typical for HD-MEA (17.5 um pitch)
                    peak_sign='neg'   # Most HD-MEA use negative peaks
                )
                
                #logging.info(f"Estimated sparsity: avg {sparsity.num_active_channels_per_unit().mean():.1f} channels/unit")
                
                sorting_analyzer = si.create_sorting_analyzer(
                    sorting_obj,
                    recording_chunk,
                    format="binary_folder",
                    sparsity=sparsity,  # Use pre-computed ChannelSparsity object
                    return_in_uV=True,
                    folder=analyzer_folder
                )

            
            # Check which extensions are already computed
            existing_extensions = sorting_analyzer.get_loaded_extension_names()
            logger.info(f"Existing extensions: {existing_extensions}")
            
            # Get all available extensions
            all_available = sorting_analyzer.get_computable_extensions()
            logger.info(f"All available extensions: {all_available}")

            # Determine which extensions need to be computed
            missing_extensions = [ext for ext in all_available if ext not in existing_extensions]
            
            if missing_extensions:
                logger.info(f"Computing {len(missing_extensions)} missing extensions: {missing_extensions}")
                
                if args.debug:
                    log_resource_usage()
                
                # Compute missing extensions with error handling
                successfully_computed = []
                failed_extensions = []
                
                for ext_name in missing_extensions:
                    try:
                        logger.info(f"  Computing {ext_name}...")
                        ext_start = timer()
                        
                        # Compute individual extension
                        sorting_analyzer.compute(ext_name, save=True, **job_kwargs)
                        
                        duration = timer() - ext_start
                        successfully_computed.append(ext_name)
                        logger.info(f"  ✓ {ext_name} completed in {duration:.2f}s")
                        
                    except Exception as e:
                        logger.error(f"  ✗ Failed to compute {ext_name}: {e}")
                        failed_extensions.append(ext_name)
                        
                        # Only fail if it's a critical extension
                        critical_extensions = ['waveforms', 'templates', 'noise_levels']
                        if ext_name in critical_extensions:
                            raise RuntimeError(f"Critical extension {ext_name} failed: {e}")
                
                if args.debug:
                    log_resource_usage()

                logger.info(f"Successfully computed: {successfully_computed}")
                if failed_extensions:
                    logger.warning(f"Failed extensions (non-critical): {failed_extensions}")
            else:
                logger.info("All extensions already computed!")

            # Verify critical extensions exist
            critical_extensions = ['waveforms', 'templates', 'quality_metrics', 'template_metrics']
            missing_critical = [ext for ext in critical_extensions 
                            if not sorting_analyzer.has_extension(ext)]
            
            if missing_critical:
                raise RuntimeError(f"Critical extensions missing: {missing_critical}")
            
            total_duration = timer() - start
            logger.debug(f"Analyzer stage complete in: {total_duration:.2f}s")


            #sorting_analyzer = si.remove_redundant_units(sorting_analyzer,duplicate_threshold=0.9,remove_strategy="minimum_shift")
            
            # Save checkpoint after analyzer
            checkpoint.save_checkpoint(
                ProcessingStage.ANALYZER_COMPLETE,
                analyzer_folder=analyzer_folder,
                extensions_computed=sorting_analyzer.get_loaded_extension_names()
            )

        else:
            # Resume from checkpoint - load existing analyzer
            logger.info(f"[STAGE 2] Resuming - loading existing analyzer...")
            analyzer_folder = checkpoint.state['analyzer_folder']
            
        
            
            # Load analyzer
            sorting_analyzer = si.load(
                analyzer_folder
            )

            logger.info(f"Loaded analyzer with extensions: {sorting_analyzer.get_loaded_extension_names()}")

            #sorting_analyzer = si.remove_redundant_units(sorting_analyzer,duplicate_threshold=0.9,remove_strategy="minimum_shift")   #this and clean sorting should go inside a curation funciton

            
            # Check if any extensions are missing (in case of partial completion)
            existing_extensions = sorting_analyzer.get_loaded_extension_names()
            all_available = sorting_analyzer.get_computable_extensions()
            missing_extensions = [ext for ext in all_available if ext not in existing_extensions]
            
            if missing_extensions:
                logger.info(f"Found {len(missing_extensions)} missing extensions, computing them...")
                
                for ext_name in missing_extensions:
                    try:
                        logger.info(f"  Computing missing extension: {ext_name}...")
                        sorting_analyzer.compute(ext_name, save=True, **job_kwargs)
                        logger.info(f"  ✓ {ext_name} completed")
                    except Exception as e:
                        logger.warning(f"  ✗ Failed to compute {ext_name}: {e}")

                        # Only fail if critical
                        critical_extensions = ['waveforms', 'templates', 'quality_metrics']
                        if ext_name in critical_extensions:
                            raise RuntimeError(f"Critical extension {ext_name} failed: {e}")
                
                # Update checkpoint with newly computed extensions
                checkpoint.save_checkpoint(
                    ProcessingStage.ANALYZER_COMPLETE,
                    analyzer_folder=analyzer_folder,
                    extensions_computed=sorting_analyzer.get_loaded_extension_names()
                )
            else:
                logger.info("All extensions already present")


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
            
            # Extract metrics
            qual_metrics = sorting_analyzer.get_extension('quality_metrics').get_data()
            template_metrics = sorting_analyzer.get_extension('template_metrics').get_data()
            
            # Get output pattern
            pattern = r"/(\d+)/data.raw.h5"
            run_id = int(re.search(pattern, file_path).group(1))
            file_pattern = os.path.dirname(file_path)
            parts = file_pattern.split('/')
            desired_pattern = '/'.join(parts[-6:])
            
            desired_output = os.path.join(output_root, desired_pattern, stream_id)
            
            # Save unfiltered metrics
            qual_metrics.to_excel(f"{desired_output}/quality_metrics_unfiltered.xlsx")
            template_metrics.to_excel(f"{desired_output}/template_metrics_unfiltered.xlsx")
            
            # Filter units
            #this should be with and without automatic curation TODO" may be a flag??
            if args.no_curation:
                logger.info("Skipping automatic curation as per user request (--no_curation)")
                non_violated_units = qual_metrics.index.values
                numunits = len(non_violated_units)
                update_qual_metrics = qual_metrics
                #rejected_record = pd.DataFrame(columns=["unit_id", "reasons"])
            else:
                logger.info("Applying automatic curation based on quality metrics")
                update_qual_metrics, rejected_record = automatic_curation(qual_metrics, thresholds)
                rejected_record.to_excel(f"{desired_output}/automatic_curation_rejection_log.xlsx", index=False)
                non_violated_units = update_qual_metrics.index.values
                numunits = len(non_violated_units)
            
                if numunits == 0:
                    logger.warning(f"No units passed quality criteria for {stream_id}/{rec_name}")
                    checkpoint.save_checkpoint(
                        ProcessingStage.REPORTS_COMPLETE,
                        num_units_filtered=0,
                        warning="No units passed quality criteria"
                    )
                    return None, 0
            
            curated_sorting = sorting_obj.select_units(non_violated_units)
            curated_analyzer = sorting_analyzer.select_units(non_violated_units)
            #curated_analyzer.save(folder=f"{output_dir}/analyzer_curated", overwrite=True)

            # Save filtered metrics
            template_metrics_filtered = template_metrics.loc[non_violated_units]
            template_metrics_filtered.to_excel(f"{desired_output}/template_metrics.xlsx")
            
            # Get locations
            if sorting_analyzer.has_extension('unit_locations'):
                locations = sorting_analyzer.get_extension('unit_locations').get_data()
            else:
                sorting_analyzer.compute('unit_locations', method='monopolar_triangulation')
                locations = sorting_analyzer.get_extension('unit_locations').get_data()
            
            qual_metrics_filtered = qual_metrics.loc[non_violated_units]
            qual_metrics_filtered['location_X'] = locations[non_violated_units, 0]
            qual_metrics_filtered['location_Y'] = locations[non_violated_units, 1]
            qual_metrics_filtered['location_Z'] = locations[non_violated_units, 2]
            qual_metrics_filtered.to_excel(f"{desired_output}/quality_metrics.xlsx")
            
            # Generate plots
            unit_ids = sorting_analyzer.unit_ids
            unit_locations = dict(zip(unit_ids, locations))
            
            # Plot unfiltered units
            fig1, ax1 = plt.subplots(figsize=(10.5, 6.5))
            si.plot_probe_map(recording_chunk, ax=ax1, with_channel_ids=False)
            for unit_id, (x, y, z) in unit_locations.items():
                ax1.scatter(x, y, s=10, c='blue', alpha=0.6)
            ax1.invert_yaxis()
            fig1.savefig(f"{desired_output}/locations_unfiltered_units.pdf")
            plt.close(fig1)
            
            # Plot filtered units
            fig2, ax2 = plt.subplots(figsize=(10.5, 6.5))
            si.plot_probe_map(recording_chunk, ax=ax2, with_channel_ids=False)
            for unit_id, (x, y, z) in unit_locations.items():
                if unit_id in non_violated_units:
                    ax2.scatter(x, y, s=10, c='blue', alpha=0.6)
            ax2.invert_yaxis()
            fig2.savefig(f"{desired_output}/locations_{numunits}_units.pdf")
            plt.close(fig2)
            
            # Generate waveform plots
            unit_extremum_channel = si.get_template_extremum_channel(
                sorting_analyzer,
                peak_sign='neg'
            )
            unit_extremum_channel = {
                key: value for key, value in unit_extremum_channel.items()
                if key in non_violated_units
            }
            
            os.makedirs(f"{output_dir}/waveforms/", mode=0o777, exist_ok=True)
            import matplotlib.backends.backend_pdf as pdf

            pdf_file = os.path.join(desired_output, "waveforms", "waveforms_subplots.pdf")

            with pdf.PdfPages(pdf_file) as pdf_pages:
                num_cols = 4
                num_rows = 3
                
                # Get waveforms extension and sparsity information
                waveform_ext = sorting_analyzer.get_extension('waveforms')
                sparsity = waveform_ext.sparsity
                
                # Get templates extension for fallback channel selection
                templates_ext = sorting_analyzer.get_extension('templates')
                
                for i in range(0, len(non_violated_units), num_rows * num_cols):
                    units_to_plot = non_violated_units[i:i + num_rows * num_cols]
                    
                    fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 8), squeeze=False)
                    plt.subplots_adjust(hspace=1, wspace=0.5)
                    
                    for j, unit_id in enumerate(units_to_plot):
                        row = j // num_cols
                        col = j % num_cols
                        ax = axes[row, col]

                        try:
                            # Get waveforms for this unit
                            wf = waveform_ext.get_waveforms_one_unit(unit_id)
                            
                            # Get extremum channel ID (global)
                            extremum_channel_id = unit_extremum_channel[unit_id]
                            channel_id_str = str(int(extremum_channel_id))
                            
                            # Handle sparsity correctly
                            if sparsity is not None:
                                # Get the sparse channel IDs for this unit
                                unit_channel_ids = sparsity.unit_id_to_channel_ids[unit_id]
                                unit_channel_ids_str = [str(int(ch)) for ch in unit_channel_ids]
                                
                                if channel_id_str in unit_channel_ids_str:
                                    # Extremum channel is in sparse set
                                    sparse_idx = unit_channel_ids_str.index(channel_id_str)
                                    channel_label = f"extremum ch {channel_id_str}"
                                else:
                                    # Extremum channel not in sparse set - find channel with maximum amplitude
                                    template = templates_ext.get_unit_template(unit_id=unit_id)
                                    max_amp_sparse_idx = np.argmax(np.abs(template).max(axis=0))
                                    sparse_idx = max_amp_sparse_idx
                                    actual_channel = unit_channel_ids[sparse_idx]

                                    logger.warning(
                                        f"Unit {unit_id}: extremum channel {channel_id_str} not in sparse set. "
                                        f"Using max amplitude channel {actual_channel} instead."
                                    )
                                    channel_label = f"ch {actual_channel} (max amp)"
                            else:
                                # Dense waveforms - use global channel index
                                sparse_idx = sorting_analyzer.channel_ids_to_indices([channel_id_str])[0]
                                channel_label = f"extremum ch {channel_id_str}"
                            
                            # Plot waveforms
                            ax.plot(wf[:, :, sparse_idx].T, lw=1, color='black', alpha=0.1, 
                                linestyle='-', marker='', markersize=0)
                            ax.set_title(f"Unit {unit_id}\n{channel_label}", fontsize=9)
                            ax.set_ylabel("Amplitude (µV)")
                            ax.set_ylim(-700, 400)
                            ax.set_xlabel("Sampled timepoints (5e-2 ms)")
                            ax.set_facecolor('white')
                            ax.tick_params(axis='x', colors='black')
                            ax.tick_params(axis='y', colors='black')
                            ax.spines['top'].set_color('white')
                            ax.spines['right'].set_color('white')

                            # Save individual subplot
                            fig_single = plt.figure()
                            ax_single = fig_single.add_subplot(111)
                            ax_single.plot(wf[:, :, sparse_idx].T, lw=1, color='black', alpha=0.1, 
                                        linestyle='-', marker='', markersize=0)
                            ax_single.set_title(f"Waveforms of Unit {unit_id}\n{channel_label}")
                            ax_single.set_ylabel("Amplitude (µV)")
                            #ax_single.set_ylim(-700, 400)
                            ax_single.set_xlabel("Sampled timepoints (5e-2 ms)")
                            ax_single.set_facecolor('white')
                            ax_single.tick_params(axis='x', colors='black')
                            ax_single.tick_params(axis='y', colors='black')
                            ax_single.spines['top'].set_color('white')
                            ax_single.spines['right'].set_color('white')
                            #fig_single.savefig(
                            #   f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}/waveforms/{unit_id}.svg", 
                            #    format="svg"
                           #)
                            plt.close(fig_single)
                            
                        except Exception as e:
                            # Fallback for plotting errors
                            logger.error(f"Failed to plot waveforms for unit {unit_id}: {e}")
                            ax.text(0.5, 0.5, f"Unit {unit_id}\nPlot failed", 
                                ha='center', va='center', transform=ax.transAxes, fontsize=8)
                            ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

                    # Remove empty subplots
                    num_units_to_plot = len(units_to_plot)
                    for k in range(num_units_to_plot, num_rows * num_cols):
                        row = k // num_cols
                        col = k % num_cols
                        ax = axes[row, col]
                        ax.axis('off')

                    # Save the subplots to the PDF file
                    pdf_pages.savefig(fig)
                    plt.close(fig)

            logger.info(f"Waveform subplots saved to: {pdf_file}")

    
            # Generate spike train analysis
            fs = recording_chunk.get_sampling_frequency()
            frame_start = 0
            if time_in_s is None:
                time_in_s = (recording_chunk.get_num_frames() / fs ) -1
            frame_end = recording_chunk.get_num_frames()                                            ##for the resume option it is required to reload the recording chunk
            spike_times = {}
            
            for idx, unit_id in enumerate(non_violated_units):
                spike_train = sorting_obj.get_unit_spike_train(
                    unit_id,
                    start_frame=frame_start,
                    end_frame=frame_end
                )
                if len(spike_train) > 0:
                    spike_times[idx] = spike_train / float(fs)
            
            np.save(f"{desired_output}/spikesorted_spike_times_dict.npy", spike_times)

            ##BURST ANALYSIS CODE
       
            fig, axs = plt.subplots(2, 1, figsize=(8, 8),sharex=True)
            # Define the ISI threshold for burst detection (e.g., 0.1 seconds)
            isi_threshold = 0.1
            # Detect bursts for each unit
            burst_statistics = helper.detect_bursts_statistics(spike_times, isi_threshold)
            bursts = [unit_stats['bursts'] for unit_stats in burst_statistics.values()]
            # Calculate spike counts for each unit
            spike_counts = {unit: len(times) for unit, times in spike_times.items()}

            # Sort units by ascending spike counts
            sorted_units = sorted(spike_counts, key=spike_counts.get)

            axs[0]= helper.plot_raster_with_bursts(axs[0],spike_times, bursts,sorted_units=sorted_units, title_suffix="(Sorted Raster Order)")
            network_data = compute_network_bursts(ax_raster=None,ax_macro=axs[1],SpikeTimes=spike_times,plot=True)

            network_data['NumUnits'] = len(non_violated_units)
            network_data["fileName"]=f"{desired_pattern}/{stream_id}"

        


            plt.tight_layout()
            plt.xlim(0, 60)
            plt.savefig(f"{desired_output}/spike_sorted_raster_plot.svg", format="svg")
        
            #fig2, axs2 = plt.subplots(2, 1, figsize=(8, 8),sharex=True)
            #axs2[0] = helper.plot_raster_with_bursts(axs[0],spike_times, bursts,sorted_units=None, title_suffix="(Origininal Raster Order)")
            # Copy the second plot to the new figure
            #fig2._axstack.add(fig2._make_key(axs2[1]), axs[1])
            #plt.tight_layout()
            #plt.xlim(0, 60)
            #plt.savefig(f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/{stream_id}/original_raster_plot.eps", format="eps")
            # Save the network data to a JSON file
            helper.save_json(f"{desired_output}/network_data.json", network_data)
            
            # compiledNetworkData =f"{BASE_FILE_PATH}/../AnalyzedData/{desired_pattern}/../../../../compiledNetworkData.csv"
            # file_exists = os.path.isfile(compiledNetworkData)
            # with open(compiledNetworkData, 'a' if file_exists else 'w',newline='') as csvfile:
            #     fieldnames = network_data.keys()
            #     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            #     if not file_exists:
            #         writer.writeheader()
            #     writer.writerow(network_data)

            if args.export_to_phy:
                # Export to Phy format
                phy_folder = os.path.join(desired_output, "phy_output")
                if not os.path.exists(phy_folder):
                    os.makedirs(phy_folder, exist_ok=True)
                si.export_to_phy(
                    sorting_analyzer=curated_analyzer,
                    output_folder=phy_folder,
                    remove_if_exists=True,
                    **job_kwargs
                )
                logger.info(f"Exported sorting to Phy format at: {phy_folder}")
            # Save final checkpoint
            checkpoint.save_checkpoint(
                ProcessingStage.REPORTS_COMPLETE,
                num_units_filtered=numunits
            )

            logger.info(f"[STAGE 3] Reports complete - {numunits} units detected")


            electrodes = None


        # Cleanup folders
        if clear_temp_files and checkpoint.is_completed():
            logger.info("Clearing stored files...")
            binary_folder = os.path.join(desired_output, "binary")
            if os.path.exists(binary_folder):
                shutil.rmtree(binary_folder)
                logger.info(f"Removed binary folder: {binary_folder}")
            sorting_folder = os.path.join(desired_output, "analyzer_output")
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
    parser.add_argument('--no-curation', action='store_true', help='Skip automatic curation')
    parser.add_argument("--output-dir", type=str, help="Output directory (overrides default)")

    args = parser.parse_args()
    path = args.file_dir_path
    well = args.well
    output_root = args.output_dir or os.path.join(BASE_FILE_PATH, "../AnalyzedData")
    output_root = os.path.normpath(os.path.abspath(output_root))

    # Ensure output directory exists
    os.makedirs(output_root, exist_ok=True)

    # Override checkpoint_dir if not provided
    if args.checkpoint_dir:
        checkpoint_root = os.path.normpath(os.path.abspath(args.checkpoint_dir))
    else:
        checkpoint_root = os.path.join(output_root, "checkpoints")
    #parse file path
    run_id,project_name,date,chip_id,desired_pattern = parse_h5_path(args.file_dir_path)
    logger = setup_logger(log_file=os.path.join(output_root, desired_pattern, "logs", f"{run_id}_{chip_id}_{date}_{well}.log"))
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
        process_block(path, stream_id=well, thresholds=thresholds, clear_temp_files=clear_temp_files, checkpoint_root=checkpoint_root, output_root=output_root)
    except Exception as e:
        logger.error(f"Processing failed for {well}: {e}")
        logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
    


