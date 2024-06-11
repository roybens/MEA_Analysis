# This is an edited version of the template extraction script from the axon_tracking package by Philipp Hornauer.
# All edits tracked on github repo: roybens/MEA_analysis
# Source: https://github.com/hornauerp/axon_tracking/blob/main/axon_tracking/

import os, h5py
import spikeinterface.full as si
import numpy as np
import scipy as sp
import sklearn as sk
from tqdm import tqdm
import shutil
import logging
import sys

#local
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from MEAProcessingLibrary import mea_processing_library as MPL
#from merge_templates import merge_templates
#from old_merg_func import merge_templates
from old_merg_func_multi import merge_templates

#from axon_tracking import spike_sorting as ss
import spike_sorting as ss

#Logger Setup
#Create a logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Create handlers
stream_handler = logging.StreamHandler()  # logs to console
#file_handler = logging.FileHandler('file.log')  # logs to a file

# Set level of handlers
stream_handler.setLevel(logging.DEBUG)
#file_handler.setLevel(logging.ERROR)

# Add handlers to the logger
logger.addHandler(stream_handler)
#logger.addHandler(file_handler)

# Create formatters and add it to handlers
#formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - %(module)s.%(funcName)s')
stream_handler.setFormatter(formatter)

def extract_templates_from_sorting_dict(sorting_list, h5_file_path, qc_params={}, te_params={}, stream_select = None, just_load_waveforms = False):
    #rec_list = list(sorting_dict.keys())
    if isinstance(h5_file_path, str):
        rec_list = [h5_file_path]    

    for rec_path in rec_list:
        #sorting_list = sorting_dict[rec_path]
        #sorting_list = sorting_dict
        #sorting_list.sort()

        for sorting in sorting_list:
        #for sorting_path in sorting_list:
            #sorting = si.KiloSortSortingExtractor(sorting_path)
            sorting_path = sorting._kwargs['folder_path']
            #Temporary debug
            #if 'Network' in stream_id, replace with 'AxonTracking'
            if 'Network' in sorting_path:
                sorting_path = sorting_path.replace('Network', 'AxonTracking')
            stream_id = [p for p in sorting_path.split('/') if p.startswith('well')][0] #Find out which well this belongs to
            if stream_select is not None:
                #stream_select = 3  # for example
                stream_select_str = f'well{stream_select:03}'
                #stream_ids = stream_ids[stream_select]
                if stream_select_str not in stream_id:
                    continue
            
            # print(stream_id)
            # #debug:
            # stream_number = int(stream_id[-1])
            # if stream_number > 0:
            #     continue
            # #debug end
            
            #rec_names, common_el, pos = ss.find_common_electrodes(rec_path, stream_id)
            multirecording, pos = ss.concatenate_recording_slices(rec_path, stream_id)

            #duration = int(h5['assay']['inputs']['record_time'][0].decode('UTF-8')) * n_recs #In case we want to use firing rate as criterion
            cleaned_sorting = select_units(sorting, **qc_params)
            cleaned_sorting = si.remove_excess_spikes(cleaned_sorting, multirecording) #Relevant if last spike time == recording_length
            cleaned_sorting.register_recording(multirecording)
            segment_sorting = si.SplitSegmentSorting(cleaned_sorting, multirecording)
            extract_all_templates(h5_file_path, stream_id, segment_sorting, pos, te_params, save_root=None, just_load_waveforms = just_load_waveforms)          
def get_assay_information(rec_path):
    h5 = h5py.File(rec_path)
    pre, post, well_id = -1, -1, 0
    while pre <= 0 or post <= 0: #some failed axon trackings give negative trigger_post values, so we try different wells
        well_name = list(h5['wells'].keys())[well_id]
        rec_name = list(h5['wells'][well_name].keys())[well_id]
        pre = h5['wells'][well_name][rec_name]['groups']['routed']['trigger_pre'][0]
        post = h5['wells'][well_name][rec_name]['groups']['routed']['trigger_post'][0]
        well_id += 1
        
    return [pre, post]
def find_successful_sortings(path_list, save_path_changes):

    sorting_dict = dict()
    for rec_path in path_list:
        save_root = ss.convert_rec_path_to_save_path(rec_path, save_path_changes)
        
        #Takes into account different sorting folder names, subfolder depth, well IDs etc.
        sorting_files = [root
                         for root, dirs, files in os.walk(save_root)
                         for name in files
                         if name == "templates.npy"]
        sorting_dict[rec_path] = sorting_files
        
    return sorting_dict         
def postprocess_sorting():
    #Maybe we will do some postprocessing before we use them
    return
def select_units(sorting, min_n_spikes=500, exclude_mua=True):
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
    
    return cleaned_sorting
def extract_waveforms(h5_file_path, segment_sorting, stream_id, save_root, n_jobs, overwrite_wf):
    try:
        save_root = segment_sorting._kwargs['recording_or_recording_list'][0]._kwargs['parent_recording']._kwargs['recording']._kwargs['file_path']
    except:
        save_root = segment_sorting._kwargs['recording_or_recording_list'][0]._parent_recording._kwargs['folder_path']
        #segment_sorting._kwargs['recording_list'][0]._kwargs['parent_recording']._kwargs['folder_path']
        #get parent folder instead of file path
        save_root = os.path.dirname(save_root)
        #replace 'recordings' with waveforms
        save_root = save_root.replace('recordings', 'waveforms')
    full_path = h5_file_path
    cutout = [x / (segment_sorting.get_sampling_frequency()/1000) for x in get_assay_information(full_path)] #convert cutout to ms
    h5 = h5py.File(full_path)
    rec_names = list(h5['wells'][stream_id].keys())
    
    for sel_idx, rec_name in enumerate(rec_names):
        wf_path = os.path.join(save_root, 'waveforms', 'seg' + str(sel_idx))        
        rec = si.MaxwellRecordingExtractor(full_path,stream_id=stream_id,rec_name=rec_name)
        chunk_size = np.min([10000, rec.get_num_samples()]) - 100 #Fallback for ultra short recordings (too little activity)
        rec_centered = si.center(rec, chunk_size=chunk_size)            
        seg_sort = si.SelectSegmentSorting(segment_sorting, sel_idx)
        seg_sort = si.remove_excess_spikes(seg_sort, rec_centered)
        seg_sort.register_recording(rec_centered) 
        #Try loading waveforms, if they don't exist or overwrite_wf is True, extract them
        if not overwrite_wf:
            try: 
                seg_we = si.WaveformExtractor.load(wf_path,
                                                with_recording=True, 
                                                sorting = seg_sort)                
                # seg_we = si.load_waveforms(wf_path,
                #                            with_recording = True, 
                #                            sorting = seg_sort)
                if len(seg_sort.get_unit_ids()) == len(seg_we.unit_ids):
                    logger.info(f'Waveforms loaded from {wf_path} already exist')
                    continue
                    #break
                else:
                    raise ValueError('Unit IDs do not match')    
            except:
                if os.path.exists(wf_path):
                    shutil.rmtree(wf_path)
                pass

        if not os.path.exists(wf_path) or overwrite_wf:                           
            logger.info(f'Extracting waveforms to {wf_path}')
            os.makedirs(wf_path, exist_ok=True)
            seg_we = si.WaveformExtractor.create(rec_centered, seg_sort,
                                                 wf_path, 
                                                 allow_unfiltered=True,
                                                 remove_if_exists=True)
            seg_we.set_params(ms_before=cutout[0], ms_after=cutout[1], return_scaled = True)
            seg_we.run_extract_waveforms(n_jobs=n_jobs)
def align_waveforms(seg_we, sel_unit_id, cutout, ms_peak_cutout, upsample, align_cutout, rm_outliers, n_jobs, n_neighbors):
    
    sample_peak_cutout = ms_peak_cutout * upsample
    peak_idx = cutout[0] * upsample
    peak_cutout = range(np.int16(peak_idx - sample_peak_cutout), np.int16(peak_idx + sample_peak_cutout))
    wfs = seg_we.get_waveforms(sel_unit_id)
    interp_wfs = sp.interpolate.pchip_interpolate(list(range(wfs.shape[1])), wfs, np.linspace(0,wfs.shape[1], num = wfs.shape[1]*upsample), axis=1)
    interp_wfs = interp_wfs - np.median(interp_wfs, axis=1)[:,np.newaxis,:]

    if align_cutout:
        peak_el = [np.where(interp_wfs[w,peak_cutout,:] == np.nanmin(interp_wfs[w,peak_cutout,:]))[1][0] for w in range(interp_wfs.shape[0])]
        ref_el, count = sp.stats.mode(peak_el, keepdims=False)
        peak_shift = [np.where(interp_wfs[w,peak_cutout,ref_el] == np.nanmin(interp_wfs[w,peak_cutout,ref_el]))[0][0] for w in range(interp_wfs.shape[0])]
        aligned_length = interp_wfs.shape[1] - 2*sample_peak_cutout
        aligned_wfs = np.full([interp_wfs.shape[0], np.int16(aligned_length), interp_wfs.shape[2]], np.nan)
        for w in range(interp_wfs.shape[0]):
            aligned_wfs[w,:,:] = interp_wfs[w,peak_shift[w]:np.int16(peak_shift[w]+aligned_length),:]
    else:
        aligned_wfs = interp_wfs
    

    if rm_outliers:
        peak_el = [np.where(interp_wfs[w,peak_cutout,:] == np.nanmin(interp_wfs[w,peak_cutout,:]))[1][0] for w in range(interp_wfs.shape[0])]
        ref_el, count = sp.stats.mode(peak_el, keepdims=False)
        aligned_wfs = remove_wf_outliers(aligned_wfs, ref_el, n_jobs, n_neighbors)

    aligned_template = np.median(aligned_wfs, axis=0)
    
    return aligned_template
def remove_wf_outliers(aligned_wfs, ref_el, n_jobs, n_neighbors):
    clf = sk.neighbors.LocalOutlierFactor(n_jobs=n_jobs, n_neighbors=n_neighbors)
    outlier_idx = clf.fit_predict(aligned_wfs[:,:,ref_el])
    #print(f'Detected {sum(outlier_idx==-1)} outliers')
    outlier_rm = np.delete(aligned_wfs, outlier_idx==-1, axis=0)
    
    return outlier_rm

import reconstruct_axons as ra
def combine_templates(h5_file_path, stream_id, segment_sorting, sel_unit_id, save_root, peak_cutout=2, align_cutout=True, upsample=2, rm_outliers=True, n_jobs=16, n_neighbors=10, overwrite_wf=False, overwrite_tmp = True, just_load_waveforms = False):
    #full_path = segment_sorting._kwargs['recording_list'][0]._parent_recording._kwargs['recording']._kwargs['file_path']
    #full_path = segment_sorting._kwargs['recording_or_recording_list'][0]._kwargs['parent_recording']._kwargs['recording']._kwargs['file_path']
    full_path = h5_file_path
    # cutout = get_assay_information(full_path)
    # if align_cutout:
    #     wf_length = np.int16((sum(cutout) - 2*peak_cutout) * upsample) #length of waveforms after adjusting for potential peak alignments
    #     template_matrix = np.full([wf_length, 26400], np.nan)
    # else:
    #     wf_length = np.int16(sum(cutout) * upsample)
    #     template_matrix = np.full([wf_length, 26400], np.nan)
        

        
    
    if not just_load_waveforms:
        extract_waveforms(h5_file_path, segment_sorting, stream_id, save_root, n_jobs, overwrite_wf)

    h5 = h5py.File(full_path)
    rec_names = list(h5['wells'][stream_id].keys())
    #Find a way to average common electrodes
    # Initialize count_matrix with the same shape as template_matrix, filled with zeros
    #count_matrix = np.zeros(template_matrix.shape)
    template_list = []
    channel_locations_list = []
    #replace templates for waveforms in save_root
    if 'templates' in save_root:
        wf_load_dir = save_root.replace('templates', 'waveforms')
        
    # aw = True
    # if aw:
    for sel_idx, rec_name in enumerate(rec_names): 
        try: rec = si.MaxwellRecordingExtractor(full_path,stream_id=stream_id,rec_name=rec_name)
        except Exception as e: 
            print(f'Error loading recording {rec_name} in segment {sel_idx}:\n{e}')
            continue
        els = rec.get_property("contact_vector")["electrode"]
        seg_sort = si.SelectSegmentSorting(segment_sorting, sel_idx)
        seg_we = si.load_waveforms(os.path.join(wf_load_dir, 'waveforms', 'seg' + str(sel_idx)), sorting = seg_sort)
        template = seg_we.get_template(sel_unit_id)
        if np.isnan(template).all():
            logger.info(f'Unit {sel_unit_id} in segment {sel_idx} is empty')
            continue
        channel_locations = seg_we.get_channel_locations()
        # Define the directory and filename
        dir_path = os.path.join(save_root, 'partial_templates')
        file_name = f'seg{sel_idx}_unit{sel_unit_id}.npy'
        channel_loc_file_name = f'seg{sel_idx}_unit{sel_unit_id}_channels.npy'
        # Create the directory if it doesn't exist
        os.makedirs(dir_path, exist_ok=True)
        # Combine the directory and filename to get the full path
        channel_loc_save_file = os.path.join(dir_path, channel_loc_file_name)
        template_save_file = os.path.join(dir_path, file_name)
        
        if not np.isnan(template).all() and np.isnan(template).any():
            logger.info(f'Unit {sel_unit_id} in segment {sel_idx} has NaN values')
            #template = np.nan_to_num(template)    
        np.save(template_save_file, template)
        np.save(channel_loc_save_file, channel_locations)
        logger.info(f'Partial template saved to {template_save_file}')
        channel_locations_list.append(channel_locations)
        template_list.append(template)
    
    logger.info(f'Merging partial templates')
    merged_template, merged_channel_loc, merged_template_filled, merged_channel_loc_filled = merge_templates(template_list, 
                                                         channel_locations_list, 
                                                         plot_dir = save_root+'/merged_templates')
    
    return merged_template, merged_channel_loc, merged_template_filled, merged_channel_loc_filled       

    
    # ph = False
    # if ph:
    #     for sel_idx, rec_name in enumerate(rec_names): 
    #         rec = si.MaxwellRecordingExtractor(full_path,stream_id=stream_id,rec_name=rec_name)
    #         els = rec.get_property("contact_vector")["electrode"]
    #         seg_sort = si.SelectSegmentSorting(segment_sorting, sel_idx)
    #         seg_we = si.load_waveforms(os.path.join(save_root, 'waveforms', 'seg' + str(sel_idx)), sorting = seg_sort)
    #         aligned_wfs = align_waveforms(seg_we, 
    #                                     sel_unit_id, 
    #                                     cutout, 
    #                                     peak_cutout, 
    #                                     upsample, 
    #                                     align_cutout, 
    #                                     rm_outliers, 
    #                                     n_jobs, 
    #                                     n_neighbors)
    #         template_matrix[:,els] = aligned_wfs #find way to average common electrodes
    
    
    #     return template_matrix

def convert_to_grid(template_matrix, pos):

    clean_template = np.delete(template_matrix, np.isnan(pos['x']), axis = 1)
    clean_x = pos['x'][~np.isnan(pos['x'])]
    clean_y = pos['y'][~np.isnan(pos['y'])]
    x_idx = np.int16(clean_x / 17.5)
    y_idx = np.int16(clean_y / 17.5)
    grid = np.full([np.max(x_idx) + 1, np.max(y_idx) + 1, clean_template.shape[0]],0).astype('float32')
    for i in range(len(y_idx)):
        grid[x_idx[i],y_idx[i],:] = clean_template[:,i]    
    return grid

def get_all_templates(h5_file_path, stream_id, segment_sorting, pos, te_params, save_root=None, just_load_waveforms = False):
    sel_unit_ids = segment_sorting.get_unit_ids()
    if save_root is None:
        recording_details = MPL.extract_recording_details(h5_file_path)
        date = recording_details[0]['date']
        chip_id = recording_details[0]['chipID']
        scanType = recording_details[0]['scanType']
        run_id = recording_details[0]['runID']
        save_root = f'./AxonReconPipeline/data/temp_data/templates/{date}/{chip_id}/{scanType}/{run_id}/{stream_id}'
    template_save_path = os.path.join(save_root, 'templates')
    if not os.path.exists(template_save_path):
        os.makedirs(template_save_path)
        
    for sel_unit_id in tqdm(sel_unit_ids):

        template_save_file = os.path.join(template_save_path, str(sel_unit_id) + '.npy')
        channel_loc_save_file = os.path.join(template_save_path, str(sel_unit_id) + '_channels.npy')
        template_save_file_fill = os.path.join(template_save_path, str(sel_unit_id) + '_filled.npy')
        channel_loc_save_file_fill = os.path.join(template_save_path, str(sel_unit_id) + '_channels_filled.npy')
        
        if not os.path.isfile(template_save_file) or te_params['overwrite_tmp']:
            try:
                merged_template, merged_channel_loc, merged_template_filled, merged_channel_locs_filled = combine_templates(h5_file_path, stream_id, segment_sorting, sel_unit_id, save_root, just_load_waveforms = just_load_waveforms, **te_params)
                #grid = convert_to_grid(template_matrix, pos)
                np.save(channel_loc_save_file, merged_channel_loc)
                np.save(template_save_file, merged_template)
                np.save(channel_loc_save_file_fill, merged_channel_locs_filled)
                np.save(template_save_file_fill, merged_template_filled)
                logger.info(f'Merged template saved to {save_root}/merged_templates')
            except Exception as e:
                print(f'Unit {sel_unit_id} encountered the following error:\n {e}')

