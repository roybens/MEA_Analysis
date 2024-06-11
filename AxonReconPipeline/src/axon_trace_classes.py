#axon_trace_object_classes

import numpy as np
import matplotlib.pyplot as plt
from tsmoothie.smoother import GaussianSmoother
import spikeinterface
import spikeinterface.full as si
import spikeinterface.extractors as se
import spikeinterface.sorters as ss
import spikeinterface.comparison as sc
import spikeinterface.widgets as sw
import spikeinterface.postprocessing as sp
import spikeinterface.preprocessing as spre
import spikeinterface.qualitymetrics as qm
import lib_helper_functions as helper
from pprint import pprint
import os
import re
from spikeinterface import compute_sparsity
import numpy as np
from scipy.signal import correlate
from sklearn.metrics.pairwise import cosine_similarity

class axon_trace_objects:
    """Class for encapsulating data for a recording"""
    
    def __init__(self, base_path, dirs):
        self.base_path = base_path
        self.sort_subdir = os.path.basename(os.path.normpath(base_path))  # Load subdirectory name

        # Initialize data containers as empty dictionaries
        self.waveforms = {}   
        self.templates = {}
        self.unit_ids = {}
        self.unit_locations = {}
        self.unit_amplitudes = {}
        self.unit_extremum_channel_ID = {}
        self.unit_extremum_channel_index = {}
        #self.extremum_channel_location = {}
        self.sparse_unit_templates = {}
        self.unit_templates = {}
        self.channel_locations = {}
        self.sampling_frequency = {}
        self.sparse_channel_ids = {}
        self.sparse_channel_locations = {}
        self.sparse_channel_indices = {}  
        #self.sparse_templates = {}

        self._load_data(dirs)  # Load data


    def _load_data(self, dirs):
        """Helper method to load all data for this recording"""
        
        # Extract scan types and run IDs from dirs
        scan_types_and_run_ids = [(re.search(r'(ActivityScan|AxonTracking|Network)/(\d{1,6})', dir).group(1), 
                                   re.search(r'(ActivityScan|AxonTracking|Network)/(\d{1,6})', dir).group(2)) 
                                  for dir in dirs if re.search(r'(ActivityScan|AxonTracking|Network)/(\d{1,6})', dir)]

        # Iterate over all subdirectories in base_path
        for sub_dir in sorted(os.listdir(self.base_path)):
            # Check if 'kilosort2' is in sub_dir
            if 'kilosort2' in sub_dir:
                # Extract scan type and run ID from sub_dir
                match = re.search(r'(ActivityScan|AxonTracking|Network)(\d{1,6})', sub_dir)
                if match:
                    scan_type_dir = match.group(1)
                    run_id_dir = match.group(2)

                    # Check if scan type and run ID match any in scan_types_and_run_ids
                    if (scan_type_dir, run_id_dir.zfill(6)) in scan_types_and_run_ids:
                        # Define run_key and dict_key
                        run_key = scan_type_dir + run_id_dir
                        dict_key = re.search(r'rec\d{4}', sub_dir).group(0) if re.search(r'rec\d{4}', sub_dir) else None

                        # Continue processing this sub_dir...                                        

                        # Initialize run_key in each dictionary if it doesn't exist
                        for dict_name in [
                            #'self.sparse_templates',
                            'self.sparse_channel_locations',
                            'self.sparse_channel_ids',
                            'self.sparse_channel_indices',
                            'self.unit_extremum_channel_index',
                            'self.sampling_frequency', 
                            'self.channel_locations', 
                            'self.waveforms', 
                            'self.unit_ids', 
                            'self.templates', 
                            'self.sparse_unit_templates', 
                            'self.unit_templates', 
                            'self.unit_amplitudes', 
                            'self.unit_locations', 
                            'self.unit_extremum_channel_ID']:
                            if run_key not in eval(dict_name) or not isinstance(eval(dict_name)[run_key], dict):
                                eval(dict_name)[run_key] = {}

                        try:
                            # Load waveforms
                            wf_subfolder = sub_dir.replace('kilosort2', 'waveforms_good')
                            wf_subfolder = os.path.join(self.base_path, wf_subfolder)
                            self.waveforms[run_key][dict_key] = si.load_waveforms(wf_subfolder)
                            
                            # Get unit ids
                            waveform_folder = self.waveforms[run_key][dict_key].folder
                            waveform_folder_str = str(waveform_folder)
                            #sorter_folder_str = waveform_folder_str.replace("waveforms_good", "kilosort2")
                            #sorting_KS3 = ss.Kilosort2Sorter._get_result_from_folder(sorter_folder_str + "/sorter_output")
                            wf_unit_ids = self.waveforms[run_key][dict_key].unit_ids

                            #Check if sorted_unit_ids is empty, if so, continue
                            #if not sorted_unit_ids:
                            if wf_unit_ids.size == 0:
                                print(f"Warning: sorted_unit_ids is empty for {run_key} {dict_key}")
                                self.unit_ids[run_key][dict_key] = None
                                continue
                            self.unit_ids[run_key][dict_key] = wf_unit_ids
                            
                            # Extract sparse channel IDs and indices
                            sparse_channel_ids = {}
                            sparse_channel_indices = {}
                            i = 0
                            for unit_id in self.unit_ids[run_key][dict_key]:
                                if run_key not in sparse_channel_ids:
                                    sparse_channel_ids[run_key] = {}
                                    sparse_channel_indices[run_key] = {}
                                if dict_key not in sparse_channel_ids[run_key]:
                                    sparse_channel_ids[run_key][dict_key] = {}
                                    sparse_channel_indices[run_key][dict_key] = {}
                                sparse_channel_ids[run_key][dict_key][unit_id] = self.waveforms[run_key][dict_key].sparsity.unit_id_to_channel_ids[wf_unit_ids[i]]
                                sparse_channel_indices[run_key][dict_key][unit_id] = self.waveforms[run_key][dict_key].sparsity.unit_id_to_channel_indices[wf_unit_ids[i]]
                                i += 1
                            self.sparse_channel_ids[run_key][dict_key] = sparse_channel_ids[run_key][dict_key]
                            self.sparse_channel_indices[run_key][dict_key] = sparse_channel_indices[run_key][dict_key]                       

                            # Compute templates
                            template_means = self.waveforms[run_key][dict_key].get_all_templates(mode="average")
                            self.templates[run_key][dict_key] = template_means                            
                            #print(self.waveforms[run_key][dict_key])

                           # Replace template dict indexes with keys matching unit_ids
                            unit_id_to_template = {}
                            unit_id_to_template_sparse = {}
                            i = 0
                            for unit_id in self.unit_ids[run_key][dict_key]:
                                if run_key not in unit_id_to_template:
                                    unit_id_to_template[run_key] = {}
                                    unit_id_to_template_sparse[run_key] = {}
                                if dict_key not in unit_id_to_template[run_key]:
                                    unit_id_to_template[run_key][dict_key] = {}
                                    unit_id_to_template_sparse[run_key][dict_key] = {}
                                # Implement sparsity mask on templates
                                mask = self.waveforms[run_key][dict_key].sparsity.mask[i]
                                unit_id_to_template[run_key][dict_key][unit_id] = self.templates[run_key][dict_key][i]
                                unit_id_to_template_sparse[run_key][dict_key][unit_id] = self.templates[run_key][dict_key][i][:, mask]
                                i += 1
                            self.sparse_unit_templates[run_key][dict_key] = unit_id_to_template_sparse[run_key][dict_key]  
                            self.unit_templates[run_key][dict_key] = unit_id_to_template[run_key][dict_key]          

                            # Get amplitudes
                            self.unit_amplitudes[run_key][dict_key] = si.compute_spike_amplitudes(self.waveforms[run_key][dict_key], load_if_exists=True, outputs="by_unit")
                            
                            # Get locations
                            self.unit_locations[run_key][dict_key] = sp.compute_unit_locations(self.waveforms[run_key][dict_key], load_if_exists=True, outputs="by_unit")

                            # Get channel locations
                            self.channel_locations[run_key][dict_key] = self.waveforms[run_key][dict_key].get_channel_locations() 

                            #Get sparse channel locations, based on sparse_channel_indices (index)
                            sparse_channel_locations = {}
                            i = 0
                            for unit_id in self.unit_ids[run_key][dict_key]:
                                if run_key not in sparse_channel_locations:
                                    sparse_channel_locations[run_key] = {}
                                if dict_key not in sparse_channel_locations[run_key]:
                                    sparse_channel_locations[run_key][dict_key] = {}
                                #sparse_channel_indices_int = np.array(self.sparse_channel_indices[run_key][dict_key][unit_id]).astype(int)
                                sparse_channel_indices_int = self.sparse_channel_indices[run_key][dict_key][unit_id]
                                sparse_channel_locations[run_key][dict_key][unit_id] = self.channel_locations[run_key][dict_key][sparse_channel_indices_int]
                                #sparse_channel_locations[run_key][dict_key][unit_id] = self.channel_locations[run_key][dict_key][self.sparse_channel_ids[run_key][dict_key][unit_id]]
                                i += 1
                            self.sparse_channel_locations[run_key][dict_key] = sparse_channel_locations[run_key][dict_key]
                            
                            # Get extremum electrode for each unit
                            self.unit_extremum_channel_ID[run_key][dict_key] = spikeinterface.full.get_template_extremum_channel(self.waveforms[run_key][dict_key], peak_sign='neg', mode='extremum')
                            self.unit_extremum_channel_index[run_key][dict_key] = spikeinterface.full.get_template_extremum_channel(self.waveforms[run_key][dict_key], peak_sign='neg', mode='extremum', outputs='index')
                            #extremum_channel_index = spikeinterface.full.get_template_extremum_channel(self.waveforms[run_key][dict_key], peak_sign='neg', mode='extremum')
                            #self.extremum_channel_location[run_key][dict_key] = self.channel_locations[run_key][dict_key][extremum_channel_index]

                            # Get channel locations
                            self.sampling_frequency[run_key][dict_key] = self.waveforms[run_key][dict_key].sampling_frequency

                                           
                                                                        
                        except Exception as e:
                            print(f"Error loading data for {run_key} {dict_key}")
                            print(e)

                            # If error, set everything to None
                            self.waveforms[run_key][dict_key] = None
                            self.templates[run_key][dict_key] = None
                            self.unit_ids[run_key][dict_key] = None
                            self.unit_locations[run_key][dict_key] = None
                            self.unit_amplitudes[run_key][dict_key] = None
                            self.unit_extremum_channel_ID[run_key][dict_key] = None
                            self.unit_extremum_channel_index[run_key][dict_key] = None
                            self.sparse_unit_templates[run_key][dict_key] = None
                            self.channel_locations[run_key][dict_key] = None
                            self.sampling_frequency[run_key][dict_key] = None
            
            #debug
            #break
        
        self._collect_information_by_unit()

 
    def _collect_information_by_unit(self):
        
        ##Collect detailed info on each unit
        i = 0
        self.unit_templates_informed = {}
        for run_key in self.waveforms:
            for dict_key in self.waveforms[run_key]:
                ##
                # Check if sorted_unit_ids is None, if so, continue
                if self.unit_ids[run_key][dict_key] is None:
                    print(f"Warning: sorted_unit_ids is None for {run_key} {dict_key}")
                    continue                        
                for unit_id in self.unit_ids[run_key][dict_key]:
                    self.unit_templates_informed[i] = {}
                    
                    #Get unit info
                    unit_channel_locs = self.channel_locations[run_key][dict_key]
                    unit_sparse_channel_ids = self.sparse_channel_ids[run_key][dict_key][unit_id]
                    unit_sparse_channel_locs = self.sparse_channel_locations[run_key][dict_key][unit_id]            
                    unit_extremum_channel_ID = self.unit_extremum_channel_ID[run_key][dict_key][unit_id]
                    unit_extremum_channel_index = self.unit_extremum_channel_index[run_key][dict_key][unit_id]
                    #unit_extremum_channel_location = self.extremum_channel_location[run_key][dict_key] 
                    unit_sampling_frequency = self.sampling_frequency[run_key][dict_key]
                    sparse_unit_template = self.sparse_unit_templates[run_key][dict_key][unit_id]
                    unit_template = self.unit_templates[run_key][dict_key][unit_id]
                    unit_amps = self.unit_amplitudes[run_key][dict_key][0][unit_id]
                    unit_loc = self.unit_locations[run_key][dict_key][unit_id]
                    unit_source_path = self.base_path

                    #Assign values to unit_templates_informed
                    self.unit_templates_informed[i]["source_path"] = unit_source_path
                    self.unit_templates_informed[i]["run_key"] = run_key
                    self.unit_templates_informed[i]["rec_num"] = dict_key
                    self.unit_templates_informed[i]["unit_id"] = unit_id
                    self.unit_templates_informed[i]["extremum_channel_ID"] = unit_extremum_channel_ID
                    self.unit_templates_informed[i]["extremum_channel_index"] = unit_extremum_channel_index
                    self.unit_templates_informed[i]["channel_locations"] = unit_channel_locs
                    self.unit_templates_informed[i]["sparse_channel_ids"] = unit_sparse_channel_ids
                    self.unit_templates_informed[i]["sparse_channel_locations"] = unit_sparse_channel_locs
                    self.unit_templates_informed[i]["extremum_channel_location"] = unit_channel_locs[int(unit_extremum_channel_index)]
                    self.unit_templates_informed[i]["sampling_frequency"] = unit_sampling_frequency
                    self.unit_templates_informed[i]["sparse_unit_template"] = sparse_unit_template
                    self.unit_templates_informed[i]["template"] = unit_template
                    self.unit_templates_informed[i]["amplitudes"] = unit_amps
                    self.unit_templates_informed[i]["location"] = unit_loc
                    i+=1

        self._group_templates()

    def _group_templates(self):
        """
        This method groups templates based on various experimental methods.
        1. Collection of units by extremum and proximity.
        2. Grouping based on cross-correlation.
        3. Grouping based on cosine similarity.
        """
        self = self._collect_units_by_extremum_and_proximity()
        
        #under construction
        #self = self._cross_correlation_grouping()
        #self = self._cosine_similarity_grouping()
        return self

    def _collect_units_by_extremum_and_proximity(self):
        """ AW 03Jan2024
        This method collects neuronal units to be merged based on their extremum and proximity.
        It's important to note that this is done purely as a proof of concept. There are many 
        open questions around the 'right' way to determine that two units are indeed the same neuron. 
        Additionally, there may be situations where one unit should be split into two. 
        Future work should aim to refine this approach and address these considerations.
        """
        self.units_by_extremum = {}

        for unit_key, unit_info in self.unit_templates_informed.items():
            # Get the extremum_channel_location from the unit_info dictionary
            extremum_channel_location = unit_info['extremum_channel_location']

            # Convert the tuple to a string and format it to ensure each coordinate has one decimal place
            extremum_channel_location_key = f'[{extremum_channel_location[0]:.1f} {extremum_channel_location[1]:.1f}]'

            # Use a regular expression to remove any spaces after the opening bracket or before the closing bracket
            extremum_channel_location_key = re.sub(r'\[\s+', '[', extremum_channel_location_key)
            extremum_channel_location_key = re.sub(r'\s+\]', ']', extremum_channel_location_key)
            
            
            if extremum_channel_location_key not in self.units_by_extremum:
                self.units_by_extremum[extremum_channel_location_key] = {
                    'source_path': [], 
                    'run_key': [], 
                    'rec_nums': [], 
                    'unit_ids': [], 
                    'channel_locations': [],
                    'sparse_channel_ids': [],
                    'sparse_channel_locations': [],
                    'sampling_frequency': [], 
                    'sparse_unit_templates': [],
                    'templates': [],
                    'amplitudes': [],
                    'locations': [],
                    'extremum_channel_location': [], # 'extremum_channel_location' is the location of the extremum channel, e.g. '(0, 0)
                    'extemum_channel_ID': [] # 'extremum_channel_ID' is the ID of the extremum channel, e.g. '0
                }

            source_path = unit_info['source_path']
            if source_path not in self.units_by_extremum[extremum_channel_location_key]['source_path']:
                self.units_by_extremum[extremum_channel_location_key]['source_path'].append(source_path)
            self.units_by_extremum[extremum_channel_location_key]['extremum_channel_location'].append(unit_info['extremum_channel_location'])
            self.units_by_extremum[extremum_channel_location_key]['extemum_channel_ID'].append(unit_info['extremum_channel_ID'])
            self.units_by_extremum[extremum_channel_location_key]['run_key'].append(unit_info['run_key'])    
            self.units_by_extremum[extremum_channel_location_key]['rec_nums'].append(unit_info['rec_num'])
            self.units_by_extremum[extremum_channel_location_key]['unit_ids'].append(unit_info['unit_id']) 
            self.units_by_extremum[extremum_channel_location_key]['channel_locations'].append(unit_info['channel_locations'])
            self.units_by_extremum[extremum_channel_location_key]['sparse_channel_ids'].append(unit_info['sparse_channel_ids'])
            self.units_by_extremum[extremum_channel_location_key]['sparse_channel_locations'].append(unit_info['sparse_channel_locations'])
            self.units_by_extremum[extremum_channel_location_key]['sampling_frequency'].append(unit_info['sampling_frequency'])
            self.units_by_extremum[extremum_channel_location_key]['sparse_unit_templates'].append(unit_info['sparse_unit_template'])
            self.units_by_extremum[extremum_channel_location_key]['templates'].append(unit_info['template'])
            self.units_by_extremum[extremum_channel_location_key]['amplitudes'].append(unit_info['amplitudes'])
            self.units_by_extremum[extremum_channel_location_key]['locations'].append(unit_info['location'])

        # Convert keys back to tuples for distance calculation
        units_by_extremum_np = {tuple(map(float, k.strip('[]').split())): v for k, v in self.units_by_extremum.items()}

        # Calculate Euclidean distance for each key
        distances = {k: np.linalg.norm(np.array(k)) for k in units_by_extremum_np.keys()}

        # Sort the dictionary by distance to the origin (0, 0)
        sorted_units_by_extremum = dict(sorted(self.units_by_extremum.items(), key=lambda item: distances[tuple(map(float, item[0].strip('[]').split()))]))
        self.units_by_extremum = sorted_units_by_extremum

        # Calculate Euclidean distance from each key to every other key
        distance_matrix = {k1: {k2: np.linalg.norm(np.array(k1) - np.array(k2)) for k2 in units_by_extremum_np.keys()} for k1 in units_by_extremum_np.keys()}

        # Create a new distance matrix with distances <= 60 microns (60 is arbitrary)
        #euclidian_distance_val = 120
        euclidian_distance_val = 0 #Use this value to turn off euclidian distance filter
        distance_matrix = {k1: {k2: v for k2, v in k1_dict.items() if v <= euclidian_distance_val} for k1, k1_dict in distance_matrix.items()}
   
        # Create a list of all distances and sort it
        unit_distances = [(k1, k2, v) for k1, inner_dict in distance_matrix.items() for k2, v in inner_dict.items()]
        unit_distances.sort(key=lambda x: x[2])

        # Initialize a new dictionary to store units to be merged
        self.units_to_be_merged = {}

        # Initialize a set to store units that have already been added to a group
        added_units = set()

        # Initialize a list to store units to be merged
        units = []

        # Initialize a counter for the group index
        i = 0

        # First pass: handle all cases where distance > 0
        for k1, k2, distance in unit_distances:
            if distance > 0:
                k1_str = '[' + ' '.join(map(str, k1)) + ']'
                k2_str = '[' + ' '.join(map(str, k2)) + ']'

                if k1_str not in added_units:
                    units.append(self.units_by_extremum[k1_str])
                    added_units.add(k1_str)
                if k2_str not in added_units:
                    units.append(self.units_by_extremum[k2_str])
                    added_units.add(k2_str)

                if (k1_str in added_units and k2_str in added_units) or (k1, k2, distance) == unit_distances[-1]:
                    self.units_to_be_merged[i] = units
                    i += 1
                    units = []

        # Second pass: handle all remaining cases where distance = 0 and they haven't been grouped elsewhere
        for k1, k2, distance in unit_distances:
            if distance == 0:
                k1_str = '[' + ' '.join(map(str, k1)) + ']'
                k2_str = '[' + ' '.join(map(str, k2)) + ']'

                if k1_str not in added_units:
                    units.append(self.units_by_extremum[k1_str])
                    added_units.add(k1_str)
                if k2_str not in added_units:
                    units.append(self.units_by_extremum[k2_str])
                    added_units.add(k2_str)

                if (k1_str in added_units and k2_str in added_units) or (k1, k2, distance) == unit_distances[-1]:
                    self.units_to_be_merged[i] = units
                    i += 1
                    units = []
        
        #remove any empty keys from self.units_to_be_merged
        self.units_to_be_merged = {k: v for k, v in self.units_to_be_merged.items() if v}
        
        # Flatten the dictionary of units to be merged, combine all fields into a single list
        self.units_to_be_merged = {k: {k2: [item for sublist in v for item in sublist[k2]] for k2 in v[0].keys()} for k, v in self.units_to_be_merged.items()}
        
        # Remove duplicate values from "source_path" and "sampling_frequency" in each sublist
        self.units_to_be_merged = {
            k: {
                k2: list(set(v)) if k2 in ["source_path", "sampling_frequency"] else v 
                for k2, v in v.items()
            } 
            for k, v in self.units_to_be_merged.items()
        }

        # Count the total number of 'templates' in self.units_by_extremum
        total_templates_extremum = sum(len(v['sparse_unit_templates']) for v in self.units_by_extremum.values())

        # Count the total number of 'templates' in self.units_to_be_merged
        total_templates_merged = sum(len(v['sparse_unit_templates']) for v in self.units_to_be_merged.values())

        # Quality check, make sure there are equal number of templates in self.units_by_extremum and self.units_to_be_merged
        assert total_templates_extremum == total_templates_merged, "Error: total number of templates in self.units_by_extremum and self.units_to_be_merged are not equal."

        return self
    
    def _cross_correlation_grouping(self, threshold = 0.95):
        """
        Perform cross-correlation between all pairs of templates and group them based on a threshold value.

        Parameters:
        templates (list of np.array): List of templates to be cross-correlated.
        threshold (float): Threshold value for grouping templates.

        Returns:
        dict: Dictionary of grouped templates.
        """
        templates = self.unit_templates_informed 

        # Initialize dictionary for grouped templates
        grouped_templates = {}

        # Calculate cross-correlation for each pair of templates
        for i in range(len(templates)):
            for j in range(i+1, len(templates)):
                # Calculate cross-correlation
                cross_corr = correlate(templates[i]["sparse_unit_template"], templates[j]["sparse_unit_template"], mode='full')            

                # Normalize the cross-correlation
                # Calculate the norm factor
                norm_factor = np.sqrt(np.sum(templates[i]["sparse_unit_template"]**2) * np.sum(templates[j]["sparse_unit_template"]**2))
                normalized_cross_corr = cross_corr / norm_factor

                # Find maximum normalized cross-correlation value
                max_normalized_cross_corr = np.max(normalized_cross_corr)

                #assert that max_normalized_cross_corr is between -1 and 1
                assert max_normalized_cross_corr >= -1 and max_normalized_cross_corr <= 1, "Error: max_normalized_cross_corr is not between -1 and 1."

                # If maximum cross-correlation value is above threshold, group templates
                if max_normalized_cross_corr >= threshold:
                    # Create a key for the pair of templates
                    key = f'template_{i}_and_template_{j}'

                    # Add templates to the dictionary
                    grouped_templates[key] = (templates[i], templates[j])

                    #print extreemum channel location for each template
                    print(f"template_{i} extremum channel location: {templates[i]['extremum_channel_location']}")
                    print(f"template_{j} extremum channel location: {templates[j]['extremum_channel_location']}")
                    
                    #print euclidian distance between extreemum channel locations
                    print(f"euclidian distance between template_{i} and template_{j}: {np.linalg.norm(np.array(templates[i]['extremum_channel_location']) - np.array(templates[j]['extremum_channel_location']))}")

        return self

    def _cosine_similarity_grouping(self, threshold = 0.8):
        """
        Perform cosine similarity between all pairs of templates and group them based on a threshold value.

        Parameters:
        templates (list of np.array): List of templates to be cross-correlated.
        threshold (float): Threshold value for grouping templates.

        Returns:
        dict: Dictionary of grouped templates.
        """

        templates = self.unit_templates_informed
        
        # Initialize dictionary for grouped templates
        grouped_templates = {}

        # Calculate cosine similarity for each pair of templates
        for i in range(len(templates)):
            for j in range(i+1, len(templates)):
                # Calculate cosine similarity
                cos_sim = cosine_similarity(templates[i]["sparse_unit_template"].reshape(1, -1), templates[j]["sparse_unit_template"].reshape(1, -1))[0][0]

                # If cosine similarity is above threshold, group templates
                if cos_sim >= threshold:
                    # Create a key for the pair of templates
                    key = f'template_{i}_and_template_{j}'

                    # Add templates to the dictionary
                    grouped_templates[key] = (templates[i], templates[j])

        return self
