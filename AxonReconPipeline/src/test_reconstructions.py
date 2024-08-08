from AxonReconPipeline.axon_velocity import axon_velocity as av
import numpy as np
import os
## define initial params, see run_pipeline.py av_params for more details

### import needed things from lib_axon_velocity_functions.py and func_analyze_and_reconstruct

# load templates (time and amplitude info)
template_dir = '' #put the dir where the templates are stored
dirs = os.path.listdir(template_dir) #list of all the directories in the template_dir
templates = {} # makes a dictionary with unit_id as key and template as value
for dir in dirs: 
    unit_id = dir.split('_')[-1] #this might not work correctly, make sure this gets unit_id
    templates[unit_id] = np.load(dir + "templates.npy")

# load channel locations (xy of channels)
channel_loc_dir = ''    
dirs = os.listdir(channel_loc_dir)
channel_locs = {}

#transform templates and locations
#import transform template from lib_axonb_velocity_functions.py
#use it, get transformed temps and chan locs
#this just prepares data for graphaxontracking class

#modify params as needed

# pass to axon_velocity (e.g transformed_template = templates[some_unit_id].T)
# gtr = GraphAxonTracking(transformed_template, trans_loc, 10000, **params)
# gtr.track_axons()


