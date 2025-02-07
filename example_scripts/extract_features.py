# this test script is for testing typical use case of extract_features.py
import os
import glob
#from RBS_network_models import extract_features as ef
from . import extract_features_helper as ef

#load in conv_params and mega_params
#from RBS_network_models.CDKL5.DIV21.src.conv_params import conv_params, mega_params
# from RBS_network_models.Organoid_RTT_R270X.DIV112_WT.src.conv_params import conv_params, mega_params
  #NOTE: You will need to import this from the appropriate directory for your this repo.
  # "." => allows you to import from the current directory
from .conv_params import conv_params, mega_params
# =============================================================================
#script_path = os.path.abspath(__file__)
#os.chdir(os.path.dirname(script_path))
raw_data_paths = [ 
    # NOTE: this is a list of paths to raw data files that you want to extract features from
    #      this is useful for batch processing.
    #      Also, NOTE: if parent dirs are provided, each path will be searched recursively for .h5 files
    #'/pscratch/sd/a/adammwea/workspace/_raw_data/CDKL5-E6D_T2_C1_05212024/240611/M08029/Network/000091/data.raw.h5',
    '/pscratch/sd/a/adammwea/workspace/_raw_data/Organoid_RTT_R270X_pA_pD_B1_d91/250107/M07297/Network/000028/data.raw.h5',
        # NOTE: change to directory of raw data
]
sorted_data_dir = (
    # this should be the parent directory of all sorted data files, ideally following proper data structure, naming conventions
    #'../data/CDKL5/DIV21/sorted'
    #'**/data/CDKL5/DIV21/sorted'   #syntax for glob.glob
    '**/data/Organoid_RTT_R270X/DIV112_WT/sorted'
    # NOTE: change to directory where you want sorted data to be saved
    #'**'
)
sorted_data_dir = glob.glob(sorted_data_dir, recursive=True)[0]
# =============================================================================
'''main'''
#conv_params = CDKL5.DIV21.src.conv_params
network_metrics_output_dir = os.path.join(os.path.dirname(sorted_data_dir), 'network_metrics')
network_metrics_output_dir = os.path.abspath(network_metrics_output_dir)
sorted_data_dir = os.path.abspath(sorted_data_dir)
feature_data = ef.extract_network_features(
    raw_data_paths,
    sorted_data_dir=sorted_data_dir,
    output_dir = network_metrics_output_dir,
    stream_select=None,
    conv_params=conv_params,
    mega_params=mega_params,
    # plot=True,
    )
print("Network metrics saved.")

'''
module load conda
conda activate netsims_env
python /pscratch/sd/a/adammwea/workspace/RBS_network_models/scripts/Organoid_RTT_R270X/DIV112_WT/extract_features.py
'''

# NOTE: use this script to calibrate conv params and mega params to experimental data. Rerun to update network metric plots.