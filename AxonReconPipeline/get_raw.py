path = '/mnt/disk15tb/adam/git_workspace/MEA_Analysis/AxonReconPipeline/data/temp_data/reconstructions/240503/M08034/AxonTracking/000082/well000'
#get all .pngs with 'raw' in the name and copy them to a folder inside AxonReconPipleine
# use walk
import os
import shutil
for root, dirs, files in os.walk(path):
    for file in files:
        if 'raw' in file:
            destination_in_axonrecon = '/mnt/disk15tb/adam/git_workspace/MEA_Analysis/AxonReconPipeline/data/reconsraw'
            unit_id = os.path.basename(os.path.dirname(root))
            if not os.path.exists(destination_in_axonrecon):
                os.makedirs(destination_in_axonrecon)
            shutil.copyfile(root + '/' + file, destination_in_axonrecon + '/' + f'{unit_id}_{file}')

