import os
import shutil

def extract_csv_details(csv_dirs):
    
    if isinstance(csv_dirs, str):
        csv_dirs = [csv_dirs]

    records = []
    for csv_dir in csv_dirs:

        parent_dir = os.path.dirname(csv_dir)
        runID = os.path.basename(parent_dir)

        grandparent_dir = os.path.dirname(parent_dir)
        scan_type = os.path.basename(grandparent_dir)

        great_grandparent_dir = os.path.dirname(grandparent_dir)
        chipID = os.path.basename(great_grandparent_dir)

        ggg_dir = os.path.dirname(great_grandparent_dir)
        date = os.path.basename(ggg_dir)

        record = {'csv_file_path': csv_dir, 
                  'runID': runID, 
                  'scanType': scan_type, 
                  'chipID': chipID,
                  'date': date}
        records.append(record)

    return records

# walk through data directory and copy all csv files to a new directory
data_dir = '/home/adam/workspace/git_workspace/MEA_Analysis_fork/AxonReconPipeline/data'
destination_dir = '/home/adam/workspace/git_workspace/MEA_Analysis_fork/AxonReconPipeline/data_csvs'

for root, dirs, files in os.walk(data_dir):
    for file in files:
        if file.endswith('.csv'):
            #turn root into prefixes in file name
            recon_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(root))))
            recon_folder = os.path.basename(recon_dir)
            assert 'reconstructions'in recon_folder, 'Some pathing error'
            csv_details = extract_csv_details([os.path.join(root, file)])
            date = csv_details[0]['date']
            chipID = csv_details[0]['chipID']
            runID = csv_details[0]['runID']
            scanType = csv_details[0]['scanType']
            #prefixes = '_'.join([recon_folder, date, chipID, runID, scanType])
            prefixes = '_'.join([date, chipID, runID, scanType])
            new_file_name = prefixes + '_' + file
            final_destination_dir = os.path.join(destination_dir, recon_folder, new_file_name)           
            final_parent_dir = os.path.dirname(final_destination_dir)
            if not os.path.exists(final_parent_dir): os.makedirs(final_parent_dir)
            shutil.copy(os.path.join(root, file), final_destination_dir)
            print(f'Copied {file} to {final_destination_dir}')