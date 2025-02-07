import os
import shutil

project_name = 'B6J_DensityTest_10012024_AR'

# Define the source directory
source_dir = '/pscratch/sd/a/adammwea/zRBS_axon_reconstruction_output'
source_dir = os.path.join(source_dir, project_name)

# Define the destination directories
csv_destination_dir = '/pscratch/sd/a/adammwea/zRBS_axon_recon_output_CSVs'
csv_destination_dir = os.path.join(csv_destination_dir, project_name)
dill_destination_dir = '/pscratch/sd/a/adammwea/zRBS_axon_recon_output_DILLs'
dill_destination_dir = os.path.join(dill_destination_dir, project_name)

# Ensure the destination directories exist
os.makedirs(csv_destination_dir, exist_ok=True)
os.makedirs(dill_destination_dir, exist_ok=True)

# Iterate over the files in the source directory
for root, _, files in os.walk(source_dir):
    for filename in files:
        # Check if the file contains 'gtr' in the name
        # Construct the full file path
        source_file = os.path.join(root, filename)
        # Create a new filename excluding all path parts above the source_dir
        relative_path = os.path.relpath(source_file, source_dir)
        new_filename = relative_path.replace('/', '_').replace(':', '_')
        
        # Determine the destination directory based on file extension
        if filename.endswith('.csv'):
            destination_file = os.path.join(csv_destination_dir, new_filename)
        elif filename.endswith('.dill'):
            if 'gtr' in filename:
                destination_file = os.path.join(dill_destination_dir, new_filename)
        else:
            continue  # Skip files that are not .csv or .dill
        
        # Copy the file to the destination directory
        shutil.copy2(source_file, destination_file)
        print(f"Copied {filename} to {destination_file}")