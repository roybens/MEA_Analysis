import os
import shutil
import re

def reorganize_files(file_list, base_dir, new_output_dir):
    """
    Reorganize files from an old directory structure to a new one.
    
    Args:
        file_list (str): Path to the text file containing the list of files.
        base_dir (str): Path to the base directory of the old files.
        new_output_dir (str): Path to the base directory of the new file structure.
    """
    # Read the file list
    with open(file_list, 'r') as f:
        paths = f.read().splitlines()

    for old_path in paths:
        # Ensure the file exists
        if not os.path.exists(old_path):
            print(f"File not found: {old_path}")
            continue
        
        # Extract components from the old path
        match = re.search(r'(?P<projectName>[^/]+)/(?P<date>\d{6})_(?P<chipID>\w+)_(?P<runID>\d{6})_well(?P<wellID>\d+)_axon_reconstruction(?P<error>_error)?.log', old_path)
        if match:
            components = match.groupdict()
            projectName = components["projectName"]
            date = components["date"]
            chipID = components["chipID"]
            runID = components["runID"]
            wellID = f"well{components['wellID']}"
            reconstructorID = "axon_reconstruction" + (components["error"] or "")
            scanType = "AxonTracking"  # Assuming scanType is fixed as per example

            # Build the new path
            new_path = os.path.join(new_output_dir, projectName, date, chipID, scanType, runID, wellID, f"{reconstructorID}.log")

            # Ensure the directory exists
            os.makedirs(os.path.dirname(new_path), exist_ok=True)
            
            # Move the file
            print(f"Moving \n{old_path} -> \n{new_path}")
            shutil.move(old_path, new_path)
            print(f"Moved {old_path} -> {new_path}")
        else:
            print(f"Failed to parse path: {old_path}")

# Usage example:
# reorganize_files('/path/to/file_list.txt', '/old/base/dir', '/new/base/dir')

reorganize_files('/pscratch/sd/a/adammwea/file_list.txt', '/pscratch/sd/a/adammwea/zRBS_axon_reconstruction_output', '/pscratch/sd/a/adammwea/zRBS_axon_reconstruction_output')