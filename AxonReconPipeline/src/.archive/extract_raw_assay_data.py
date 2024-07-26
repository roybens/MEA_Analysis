import os
import re
from datetime import datetime

def find_chip_ids(directories):
    include_dirs = []
    exclude_dirs = []
    date_dirs = []
    chip_dirs = []
    runID_dirs = []

    def process_chip_dir(chip_dir):
        # Check if the ActivityScan, Network, and AxonTracking directories exist
        activity_scan_dir = os.path.join(chip_dir, "ActivityScan")
        network_dir = os.path.join(chip_dir, "Network")
        axon_tracking_dir = os.path.join(chip_dir, "AxonTracking")

        # Flags to check if at least one of Network or AxonTracking scans exist
        network_or_axon_exists = os.path.exists(network_dir) or os.path.exists(axon_tracking_dir)

        if os.path.exists(activity_scan_dir) and network_or_axon_exists:
            # Include this directory
            for run_dir in os.listdir(activity_scan_dir):
                for file in os.listdir(os.path.join(activity_scan_dir, run_dir)):
                    if os.path.splitext(file)[1] == ".h5":
                        include_dirs.append(os.path.join(activity_scan_dir, run_dir, file))

            if os.path.exists(network_dir):
                for run_dir in os.listdir(network_dir):
                    for file in os.listdir(os.path.join(network_dir, run_dir)):
                        if os.path.splitext(file)[1] == ".h5":
                            include_dirs.append(os.path.join(network_dir, run_dir, file))

            if os.path.exists(axon_tracking_dir):
                for run_dir in os.listdir(axon_tracking_dir):
                    for file in os.listdir(os.path.join(axon_tracking_dir, run_dir)):
                        if os.path.splitext(file)[1] == ".h5":
                            include_dirs.append(os.path.join(axon_tracking_dir, run_dir, file))
        else:
            # Exclude this directory
            exclude_dirs.append(chip_dir)

    def process_runID_dir(runID_dir):
        # Include this directory
        for file in os.listdir(runID_dir):
            if os.path.splitext(file)[1] == ".h5":
                include_dirs.append(os.path.join(runID_dir, file))
        else:
            # Exclude this directory
            exclude_dirs.append(runID_dir)

    # If directories is a string, convert it to a list
    if isinstance(directories, str):
        directories = [directories]

    # Loop through the directories
    for directory in directories:
        dir_name = os.path.basename(directory)
        # Check if the directory name matches the format YYMMDD
        if re.match(r'\d{6}$', dir_name):
            try:
                # Try to convert the directory name to a date
                datetime.strptime(dir_name, '%y%m%d')
                date_dirs.append(directory)
            except ValueError:
                # The directory name could not be converted to a date
                try: 
                    # Check if the directory name matches the format of a runID (6 digit numbers)
                    if re.match(r'\d{6}$', dir_name):
                        parent_dir_name = os.path.basename(os.path.dirname(directory))
                        if parent_dir_name in ["ActivityScan", "Network", "AxonTracking"]:
                            runID_dirs.append(directory)
                except:
                    # The directory name could not be converted to a runID
                    continue
        # Check if the directory name matches the format MXXXXX or XXXXX
        elif re.match(r'M?\d{5}$', dir_name):
            chip_dirs.append(directory)

        for root, dirs, files in os.walk(directory):
            for dir in dirs:
                # Check if the directory name matches the format YYMMDD
                if re.match(r'\d{6}$', dir):
                    try:
                        # Try to convert the directory name to a date
                        datetime.strptime(dir, '%y%m%d')
                        date_dirs.append(os.path.join(root, dir))
                    except ValueError:
                         # The directory name could not be converted to a date
                        try: 
                            # Check if the directory name matches the format of a runID (6 digit numbers)
                            re.match(r'\d{6}$', dir_name)
                            runID_dirs.append(directory)
                        except:
                            # The directory name could not be converted to a runID
                            continue
                # Check if the directory name matches the format MXXXXX or XXXXX
                elif re.match(r'M?\d{5}$', dir):
                    chip_dirs.append(os.path.join(root, dir))
                    
    for chip_dir in chip_dirs:
        #if chip_dir.startswith(date_dir):
        process_chip_dir(chip_dir)
    for runID_dir in runID_dirs:
        #if runID_dir.startswith(chip_dir):
        process_runID_dir(runID_dir)

    return sorted(set(include_dirs)), sorted(set(exclude_dirs))

def main(directories):
    #directories = select_folders()
    include_dirs, exclude_dirs = find_chip_ids(directories)
    maxone_dirs = []
    maxtwo_dirs = []
    for dir in include_dirs:
        # Extract the chip ID from the path to the .h5 file
        chip_id = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(dir))))
        # Check if the chip ID starts with "M"
        if chip_id.startswith("M"):
            maxtwo_dirs.append(dir)
        else:
            maxone_dirs.append(dir)

    # Print the MaxOne directories or None if the list is empty
    print("\nMaxOne directories:")
    if maxone_dirs:
        for dir in maxone_dirs:
            print(dir)
    else:
        print(None)

    # Print the MaxTwo directories or None if the list is empty
    print("\nMaxTwo directories:")
    if maxtwo_dirs:
        for dir in maxtwo_dirs:
            print(dir)
    else:
        print(None)

    # Print the excluded directories or None if the list is empty
    print("\nExcluded directories:")
    if exclude_dirs:
        for dir in exclude_dirs:
            print(dir)
    else:
        print(None)

    return maxone_dirs, maxtwo_dirs, exclude_dirs, include_dirs

if __name__ == "__main__":
    #for development only - make sure various levels of directories can be found
    mock_directories = ["/mnt/disk20tb/PrimaryNeuronData/Maxtwo/FolicAcid/FolicAcid/230921", 
                            "/mnt/disk20tb/PrimaryNeuronData/Maxtwo/FolicAcid/FolicAcid/231002", 
                            "/mnt/harddrive-2/ADNP/ADNP/230510/16821",
                            "/mnt/harddrive-2/CDKL5"]
    main(mock_directories)