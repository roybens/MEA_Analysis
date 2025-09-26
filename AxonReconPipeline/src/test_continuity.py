import spikeinterface.extractors as se
import os
import re
from termcolor import colored

def check_MEA_data_continuity(file_path, verbose = False):
    # Initialize an array to store the check results
    check_array = []

    # Function to try reading the maxwell file and return True if successful, False otherwise
    def try_read_maxwell(rec_num, well_num=None, verbose=False):
        rec_name = 'rec' + rec_num
        well_name = 'well' + well_num if well_num else None
        try:
            if well_name:
                print(f"Checking {well_name}")
                print(f"\tChecking {rec_name}...", end="")
            else: print(f"Checking {rec_name}...", end="")
            recording = se.read_maxwell(file_path,rec_name=rec_name, stream_id=well_name)
            check_array.append(True)
            print("continuous")
            return True
        except Exception as e:
            # Handle different types of exceptions
            if "Unable to open object (object 'rec" in str(e) and "doesn't exist)" in str(e):
                print("doesn't exist")
                if verbose:
                    print(colored(e, 'yellow'))
                    print(colored(f"This error is expected. All recordings found. rec_name:{rec_name} does not exist.", 'yellow'))
            elif "Unable to open object (object 'routed' doesn't exist)" in str(e):
                print("non-continuous")
                check_array.append(False)
                if verbose:
                    print(colored(e, 'yellow'))
                    print(colored("This error indicates that 'RecordSpikesOnly' was active during this recording. Data are not continuous.", 'yellow'))
            elif "stream_id well006 is not in ['well000', 'well001', 'well002', 'well003', 'well004', 'well005']" in str(e):
                print("well006 doesn't exist")
                if verbose:
                    print(colored(e, 'yellow'))
                    print(colored("This error is expected. MaxTwo does not accomodate a 7th well. Stream ID well006 is not found in the list.", 'yellow'))
            else:
                print("unknown error")
                check_array.append(False)
                if verbose:
                    print(colored(e, 'red'))
                    print(colored("This error is unexpected. Please check the error message and try to resolve it.", 'red'))
            return False

    # Get the scan type and chip id
    scan_type = get_scan_type(file_path)
    chip_id = get_chip_id(file_path)

    print(f"Chip: {chip_id}")
    print(f"Scan type: {scan_type}")

    # Initialize the recording and well numbers
    recnumber = 0
    wellnumber = 0 if chip_id[0] == 'M' else None

    # Check the data continuity based on the scan type
    if scan_type == 'Network':
        # For Network scans, increment the well number until data is no longer readable
        while try_read_maxwell(str(recnumber).zfill(4), str(wellnumber).zfill(3) if wellnumber is not None else None, verbose=verbose):
            if wellnumber is not None:
                wellnumber += 1  # Increment the well number
            else:
                break  # If well number is None, break the loop
    elif scan_type == 'ActivityScan' or scan_type == 'AxonTracking':
        # For ActivityScan and AxonTracking scans, increment the well number first, then the recording number
        while try_read_maxwell(str(recnumber).zfill(4), str(wellnumber).zfill(3) if wellnumber is not None else None, verbose=verbose):
            if wellnumber is not None:
                wellnumber += 1  # Increment the well number
            else:
                recnumber += 1  # If well number is None, increment the recording number

    return check_array, len(check_array)

# Function to get the scan type from the directory name
def get_scan_type(directory):
    dir_name = os.path.basename(directory)
    if dir_name in ['Network', 'ActivityScan', 'AxonTracking']:
        return dir_name
    else:
        parent_dir = os.path.dirname(directory)
        if parent_dir == directory:  # We've reached the root directory
            return None
        return get_scan_type(parent_dir)

# Function to get the chip id from the directory name
def get_chip_id(directory):
    # Get the base name of the directory
    dir_name = os.path.basename(directory)
    # If the directory name matches the pattern 'M?\d{5}$', return the directory name
    if re.match(r'M?\d{5}$', dir_name):
        return dir_name
    else:
        # Get the parent directory
        parent_dir = os.path.dirname(directory)
        # If we've reached the root directory, return None
        if parent_dir == directory:
            return None
        # Recursively call the function with the parent directory
        return get_chip_id(parent_dir)
    
# Main function to check the data continuity for multiple file paths
def main(file_paths):
    # Function to get the chip id and date from the file path
    def get_chip_id_and_date(file_path):
        # Assuming the chip_id and date are part of the file_path
        # Modify the regex patterns as per your requirements
        chip_id = get_chip_id(file_path)
        date = re.search(r'\d{6}', file_path)
        return chip_id if chip_id else None, date.group(0) if date else None

    # Loop through each file path
    for i, file_path in enumerate(file_paths, start=1):
        # Print the file path being checked
        print(colored(f"\nChecking {file_path} for continuity...", 'yellow'))
        # Get the chip id and date from the file path
        chip_id, date = get_chip_id_and_date(file_path)
        print(f"Date: {date}")
        # Check the data continuity for the file path
        check_array, num_recordings = check_MEA_data_continuity(file_path, verbose = True)
        # If false is in check_array, print that data is not continuous and remove the file path from file_paths
        if False in check_array:
            print(f"Data for chip {chip_id} on date {date} is not continuous.")
            file_paths.remove(file_path)
        else:
            print(f"Data for chip {chip_id} on date {date} is continuous.")

    # Return the file paths that have continuous data
    # Print the directories with continuous data
    print("\nContinuous directories:")
    file_paths.sort()
    for dir in file_paths:
        print(dir)
    return file_paths

# Main execution
if __name__ == "__main__":
    #development only
    # Define the directories for MaxOne and MaxTwo
    MaxOne_dir = ["/mnt/ben-shalom_nas/irc/media/harddrive8tb/Adam/MEA_AxonTraceDevScans/230828/16862/ActivityScan/000014/data.raw.h5",
                  "/mnt/ben-shalom_nas/irc/media/harddrive8tb/Adam/MEA_AxonTraceDevScans/230828/16862/Network/000015/data.raw.h5"]
    MaxTwo_dir = ["/mnt/disk20tb/PrimaryNeuronData/Maxtwo/FolicAcid/FolicAcid/230921/M05506/ActivityScan/000007/data.raw.h5",
                    "/mnt/disk20tb/PrimaryNeuronData/Maxtwo/FolicAcid/FolicAcid/230921/M05506/Network/000008/data.raw.h5"]    
    
    # Check the data continuity for the directories and get the directories with continuous data
    continuous_MaxOne_dir = main(MaxOne_dir)
    continuous_MaxTwo_dir = main(MaxTwo_dir)

    # # Print the directories with continuous data
    # print("\nContinuous MaxOne directories:")
    # continuous_MaxOne_dir.sort()
    # for dir in continuous_MaxOne_dir:
    #     print(dir)
    # print("\nContinuous MaxTwo directories:")
    # continuous_MaxTwo_dir.sort()
    # for dir in continuous_MaxTwo_dir:
    #     print(dir)