import os
import re
from collections import defaultdict
from termcolor import colored


def group_directories(directories):
    # Initialize a dictionary to store the grouped directories
    grouped_dirs = defaultdict(list)

    # Function to extract chip_id and date from a directory
    def get_chip_id_and_date(directory):
        chip_id = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(directory))))
        date = re.search(r'\d{6}', directory)
        return chip_id, date.group(0) if date else None

    # Group directories
    for dir in directories:
        chip_id, date = get_chip_id_and_date(dir)
        grouped_dirs[(chip_id, date)].append(dir)

    return grouped_dirs

#check each directory at each index of grouped_dirs. Validate the following:
# Rule 1: At least one network scan or axon tracking scan
# Rule 2: No lone activity scans
def validate_directories(grouped_dirs):
    # Initialize a dictionary to store the validated directories
    validated_dirs = {}

    # Loop through each group of directories
    for (chip_id, date), dirs in grouped_dirs.items():
        # Initialize counters for network scans, activity scans, and axon tracking scans
        network_scans = 0
        activity_scans = 0
        axon_tracking_scans = 0

        # Count the number of network scans, activity scans, and axon tracking scans in the group
        for dir in dirs:
            if 'Network' in dir:
                network_scans += 1
            elif 'ActivityScan' in dir:
                activity_scans += 1
            elif 'AxonTracking' in dir:
                axon_tracking_scans += 1

        # Validate the group according to the rules
        # Rule 1: At least one network scan or axon tracking scan
        # Rule 2: No lone activity scans
        if (network_scans > 0 or axon_tracking_scans > 0) and not (activity_scans > 0 and network_scans == 0 and axon_tracking_scans == 0):
            validated_dirs[(chip_id, date)] = dirs

    return validated_dirs

def check_grouped_dirs(grouped_dirs):
    for key in grouped_dirs:
        if len(grouped_dirs[key]) == 1:
            if 'Network' in grouped_dirs[key][0]:
                print(colored(f"Warning: {key} has a lone Network scan.", 'yellow'))
            elif 'AxonTracking' in grouped_dirs[key][0]:
                print(colored(f"Warning: {key} has a lone AxonTracking scan.", 'yellow'))
            elif 'ActivityScan' in grouped_dirs[key][0]:
                print(colored(f"Warning: {key} has a lone Activity scan. These data should have been excluded by this point. Revise group_dirs.py as needed.", 'red'))
        else:
            print(f"{key} has {len(grouped_dirs[key])} scans.")

def main(directories):
    grouped_dirs = group_directories(directories)
    validated_dirs = validate_directories(grouped_dirs)
    check_grouped_dirs(validated_dirs)

    return validated_dirs

if __name__ == "__main__":
    directories = "something" # Replace with your list of directories
    main