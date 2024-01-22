#library imports
import os
import platform
import re

#local imports
from file_selector import main as file_selector_main
from extract_raw_assay_data import main as extract_raw_assay_data_main
from test_continuity import main as test_continuity_main
from mea_analysis_pipeline_aw import process_block as temp_mea_analysis_pipeline
from group_dirs import main as group_and_validate_directories
from test_continuity import get_scan_type
from test_continuity import get_chip_id
from axon_trace_classes import axon_trace_objects
from count_recordings_and_wells import count_recordings_and_wells
from reconstruct_axons import reconstruct_axon_trace_objects

def main(debug_mode, selected_folders=None):
     # Clear the terminal
    if platform.system() == "Windows":
        os.system('cls')
    else:
        os.system('clear')

    # Select the folders to process
    if debug_mode:
        selected_folders = file_selector_main(selected_folders, debug_mode=True)
    else:
        selected_folders = file_selector_main(selected_folders=None, debug_mode=False)

    #get the maxone and maxtwo raw assay data paths
    maxone_dirs, maxtwo_dirs, exclude_dirs, include_dirs = extract_raw_assay_data_main(selected_folders)

    #Test each raw data file for continuous data. Spike-only data will be excluded.
    cont_maxone_dirs = test_continuity_main(maxone_dirs)
    cont_maxtwo_dirs = test_continuity_main(maxtwo_dirs)

    #Pair up data. For every activity scan, there must be a network scan. Network scans without activity scans will be included, but flagged.
    grouped_dirs_maxone = group_and_validate_directories(cont_maxone_dirs)
    grouped_dirs_maxtwo = group_and_validate_directories(cont_maxtwo_dirs)
    
    #if not empty, run mea_analysis_pipeline_aw on each pair of activity and network scans
    
    #MaxOne Pipeline
    #Turn this off for now, since it's not working yet.
    run_maxone_analysis = False
    if grouped_dirs_maxone and run_maxone_analysis:
        print("Running MEA Analysis Pipeline on MaxOne data...")
        for key, dirs in grouped_dirs_maxone.items():
            # Now 'key' is the dictionary key and 'dirs' is the list of directories
            # You can process each group of directories here
            for dir in dirs:
                print(f"\tAnalyzing MEA, chip_id: {key[0]}, record date: {key[1]}")
                print(f"\tfile_path: [{dir}]")
                scan_type = get_scan_type(dir)
                _, _, recording_count  = count_recordings_and_wells(dir, scan_type)  # Ignore well_count for MaxOne
                chip_id = get_chip_id(dir)
                date = re.search(r'\d{6}', dir).group(0)
                channel_locations = {}
                waveforms = {}
                if scan_type in ["ActivityScan", "AxonTracking", "Network"]:
                    for j in range(recording_count):
                        print(f"\t\trecording {j}...")
                        recnumber = j
                        rec_name = 'rec' + str(recnumber).zfill(4)
                        file_path = dir
                        waveforms[rec_name], channel_locations[rec_name], recording_dir = temp_mea_analysis_pipeline(recnumber, scan_type, chip_id, date, file_path, clear_temp_files=False)
                    #axon_trace_precursors = axon_trace_objects(recording_dir, channel_locations)
                    #reconstruct_axon_trace_objects(axon_trace_precursors)
                    print()
                else:
                    print("Error: scan type not recognized.")
                    #return
    
    #MaxTwo Pipeline:
    if grouped_dirs_maxtwo:
        print("Running MEA Analysis Pipeline on MaxTwo data...")
        for key, dirs in grouped_dirs_maxtwo.items():
            # Now 'key' is the dictionary key and 'dirs' is the list of directories
            # process each group of directories here
            for dir in dirs:
                print(f"\tAnalyzing MEA, chip_id: {key[0]}, record date: {key[1]}")
                print(f"\tfile_path: [{dir}]")                
                scan_type = get_scan_type(dir)
                wells_per_chip, recordings_per_well, _ = count_recordings_and_wells(dir, scan_type)
                chip_id = get_chip_id(dir)
                date = re.search(r'\d{6}', dir).group(0)
                channel_locations = {}
                waveforms = {}
                recording_dir = {}
                for well_name in recordings_per_well.keys():
                    print(f"\t\t{well_name}...")
                    #debug
                    if well_name != "well000":
                        continue
                    #debug                    
                    channel_locations[well_name] = {}
                    waveforms[well_name] = {}
                    well_number = int(well_name[4:])  # Extract the well number from the well name
                    if scan_type in ["ActivityScan", "AxonTracking", "Network"]:
                        for j in range(recordings_per_well[well_name]):                            
                            print(f"\t\t\trecording {j}...")
                            recnumber = j
                            rec_name = 'rec' + str(recnumber).zfill(4)
                            file_path = dir
                            recording_dir[well_name] = temp_mea_analysis_pipeline(recnumber,scan_type, chip_id, date, file_path, clear_temp_files=False, well_number=well_number)                       
                    else:
                        print("Error: scan type not recognized.")
            for well_name in recordings_per_well.keys():
                #debug
                if well_name != "well000":
                    continue
                #debug                
                axon_trace_precursors = axon_trace_objects(recording_dir[well_name], dirs)
                reconstruct_axon_trace_objects(axon_trace_precursors)       
    print("done.")
    
if __name__ == "__main__":
    debug_mode = True
    #testing Full Activity Scan and Network Scan...
    selected_folders_test1 = [        
        "/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/ActivityScan",
        "/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/Network",
        ]
    #...vs. Axon Tracking Scan Alone
    selected_folders_test2 = [  
        "/mnt/ben-shalom_nas/rbs_maxtwo/rbsmaxtwo/media/rbs-maxtwo/harddisk20tb/Tests_Adam/Tests_Adam/240118/M06844/AxonTracking",
        ]  
    if debug_mode:
        main(debug_mode, selected_folders=selected_folders_test1)
        main(debug_mode, selected_folders=selected_folders_test2)
    else:
        main(debug_mode)