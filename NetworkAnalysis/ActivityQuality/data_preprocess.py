# data_processing.py

import os
import pandas as pd
import warnings

# Ignore specific warnings
warnings.filterwarnings('ignore', category=pd.errors.SettingWithCopyWarning)
warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)


def process_activity_recording_to_csv(base_path, excel_filename):
    """
    Process recording Excel sheet data to well-formatted Compiled_ActivityScan.csv files, extract the active electrode to 'Active_area' features for further processing.
    The recording Excel sheet should contain these information: 
    1. A table with general information (# of chips, Plate #, Well #, Plate ID, Genotype) 
    2. A table with active electrode data (Date, DIV, P1W1, P1W2, ..., P2W1, P2W2, ..., etc.)

    Input:
    - base_path (str): The base directory containing the Excel file.
    - excel_filename (str): The name of the Excel file to process.
    Output:
    - CSV files saved in the respective directories.
    """
    # Load the Excel file
    excel_path = os.path.join(base_path, excel_filename)
    sheets = pd.read_excel(excel_path, sheet_name=None)  # Load all sheets into a dictionary

    # Process each sheet
    for sheet_name, simple_check_df in sheets.items():
        # Identify the index where the second part starts by locating the header "P1W1"
        active_area_start_idx = simple_check_df[simple_check_df.eq("P1W1").any(axis=1)].index.min()

        # Extract the general information data and the active area data separately, making a copy to avoid SettingWithCopyWarning
        general_info_df = simple_check_df.iloc[:active_area_start_idx].copy()
        active_area_df = simple_check_df.iloc[active_area_start_idx + 1:].copy()  # skip the row with headers and make a copy
        active_area_df.columns = simple_check_df.iloc[active_area_start_idx]  # Set new header for active area data

        # Construct correct keys for the mapping
        general_info_df['Mapping_Key'] = general_info_df.apply(lambda x: f'P{x["Plate #"]}W{x["Well #"]}', axis=1)
        well_plate_genotype_mapping = general_info_df.set_index('Mapping_Key')[['Plate ID', 'Genotype']]

        # Create the new DataFrame structured as per the requirements
        new_csv_data = []
        for column in active_area_df.columns[2:]:  # skip 'Date' and 'DIV' columns
            for index, row in active_area_df.iterrows():
                if column in well_plate_genotype_mapping.index:
                    plate_id = well_plate_genotype_mapping.loc[column, 'Plate ID']
                    genotype = well_plate_genotype_mapping.loc[column, 'Genotype']
                    well_number = int(column.split('W')[1])  # Assuming well numbers are like 'W1', 'W2', etc.
                    plate_number = general_info_df[general_info_df['Mapping_Key'] == column]['Plate #'].values[0]  # Get the Plate #
                    # Standardize NeuronType for any genotype containing "WT"
                    neuron_type = "WT" if "WT" in genotype else genotype
                    # neuron_type = genotype
                    new_csv_data.append({
                        'DIV': row['DIV'],
                        'Chip_ID': plate_id,
                        'Well': well_number,
                        'Plate_ID': plate_number,  # Add Plate_ID from the general info
                        'NeuronType': neuron_type,
                        'Active_area': row[column]
                    })

        # Convert list of dictionaries into a DataFrame
        new_csv_df = pd.DataFrame(new_csv_data)

        # Sort the DataFrame by 'DIV', 'Plate_ID', and 'Well'
        new_csv_df['Well'] = pd.to_numeric(new_csv_df['Well'])  # Ensure 'Well' is an integer
        new_csv_df.sort_values(by=['DIV', 'Plate_ID', 'Well'], ascending=[True, True, True], inplace=True)

        # Define the full path for saving the file
        full_path = os.path.join(base_path, sheet_name, 'Activity', 'Compiled_ActivityScan.csv')
        os.makedirs(os.path.dirname(full_path), exist_ok=True)  # Create the directory if it doesn't exist

        # Save the DataFrame to a new CSV file
        new_csv_df.to_csv(full_path, index=False)

    print("Data processing complete. Files have been saved in their respective directories.")

def adjust_and_copy_csvs(source_dir, destination_dir):
    """
    Run for global WT quality check. Run this function before merging activity data and network data.
    Adjust the multiple version activity/network compiled CSV files (from different version of analysis pipeline) in the source directory and copy them to the destination directory.
    The adjustments include standardizing 'NeuronType' values to general 'WT', activity/network feature columns, and saving the modified CSVs.
    """
    # Ensure the destination directory exists
    os.makedirs(destination_dir, exist_ok=True)
    
    # Iterate over folders in the source directory
    for folder in os.listdir(source_dir):
        folder_path = os.path.join(source_dir, folder)
        if os.path.isdir(folder_path) and not folder.startswith('.'):
            # Define paths for activity and network CSVs
            activity_csv_path = os.path.join(folder_path, 'Compiled_ActivityScan.csv')
            network_csv_path = os.path.join(folder_path, 'Compiled_Networks.csv')
            output_folder = os.path.join(destination_dir, folder)
            os.makedirs(output_folder, exist_ok=True)  # Create corresponding folder in destination
            
            # Process Network CSV
            if os.path.exists(network_csv_path):
                df_network = pd.read_csv(network_csv_path)
                # Standardize 'NeuronType' values and strip spaces
                if 'NeuronType' in df_network.columns:
                    df_network['NeuronType'] = df_network['NeuronType'].str.strip().replace(regex={r'^.*WT.*$': 'WT'})
                # Rename 'IBI' column to 'mean_IBI' if it exists
                if 'IBI' in df_network.columns:
                    df_network.rename(columns={'IBI': 'mean_IBI'}, inplace=True)
                # Rename 'Burst_Peak' column to 'mean_Spike_per_Burst' if it exists
                if 'Burst_Peak' in df_network.columns:
                    df_network.rename(columns={'Burst_Peak': 'mean_Burst_Peak'}, inplace=True)
                # Rename 'Spike_per_Burst' column to 'mean_Spike_per_Burst' if it exists
                if 'Spike_per_Burst' in df_network.columns:
                    df_network.rename(columns={'Spike_per_Burst': 'mean_Spike_per_Burst'}, inplace=True)
                if 'BurstDuration' in df_network.columns:
                    df_network.rename(columns={'BurstDuration': 'mean_BurstDuration'}, inplace=True)
                # Save the modified CSV to the destination path
                df_network.to_csv(os.path.join(output_folder, 'Compiled_Networks.csv'), index=False)
            
            # Process Activity CSV
            if os.path.exists(activity_csv_path):
                df_activity = pd.read_csv(activity_csv_path)
                # Standardize 'NeuronType' values and strip spaces
                if 'NeuronType' in df_activity.columns:
                    df_activity['NeuronType'] = df_activity['NeuronType'].str.strip().replace(regex={r'^.*WT.*$': 'WT'})
                # Rename 'Active_Electrodes' to 'Active_area' if necessary
                if 'Active_area' not in df_activity.columns and 'Active_Electrodes' in df_activity.columns:
                    df_activity.rename(columns={'Active_Electrodes': 'Active_area'}, inplace=True)
                # Save the modified CSV to the destination path
                df_activity.to_csv(os.path.join(output_folder, 'Compiled_ActivityScan.csv'), index=False)

def merge_activity_data_and_update_networks(homocheck_dir, quickcheck_dir):
    homo_folders = os.listdir(homocheck_dir)
    quick_folders = os.listdir(quickcheck_dir)

    for folder in homo_folders:
        homo_path = os.path.join(homocheck_dir, folder)
        activity_csv_path = os.path.join(homo_path, 'Compiled_ActivityScan.csv')
        network_csv_path = os.path.join(homo_path, 'Compiled_Networks.csv')
        
        # Update Compiled_Networks.csv if it exists
        if os.path.exists(network_csv_path):
            df_network = pd.read_csv(network_csv_path)
            if 'NeuronType' in df_network.columns:
                df_network['NeuronType'] = df_network['NeuronType'].str.strip().replace(regex={r'^.*WT.*$': 'WT'})
            df_network.to_csv(network_csv_path, index=False)
        # 
        if os.path.exists(activity_csv_path):
            df_homo = pd.read_csv(activity_csv_path)
            if 'Active_area' not in df_homo.columns and folder in quick_folders:
                quick_path = os.path.join(quickcheck_dir, folder, 'Activity', 'Compiled_ActivityScan.csv')
                if os.path.exists(quick_path):
                    df_quick = pd.read_csv(quick_path)
                    if 'Active_area' in df_quick.columns:
                        merged_df = pd.merge(df_homo, df_quick[['Well', 'DIV', 'Chip_ID', 'Active_area']],
                                             on=['Well', 'DIV', 'Chip_ID'], how='left')
                        merged_df.to_csv(activity_csv_path, index=False)

def process_datasets(source_dir, destination_dir, quickcheck_dir):
    # Adjust and copy datasets
    adjust_and_copy_csvs(source_dir, destination_dir)

    # Merge additional data into the datasets
    merge_activity_data_and_update_networks(destination_dir, quickcheck_dir)



def make_reffile(file_path, output_path):
    """
    Run for global WT quality check. Run this function before merging activity data and network data.
    Process recording Excel sheet (non-formatted reference files) to well-formatted reference files.
    The recording Excel sheet should contain these information: 
    A table with recording information (typical reference sheet in pipeline) inclduing columns: 'Date', 'Div', 'Assay', 'Run #', 'Wells_Recorded', 'ID', 'Neuron Source'

    Input:
    - file_path (str): The path to the Excel file to process.
    - output_path (str): The path to save the processed reference file.
    Output:
    - Excel file saved with processed data.
    """

    # Load the Excel file with all sheets
    xls = pd.ExcelFile(file_path)

    # Create a writer object to write multiple sheets
    with pd.ExcelWriter(output_path) as writer:
        for sheet_name in xls.sheet_names:
            # Read each sheet, initially loading all columns normally
            df = pd.read_excel(xls, sheet_name=sheet_name)

            # Columns to convert to string
            str_columns = ['Div', 'Assay', 'Run #', 'Wells_Recorded', 'ID', 'Neuron Source']

            # Process each column that needs to be string
            for col in str_columns:
                if col in df.columns:  # Check if column exists in DataFrame
                    # Convert numeric values to string, and ensure no floating point representation
                    df[col] = df[col].apply(lambda x: f'{int(x):d}' if pd.notnull(x) and isinstance(x, (int, float)) else x)
                    # Strip any leading/trailing white space
                    df[col] = df[col].astype(str).str.strip()

            # Rename 'Div' column to 'DIV' if it exists
            if 'Div' in df.columns:
                df.rename(columns={'Div': 'DIV'}, inplace=True)

            # Save each processed DataFrame to a separate sheet in the same Excel file
            df.to_excel(writer, sheet_name=sheet_name, index=False)

########################################################################################################################################
            
def merge_csvs(base_path, ref_file_path):
    """
    Run for global WT quality check. 
    Merge activity and network CSV files for each directory in base_path,
    based on the reference Excel file at ref_file_path.

    Parameters:
    - base_path (str): The base directory containing the folders for processing.
    - ref_file_path (str): The path to the reference Excel file.

    Returns:
    - combined_dataframes (dict): A dictionary with directory names as keys and
      combined DataFrames as values.
    """
    # Load the reference file
    ref_data = pd.read_excel(ref_file_path, sheet_name=None)

    # Initialize a dictionary to store combined DataFrames
    combined_dataframes = {}

    # Loop over directories in the base_path
    directories = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]

    for directory in directories:
        folder_path = os.path.join(base_path, directory)
        activity_file = os.path.join(folder_path, 'Compiled_ActivityScan.csv')
        network_file = os.path.join(folder_path, 'Compiled_Networks.csv')

        if os.path.exists(activity_file) and os.path.exists(network_file):
            # Process activity file
            df_activity = pd.read_csv(activity_file)


            # Process network file
            df_network = pd.read_csv(network_file)
            ref_sheet_name = directory  # Assuming the sheet name corresponds to the directory

            if ref_sheet_name in ref_data:
                ref_df = ref_data[ref_sheet_name]
                # Filter rows based on 'Assay' for nerwork
                valid_ref_df = ref_df[ref_df['Assay'].str.lower().isin(['network today', 'network', 'network today/best', 'sparse 7x', 'sparse7x'])]

                # Convert "Wells Recorded" to string, split into individual wells, and explode into separate rows
                valid_ref_df['Wells_Recorded'] = valid_ref_df['Wells_Recorded'].astype(str).str.split(',')
                valid_ref_df = valid_ref_df.explode('Wells_Recorded').reset_index(drop=True)

                # Convert all relevant columns to strings and strip whitespace
                valid_ref_df[['Run #', 'DIV', 'Wells_Recorded', 'ID']] = (
                    valid_ref_df[['Run #', 'DIV', 'Wells_Recorded', 'ID']]
                    .astype(str)
                    .apply(lambda x: x.str.strip())
                )
                df_network[['Run_ID', 'DIV', 'Well', 'Chip_ID']] = (
                    df_network[['Run_ID', 'DIV', 'Well', 'Chip_ID']]
                    .astype(str)
                    .apply(lambda x: x.str.strip())
                )

                # Filter df_network using the adjusted mask
                df_network = df_network[
                    df_network.set_index(['Run_ID', 'DIV', 'Well', 'Chip_ID']).index.isin(
                        valid_ref_df.set_index(['Run #', 'DIV', 'Wells_Recorded', 'ID']).index
                    )
                ]

                df_activity[['Run_ID', 'DIV', 'Well', 'Chip_ID']] = (
                    df_activity[['Run_ID', 'DIV', 'Well', 'Chip_ID']]
                    .astype(str)
                    .apply(lambda x: x.str.strip())
                )
                
               # Filter df_activity using the adjusted mask
                df_activity = df_activity[
                    df_activity.set_index(['Run_ID', 'DIV', 'Well', 'Chip_ID']).index.isin(
                        valid_ref_df.set_index(['Run #', 'DIV', 'Wells_Recorded', 'ID']).index
                    )
                ]

            df_activity = df_activity.drop(['Run_ID', 'Time'], axis=1, errors='ignore')

            df_network = df_network.drop(['Run_ID', 'Time'], axis=1, errors='ignore')

            # Merge df_activity and df_network directly
            combined_df = pd.merge(
                df_activity, df_network, on=['DIV', 'Chip_ID', 'Well', 'NeuronType'], how='inner'
            )  # Inner join to keep only matching rows

            # Remove duplicate rows
            combined_df = combined_df.drop_duplicates()
            combined_df.to_csv(os.path.join(folder_path, 'combined_data.csv'), index=False)

            # Store the combined DataFrame in the dictionary
            combined_dataframes[directory] = combined_df
        else:
            print(f"Missing activity or network file in {directory}")

    print("Data processing completed for all directories")
    return combined_dataframes

def combine_all_data(combined_dataframes):
    """
    Combine all DataFrames from multiple directories into a single DataFrame.

    Parameters:
    - combined_dataframes (dict): A dictionary with directory names as keys and
      combined DataFrames as values.

    Returns:
    - all_data (DataFrame): The combined DataFrame.
    """
    # Initialize an empty DataFrame to store all combined data
    all_data = pd.DataFrame()

    for directory, df in combined_dataframes.items():
        # Add a new column 'Trial' with the folder name
        # Ensure 'Trial' is inserted right after 'Chip_ID'
        if 'Chip_ID' in df.columns:
            # Find index of 'Chip_ID' column
            loc = df.columns.get_loc('Chip_ID') + 1
            # Insert 'Trial' column right after 'Chip_ID'
            df.insert(loc, 'Trial', directory)
        else:
            df['Trial'] = directory  # Fallback if 'Chip_ID' is not in the columns

        # Append the DataFrame to the all_data DataFrame
        all_data = pd.concat([all_data, df], ignore_index=True)

    # Print the shape of the final DataFrame to confirm the number of rows and columns
    print(f"Final combined data shape: {all_data.shape}")
    return all_data
