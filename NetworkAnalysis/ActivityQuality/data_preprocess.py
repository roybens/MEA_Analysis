import os
import pandas as pd
import warnings

# Ignore specific warnings
warnings.filterwarnings('ignore', category=pd.errors.SettingWithCopyWarning)
warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)


def process_activity_recording_to_csv(base_path):
    """
    Process 'activity_record.xlsx' in the base_path to generate Compiled_ActivityScan.csv files.
    """
    excel_filename = 'activity_record.xlsx'
    excel_path = os.path.join(base_path, excel_filename)
    sheets = pd.read_excel(excel_path, sheet_name=None)

    for sheet_name, simple_check_df in sheets.items():
        active_area_start_idx = simple_check_df[simple_check_df.eq("P1W1").any(axis=1)].index.min()
        general_info_df = simple_check_df.iloc[:active_area_start_idx].copy()
        active_area_df = simple_check_df.iloc[active_area_start_idx + 1:].copy()
        active_area_df.columns = simple_check_df.iloc[active_area_start_idx]

        general_info_df['Mapping_Key'] = general_info_df.apply(lambda x: f'P{x["Plate #"]}W{x["Well #"]}', axis=1)
        well_plate_genotype_mapping = general_info_df.set_index('Mapping_Key')[['Plate ID', 'Genotype']]

        new_csv_data = []
        for column in active_area_df.columns[2:]:
            for index, row in active_area_df.iterrows():
                if column in well_plate_genotype_mapping.index:
                    plate_id = well_plate_genotype_mapping.loc[column, 'Plate ID']
                    genotype = well_plate_genotype_mapping.loc[column, 'Genotype']
                    well_number = int(column.split('W')[1])
                    plate_number = general_info_df[general_info_df['Mapping_Key'] == column]['Plate #'].values[0]
                    neuron_type = "WT" if "WT" in genotype else genotype
                    new_csv_data.append({
                        'DIV': row['DIV'],
                        'Chip_ID': plate_id,
                        'Well': well_number,
                        'Plate_ID': plate_number,
                        'NeuronType': neuron_type,
                        'Active_area': row[column]
                    })

        new_csv_df = pd.DataFrame(new_csv_data)
        new_csv_df['Well'] = pd.to_numeric(new_csv_df['Well'])
        new_csv_df.sort_values(by=['DIV', 'Plate_ID', 'Well'], ascending=[True, True, True], inplace=True)

        full_path = os.path.join(base_path, 'ActivityElectrode', sheet_name, 'Activity', 'Compiled_ActivityScan.csv')
        os.makedirs(os.path.dirname(full_path), exist_ok=True)
        new_csv_df.to_csv(full_path, index=False)

    print("Activity recording processing complete.")


def adjust_and_copy_csvs(base_path):
    """
    Adjust and copy CSV files from 'CSVs' to 'Final_CSVs' in the base_path.
    """
    source_dir = os.path.join(base_path, 'CSVs')
    destination_dir = os.path.join(base_path, 'Final_CSVs')
    os.makedirs(destination_dir, exist_ok=True)
    
    for folder in os.listdir(source_dir):
        folder_path = os.path.join(source_dir, folder)
        if os.path.isdir(folder_path) and not folder.startswith('.'):
            activity_csv_path = os.path.join(folder_path, 'Compiled_ActivityScan.csv')
            network_csv_path = os.path.join(folder_path, 'Compiled_Networks.csv')
            output_folder = os.path.join(destination_dir, folder)
            os.makedirs(output_folder, exist_ok=True)
            
            if os.path.exists(network_csv_path):
                df_network = pd.read_csv(network_csv_path)
                if 'NeuronType' in df_network.columns:
                    df_network['NeuronType'] = df_network['NeuronType'].str.strip().replace(regex={r'^.*WT.*$': 'WT'})
                if 'IBI' in df_network.columns:
                    df_network.rename(columns={'IBI': 'mean_IBI'}, inplace=True)
                if 'Burst_Peak' in df_network.columns:
                    df_network.rename(columns={'Burst_Peak': 'mean_Burst_Peak'}, inplace=True)
                if 'Spike_per_Burst' in df_network.columns:
                    df_network.rename(columns={'Spike_per_Burst': 'mean_Spike_per_Burst'}, inplace=True)
                if 'BurstDuration' in df_network.columns:
                    df_network.rename(columns={'BurstDuration': 'mean_BurstDuration'}, inplace=True)
                df_network.to_csv(os.path.join(output_folder, 'Compiled_Networks.csv'), index=False)
            
            if os.path.exists(activity_csv_path):
                df_activity = pd.read_csv(activity_csv_path)
                if 'NeuronType' in df_activity.columns:
                    df_activity['NeuronType'] = df_activity['NeuronType'].str.strip().replace(regex={r'^.*WT.*$': 'WT'})
                if 'Active_area' not in df_activity.columns and 'Active_Electrodes' in df_activity.columns:
                    df_activity.rename(columns={'Active_Electrodes': 'Active_area'}, inplace=True)
                df_activity.to_csv(os.path.join(output_folder, 'Compiled_ActivityScan.csv'), index=False)


def merge_activity_data_and_update_networks(base_path):
    homocheck_dir = os.path.join(base_path, 'Final_CSVs')
    quickcheck_dir = os.path.join(base_path, 'ActivityElectrode')
    homo_folders = os.listdir(homocheck_dir)
    quick_folders = os.listdir(os.path.join(quickcheck_dir))

    for folder in homo_folders:
        homo_path = os.path.join(homocheck_dir, folder)
        activity_csv_path = os.path.join(homo_path, 'Compiled_ActivityScan.csv')
        network_csv_path = os.path.join(homo_path, 'Compiled_Networks.csv')
        
        if os.path.exists(network_csv_path):
            df_network = pd.read_csv(network_csv_path)
            if 'NeuronType' in df_network.columns:
                df_network['NeuronType'] = df_network['NeuronType'].str.strip().replace(regex={r'^.*WT.*$': 'WT'})
            df_network.to_csv(network_csv_path, index=False)
        
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


def make_reffile(base_path):
    """
    Process 'reference_file.xlsx' in the base_path to generate 'Reffile.xlsx'.
    """
    source_excel_path = os.path.join(base_path, 'reference_file.xlsx')
    output_path = os.path.join(base_path, 'Reffile.xlsx')
    xls = pd.ExcelFile(source_excel_path)

    with pd.ExcelWriter(output_path) as writer:
        for sheet_name in xls.sheet_names:
            df = pd.read_excel(xls, sheet_name=sheet_name)
            str_columns = ['Div', 'Assay', 'Run #', 'Wells_Recorded', 'ID', 'Neuron Source']

            for col in str_columns:
                if col in df.columns:
                    df[col] = df[col].apply(lambda x: f'{int(x):d}' if pd.notnull(x) and isinstance(x, (int, float)) else x)
                    df[col] = df[col].astype(str).str.strip()

            if 'Div' in df.columns:
                df.rename(columns={'Div': 'DIV'}, inplace=True)

            df.to_excel(writer, sheet_name=sheet_name, index=False)


def merge_csvs_and_combine(base_path):
    """
    Merge CSV files based on 'Reffile.xlsx' and combine all data into 'final_combined_data.csv'.
    """
    ref_file_path = os.path.join(base_path, 'Reffile.xlsx')
    ref_data = pd.read_excel(ref_file_path, sheet_name=None)
    homocheck_dir = os.path.join(base_path, 'Final_CSVs')
    combined_dataframes = {}

    directories = [d for d in os.listdir(homocheck_dir) if os.path.isdir(os.path.join(homocheck_dir, d))]

    for directory in directories:
        folder_path = os.path.join(homocheck_dir, directory)
        activity_file = os.path.join(folder_path, 'Compiled_ActivityScan.csv')
        network_file = os.path.join(folder_path, 'Compiled_Networks.csv')

        if os.path.exists(activity_file) and os.path.exists(network_file):
            df_activity = pd.read_csv(activity_file)
            df_network = pd.read_csv(network_file)
            ref_sheet_name = directory

            if ref_sheet_name in ref_data:
                ref_df = ref_data[ref_sheet_name]
                valid_ref_df = ref_df[ref_df['Assay'].str.lower().isin(['network today', 'network', 'network today/best', 'sparse 7x', 'sparse7x'])]
                valid_ref_df['Wells_Recorded'] = valid_ref_df['Wells_Recorded'].astype(str).str.split(',')
                valid_ref_df = valid_ref_df.explode('Wells_Recorded').reset_index(drop=True)
                valid_ref_df[['Run #', 'DIV', 'Wells_Recorded', 'ID']] = valid_ref_df[['Run #', 'DIV', 'Wells_Recorded', 'ID']].astype(str).apply(lambda x: x.str.strip())
                df_network[['Run_ID', 'DIV', 'Well', 'Chip_ID']] = df_network[['Run_ID', 'DIV', 'Well', 'Chip_ID']].astype(str).apply(lambda x: x.str.strip())

                df_network = df_network[
                    df_network.set_index(['Run_ID', 'DIV', 'Well', 'Chip_ID']).index.isin(
                        valid_ref_df.set_index(['Run #', 'DIV', 'Wells_Recorded', 'ID']).index
                    )
                ]

                df_activity[['Run_ID', 'DIV', 'Well', 'Chip_ID']] = df_activity[['Run_ID', 'DIV', 'Well', 'Chip_ID']].astype(str).apply(lambda x: x.str.strip())
                df_activity = df_activity[
                    df_activity.set_index(['Run_ID', 'DIV', 'Well', 'Chip_ID']).index.isin(
                        valid_ref_df.set_index(['Run #', 'DIV', 'Wells_Recorded', 'ID']).index
                    )
                ]

            df_activity = df_activity.drop(['Run_ID', 'Time'], axis=1, errors='ignore')
            df_network = df_network.drop(['Run_ID', 'Time'], axis=1, errors='ignore')

            combined_df = pd.merge(
                df_activity, df_network, on=['DIV', 'Chip_ID', 'Well', 'NeuronType'], how='inner'
            ).drop_duplicates()
            combined_df.to_csv(os.path.join(folder_path, 'combined_data.csv'), index=False)
            combined_dataframes[directory] = combined_df
        else:
            print(f"Missing activity or network file in {directory}")

    all_data = pd.DataFrame()
    for directory, df in combined_dataframes.items():
        if 'Chip_ID' in df.columns:
            loc = df.columns.get_loc('Chip_ID') + 1
            df.insert(loc, 'Trial', directory)
        else:
            df['Trial'] = directory
        all_data = pd.concat([all_data, df], ignore_index=True)

    final_output_path = os.path.join(base_path, 'final_combined_data.csv')
    all_data.to_csv(final_output_path, index=False)
    print(f"Final combined data saved to {final_output_path}")
    print(f"Final combined data shape: {all_data.shape}")


def process_all_data(base_path):
    """
    Main function to process all data in the base_path.
    """
    process_activity_recording_to_csv(base_path)
    adjust_and_copy_csvs(base_path)
    merge_activity_data_and_update_networks(base_path)
    make_reffile(base_path)
    merge_csvs_and_combine(base_path)
