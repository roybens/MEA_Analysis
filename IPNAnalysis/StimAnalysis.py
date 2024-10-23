import spikeinterface
import spikeinterface.full as si
import spikeinterface.extractors as se
import h5py
import matplotlib.pyplot as plt

from helper_functions import detect_peaks
import pandas as pd
import numpy as np
import time
from multiprocessing import Pool

import helper_functions as helper

class StimulationAnalysis:
    def __init__(self, file_path, recording_electrode, stim_electrode, artifact_electrode=None):
        self.visible_artifact = False
        self.file_path = file_path
        self.rec_electrode = recording_electrode
        self.stim_electrode = stim_electrode
        self.artifact_electrode = artifact_electrode
        self.rec_channel = self.get_channel_id(recording_electrode)
        self.stim_channel = self.get_channel_id(stim_electrode)

        self.channel_dict = {'Stim Channel': self.stim_channel,
                             'Recording Channel': self.rec_channel
                             }
        self.electrode_dict = {'Stim Electrode': self.stim_electrode,
                               'Recording Electrode': self.rec_electrode}
        
        if self.artifact_electrode is not None:
            self.artifact_channel = self.get_channel_id(artifact_electrode)
            self.channel_dict['Artifact Channel'] = self.artifact_channel
            self.electrode_dict['Artifact Electrode'] = self.artifact_electrode
        
        
        self.recording = se.read_maxwell(self.file_path)
        self.recording_bp = si.bandpass_filter(self.recording, freq_min=300, freq_max=3000)

        self.all_channel_ids = self.recording.get_channel_ids()
        self.fs = self.recording.get_sampling_frequency()
        self.num_chan = self.recording.get_num_channels()
        self.num_seg = self.recording.get_num_segments()
        self.num_samples = self.recording.get_num_samples(segment_index=0)
        self.total_recording = self.recording.get_total_duration()
    

    def print_file_structure(self):
        with h5py.File(self.file_path, 'r') as h5file:
            def print_hdf5_structure(name, obj):
                print(name)

            h5file.visititems(print_hdf5_structure)

    def plot_neuron_print(self):
        fig, ax = plt.subplots(figsize=(8,8))
        si.plot_probe_map(self.recording_bp, ax=ax, with_channel_ids=False)
        ax.invert_yaxis()

        return

    def get_channel_id(self, electrode_id):
        # get the corresponding channel ID of an electrode 
    
        with h5py.File(self.file_path, 'r') as h5file:

            mapping = h5file['data_store/data0000/settings/mapping']

            return [entry for entry in mapping if entry[1] == electrode_id][0][0]
        
    def process_batch(self, time_range):
            start_frame, end_frame = time_range

            try:
                #print(f"Processing frames {start_frame} to {end_frame}")

                # Retrieve batch of data
                traces = self.recording_bp.get_traces(start_frame=start_frame, end_frame=end_frame, segment_index=0, return_scaled=False)

                # Collect peaks for each channel in this batch
                batch_peaks = {'Time Range': f"{start_frame} to {end_frame}"}
                batch_peak_counts = {'Time Range': f"{start_frame} to {end_frame}"}
                for channel_type,channel in self.channel_dict.items():
                    channel_indices = np.where(self.all_channel_ids == str(channel))[0]

                    channel_index = channel_indices[0]
                    if channel_type == 'Recording Channel':
                        threshold = 9
                    else:
                        threshold = 200

                    x, y = detect_peaks(traces[:, channel_index], peak_sign="neg", abs_threshold=threshold)
                    batch_peaks[f'Channel {channel}'] = list(x)

                return batch_peaks

            except ValueError as e:
                print(f"Skipping range {start_frame} to {end_frame} due to error: {e}")
                return None # Skip this batch if an error occurs


    def get_spike_counts(self):
        # Returns dataframe with spike counts for stim and recording electrodes

        start_time = 0
        end_time = self.num_samples / self.fs

        batch_size = int(self.fs * 10)
        total_frames = int((end_time - start_time)*self.fs)

        # Prepare time ranges 
        time_ranges = [(start, min(start + batch_size, total_frames)) for start in range(0, total_frames, batch_size) if start + batch_size <= total_frames] 
        
        # Use multiprocessing to process batches in parallel 
        with Pool(16) as pool: 
            results = pool.map(self.process_batch, time_ranges) 

        # Filter out None results and extract valid results 
        valid_results = [result for result in results if result is not None] 

        self.peak_counts_df = pd.DataFrame(valid_results) 

        return self.peak_counts_df
    
    def get_spike_counts_in_range(self, start_time=0, end_time=None):
        """
        Returns a DataFrame with spike counts for stim and recording electrodes within a specific time range.

        Parameters:
            start_time (float): The start time in seconds (default 0).
            end_time (float): The end time in seconds (default None, which means the entire duration).

        Returns:
            pd.DataFrame: A DataFrame with spike counts for the specified time range.
        """
        # If end_time is not provided, use the total recording duration
        if end_time is None:
            end_time = self.num_samples / self.fs

        # Convert time to frames
        start_frame = int(start_time * self.fs)
        end_frame = int(end_time * self.fs)

        # Define the batch size (10 seconds of data at a time)
        batch_size = int(self.fs * 10)
        total_frames = end_frame - start_frame

        # Prepare time ranges for multiprocessing
        time_ranges = [(start, min(start + batch_size, total_frames + start_frame)) 
                    for start in range(start_frame, end_frame, batch_size)]

        # Use multiprocessing to process batches in parallel
        with Pool(16) as pool:
            results = pool.map(self.process_batch, time_ranges)

        # Filter out None results and extract valid results
        valid_results = [result for result in results if result is not None]

        # Convert the valid results to a DataFrame
        self.peak_counts_df = pd.DataFrame(valid_results)

        return self.peak_counts_df



    def plot_spike_counts_bar_graph(self, electrode_type, trial_no):
        # electrode_type: 'stim', 'recording', 'artifact'

        # Convert 'Time Range' from samples to seconds by dividing by the sample frequency
        time_values = pd.Series([int(time.split(' ')[0]) / self.fs for time in self.peak_counts_df['Time Range']])

        total_time_range = max(time_values) - min(time_values)
        first_third = min(time_values) + total_time_range / 3
        second_third = min(time_values) + 2 * total_time_range / 3

        # Select the spike counts based on the electrode type
        if electrode_type == 'recording':
            spike_counts = self.peak_counts_df[f'Channel {self.rec_channel}']
        elif electrode_type == 'stim':
            spike_counts = self.peak_counts_df[f'Channel {self.stim_channel}']
        elif electrode_type == 'artifact':
            spike_counts = self.peak_counts_df[f'Channel {self.artifact_channel}']
        else: 
            raise ValueError("Invalid input parameter for plot_spike_counts(). electrode_type must be 'stim', 'recording', or 'artifact'" )        

        # Plot the bar graph for each batch (each bar corresponds to a batch)
        plt.figure(figsize=(10, 6))

        # Plot the bars with time_values as the x-axis and spike_counts as the y-axis
        plt.bar(time_values, spike_counts, width=(time_values.iloc[1] - time_values.iloc[0]), color='b')

        # Add labels and title
        plt.title(f'Spike Counts per Batch for Trial {trial_no} - {electrode_type.capitalize()} Electrode')
        plt.xlabel('Time (seconds)')
        plt.ylabel('Spike Count per Batch')
        plt.axvline(x=first_third, color='g', linestyle='--')
        plt.axvline(x=second_third, color='g', linestyle='--')

        # Display the plot
        plt.show()


    def plot_spike_counts(self, electrode_type, trial_no):
        # electrode_type: 'stim', 'recording', 'artifact'

        # convert time values to seconds by dividing by sample frequency
        time_values = pd.Series([int(time.split(' ')[0]) / self.fs for time in self.peak_counts_df['Time Range']]) 
            
        total_time_range = max(time_values) - min(time_values)
        first_third = min(time_values) + total_time_range / 3
        second_third = min(time_values) + 2 * total_time_range / 3

        if electrode_type == 'recording':
            spike_counts = self.peak_counts_df[f'Channel {self.rec_channel}']
        elif electrode_type == 'stim':
            spike_counts = self.peak_counts_df[f'Channel {self.stim_channel}']
        elif electrode_type == 'artifact':
            spike_counts = self.peak_counts_df[f'Channel {self.artifact_channel}']
        else:
            raise ValueError("Invalid input parameter for plot_spike_counts(). electrode_type must be 'stim', 'recording', or 'artifact'" )

        spike_counts = spike_counts.apply(len)

        # Calculate spike counts for each third 
        first_third_spike_count = spike_counts[time_values <= first_third].sum() 
        second_third_spike_count = spike_counts[(time_values > first_third) & (time_values <= second_third)].sum() 
        third_third_spike_count = spike_counts[time_values > second_third].sum() 
        
        
        
        print(f"Pre-stim total spike count: {first_third_spike_count}") 
        print(f"Stim total spike count: {second_third_spike_count}") 
        print(f"Post-stim spike count: {third_third_spike_count}") 


        plt.figure(figsize=(10,6))
        plt.plot(time_values, spike_counts, marker='o', linestyle='-', color='b')

        plt.axvline(x=first_third, color='g', linestyle='--')
        plt.axvline(x=second_third, color='g', linestyle='--')


        plt.title(f'Spike Counts Over Time - Trial {trial_no} {electrode_type.capitalize()} Electrode')
        plt.xlabel('Time (in samples)')
        plt.ylabel('Spike Count')
        plt.grid(True)

        plt.show()

    def plot_individual_traces(self, electrode_type, trial_no, bp_filter=True, time_range=None, start_at=0):
        # electrode_type: 'Stim', 'recording', 'artifact'
        # time_range: number of seconds to show on x-axis (smaller value is more zoomed in)
        #             --> default is entire stim phase (60 seconds)

        if time_range is None:
            time_range = self.total_recording / 3

        if electrode_type.lower() == 'stim':
            channel = str(self.stim_channel)
        elif electrode_type.lower() == 'recording':
            channel = str(self.rec_channel)
        elif electrode_type.lower() == 'artifact':
            channel = str(self.artifact_channel)
        else:
            raise ValueError()
    
        if bp_filter:
            recording_data = self.recording_bp
        else:
            recording_data = self.recording

        chunk_duration = self.total_recording / 3

        samples_per_chunk = int(chunk_duration * self.fs)

        pre_stim_start = 0
        pre_stim_end = samples_per_chunk
        pre_stim_data = recording_data.get_traces(channel_ids=[channel], start_frame=pre_stim_start, end_frame=pre_stim_end)

        during_stim_start = pre_stim_end 
        during_stim_end = pre_stim_end + samples_per_chunk 
        during_stim_data = recording_data.get_traces(channel_ids=[channel], start_frame=during_stim_start, end_frame=during_stim_end) 

        post_stim_start = during_stim_end 
        post_stim_end = self.num_samples 
        post_stim_data = recording_data.get_traces(channel_ids=[channel], start_frame=post_stim_start, end_frame=post_stim_end) 

        time_pre_stim = np.arange(0, len(pre_stim_data)) / self.fs
        time_during_stim = np.arange(0, len(during_stim_data)) / self.fs
        time_post_stim = np.arange(0, len(post_stim_data)) / self.fs

        plt.figure()
        plt.plot(time_pre_stim, pre_stim_data)
        plt.xlabel('Time (S)')
        plt.ylabel('Amplitude')
        plt.title(f'{electrode_type.capitalize()} Electrode Pre-Stim Trace Trial {trial_no}')
        plt.xlim(start_at, start_at + time_range)
        plt.show()

        plt.figure()
        plt.plot(time_during_stim, during_stim_data)
        plt.xlabel('Time (S)')
        plt.ylabel('Amplitude')
        plt.title(f'{electrode_type.capitalize()} Electrode During-Stim Trace Trial {trial_no}')
        plt.xlim(start_at, start_at + time_range)
        
        # Add stim lines
        for x in np.arange(start_at, start_at + time_range, 0.25):
            plt.axvline(x=x, color='r', linestyle=':', linewidth = 1)

        plt.show()

        plt.figure()
        plt.plot(time_post_stim, post_stim_data)
        plt.xlabel('Time (S)')
        plt.ylabel('Amplitude')
        plt.title(f'{electrode_type.capitalize()} Electrode Post-Stim Trace Trial {trial_no}')
        plt.xlim(start_at, start_at + time_range)
        plt.show()
        

    def plot_stim_traces(self, trial_no, bp_filter=True, time_range=None, start_at=0):
        if time_range is None:
            time_range = self.total_recording / 3

        if bp_filter:
            recording_data = self.recording_bp
        else:
            recording_data = self.recording

        chunk_duration = self.total_recording / 3

        samples_per_chunk = int(chunk_duration * self.fs)

        pre_stim_start = 0
        pre_stim_end = samples_per_chunk

        stim_start = samples_per_chunk 
        stim_end = stim_start + samples_per_chunk 
        channels = [str(self.rec_channel), str(self.stim_channel)]
        if self.artifact_electrode is not None:
            channels.append(str(self.artifact_channel))
        stim_data = recording_data.get_traces(
            channel_ids=channels, 
            start_frame=stim_start, 
            end_frame=stim_end
        ) 

        time_during_stim = np.arange(0, len(stim_data)) / self.fs

        # start and end times of the graph
        start_time = stim_start / self.fs + start_at
        end_time = start_time + time_range

        peaks = self.get_spike_counts_in_range(start_time, end_time)

        print(type(peaks[f'Channel {self.stim_channel}']))
        print(peaks[f'Channel {self.stim_channel}'])

        stim_peaks = np.hstack(peaks[f'Channel {self.stim_channel}']) / self.fs + start_at
        rec_peaks = np.hstack(peaks[f'Channel {self.rec_channel}']) / self.fs + start_at
        


        # plot artifact electrode trace if it exists
        if self.artifact_electrode is not None:
            artifact_peaks = np.hstack(peaks[f'Channel {self.artifact_channel}']) / self.fs + start_at
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(20,4))
        else:
            fig, (ax1, ax3) = plt.subplots(2, 1, sharex=True, figsize=(20,4))


        # plot stim electrode trace
        ax1.plot(time_during_stim, stim_data[:,1], 'b-')
        ax1.set_ylabel('Stim ')
        ax1.set_title(f'Trial {trial_no} Traces')

        # plot recording electrode trace
        ax3.plot(time_during_stim, stim_data[:,0], 'b-')
        ax3.set_ylabel('Recording')
        ax3.set_xlabel('Time (s)')

        # add stim lines
        for t in rec_peaks:
            ax3.axvline(x=t, color='r', linestyle='--', linewidth=1)

        if self.artifact_electrode is not None:
            ax2.plot(time_during_stim, stim_data[:,2], 'b-')
            ax2.set_ylabel('Artifact ')

        plt.xlim(start_at, start_at + time_range)

        plt.show()

        print(f"Start time: {start_time}")
        print(f"End Time: {end_time}")
        return start_time, end_time
    
