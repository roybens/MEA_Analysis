import spikeinterface
import spikeinterface.full as si
import spikeinterface.extractors as se
import h5py
import matplotlib.pyplot as plt
import math

from helper_functions import detect_peaks
import pandas as pd
import numpy as np
import time
from multiprocessing import Pool

import helper_functions as helper

class StimulationAnalysis:
    def __init__(self, file_path, stim_frequency, recording_electrode, stim_electrode, artifact_electrode=None, spike_threshold=9, peak_sign="neg"):
        self.visible_artifact = False
        self.spike_threshold = spike_threshold
        self.file_path = file_path
        self.rec_electrode = recording_electrode
        self.stim_electrode = stim_electrode
        self.artifact_electrode = artifact_electrode
        self.rec_channel = self.get_channel_id(recording_electrode)
        self.stim_channel = self.get_channel_id(stim_electrode)
        self.peak_sign = peak_sign

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
        self.pre_stim_length = None
        self.stim_length = None
        self.post_stim_length = None
        self.stim_freq = stim_frequency 
    

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
                        threshold = self.spike_threshold
                    else:
                        threshold = 200

                    x, y = detect_peaks(traces[:, channel_index], peak_sign=self.peak_sign, abs_threshold=threshold)
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

        # time ranges
        time_ranges = [(start, min(start + batch_size, total_frames)) for start in range(0, total_frames, batch_size) if start + batch_size <= total_frames] 
        
        # multiprocessing
        with Pool(16) as pool: 
            results = pool.map(self.process_batch, time_ranges) 

        valid_results = [result for result in results if result is not None] 

        self.peak_counts_df = pd.DataFrame(valid_results) 

        return self.peak_counts_df
    
    def get_spike_counts_in_range(self, start_time=0, end_time=None):
        # returns df with spike counts for stim and rec electrodes within specific time range
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
        peaks_df = pd.DataFrame(valid_results)

        return peaks_df
    

    def get_spike_counts_in_range2(self, start_time=0, end_time=None):
        if end_time is None:
            end_time = self.num_samples / self.fs

        if self.visible_artifact == True:
            _, peaks = self.filter_artifacts()
        else:
            peaks = self.peak_counts_df

        row_duration_samples = int(self.fs * 10)

        spikes = []

        for index, row in peaks.iterrows():
            row_offset = index * row_duration_samples

            curr_spikes = [spike + row_offset for spike in row[f'Channel {self.rec_channel}']]

            spikes.extend(curr_spikes)

        return [spike for spike in spikes if ((spike >= start_time * self.fs) and (spike <= end_time * self.fs))]
        

    def filter_artifacts(self):
        total_artifacts = 0
        peak_counts = self.peak_counts_df.copy()
        for index, row in peak_counts.iterrows():
            rec_electrode_channel = f"Channel {self.channel_dict['Recording Channel']}"
            rec_electrode_peaks = row[rec_electrode_channel]

            stim_electrode_channel = f"Channel {self.channel_dict['Stim Channel']}"
            stim_electrode_peaks = row[stim_electrode_channel]

            rec_electrode_peaks = [int(value) for value in rec_electrode_peaks]
            stim_electrode_peaks = [int(value) for value in stim_electrode_peaks]
            
            artifact_indices = [
                index for index, value in enumerate(rec_electrode_peaks)
                if any(abs(value - artifact) <= 2 for artifact in stim_electrode_peaks)
            ]

            total_artifacts += len(artifact_indices)

            filtered_peaks = [
                value for idx, value in enumerate(rec_electrode_peaks) if idx not in artifact_indices
            ]

            peak_counts.at[index, rec_electrode_channel] = filtered_peaks

        
        return total_artifacts, peak_counts
    
    def avg_inter_stim_sc(self):
        # calculates the average no. of spikes between each stim, and average latency between stim

        if self.visible_artifact == True:
            _, peak_counts = self.filter_artifacts()
        else:
            peak_counts = self.peak_counts_df

        # print(peak_counts)

        num_stim_rows = 0
        num_stims = 0

        evoked_peaks = 0

        total_latency = 0
        n = 0

        for index, row in peak_counts.iterrows():
            rec_electrode_channel = f"Channel {self.channel_dict['Recording Channel']}"
            rec_electrode_peaks = row[rec_electrode_channel]

            stim_electrode_channel = f"Channel {self.channel_dict['Stim Channel']}"
            stim_electrode_peaks = row[stim_electrode_channel]

            stim_row = False
            last_stim_row = False
            
            # check if this row (time frame) has stims
            if len(stim_electrode_peaks) > 0:
                num_stim_rows += 1
                stim_row = True

            # check if this the last row with stims
            if index + 1 < len(peak_counts):
                next_row = peak_counts.iloc[index + 1]
                next_stim_electrode_peaks = next_row[stim_electrode_channel]
                if stim_row and len(next_stim_electrode_peaks) == 0:
                    last_stim_row = True

            if num_stim_rows == 1:
                relevant_peaks = [peak for peak in rec_electrode_peaks if peak > stim_electrode_peaks[0]]
            elif last_stim_row:
                relevant_peaks = [peak for peak in rec_electrode_peaks if peak < stim_electrode_peaks[len(stim_electrode_peaks) - 1]]
            elif stim_row:
                relevant_peaks = rec_electrode_peaks
            else:
                continue

            evoked_peaks += len(relevant_peaks)
            no_inter_stim_spikes = 0

            i = 0 # stim peaks index
            while i < (len(stim_electrode_peaks) - 2):
                inter_stim_spikes = False
                for j in range(len(rec_electrode_peaks)):
                    if rec_electrode_peaks[j] > stim_electrode_peaks[i] and rec_electrode_peaks[j] < stim_electrode_peaks[i+1]:
                        # finds the first action potential after a stim
                        latency = (rec_electrode_peaks[j] - stim_electrode_peaks[i]) / self.fs
                        total_latency += latency
                        n += 1
                        inter_stim_spikes = True
                        break
                
                if(not inter_stim_spikes):
                    no_inter_stim_spikes += 1

                i += 1

        print(f"Number of stims without spikes: {no_inter_stim_spikes}")
        avg_latency = total_latency / n
        print(f'evoked_peaks: {evoked_peaks}')
        print(f"Num stims {((self.stim_length / self.stim_freq) - 1)}")
        avg_spikes_per_stim = evoked_peaks / ((self.stim_length / self.stim_freq) - 1)

        return avg_latency, avg_spikes_per_stim


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

        # convert time values to seconds 
        time_values = pd.Series([int(time.split(' ')[0]) / self.fs for time in self.peak_counts_df['Time Range']]) 
        
        if electrode_type == 'recording' and self.visible_artifact == True:
            num_artifacts, peak_counts = self.filter_artifacts()
        else:
            peak_counts = self.peak_counts_df

        if electrode_type == 'recording':
            spike_counts = peak_counts[f'Channel {self.rec_channel}']
        elif electrode_type == 'stim':
            spike_counts = peak_counts[f'Channel {self.stim_channel}']
        elif electrode_type == 'artifact':
            spike_counts = peak_counts[f'Channel {self.artifact_channel}']
        else:
            raise ValueError("Invalid input parameter for plot_spike_counts(). electrode_type must be 'stim', 'recording', or 'artifact'" )


        spike_counts = spike_counts.apply(len)

        # Calculate spike counts for each phase
        pre_stim_sc = spike_counts[time_values < self.pre_stim_length].sum() 
        stim_sc = spike_counts[(time_values >= self.pre_stim_length) & (time_values <= (self.pre_stim_length + self.stim_length))].sum()
        post_stim_sc = spike_counts[time_values > (self.pre_stim_length + self.stim_length)].sum() 
        
        
        print(f"Pre-stim total spike count: {pre_stim_sc}") 
        print(f"Stim total spike count: {stim_sc}") 
        print(f"Post-stim spike count: {post_stim_sc}") 


        plt.figure(figsize=(10,6))
        plt.plot(time_values, spike_counts, marker='o', linestyle='-', color='b')

        plt.axvline(x=self.pre_stim_length, color='g', linestyle='--')
        plt.axvline(x=(self.pre_stim_length + self.stim_length), color='g', linestyle='--')


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

        plt.show()

        plt.figure()
        plt.plot(time_post_stim, post_stim_data)
        plt.xlabel('Time (S)')
        plt.ylabel('Amplitude')
        plt.title(f'{electrode_type.capitalize()} Electrode Post-Stim Trace Trial {trial_no}')
        plt.xlim(start_at, start_at + time_range)
        plt.show()
        

    def plot_stim_traces(self, trial_no, bp_filter=True, time_range=None, start_at=0, draw_stim_lines=False):
        if time_range is None:
            time_range = self.total_recording / 3

        if bp_filter:
            recording_data = self.recording_bp
        else:
            recording_data = self.recording

        start_time = start_at
        end_time = start_at + time_range

        if start_time < 0: 
            raise ValueError("Start_at must be greater than 0")
        if end_time > self.total_recording: 
            raise ValueError("End time out of range of recording length")

        channels = [str(self.rec_channel), str(self.stim_channel)]
        if self.artifact_electrode is not None:
            channels.append(str(self.artifact_channel))
        stim_data = recording_data.get_traces(
            channel_ids=channels, 
            start_frame=(start_time * self.fs), 
            end_frame=(end_time * self.fs)
        ) 

        trace_length = (np.arange(0, len(stim_data)) / self.fs) + start_at

        rec_peaks2 = self.get_spike_counts_in_range2(start_time, end_time)
        rec_peaks2_times = [(peak / self.fs) for peak in rec_peaks2]

        # plot artifact electrode trace if it exists
        if self.artifact_electrode is not None:
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(20,4))
        else:
            fig, (ax1, ax3) = plt.subplots(2, 1, sharex=True, figsize=(20,4))


        # plot stim electrode trace
        ax1.plot(trace_length, stim_data[:,1], 'b-')
        ax1.set_ylabel('Stim ')
        ax1.set_title(f'Trial {trial_no} Traces')

        # plot recording electrode trace
        ax3.plot(trace_length, stim_data[:,0], 'b-')
        ax3.set_ylabel('Recording')
        ax3.set_xlabel('Time (s)')

        # add stim lines
        for t in rec_peaks2_times:
            ax3.axvline(x=t, color='r', linestyle='--', linewidth=1)

        if self.artifact_electrode is not None:
            ax2.plot(trace_length, stim_data[:,2], 'b-')
            ax2.set_ylabel('Artifact ')

        plt.xlim(start_at, start_at + time_range)

        plt.show()

        print(f"Start time: {start_time}")
        print(f"End Time: {end_time}")
        return start_time, end_time


    def isi(self):
        # calculates the inter spike interval (ISI) for each phase of the Stim Assay
        # returns dictionary containg mean ISI and std for each phase
    
        if self.visible_artifact == True: 
            _, peaks = self.filter_artifacts()
        else:
            peaks = self.peak_counts_df


        pre_stim_spikes = []
        stim_spikes = []
        post_stim_spikes = []

        row_duration_samples = int(self.fs * 10)

        for index, row in peaks.iterrows():
            row_offset = index * row_duration_samples

            spikes = [spike + row_offset for spike in row[f'Channel {self.rec_channel}']]

            if 0 <= index < 12: # pre stim rows
                pre_stim_spikes.extend(spikes)
            elif 12 <= index < 24:    # Stim rows
                stim_spikes.extend(spikes)
            elif 24 <= index <= 35:    # Post-stim rows
                post_stim_spikes.extend(spikes)

        pre_stim_isis = np.diff(sorted(pre_stim_spikes)) / self.fs if len(pre_stim_spikes) > 1 else []
        stim_isis = np.diff(sorted(stim_spikes)) / self.fs if len(stim_spikes) > 1 else []
        post_stim_isis = np.diff(sorted(post_stim_spikes)) / self.fs if len(post_stim_spikes) > 1 else []

        self.plot_isi_distribution(pre_stim_isis, "Pre-Stim")
        self.plot_isi_distribution(stim_isis, "Stim")
        self.plot_isi_distribution(post_stim_isis, "Post-Stim")
        
        return (pre_stim_isis, stim_isis, post_stim_isis)
    
    def plot_isi_distribution(self, isis, stim_phase):
        # num_bins = round(math.sqrt(len(isis)))
        num_bins = 30
        plt.figure(figsize=(6, 4))
        plt.hist(isis, bins=num_bins, color='blue', alpha=0.7, edgecolor='black')
        plt.title(f"{stim_phase} ISI Distribution")
        plt.xlabel("ISI (s)")
        plt.ylabel("Frequency")
        plt.xlim(0, 6)
        plt.ylim(0, 400)
        plt.show()



    def calc_mean_isi(self, isis):

        mean_pre_isi = np.mean(isis[0]) if len(isis[0]) > 0 else np.nan
        std_pre_isi = np.std(isis[0]) if len(isis[0]) > 0 else np.nan

        mean_stim_isi = np.mean(isis[1]) if len(isis[1]) > 0 else np.nan
        std_stim_isi = np.std(isis[1]) if len(isis[1]) > 0 else np.nan

        mean_post_isi = np.mean(isis[2]) if len(isis[2]) > 0 else np.nan
        std_post_isi = np.std(isis[2]) if len(isis[2]) > 0 else np.nan


        return {
            "Pre-stim ISI": (mean_pre_isi, std_pre_isi),
            "Stim ISI": (mean_stim_isi, std_stim_isi),
            "Post-stim ISI": (mean_post_isi, std_post_isi)
        }
    
    def calculate_fano_factor(self, isi_stats):
        # calculates fanofactor of ISI for each phase of stim assay
        # parameter - isi_stats: df returned from self.isi()
            # FANOFACTOR closer to 1 indicates more random spiking (poission distr.)
            # FANOFACTOR less than 1 indicates more regular and consistent firing
            # FANOFACTOR greater than 1 indicates super-Poisson variabilty - more bursts / irregular firing
        # returns dictionary of fanofactor for each phase of Stim Assay
        fano_factors = {}

        for phase, (mean_isi, std_isi) in isi_stats.items():
            if mean_isi > 0:  # Ensure the mean is non-zero to avoid division by zero
                fano_factor = (std_isi ** 2) / mean_isi
            else:
                fano_factor = float('nan')  # Assign NaN if mean ISI is zero

            fano_factors[phase] = fano_factor

        return fano_factors
    
    def run_full_analysis(self):
        if self.stim_length == None:
            length = input("Please input length of stim period in seconds: ")
            self.stim_length = float(length)
            self.pre_stim_length = float(length)
            self.post_stim_length = float(length)
        
        self.get_spike_counts()


        self.plot_stim_traces(1, time_range=(8 * self.stim_freq), start_at=(1.5*self.pre_stim_length))

        vis_artifact = input("Is the artifact being detected? (y/n): ").strip().lower()

        if vis_artifact == 'y':
            self.visible_artifact = True
            print("Artifact detection set to visible (True).")
            self.plot_stim_traces(1, time_range=(8 * self.stim_freq), start_at=(1.5*self.pre_stim_length))
        elif vis_artifact == 'n':
            self.visible_artifact = False
        else:
            print("Invalid input! Please enter 'y' for yes or 'n' for no.")

        
        self.plot_spike_counts('recording', 1)

        isi = self.isi()
        isi_mean_std = self.calc_mean_isi(isi)
        fanofactors = self.calculate_fano_factor(isi_mean_std)
        print(f"Mean and std ISI: {isi_mean_std}")
        print(f"Fanofactors: {fanofactors}")
    

    