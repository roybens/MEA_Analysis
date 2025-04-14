import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import spikeinterface.full as si
import spikeinterface.extractors as se
import h5py
import matplotlib.pyplot as plt
from helper_functions import detect_peaks
import pandas as pd
import numpy as np
from multiprocessing import Pool
from scipy.interpolate import CubicSpline
from stim_helper_functions import filter_spikes_by_isi
from sklearn.manifold import TSNE


class StimulationAnalysis:
    """
    A class to analyze neuronal responses to electrical stimulation in MEA recordings.

    Attributes
    ----------
    file_path : str
        Path to the raw .h5 file from the MEA recording.
    stim_frequency : float
        Frequency (Hz) of stimulation during the experiment.
    recording_electrode : int
        Electrode ID used for recording spikes.
    stim_electrode : int
        Electrode ID used for stimulation.
    artifact_electrode : int, optional
        Electrode ID used for measuring visible artifacts (if applicable).
    spike_threshold : float, optional
        Amplitude threshold for spike detection on the recording electrode (default: 9).
    peak_sign : str, optional
        Sign of the peaks to detect: "neg", "pos", or "both" (default: "neg").

    Usage
    -----
    analysis = StimulationAnalysis(
        file_path='path/to/file.raw.h5',
        stim_frequency=4,
        recording_electrode=73,
        stim_electrode=141,
        spike_threshold=9
    )

    analysis.run_full_analysis(stim_start=60, stim_length=120)

    Notes
    -----
    After initializing, call `run_full_analysis()` to generate summary plots
    and perform preliminary spike analysis.
    Additional methods support artifact filtering, waveform visualization,
    and ISI/fano factor analysis.
    """

    def __init__(self, file_path, stim_frequency, recording_electrode, stim_electrode, artifact_electrode=None, spike_threshold=9, peak_sign="neg"):
        self.visible_artifact = True # whether artifact is visible on graph
        self.spike_threshold = spike_threshold # threshold for spike detection
        self.file_path = file_path
        self.rec_electrode = recording_electrode
        self.stim_electrode = stim_electrode
        self.artifact_electrode = artifact_electrode
        self.rec_channel = self.get_channel_id(recording_electrode)
        self.stim_channel = self.get_channel_id(stim_electrode)
        self.peak_sign = peak_sign # "neg" or "pos" or "both" for spike detection
        self.stim_start = None

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

        self.peak_counts_df = self.get_spike_counts()


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
                    
                    # convert to seconds
                    spike_times_sec = x / self.fs
            
                    # apply ISI filtering
                    filtered_times = filter_spikes_by_isi(spike_times_sec, isi_threshold=0.01)
                    
                    # convert back to samples
                    filtered_samples = (filtered_times * self.fs).astype(int)
                    
                    batch_peaks[f'Channel {channel}'] = list(filtered_samples)

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
    

    def get_spike_counts_in_range2(self, channel_type, start_time=0, end_time=None):
        # gets the spikes in certain time range, returned as a list of sample indeces
        # channel_type: 'recording' or 'stim'
        if end_time is None:
            end_time = self.num_samples / self.fs

        if self.visible_artifact == True:
            _, peaks = self.filter_artifacts()
        else:
            peaks = self.peak_counts_df

        if channel_type == 'recording':
            channel = f'Channel {self.rec_channel}'
        elif channel_type == 'stim':
            channel = f'Channel {self.stim_channel}'
        else:
            raise ValueError(f"Parameter channel_type must either be \'recording\' or \'stim\'. It's currently set as {channel_type}")

        row_duration_samples = int(self.fs * 10)

        spikes = []
    
        for index, row in peaks.iterrows():
            row_offset = index * row_duration_samples

            curr_spikes = [spike + row_offset for spike in row[channel]]

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
                if any(abs(value - artifact) <= 10 for artifact in stim_electrode_peaks)
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
        print(f"Num stims {(self.stim_length * self.stim_freq - 1)}")  
        avg_spikes_per_stim = evoked_peaks / (self.stim_length * self.stim_freq - 1)  

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

        row_duration_samples = int(self.fs * 10)

        spikes = []
        for index, row in peak_counts.iterrows():
            row_offset = index * row_duration_samples

            row_spikes = [spike + row_offset for spike in row[f'Channel {self.rec_channel}']]
        
            spikes.extend(row_spikes)
        print(len(spikes))
        # Calculate spike counts for each phase
        pre_stim_spikes = [spike for spike in spikes if spike < (self.fs * self.stim_start)]
        stim_spikes = [spike for spike in spikes if spike >= (self.fs * self.stim_start) and spike < (self.fs * (self.stim_start + self.stim_length))]
        post_stim_spikes = [spike for spike in spikes if spike >= (self.fs * (self.stim_start + self.stim_length))]
        
        
        print(f"Pre-stim total spike count: {len(pre_stim_spikes)}") 
        print(f"Stim total spike count: {len(stim_spikes)}") 
        print(f"Post-stim spike count: {len(post_stim_spikes)}") 

        plt.figure(figsize=(10,6))
        plt.plot(time_values, spike_counts, marker='o', linestyle='-', color='b')

        plt.axvline(x=self.stim_start, color='g', linestyle='--')
        plt.axvline(x=(self.stim_start + self.stim_length), color='g', linestyle='--')


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
            time_range = self.total_recording

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

        rec_peaks = self.get_spike_counts_in_range2('recording', start_time, end_time)
        rec_peaks_times = [(peak / self.fs) for peak in rec_peaks]

        stim_peaks = self.get_spike_counts_in_range2('stim', start_time, end_time)
        stim_peaks_times = [(peak / self.fs) for peak in stim_peaks]

        # plot artifact electrode trace if it exists
        if self.artifact_electrode is not None:
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(20,4))
        else:
            fig, (ax1, ax3) = plt.subplots(2, 1, sharex=True, figsize=(20,4))


        # plot stim electrode trace
        ax1.plot(trace_length, stim_data[:,1], 'b-')
        ax1.set_ylabel('Stim ')
        ax1.set_title(f'Trial {trial_no} Traces')

        # add stim lines to stim electrode trace
        for t in stim_peaks_times:
            ax1.axvline(x=t, color='r', linestyle='--', linewidth=1)

        # plot recording electrode trace
        ax3.plot(trace_length, stim_data[:,0], 'b-')
        ax3.set_ylabel('Recording')
        ax3.set_xlabel('Time (s)')

        # add stim lines to recording trace (identifies artifacts)
        for t in rec_peaks_times:
            ax3.axvline(x=t, color='r', linestyle='--', linewidth=1)

        if self.artifact_electrode is not None:
            ax2.plot(trace_length, stim_data[:,2], 'b-')
            ax2.set_ylabel('Artifact ')

        plt.xlim(start_at, start_at + time_range)

        plt.show()

        print(f"Start time: {start_time}")
        print(f"End Time: {end_time}")
        return start_time, end_time

    def filter_artifact(self, trace, artifact_start_idx, artifact_end_idx):
        # Filter artifact region in given trace w/ interpolated values
        # trace - np.array
        # returns filtered trace in np.array

        artifact_start_idx = max(1, artifact_start_idx)
        artifact_end_idx = min(len(trace) - 2, artifact_end_idx)

        # Indices just outside the artifact region
        x = [artifact_start_idx - 1, artifact_end_idx + 1]
        y = [trace[artifact_start_idx - 1], trace[artifact_end_idx + 1]]

        artifact_indices = np.arange(artifact_start_idx, artifact_end_idx + 1)

        # fill with interpolated values
        spline = CubicSpline(x, y, bc_type='natural')
        trace[artifact_indices] = spline(artifact_indices)

        return trace

    def overlap_stim_responses(self, time_range, electrode_type='recording', filter_artifact=True):
        # overlaps immediate effect of all stims on graph
        # electrode_type: 'stim' or 'recording'
        # time_range: time in seconds before and after stim to show on x-axis
        # filter_artifact (boolean): whether to apply artifact filter 

        start = self.stim_start if self.stim_start is not None else self.pre_stim_length
        end = start + self.stim_length

        stim_spikes = self.get_spike_counts_in_range2('stim', start, end) 
        stim_times = [(spike / self.fs) for spike in stim_spikes]

        recording = self.recording_bp

        samples_per_time_range = int(time_range * self.fs)

        all_traces = []
        

        plt.figure(figsize=(10, 6))
        plt.title(f"{electrode_type} Channel Overlapped Artifact - {time_range}s Window")
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude (mV)")

        # Loop through each stim
        for stim_time in stim_times:
            # Conversion from time to samples
            stim_sample = int(stim_time * self.fs)

            channel = [str(self.rec_channel)] if electrode_type == 'recording' else [str(self.stim_channel)]

            # Immediate post-stim recording trace
            trace = recording.get_traces(
                channel_ids=channel,
                start_frame=stim_sample - samples_per_time_range,
                end_frame=stim_sample + samples_per_time_range
            ).flatten()

            # Interpolate the artifact region if filtering is enabled
            if filter_artifact:
                artifact_start_idx = int((stim_time - 0.002) * self.fs - (stim_sample - samples_per_time_range))
                artifact_end_idx = int((stim_time + 0.002) * self.fs - (stim_sample - samples_per_time_range))
                # trace = stim_helper.polynomial_interpolation(trace, artifact_start_idx, artifact_end_idx, degree=2)
                trace = stim_helper.bandpass_filter(trace, self.fs, 300, 800)
            all_traces.append(trace)

            time_axis = np.linspace(-time_range, time_range, len(trace))
            plt.plot(time_axis, trace, alpha=0.5)

        plt.show()

        all_traces = np.array(all_traces)
        mean_trace = np.mean(all_traces, axis=0)

        # Interpolate the artifact region for the mean trace
        if filter_artifact:
            artifact_start_idx = int(samples_per_time_range - 0.002 * self.fs)
            artifact_end_idx = int(samples_per_time_range + 0.002 * self.fs)
            # mean_trace = stim_helper.polynomial_interpolation(mean_trace, artifact_start_idx, artifact_end_idx, degree=2)
            mean_trace = stim_helper.bandpass_filter(mean_trace, self.fs, 300, 800)                

        # Plot the mean trace
        plt.figure(figsize=(10, 6))
        time_axis = np.linspace(-time_range, time_range, len(mean_trace))
        plt.plot(time_axis, mean_trace, label="Average Trace", color="black", linewidth=2)

        plt.axvline(x=0, color='r', linestyle='--', label="Stim Start")
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude (mV)")
        plt.title(f"{electrode_type} Channel Mean Artifact Trace - {time_range}s Window")
        plt.legend()
        plt.show()

        return mean_trace

    def isi(self):
        # calculates the inter spike interval (ISI) for each phase of the Stim Assay
        # returns dictionary containg mean ISI and std for each phase
    
        if self.visible_artifact == True: 
            _, peaks = self.filter_artifacts()
        else:
            peaks = self.peak_counts_df

        row_duration_samples = int(self.fs * 10)

        spikes = []
        for index, row in peaks.iterrows():
            row_offset = index * row_duration_samples

            row_spikes = [spike + row_offset for spike in row[f'Channel {self.rec_channel}']]
        
            spikes.extend(row_spikes)

        pre_stim_spikes = [spike for spike in spikes if spike < (self.fs * self.stim_start)]
        stim_spikes = [spike for spike in spikes if spike >= (self.fs * self.stim_start) and spike < (self.fs * (self.stim_start + self.stim_length))]
        post_stim_spikes = [spike for spike in spikes if spike >= (self.fs * (self.stim_start + self.stim_length))]

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
    
    def extract_spike_waveforms(self, window_size=50):
        """
        Extracts spike waveforms from a continuous recording trace.
        Parameter - window_size: Total number of samples per waveform.
        Returns - waveforms: dictionary that maps spike time to 1d array waveform for that spike
        """

        # get entire recording trace
        recording_trace = self.recording_bp.get_traces(channel_ids=[str(self.rec_channel)])

        # get spike indices
        spike_indices = []
        row_duration_samples = int(self.fs * 10)  

        for index, row in self.peak_counts_df.iterrows():
            row_offset = index * row_duration_samples 
            spikes = [spike + row_offset for spike in row[f'Channel {self.rec_channel}']] 
            spike_indices.extend(spikes)  


        half_window = window_size // 2
        waveforms = {}

        for idx in spike_indices:
            if idx - half_window < 0 or idx + half_window >= len(recording_trace):
                continue  # Skip spikes near the edges

            waveform = recording_trace[idx - half_window : idx + half_window].flatten()
            spike_time = idx / self.fs
            waveforms[spike_time] = waveform

        return waveforms
    
    def run_full_analysis(self, stim_start: float, stim_length: float, filter_artifact: bool = True, plot_traces: bool = True, plot_spike_counts: bool = True):
        """
        Runs the full stimulation analysis pipeline.

        Parameters:
        -----------
        stim_start : float
            Time (in seconds) when stimulation starts.
        stim_length : float
            Duration (in seconds) of the stimulation period.
        plot_traces : bool
            If True, plots stim and recording traces.
        plot_spike_counts : bool
            If True, plots spike counts across time.
        """

        self.stim_start = stim_start
        self.stim_length = stim_length
        self.pre_stim_length = stim_start
        self.post_stim_length = self.total_recording - (stim_start + stim_length)
        self.visible_artifact = filter_artifact

        self.peak_counts_df = self.get_spike_counts()

        if plot_traces:
            self.plot_stim_traces(trial_no=1, time_range=(8 / self.stim_freq), start_at=stim_start)

        if plot_spike_counts:
            self.plot_spike_counts(electrode_type='recording', trial_no=1)

        # Uncomment these for deeper metrics:
        isi = self.isi()
        isi_mean_std = self.calc_mean_isi(isi)
        fanofactors = self.calculate_fano_factor(isi_mean_std)
        print(f"Mean and std ISI: {isi_mean_std}")
        print(f"Fano factors: {fanofactors}")

    
    def plot_spike_types_on_trace(self, time_window=None):
        """
        Plot segments of the recording trace with color-coded spike types based on t-SNE clusters.
        
        Parameters:
        -----------
        time_window : tuple, optional
            (start_time, end_time) in seconds to plot. If None, plots first 5 seconds
        """
        # Set default time window if none provided
        if time_window is None:
            time_window = (0, 5)  # Plot first 5 seconds
        
        # Convert time window to samples
        start_sample = int(time_window[0] * self.fs)
        end_sample = int(time_window[1] * self.fs)
        
        # Get the trace for the recording electrode
        channel_index = np.where(self.all_channel_ids == str(self.rec_channel))[0][0]
        trace = self.recording_bp.get_traces(segment_index=0, 
                                           channel_ids=self.all_channel_ids[channel_index],
                                           start_frame=start_sample,
                                           end_frame=end_sample)
        
        # Define regions as before
        regions = {
            'Main Cluster': (self.tsne_results[:, 0] < 20) & (self.tsne_results[:, 0] > -40) & (self.tsne_results[:, 1] < 30),
            'Vertical Chain': (self.tsne_results[:, 0] > 30) & (self.tsne_results[:, 0] < 50),
            'Top Cluster': (self.tsne_results[:, 1] > 30),
            'Far Right Outliers': (self.tsne_results[:, 0] > 50)
        }
        
        # Create a color map for different spike types
        colors = {'Main Cluster': 'blue', 
                  'Vertical Chain': 'red',
                  'Top Cluster': 'green',
                  'Far Right Outliers': 'purple'}
        
        # Convert spike times from dictionary keys to array
        spike_times = np.array(list(self.spike_waveforms_dict.keys()))
        
        # Create the plot
        plt.figure(figsize=(15, 8))
        
        # Plot the raw trace
        time_axis = np.arange(start_sample, end_sample) / self.fs
        plt.plot(time_axis, trace.flatten(), 'k', alpha=0.5, label='Raw Trace')
        
        # Plot colored markers for each spike type
        for region_name, region_mask in regions.items():
            # Find spikes that belong to this cluster
            region_spikes = spike_times[region_mask]
            
            # Filter spikes within time window
            mask = (region_spikes >= time_window[0]) & (region_spikes <= time_window[1])
            region_spikes = region_spikes[mask]
            
            # Plot markers at spike times
            if len(region_spikes) > 0:
                spike_amplitudes = [trace[int((t - time_window[0]) * self.fs)] for t in region_spikes]
                plt.scatter(region_spikes, spike_amplitudes, 
                           color=colors[region_name], 
                           label=region_name,
                           marker='|',
                           s=100)
        
        plt.xlabel('Time (seconds)')
        plt.ylabel('Amplitude')
        plt.title(f'Spike Types on Raw Trace (Recording Electrode {self.rec_electrode})')
        plt.legend()
        plt.grid(True)
        plt.show()



