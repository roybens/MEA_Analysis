
# Notes =======================================================================
'''
'''
# Imports =====================================================================
import numpy as np
from scipy.stats import norm
from scipy.signal import convolve, find_peaks
import matplotlib.pyplot as plt
from MEA_Analysis.IPNAnalysis.helper_functions import detect_bursts_statistics
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor
from scipy.ndimage import gaussian_filter1d
from scipy import signal
import traceback
from scipy.spatial.distance import mahalanobis
from MEA_Analysis.NeuronClassication.classify_neurons import classify_neurons_v2
from MEA_Analysis.NetworkAnalysis.awNetworkAnalysis.dtw import dtw_burst_analysis
import time
from .helper import indent_mode_on, indent_mode_off, indent_increase, indent_decrease
import spikeinterface.postprocessing as spost

# Functions ===================================================================
'''Newer Functions'''
def plot_raster_plot_v2(ax, spiking_data_by_unit, unit_types=None):
    """Plot a raster plot for spiking data."""
    
    # Calculate the average firing rate for each unit
    firing_rates = {}
    for gid in spiking_data_by_unit:
        spike_times = spiking_data_by_unit[gid]['spike_times']
        spike_times = [spike_times] if isinstance(spike_times, (int, float)) else spike_times
        if spike_times is None: firing_rate = 0
        else: firing_rate = len(spike_times) / (max(spike_times) - min(spike_times)) if len(spike_times) > 1 else 0
        firing_rates[gid] = firing_rate
        
    # Sort the units based on their average firing rates
    sorted_units = sorted(firing_rates, key=firing_rates.get)
    
    # Create a mapping from original gid to new y-axis position
    gid_to_ypos = {gid: pos for pos, gid in enumerate(sorted_units)}
    
    # Plot the units in the sorted order
    if unit_types is None:
        for gid in sorted_units:
            spike_times = spiking_data_by_unit[gid]['spike_times']
            spike_times = [spike_times] if isinstance(spike_times, (int, float)) else spike_times
            if spike_times is not None:
                ax.plot(spike_times, [gid_to_ypos[gid]] * len(spike_times), 'b.', markersize=2)
    else:
        # Define legend markers for excitatory (E) and inhibitory (I) units
        exc_marker, = ax.plot([], [], 'b.', markersize=2, label='Excitatory (E)')
        inh_marker, = ax.plot([], [], 'r.', markersize=2, label='Inhibitory (I)')

        for gid in sorted_units:
            spike_times = spiking_data_by_unit[gid]['spike_times']
            spike_times = [spike_times] if isinstance(spike_times, (int, float)) else spike_times
            if spike_times is not None:
                if unit_types[gid] == 'E':
                    ax.plot(spike_times, [gid_to_ypos[gid]] * len(spike_times), 'b.', markersize=2)
                elif unit_types[gid] == 'I':
                    ax.plot(spike_times, [gid_to_ypos[gid]] * len(spike_times), 'r.', markersize=2)

        # Add legend with proxy markers
        ax.legend(handles=[exc_marker, inh_marker])

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Unit ID (sorted by firing rate)')
    ax.set_title('Raster Plot')
    plt.tight_layout()
    return ax

def plot_network_bursting_v2(ax, bursting_ax, mega_ax=None, mode='normal'):
    """Plot network bursting activity with clear differentiation between lines."""
    
    if mode == 'normal':
        # Copy ax features to the new ax
        ax.set_xlim(bursting_ax.get_xlim())
        ax.set_ylim(bursting_ax.get_ylim())
        ax.set_ylabel('Firing Rate (Hz)')
        ax.set_xlabel('Time (s)')
        ax.set_title(bursting_ax.get_title())

        # Plot Line 1 (blue)
        ax.plot(
            bursting_ax.get_lines()[0].get_xdata(),
            bursting_ax.get_lines()[0].get_ydata(),
            color='blue',
            #label='Bursting'
        )
        
        # plot bursting peaks
        ax.plot(
            bursting_ax.get_lines()[1].get_xdata(),
            bursting_ax.get_lines()[1].get_ydata(),
            'o', color='red', 
            label='Bursts'
        )

        # # Mark peaks with vertical dotted lines
        # for x_peak in bursting_ax.get_lines()[1].get_xdata():
        #     ax.axvline(x=x_peak, linestyle='dotted', color='red', alpha=0.8)

        # If mega_ax is not None, plot filled area and peaks
        if mega_ax is not None:
            # Plot the filled area behind other lines
            ax.fill_between(
                mega_ax.get_lines()[0].get_xdata(),
                mega_ax.get_lines()[0].get_ydata(),
                color='red',
                alpha=0.2,  # Transparency for the filled area
                label='Mega Bursts'
            )
            
            #outline the filled aread, dotted line
            ax.plot(
                mega_ax.get_lines()[0].get_xdata(),
                mega_ax.get_lines()[0].get_ydata(),
                color='red',
                linestyle='dotted'
            )

            # Mark mega peaks with vertical dotted lines
            for x_peak in mega_ax.get_lines()[1].get_xdata():
                ax.axvline(x=x_peak, linestyle='dotted', 
                        color='red', 
                        #label='Mega Bursts'
                        #alpha=0.8
                        )
                
            # remove redundant legend items
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys())

        # Add a legend for clarity
        ax.legend()

        return ax
    elif mode == 'mega':
        # this is for the case where we're plotting mega bursts only
        # Copy ax features to the new ax
        ax.set_xlim(bursting_ax.get_xlim())
        ax.set_ylim(bursting_ax.get_ylim())
        ax.set_ylabel('Firing Rate (Hz)')
        ax.set_xlabel('Time (s)')
        ax.set_title(bursting_ax.get_title())
        
        # Plot Line 1 (red)
        ax.plot(
            bursting_ax.get_lines()[0].get_xdata(),
            bursting_ax.get_lines()[0].get_ydata(),
            color='red',
            #label='Mega Bursts'
        )
        
        # plot bursting peaks
        ax.plot(
            bursting_ax.get_lines()[1].get_xdata(),
            bursting_ax.get_lines()[1].get_ydata(),
            'o', 
            color='orange',
            label='Mega Bursts'
        )
        
        # Mark peaks with vertical dotted lines
        for x_peak in bursting_ax.get_lines()[1].get_xdata():
            ax.axvline(x=x_peak, linestyle='dotted', color='red', alpha=0.8)
            
        # Add a legend for clarity
        ax.legend()
        
        return ax

def plot_network_summary_v2(network_data, **kwargs):
        
    # unpack kwargs
    #sim_data_path = kwargs.get('sim_data_path', None)
    #save_dir = sim_data_path.replace('_data.pkl', '').replace('_data.json', '')     # remove .pkl or .json from sim_data_path # HACK, just specify this ahead of time in the run script.
    limit_seconds = kwargs.get('limit_seconds', None)
    
    # unpack network data
    sim_data_path = network_data.get('sim_data_path', None)
    sim_data_basename = os.path.basename(sim_data_path)
    basename = sim_data_basename.replace('_data.pkl', '').replace('_data.json', '')
    output_dir = kwargs.get('output_dir', None)
    save_dir = os.path.join(output_dir, basename)
    #if sim_data_path is None: raise ValueError('No sim_data_path found in network_data')
    #save_dir = sim_data_path.replace('_data.pkl', '').replace('_data.json', '')
    #save_dir = kwargs.get('output_dir', None)
    bursting_ax = network_data['bursting_data']['ax']
    mega_bursting_ax = network_data['mega_bursting_data']['ax'] 
    spiking_data_by_unit = network_data['spiking_data']['spiking_metrics_by_unit'] 
    unit_types = network_data['unit_types']
    
    # prep save paths
    if not os.path.exists(save_dir): os.makedirs(save_dir, exist_ok=True)
    network_summary_mode = kwargs.get('network_summary_mode', '2p') # TODO: spec in run script    
    save_path_png_2p = os.path.join(save_dir, 'network_summary_2p.png')
    save_path_pdf_2p = os.path.join(save_dir, 'network_summary_2p.pdf')
    save_path_png_3p = os.path.join(save_dir, 'network_summary_3p.png')
    save_path_pdf_3p = os.path.join(save_dir, 'network_summary_3p.pdf')  
    
    # choose mode:
    mode = network_summary_mode
    if mode == '2p':
        
        # init plot
        fig, ax = plt.subplots(2, 1, figsize=(16, 9))
        
        # plot raster plot
        print("Generating raster plot...")

        ax[0] = plot_raster_plot_v2(ax[0], spiking_data_by_unit, unit_types)

        # plot network bursting plots
        print("Generating network bursting plot...")
        ax[1] = plot_network_bursting_v2(ax[1], bursting_ax, mega_ax=mega_bursting_ax)

        # adjust axes
        print("Adjusting axes...")         # ensure x-axes are the same for both plots
        x_min = min([ax[0].get_xlim()[0], ax[1].get_xlim()[0]])         # get the narrowes x-axis range shared by both plots
        x_max = max([ax[0].get_xlim()[1], ax[1].get_xlim()[1]])
        ax[0].set_xlim(x_min, x_max)
        ax[1].set_xlim(x_min, x_max)
        
        ## override y axis to 17.5
        ## ==================
        #ax[0].set_ylim(0, 17.5)
        ax[1].set_ylim(0, 17.5)
        #ax[2].set_ylim(0, 17.5)
        ##==================
        
        # override remove legends
        # ==================
        #ax[0].legend().remove()
        ax[1].legend().remove()
        # ==================
        
        plt.tight_layout()
        
        # set limit on time as needed
        if limit_seconds is not None:
            try:
                ax[0].set_xlim(0, limit_seconds)
                ax[1].set_xlim(0, limit_seconds)
                ax[2].set_xlim(0, limit_seconds) # just using try except to catch the error if it's a 2 panel plot
            except:
                pass
                
        #save png and pdf files
        fig.savefig(save_path_pdf_2p) #save as pdf
        print(f"Network summary plot saved to {save_path_pdf_2p}")
        fig.savefig(save_path_png_2p, dpi=600) #save as png
        print(f"Network summary plot saved to {save_path_png_2p}") 
        
    elif mode == '3p':
        # init plot
        fig, ax = plt.subplots(3, 1, figsize=(16, 9))
        
        # plot raster plot
        print("Generating raster plot...")
        ax[0] = plot_raster_plot_v2(ax[0], spiking_data_by_unit, unit_types)
        
        # plot network bursting plots
        print("Generating network bursting plot...")
        ax[1] = plot_network_bursting_v2(ax[1], bursting_ax)
        
        # plot network bursting plots
        print("Generating mega network bursting plot...")
        ax[2] = plot_network_bursting_v2(ax[2], mega_bursting_ax, mode='mega')
        
        # adjust axes
        print("Adjusting axes...")
        # ensure x-axes are the same for both plots
        # get the narrowes x-axis range shared by both plots
        x_min = min([ax[0].get_xlim()[0], ax[1].get_xlim()[0], ax[2].get_xlim()[0]])
        x_max = max([ax[0].get_xlim()[1], ax[1].get_xlim()[1], ax[2].get_xlim()[1]])
        ax[0].set_xlim(x_min, x_max)
        ax[1].set_xlim(x_min, x_max)
        ax[2].set_xlim(x_min, x_max)
        
        ## override y axis to 17.5
        ## ==================
        #ax[0].set_ylim(0, 17.5)
        ax[1].set_ylim(0, 17.5)
        ax[2].set_ylim(0, 17.5)
        ##==================

        # override remove legends
        # ==================
        #ax[0].legend().remove()
        ax[1].legend().remove()
        ax[2].legend().remove()
        # ==================
        
        plt.tight_layout()
        
        # set limit on time as needed
        if limit_seconds is not None:
            try:
                ax[0].set_xlim(0, limit_seconds)
                ax[1].set_xlim(0, limit_seconds)
                ax[2].set_xlim(0, limit_seconds) # just using try except to catch the error if it's a 2 panel plot
            except:
                pass
                
        #save png and pdf files
        fig.savefig(save_path_pdf_3p) #save as pdf
        print(f"Network summary plot saved to {save_path_pdf_3p}")
        fig.savefig(save_path_png_3p, dpi=600) #save as png
        print(f"Network summary plot saved to {save_path_png_3p}") 

def compute_bursting_metrics(network_data, kwargs):

    #
    spike_times = network_data['spiking_data']['spike_times']
    spike_times_by_unit = network_data['spiking_data']['spiking_times_by_unit']
    debug_mode = kwargs.get('debug_mode', False)
    
    #
    conv_params = kwargs.get('conv_params', None)
    try: conv_params = conv_params.conv_params if conv_params is not None else None
    except: conv_params = conv_params
    assert conv_params is not None, 'No convolution parameters provided'
    assert 'binSize' in conv_params, 'Convolution parameters must include "binSize"'
    assert 'gaussianSigma' in conv_params, 'Convolution parameters must include "gaussianSigma"'
    
    mega_conv_params = kwargs.get('mega_params', None)
    try: mega_conv_params = mega_conv_params.mega_conv_params if mega_conv_params is not None else None
    except: mega_conv_params = mega_conv_params
    assert mega_conv_params is not None, 'No mega convolution parameters provided'
    assert 'binSize' in mega_conv_params, 'Mega convolution parameters must include "binSize"'
    assert 'gaussianSigma' in mega_conv_params, 'Mega convolution parameters must include "gaussianSigma"'
    
    try:
        bursting_data = analyze_bursting_activity_v4(
                                spike_times, spike_times_by_unit,
                                binSize = conv_params['binSize'],
                                gaussianSigma = conv_params['gaussianSigma'],
                                thresholdBurst = conv_params['thresholdBurst'],
                                min_peak_distance = conv_params['min_peak_distance'],
                                prominence = conv_params['prominence'],
                                #title = 'Network Activity'
                                **kwargs
                                )
    except Exception as e:
        print(f'Error analyzing bursting activity: {e}')
        traceback.print_exc()
        bursting_data = None
    
    try:
        mega_bursting_data = analyze_bursting_activity_v4(
                                spike_times, spike_times_by_unit,
                                binSize = mega_conv_params['binSize'],
                                gaussianSigma = mega_conv_params['gaussianSigma'],
                                thresholdBurst = mega_conv_params['thresholdBurst'],
                                min_peak_distance = mega_conv_params['min_peak_distance'],
                                prominence = mega_conv_params['prominence'],
                                #title = 'Network Activity'
                                **kwargs
                                )
    except Exception as e:
        print(f'Error analyzing mega bursting activity: {e}')
        traceback.print_exc()
        mega_bursting_data = None
    
    # update network data
    network_data.update({
        'bursting_data': bursting_data,
        'mega_bursting_data': mega_bursting_data,
    })
    
    return network_data

def compute_unit_spike_metrics(unit, spike_times_by_unit, wfe, sampling_rate, plot_wfs, 
                 recording_object, sorting_object, wf_well_folder, **pkwargs):
    """Function to process a single unit."""
    try:
        print(f'Processing unit {unit}...')
        
        if wfe is None:
            try: isi_diffs = np.diff(spike_times_by_unit[unit])
            except: isi_diffs = None
            return unit, {
                'num_spikes': len(spike_times_by_unit[unit]) if unit in spike_times_by_unit else 0,
                'wf_metrics': 'Not implemented...yet',
                #'fr': get_unit_fr(recording_object, sorting_object, unit, sampling_rate),
                'fr': len(spike_times_by_unit[unit])/(spike_times_by_unit[unit][-1] - spike_times_by_unit[unit][0]) if unit in spike_times_by_unit else 0,
                'isi': {
                    'data': isi_diffs,
                    'mean': np.nanmean(isi_diffs) if isi_diffs is not None else np.nan,
                    'std': np.nanstd(isi_diffs) if isi_diffs is not None else np.nan,
                    'median': np.nanmedian(isi_diffs) if isi_diffs is not None else np.nan,
                    #'cov': np.nanstd(isi_diffs) / np.nanmean(isi_diffs) if np.nanmean(isi_diffs) > 0 else np.nan,
                    'cov': np.nanstd(isi_diffs) / np.nanmean(isi_diffs) if isi_diffs is not None and np.nanmean(isi_diffs) > 0 else np.nan,
                    'max': np.nanmax(isi_diffs) if isi_diffs is not None else np.nan,
                    'min': np.nanmin(isi_diffs) if isi_diffs is not None else np.nan,
                },
                'spike_times': spike_times_by_unit[unit] if unit in spike_times_by_unit else None,
            }        
        else:        
            # 
            unit_wf_path = wf_well_folder.replace('waveforms', 'waveform_plots') + f"/unit_{unit}_waveforms.png"
            unit_wfs = wfe.get_waveforms(unit)        
            avg_waveform = np.nanmean(unit_wfs, axis=0)  # Shape: (n_samples, n_channels)
            best_channel_idx = np.argmax(np.max(np.abs(avg_waveform), axis=0))  # Index of best channel
            best_channel_waveforms = unit_wfs[:, :, best_channel_idx]
            try: isi_diffs = np.diff(spike_times_by_unit[unit])
            except: isi_diffs = None
            
            # get waveform metrics
            wf_metrics = compute_wf_metrics(best_channel_waveforms, sampling_rate, plot_wf=plot_wfs, save_fig=True, unit=unit, fig_name=unit_wf_path)
            
            try:
                # Handle cases where isi_diffs is an empty list, None, or np.nan
                if isi_diffs is None or len(isi_diffs) == 0 or np.isnan(isi_diffs).all():
                    isi_diffs = np.array([np.nan])

                results = {
                    'num_spikes': len(spike_times_by_unit[unit]) if unit in spike_times_by_unit else 0,
                    'wf_metrics': wf_metrics,
                    'fr': len(spike_times_by_unit[unit])/(spike_times_by_unit[unit][-1] - spike_times_by_unit[unit][0]) if unit in spike_times_by_unit else 0,
                    'isi': {
                        'data': isi_diffs,
                        'mean': np.nanmean(isi_diffs),
                        'std': np.nanstd(isi_diffs),
                        'median': np.nanmedian(isi_diffs),
                        'cov': np.nanstd(isi_diffs) / np.nanmean(isi_diffs) if np.nanmean(isi_diffs) > 0 else np.nan,
                        'max': np.nanmax(isi_diffs),
                        'min': np.nanmin(isi_diffs),
                    },
                    'spike_times': spike_times_by_unit[unit] if unit in spike_times_by_unit else None,
                }
            except Exception as e:
                traceback.print_exc()
                print(f'Error processing unit {unit}: {e}')
                return unit, None
                        
                        
            # return unit, {
            #     'num_spikes': len(spike_times_by_unit[unit]),
            #     'wf_metrics': wf_metrics,
            #     'fr': get_unit_fr(recording_object, sorting_object, unit, sampling_rate),
            #     'isi': {
            #         'data': isi_diffs,
            #         'mean': np.nanmean(isi_diffs) if isi_diffs is not None else np.nan,
            #         'std': np.nanstd(isi_diffs) if isi_diffs is not None else np.nan,
            #         'median': np.nanmedian(isi_diffs) if isi_diffs is not None else np.nan,
            #         'cov': np.nanstd(isi_diffs) / np.nanmean(isi_diffs) if isi_diffs is not None and np.nanmean(isi_diffs) > 0 else np.nan,
            #         'max': np.nanmax(isi_diffs) if isi_diffs is not None else np.nan,
            #         'min': np.nanmin(isi_diffs) if isi_diffs is not None else np.nan,
            #     },
            #     'spike_times': spike_times_by_unit[unit],
            # }
            return unit, results

    except Exception as e:
        traceback.print_exc()
        print(f'Error processing unit {unit}: {e}')
        return unit, None  # Return None in case of an error

def compute_spike_metrics_by_unit(network_data, kwargs):
    """Parallelized function to extract spiking metrics from experimental data using ThreadPoolExecutor."""
    # init function ============================================================
    print("Computing spiking metrics by unit...")
    indent_increase()
    source = network_data['source'] # important to treat simulated and experimental data differently at this step.
        
    # just for testing #TODO: move to the run script later
    # kwargs['run_parallel'] = True
    # kwargs['debug_mode'] = False
    # kwargs['max_workers'] = 2
    # kwargs['plot_wfs'] = True
    
    
    # Subfunctions =============================================================
    def process_units_in_parallel(units, wfe, network_data, kwargs):
        '''Process units in parallel.'''       
        #init
        indent_increase()
        max_workers = kwargs['max_workers']
        debug_mode = kwargs.get('debug_mode', False) # Impose limit on number of units for debugging purposes
        if debug_mode:
            print(f'Debug mode: Only processing 10 units.')
            num_units = len(units)
            if num_units > 10:
                units = units[:10]
                max_workers = 10
            else:
                max_workers = num_units
                
        # Set parameters for process_unit       
        pkwargs = {
            'spike_times_by_unit': network_data['spiking_data']['spiking_times_by_unit'],
            'wfe': wfe,
            'sampling_rate': network_data['sampling_rate'],
            'plot_wfs': kwargs['plot_wfs'] if source == 'experimental' else False,
            'recording_object': kwargs.get('recording_object', None),
            'sorting_object': kwargs.get('sorting_object', None),
            'sorting_object': kwargs.get('sorting_object', None),
            'wf_well_folder': wfe.folder._str if wfe is not None else None
        }
        
        # Process units in parallel
        print(f'Processing {len(units)} units in parallel using {max_workers} workers...')
        spiking_metrics_by_unit = {}
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_unit = {
                executor.submit(compute_unit_spike_metrics, unit, **pkwargs): unit for unit in units
            }

            # Collect results as they complete
            for future in as_completed(future_to_unit):
                unit, result = future.result()
                if result is not None:
                    spiking_metrics_by_unit[unit] = result

        indent_decrease()
        network_data['spiking_data']['spiking_metrics_by_unit'] = spiking_metrics_by_unit
        return network_data
    
    def process_units_in_sequence(units, wfe, network_data, kwargs):
        raise NotImplementedError('This function is not yet implemented.')
        # for unit in units:
        #     unit_result = process_unit(unit, spike_times_by_unit, wf_extractor, sampling_rate, plot_wfs, recording_object, sorting_object, wf_well_folder)
        #     if unit_result is not None:
        #         spiking_metrics_by_unit[unit] = unit_result
        # network_data['spiking_data']['spiking_metrics_by_unit'] = spiking_metrics_by_unit
        # return network_data
    
    def get_units_ids(source, network_data, kwargs):
        '''Get unit IDs for simulated or experimental data.'''
        if source == 'simulated':
            cellData = kwargs['cellData']
            units = [int(i['gid']) for i in cellData]
            network_data['gids'] = units
            network_data['unit_ids'] = None
        elif source == 'experimental':
            sorting_object = kwargs['sorting_object']
            units = sorting_object.get_unit_ids()
            network_data['unit_ids'] = units
            network_data['gids'] = None
        return units, network_data
        
        # raise NotImplementedError('This function is not yet implemented.')
        # # Get unit IDs
        # units = sorting_object.get_unit_ids()
        # return units
    
    def define_wfe(source, network_data, kwargs):
        '''Define waveform extractor for simulated or experimental data.'''
        if 'wf_extractor' not in kwargs: 
            kwargs['wf_extractor'] = None # usual case for simulated data
        elif kwargs['wf_extractor'] is not None and 'simulated' in source:
            pass # this is an acceptable case for simulated data - just unexpected to see wf_extractor defined as None already.
        elif kwargs['wf_extractor'] is not None and 'experimental' in source:
            pass # this is the normal case for experimental data - just putting this here for clarity.
        elif kwargs['wf_extractor'] is None and 'experimental' in source:
            raise ValueError('No waveform extractor provided for experimental data.')
        wfe = kwargs['wf_extractor'] # NOTE: this is none for simulated data ^^^ 
        return wfe, kwargs
    # Main =====================================================================
    #runtime options
    debug_mode = kwargs.get('debug_mode', False) # Impose limit on number of units for debugging purposes
    run_parallel = kwargs.get('run_parallel', False) # Run in parallel or sequential mode
    
    # set wf_extractor to None by default. If the data are experimental, this should be set.
    #print('Warning: not sure if I should get wfe from network_data or kwargs.')
    wfe, kwargs = define_wfe(source, network_data, kwargs)
    
    # get units
    units, network_data = get_units_ids(source, network_data, kwargs)
    
    # process units
    if run_parallel: network_data = process_units_in_parallel(units, wfe, network_data, kwargs)
    else: network_data = process_units_in_sequence(units, wfe, network_data, kwargs)
    
    # end func =================================================================
    indent_decrease()
    print("Spiking metrics by unit computed!")
    return network_data
        # ====================
        # old code
        # ====================
        
        # # Extract objects from kwargs
        # recording_object = kwargs['recording_object']
        # sorting_object = kwargs['sorting_object']
        # wf_extractor = kwargs['wf_extractor']
        # sampling_rate = recording_object.get_sampling_frequency()
        
        # #add immediately available data to network_data
        # network_data['source'] = 'experimental'
        # network_data['timeVector'] = timeVector
        # network_data['spiking_data']['spike_times'] = spike_times
        # network_data['spiking_data']['spiking_times_by_unit'] = spike_times_by_unit

        # # Get unit IDs
        # units = sorting_object.get_unit_ids()

        # # Set directory for plot output
        # plot_wfs = kwargs['plot_wfs']
        # wf_posixpath = wf_extractor.folder
        # #wf_well_folder = str(wf_posixpath) + '_plots/unit_wf_plots'  
        # wf_well_folder = str(wf_posixpath)
        
        
        # spiking_metrics_by_unit = {}
        
        # #set max workers and threads
        # max_workers = kwargs['max_workers']
        # max_threads = max_workers//2
        
        # #if debug mode, shorten units in units to speed things up. Shorten to 10.
        # if debug_mode:
        #     print(f'Debug mode: Only processing 10 units.')
        #     units = units[:10]
        #     max_workers = 10
        #     max_threads = 5
            

        # threads = False
        # procs = True
        # # threads = True
        # # procs = False
        # if threads:
        #     # Parallel processing with at most 25 threads at a time
        #     #max_workers = 25  
        #     #max_workers = 50 # Increase number of threads to 50 - should be same amount of cpus right?
        #     #max_workers = kwargs['max_workers']
        #     #max_threads = max_workers//2
        #     with ThreadPoolExecutor(max_workers=max_threads) as executor:
        #         future_to_unit = {
        #             executor.submit(process_unit, unit, spike_times_by_unit, wf_extractor, sampling_rate, plot_wfs, 
        #                             recording_object, sorting_object, wf_well_folder): unit for unit in units
        #         }

        #         # Collect results as they complete
        #         for future in as_completed(future_to_unit):
        #             unit, result = future.result()
        #             if result is not None:
        #                 spiking_metrics_by_unit[unit] = result

        #     # Store computed spiking metrics
        #     network_data['spiking_data']['spiking_metrics_by_unit'] = spiking_metrics_by_unit
            
        #     return network_data
        # if procs:
            # # Parallel processing with at most 25 processes at a time
            # #max_workers = 25
            # #max_workers = kwargs['max_workers']  
            # with ProcessPoolExecutor(max_workers=max_workers) as executor:
            #     future_to_unit = {
            #         executor.submit(process_unit, unit, spike_times_by_unit, wf_extractor, sampling_rate, plot_wfs, 
            #                         recording_object, sorting_object, wf_well_folder): unit for unit in units
            #     }

            #     # Collect results as they complete
            #     for future in as_completed(future_to_unit):
            #         unit, result = future.result()
            #         if result is not None:
            #             spiking_metrics_by_unit[unit] = result

            # # Store computed spiking metrics
            # network_data['spiking_data']['spiking_metrics_by_unit'] = spiking_metrics_by_unit
            
            # return network_data

def compute_network_metrics(conv_params, mega_params, source, **kwargs):
    '''
    Get network metrics for simulated or experimental data.
    
    Ideally this minimizes changes to computing network metrics for either simulated or experimental data in the future by consolidating the process into a single function.
    Here, we focus on pre-processing either type of data into like formats before computing network metrics - which should be mostly identical downstream.
    
    Core philosophies here:
    - every subfunction should be able to handle either simulated or experimental data.
    - use as many subfuncs as possible to keep the main function clean and readable.
    - use as many shared functions as possible to keep the code DRY. (Don't Repeat Yourself)
    - keep the main function as simple as possible.
    - always put try and except logic in subfunctions to catch errors and print them out. Use traceback to get the full error message.
    '''
    # debug - move to run script later #HACK
    # fitness_save_path = kwargs['fitness_save_path']
    # basename = os.path.basename(fitness_save_path)
    # sa_dir = os.path.dirname(fitness_save_path)
    # # remove .json from basename
    # basename = basename.replace('.json', '')
    # dtw_dir = os.path.join(sa_dir, basename, 'dtw_temp')    
    # kwargs['dtw_temp'] = dtw_dir
    
    # Subfunctions =============================================================
    def validate_inputs(source, kwargs):
        '''Validate input data.'''
        assert source in ['simulated', 'experimental'], 'Invalid data source'
        if source == 'simulated':
            assert 'simData' in kwargs, 'No simData provided'
            assert 'popData' in kwargs, 'No popData provided'
            assert 'cellData' in kwargs, 'No cellData provided'
        elif source == 'experimental':
            assert 'sorting_object' in kwargs, 'No sorting object provided'
            assert 'recording_object' in kwargs, 'No recording object provided'
            assert 'wf_extractor' in kwargs, 'No waveform extractor provided'
        #return source, kwargs
    
    def initialize_data_dict(source, kwargs):
        '''Initialize dictionary to store network data.'''
        global network_data
        network_data = {}
        
        try:
            if source == 'simulated':
                simData = kwargs['simData']
                network_data['source'] = 'simulated'
                network_data['spiking_data'] = {}
                network_data['sim_data_path'] = kwargs.get('sim_data_path', None)
                return network_data
            elif source == 'experimental':
                sorting_object = kwargs['sorting_object']
                recording_object = kwargs['recording_object']
                network_data['source'] = 'experimental'
                assert sorting_object is not None, 'No sorting object provided'
                network_data['recording_path'] = recording_object._kwargs.get('file_path', None)
                network_data['sorting_output'] = sorting_object._kwargs.get('folder_path', None)
                assert network_data['sorting_output'] is not None, 'No sorting object source found'
                network_data['waveform_output'] = network_data['sorting_output'].replace('sorted', 'waveforms').replace('sorter_output', '')
                network_data['network_metrics_output'] = network_data['sorting_output'].replace('sorted', 'network_metrics').replace('sorter_output', '')
                network_data['spiking_data'] = {}
                #network_data['spiking_metrics'] = {}
                #network_data['bursting_metrics'] = {}
                return network_data
        except Exception as e:
            print(f'Error initializing network data: {e}')
            traceback.print_exc()
            pass
    
    def update_kwargs(conv_params, mega_params, source, kwargs):
        '''Update kwargs with additional parameters.'''
        try:
            if source == 'simulated':
                kwargs['source'] = source
                kwargs['conv_params'] = conv_params
                kwargs['mega_params'] = mega_params
            elif source == 'experimental':
                # sorting_object = kwargs['sorting_object']
                # recording_object = kwargs['recording_object']
                # wf_extractor = kwargs['wf_extractor']
                kwargs['source'] = source
                kwargs['conv_params'] = conv_params
                kwargs['mega_params'] = mega_params
                # kwargs['sorting_object'] = sorting_object
                # kwargs['recording_object'] = recording_object
                # kwargs['wf_extractor'] = wf_extractor
            return kwargs
        except Exception as e:
            print(f'Error updating kwargs: {e}')
            traceback.print_exc()
            pass
    
    def initialize_spike_data(network_data, kwargs):
        '''Initialize spike data for simulated or experimental data.'''
        try:
            source = kwargs['source']
            if source == 'simulated':
                simData = kwargs['simData']
                rasterData = simData.copy()
                if len(rasterData['spkt']) == 0:
                    #print('No spike times found in rasterData')
                    raise ValueError('No spike times found in rasterData')
                spike_times = np.array(rasterData['spkt']) / 1000
                timeVector = np.array(rasterData['t']) / 1000
                spike_times_by_unit = {int(i): spike_times[rasterData['spkid'] == i] for i in np.unique(rasterData['spkid'])}
                network_data['timeVector'] = timeVector
                # NOTE: idk, sampling_rate is kindof available in simData if I record full trace of a particular neuron. Not super important for now.
                network_data['sampling_rate'] = 'Not currently implemented for simulated data'
                network_data['spiking_data']['spike_times'] = spike_times
                network_data['spiking_data']['spiking_times_by_unit'] = spike_times_by_unit
                return spike_times, timeVector, spike_times_by_unit, network_data
            elif source == 'experimental':
                sorting_object = kwargs['sorting_object']
                recording_object = kwargs['recording_object']
                sampling_rate = recording_object.get_sampling_frequency()
                spike_times = get_spike_times(recording_object, sorting_object, sampling_rate=sampling_rate) #seconds
                timeVector = get_time_vector(recording_object, sampling_rate=sampling_rate) #seconds
                spike_times_by_unit = get_spike_times_by_unit(sorting_object, sampling_rate=sampling_rate)
                network_data['timeVector'] = timeVector
                network_data['sampling_rate'] = sampling_rate
                network_data['spiking_data']['spike_times'] = spike_times
                network_data['spiking_data']['spiking_times_by_unit'] = spike_times_by_unit
                return spike_times, timeVector, spike_times_by_unit, network_data
        except Exception as e:
            print(f'Error initializing spike data: {e}')
            traceback.print_exc()
            pass
    
    def compute_spike_metrics(network_data, kwargs):

        #extract spiking metrics from either simulated or experimental data
        #source = kwargs['source']
        try: 
            # individual unit metrics
            network_data = compute_spike_metrics_by_unit(network_data, kwargs)
        except Exception as e:
            print(f'Error extracting metrics from simulated data: {e}')
            traceback.print_exc()
            pass
        
        return network_data
        
    def compute_burst_metrics(network_data, kwargs):
        #extract bursting metrics from simulated data (but this one works for both simulated and experimental data)
        try: 
            network_data = compute_bursting_metrics(network_data, kwargs)
        except Exception as e:
            print(f'Error calculating bursting activity: {e}')
            traceback.print_exc()
            pass
        
        return network_data
        
    def classify_units(network_data, kwargs):     
        # classify neurons
        source = kwargs['source']
        try:
            if source == 'simulated':
                cellData = kwargs['cellData']
                unit_types = {int(i['gid']): i['tags']['pop'] for i in cellData}
                network_data['unit_types'] = unit_types
            elif source == 'experimental':
                network_data = classify_neurons_v2(network_data, **kwargs)
                classification_data = network_data['classification_output']
                #network_data['unit_types'] = classification_data['classified_units']
                unit_types = {}
                for unit_id, unit in classification_data['classified_units'].items():
                    desc = unit['desc']
                    if 'inhib' in desc:
                        unit_types[unit_id] = 'I'
                    elif 'excit' in desc:
                        unit_types[unit_id] = 'E'
                    else:
                        unit_types[unit_id] = 'U'
                network_data['unit_types'] = unit_types
        except Exception as e:
            print(f'Error classifying neurons: {e}')
            traceback.print_exc()
            pass
        
        return network_data
        
    def compute_dynamic_time_warping(network_data, kwargs):
        # dynamic time warping (dtw) - burst analysis
        try:
            network_data = dtw_burst_analysis(network_data, kwargs)
        except Exception as e:
            print(f'Error in dtw burst analysis: {e}')
            traceback.print_exc()
            pass  
        
        return network_data
    
    def compute_summary_metrics(network_data, kwargs):
        
        # fr
        spiking_data = network_data['spiking_data']
        frs = [i['fr'] for i in spiking_data['spiking_metrics_by_unit'].values()]
        frs = [np.nan if i == np.inf else i for i in frs] #replace inf with nan
        spiking_data['frs'] = {
            'data': frs if len(frs) > 0 else None,
            'mean': np.nanmean(frs) if len(frs) > 0 else None,
            'std': np.nanstd(frs) if len(frs) > 0 else None,
            'median': np.nanmedian(frs) if len(frs) > 0 else None,
            'cov': np.nanstd(frs) / np.nanmean(frs) if np.nanmean(frs) and len(frs) > 0 else None,
            'max': np.nanmax(frs) if len(frs) > 0 else None,
            'min': np.nanmin(frs) if len(frs) > 0 else None,
        }
        
        # E/U frs
        unit_types = network_data['unit_types']
        i_frs = [metrics['fr'] for i, metrics in spiking_data['spiking_metrics_by_unit'].items() if i in unit_types and unit_types[i] == 'I']
        e_frs = [metrics['fr'] for i, metrics in spiking_data['spiking_metrics_by_unit'].items() if i in unit_types and unit_types[i] == 'E']
        u_frs = [metrics['fr'] for i, metrics in spiking_data['spiking_metrics_by_unit'].items() if i not in unit_types]
        i_frs = [np.nan if i == np.inf else i for i in i_frs] #replace inf with nan
        e_frs = [np.nan if i == np.inf else i for i in e_frs]
        u_frs = [np.nan if i == np.inf else i for i in u_frs]
        spiking_data['i_frs'] = {
            'data': i_frs if len(i_frs) > 0 else None,
            'mean': np.nanmean(i_frs) if len(i_frs) > 0 else None,
            'std': np.nanstd(i_frs) if len(i_frs) > 0 else None,
            'median': np.nanmedian(i_frs) if len(i_frs) > 0 else None,
            'cov': np.nanstd(i_frs) / np.nanmean(i_frs) if np.nanmean(i_frs) and len(i_frs) > 0 else None,
            'max': np.nanmax(i_frs) if len(i_frs) > 0 else None,
            'min': np.nanmin(i_frs) if len(i_frs) > 0 else None,
        }
        spiking_data['e_frs'] = {
            'data': e_frs if len(e_frs) > 0 else None,
            'mean': np.nanmean(e_frs) if len(e_frs) > 0 else None,
            'std': np.nanstd(e_frs) if len(e_frs) > 0 else None,
            'median': np.nanmedian(e_frs) if len(e_frs) > 0 else None,
            'cov': np.nanstd(e_frs) / np.nanmean(e_frs) if np.nanmean(e_frs) and len(e_frs) > 0 else None,
            'max': np.nanmax(e_frs) if len(e_frs) > 0 else None,
            'min': np.nanmin(e_frs) if len(e_frs) > 0 else None,
        }
        spiking_data['u_frs'] = {
            'data': u_frs if len(u_frs) > 0 else None,
            'mean': np.nanmean(u_frs) if len(u_frs) > 0 else None,
            'std': np.nanstd(u_frs) if len(u_frs) > 0 else None,
            'median': np.nanmedian(u_frs) if len(u_frs) > 0 else None,
            'cov': np.nanstd(u_frs) / np.nanmean(u_frs) if np.nanmean(u_frs) and len(u_frs) > 0 else None,
            'max': np.nanmax(u_frs) if len(u_frs) > 0 else None,
            'min': np.nanmin(u_frs) if len(u_frs) > 0 else None,
        }
        
        # isi
        isis = [i['isi']['mean'] for i in spiking_data['spiking_metrics_by_unit'].values()]
        isis = [np.nan if i == np.inf else i for i in isis]
        spiking_data['isi'] = {
            'data': isis if len(isis) > 0 else None,
            'mean': np.nanmean(isis) if len(isis) > 0 else None,
            'std': np.nanstd(isis) if len(isis) > 0 else None,
            'median': np.nanmedian(isis) if len(isis) > 0 else None,
            'cov': np.nanstd(isis) / np.nanmean(isis) if np.nanmean(isis) and len(isis) > 0 else None,
            'max': np.nanmax(isis) if len(isis) > 0 else None,
            'min': np.nanmin(isis) if len(isis) > 0 else None,
        }
        
        # E/I isi
        i_isis = [metrics['isi']['mean'] for i, metrics in spiking_data['spiking_metrics_by_unit'].items() if i in unit_types and unit_types[i] == 'I'] 
        e_isis = [metrics['isi']['mean'] for i, metrics in spiking_data['spiking_metrics_by_unit'].items() if i in unit_types and unit_types[i] == 'E']
        u_isis = [metrics['isi']['mean'] for i, metrics in spiking_data['spiking_metrics_by_unit'].items() if i not in unit_types]  
        i_isis = [np.nan if i == np.inf else i for i in i_isis]
        e_isis = [np.nan if i == np.inf else i for i in e_isis]
        u_isis = [np.nan if i == np.inf else i for i in u_isis]
        spiking_data['i_isi'] = {
            'data': i_isis if len(i_isis) > 0 else None,
            'mean': np.nanmean(i_isis) if len(i_isis) > 0 else None,
            'std': np.nanstd(i_isis) if len(i_isis) > 0 else None,
            'median': np.nanmedian(i_isis) if len(i_isis) > 0 else None,
            'cov': np.nanstd(i_isis) / np.nanmean(i_isis) if np.nanmean(i_isis) and len(i_isis) > 0 else None,
            'max': np.nanmax(i_isis) if len(i_isis) > 0 else None,
            'min': np.nanmin(i_isis) if len(i_isis) > 0 else None,
        }
        
        spiking_data['e_isi'] = {
            'data': e_isis if len(e_isis) > 0 else None,
            'mean': np.nanmean(e_isis) if len(e_isis) > 0 else None,
            'std': np.nanstd(e_isis) if len(e_isis) > 0 else None,
            'median': np.nanmedian(e_isis) if len(e_isis) > 0 else None,
            'cov': np.nanstd(e_isis) / np.nanmean(e_isis) if np.nanmean(e_isis) and len(e_isis) > 0 else None,
            'max': np.nanmax(e_isis) if len(e_isis) > 0 else None,
            'min': np.nanmin(e_isis) if len(e_isis) > 0 else None,
        }
        
        spiking_data['u_isi'] = {
            'data': u_isis if len(u_isis) > 0 else None,
            'mean': np.nanmean(u_isis) if len(u_isis) > 0 else None,
            'std': np.nanstd(u_isis) if len(u_isis) > 0 else None,
            'median': np.nanmedian(u_isis) if len(u_isis) > 0 else None,
            'cov': np.nanstd(u_isis) / np.nanmean(u_isis) if np.nanmean(u_isis) and len(u_isis) > 0 else None,
            'max': np.nanmax(u_isis) if len(u_isis) > 0 else None,
            'min': np.nanmin(u_isis) if len(u_isis) > 0 else None,
        }
        
        # burst metrics
        #bursting_data = network_data['bursting_data']
        #burst_metrics = bursting_data['burst_metrics']
        
        
        #raise NotImplementedError('This function is not yet implemented.')
        return network_data
    
    def locate_units(network_data, kwargs):
        source = kwargs['source']
        if source == 'simulated':
            # TODO: EXTRACT from SIMDATA
            # raise NotImplementedError('This function is not yet implemented.')
            cellData = kwargs['cellData']
            unit_locations_dict = {}
            for cell in cellData:
                #cell['location'] = None
                gid = int(cell['gid'])
                x = cell['tags']['x']
                y = cell['tags']['y']
                z = cell['tags']['z']
                unit_locations_dict[gid] = (x, y, z) # TODO: extract x, y, z from cell
            
            network_data['unit_locations'] = unit_locations_dict
            
            return network_data
        elif source == 'experimental':
            we = kwargs.get('wf_extractor', None)
            
            classification_output = network_data['classification_output']
            include_unit_ids = classification_output['include_units']
            # classified_units = classification_output['classified_units']
        
            unit_locations = spost.compute_unit_locations(we)
            unit_locations_dict = {unit_id: unit_locations[i] for i, unit_id in enumerate(include_unit_ids)}
            
            network_data['unit_locations'] = unit_locations_dict
            
            return network_data
            # inhib_neuron_locs = np.array([unit_locations_dict[i] for i in include_unit_ids if classified_units[i]['desc'] == 'inhib'])
            # excit_neuron_locs = np.array([unit_locations_dict[i] for i in include_unit_ids if classified_units[i]['desc'] == 'excit'])
    
    def add_sim_specific_data(network_data, kwargs):
        # add simData, popData, cellData to network_data - stuff that should be in kwargs
        network_data['simData'] = kwargs['simData']
        network_data['popData'] = kwargs['popData']
        network_data['cellData'] = kwargs['cellData']
        return network_data
    
    # Main =====================================================================
    # check data source and validate inputs based on source
    validate_inputs(source, kwargs)
    
    #init
    network_data = initialize_data_dict(source, kwargs)
    kwargs = update_kwargs(conv_params, mega_params, source, kwargs) # NOTE: source folded back into kwargs - so no need to explicitly pass source to other functions
    _, _, _, network_data = initialize_spike_data(network_data, kwargs)
    
    # main computation steps
    network_data = compute_spike_metrics(network_data, kwargs)
    network_data = compute_burst_metrics(network_data, kwargs)
    network_data = classify_units(network_data, kwargs)
    network_data = locate_units(network_data, kwargs)
    if source == 'experimental': # for now only do with with experimental data
        network_data = compute_dynamic_time_warping(network_data, kwargs)
    if source == 'simulated': 
        network_data = add_sim_specific_data(network_data, kwargs) # for now only do with simulated data
    network_data = compute_summary_metrics(network_data, kwargs)
    
    # debug
    # unit_types = network_data['unit_types']
    # print(f'Unit types: {unit_types}')
    # import sys
    # sys.exit()
    
    return network_data

'''# aw 2025-02-26+ updates above this point'''
    
def get_simulated_network_metrics_v2(conv_params, mega_params, simData=None, popData=None, **kwargs):
    
    # Subfunctions =============================================================
    def init_network_data(simData):
        '''Initialize dictionary to store network data.'''
        global network_data
        network_data = {}
        network_data['source'] = 'simulated'
        network_data['timeVector'] = None
        network_data['spiking_data'] = {}
        
        assert simData is not None, 'No simData provided'
    
        # copy simData
        rasterData = simData.copy()
        
        #check if rasterData['spkt'] is empty
        if len(rasterData['spkt']) == 0:
            print('No spike times found in rasterData')
            return None
    
        #convert time to seconds - get initially available data
        spike_times = np.array(rasterData['spkt']) / 1000
        timeVector = np.array(rasterData['t']) / 1000
        spike_times_by_unit = {int(i): spike_times[rasterData['spkid'] == i] for i in np.unique(rasterData['spkid'])} #mea_analysis_pipeline.py expects spike_times as dictionary
        
        network_data['timeVector'] = timeVector
        network_data['spiking_data']['spike_times'] = spike_times
        network_data['spiking_data']['spiking_times_by_unit'] = spike_times_by_unit    
            
        return network_data
    
    def get_spike_metrics(spike_times, timeVector, spike_times_by_unit, **kwargs):
        '''Extract spiking metrics from simulated data.'''
        #extract spiking metrics from simulated data
        try: 
            extract_metrics_from_simulated_data(spike_times, timeVector, spike_times_by_unit, **kwargs)
        except Exception as e:
            traceback.print_exc()
            print(f'Error extracting metrics from simulated data: {e}')
            pass
        
    def get_burst_metrics(spike_times, spike_times_by_unit, **kwargs):
        '''Extract bursting metrics from simulated data.'''
        try: 
            extract_bursting_activity_data(spike_times, spike_times_by_unit, **kwargs)
        except Exception as e:
            traceback.print_exc()
            print(f'Error calculating bursting activity: {e}')
            pass
    
    def convert_single_element_arrays(data):
        if isinstance(data, dict):
            for key, value in data.items():
                data[key] = convert_single_element_arrays(value)
        elif isinstance(data, np.ndarray) and data.size == 1:
            return data.item()
        return data
    
    # Main =====================================================================
    #add conv_params to kwargs
    kwargs['conv_params'] = conv_params
    kwargs['mega_params'] = mega_params
    
    #candidate_label = kwargs.get('candidate_label', None)
    print('') #for formatting
    print('Calculating Network Activity Metrics for Simulated Data...')
    start_time = time.time()
    
    # init
    network_data = init_network_data(simData)
    spike_times = network_data['spiking_data']['spike_times']
    timeVector = network_data['timeVector']
    spike_times_by_unit = network_data['spiking_data']['spiking_times_by_unit']
        
    # main computation steps
    get_spike_metrics(spike_times, timeVector, spike_times_by_unit, **kwargs)
    get_burst_metrics(spike_times, spike_times_by_unit, **kwargs)
    network_data = convert_single_element_arrays(network_data)    
    print('Network Activity Metrics Calculated and Extracted!')
    print(f'Elapsed time: {time.time() - start_time} seconds')
    print('') #for formatting    
    
    #
    return network_data

def get_experimental_network_metrics_v3(sorting_object, recording_object, wf_extractor, conv_params, mega_params, debug_mode = False, **kwargs):
    
    # Subfunctions =============================================================
    def initialize_data_dict(sorting_object=None, recording_segment_object=None):
        '''Initialize dictionary to store network data.'''
        global network_data
        network_data = {}
        network_data['source'] = 'experimental'
        assert sorting_object is not None, 'No sorting object provided'
        network_data['recording_path'] = recording_segment_object._kwargs.get('file_path', None)
        network_data['sorting_output'] = sorting_object._kwargs.get('folder_path', None)
        assert network_data['sorting_output'] is not None, 'No sorting object source found'
        network_data['waveform_output'] = network_data['sorting_output'].replace('sorted', 'waveforms').replace('sorter_output', '')
        network_data['network_metrics_output'] = network_data['sorting_output'].replace('sorted', 'network_metrics').replace('sorter_output', '')
        network_data['spiking_data'] = {}
        #network_data['spiking_metrics'] = {}
        #network_data['bursting_metrics'] = {}
        return network_data
    
    def update_kwargs(kwargs, conv_params, mega_params, sorting_object, recording_segment_object, wf_extractor):
        '''Update kwargs with additional parameters.'''
        kwargs['conv_params'] = conv_params
        kwargs['mega_params'] = mega_params
        kwargs['sorting_object'] = sorting_object
        kwargs['recording_object'] = recording_segment_object
        kwargs['wf_extractor'] = wf_extractor
        return kwargs
    
    def compute_spike_metrics(spike_times, timeVector, spike_times_by_unit, debug_mode=False, **kwargs):

        #extract spiking metrics from simulated data
        # aw 2025-02-25 10:58:09 - plotting operations to be mindful of here:
        # 1. plot_wfs: whether to plot waveforms for individual units. Actual plotting occurrs in compute_wf_metrics
        try: 
            #extract_metrics_from_experimental_data_v2(spike_times, timeVector, spike_times_by_unit, **kwargs)
            #extract_metrics_from_experimental_data_v3(spike_times, timeVector, spike_times_by_unit, **kwargs) # aw 2025-02-21 01:12:14 - parallelized version
            extract_metrics_from_experimental_data_v3_threads(spike_times, timeVector, spike_times_by_unit, debug_mode = debug_mode, **kwargs) # aw 2025-02-21 01:12:14 - parallelized version w/ threads might be more efficient?
        except Exception as e:
            print(f'Error extracting metrics from simulated data: {e}')
            traceback.print_exc()
            pass
        
    def compute_burst_metrics(spike_times, spike_times_by_unit, debug_mode=False, **kwargs):        
        #extract bursting metrics from simulated data (but this one works for both simulated and experimental data)
        try: 
            extract_bursting_activity_data(spike_times, spike_times_by_unit, debug_mode = debug_mode, **kwargs)
        except Exception as e:
            print(f'Error calculating bursting activity: {e}')
            traceback.print_exc()
            pass
        
    def classify_putative_units(network_data, **kwargs):
        # classify neurons
        # from RBS_network_models.experimental import classify_neurons        
        try:
            network_data = classify_neurons_v2(network_data, **kwargs)
        except Exception as e:
            print(f'Error classifying neurons: {e}')
            traceback.print_exc()
            pass
        
    def compute_dynamic_time_warping(network_data, **kwargs):
        # dynamic time warping (dtw) - burst analysis
        try:
            network_data = dtw_burst_analysis(network_data, **kwargs)
        except Exception as e:
            print(f'Error in dtw burst analysis: {e}')
            traceback.print_exc()
            pass  
    # Main =====================================================================
    # init
    network_data = initialize_data_dict(sorting_object, recording_object)
    kwargs = update_kwargs(kwargs, conv_params, mega_params, sorting_object, recording_object, wf_extractor)
    
    # intial calcs
    sampling_rate = recording_object.get_sampling_frequency()
    spike_times = get_spike_times(recording_object, sorting_object, sampling_rate=sampling_rate) #seconds
    timeVector = get_time_vector(recording_object, sampling_rate=sampling_rate) #seconds
    spike_times_by_unit = get_spike_times_by_unit(sorting_object, sampling_rate=sampling_rate) 
        
    # main computation steps
    compute_spike_metrics(spike_times, timeVector, spike_times_by_unit, debug_mode=debug_mode, **kwargs)
    compute_burst_metrics(spike_times, spike_times_by_unit, debug_mode=debug_mode, **kwargs)
    classify_putative_units(network_data, **kwargs)
    compute_dynamic_time_warping(network_data, **kwargs)
    
    return network_data

def extract_metrics_from_experimental_data_v3_threads(spike_times, timeVector, spike_times_by_unit, debug_mode = False, **kwargs):
    """Parallelized function to extract spiking metrics from experimental data using ThreadPoolExecutor."""
    
    # Extract objects from kwargs
    recording_object = kwargs['recording_object']
    sorting_object = kwargs['sorting_object']
    wf_extractor = kwargs['wf_extractor']
    sampling_rate = recording_object.get_sampling_frequency()
    
    #add immediately available data to network_data
    network_data['source'] = 'experimental'
    network_data['timeVector'] = timeVector
    network_data['spiking_data']['spike_times'] = spike_times
    network_data['spiking_data']['spiking_times_by_unit'] = spike_times_by_unit

    # Get unit IDs
    units = sorting_object.get_unit_ids()

    # Set directory for plot output
    plot_wfs = kwargs['plot_wfs']
    wf_posixpath = wf_extractor.folder
    #wf_well_folder = str(wf_posixpath) + '_plots/unit_wf_plots'  
    wf_well_folder = str(wf_posixpath)
    
    
    spiking_metrics_by_unit = {}
    
    #set max workers and threads
    max_workers = kwargs['max_workers']
    max_threads = max_workers//2
    
    #if debug mode, shorten units in units to speed things up. Shorten to 10.
    if debug_mode:
        print(f'Debug mode: Only processing 10 units.')
        units = units[:10]
        max_workers = 10
        max_threads = 5
        

    threads = False
    procs = True
    # threads = True
    # procs = False
    if threads:
        # Parallel processing with at most 25 threads at a time
        #max_workers = 25  
        #max_workers = 50 # Increase number of threads to 50 - should be same amount of cpus right?
        #max_workers = kwargs['max_workers']
        #max_threads = max_workers//2
        with ThreadPoolExecutor(max_workers=max_threads) as executor:
            future_to_unit = {
                executor.submit(process_unit, unit, spike_times_by_unit, wf_extractor, sampling_rate, plot_wfs, 
                                recording_object, sorting_object, wf_well_folder): unit for unit in units
            }

            # Collect results as they complete
            for future in as_completed(future_to_unit):
                unit, result = future.result()
                if result is not None:
                    spiking_metrics_by_unit[unit] = result

        # Store computed spiking metrics
        network_data['spiking_data']['spiking_metrics_by_unit'] = spiking_metrics_by_unit
        
        return network_data
    if procs:
        # Parallel processing with at most 25 processes at a time
        #max_workers = 25
        #max_workers = kwargs['max_workers']  
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_unit = {
                executor.submit(process_unit, unit, spike_times_by_unit, wf_extractor, sampling_rate, plot_wfs, 
                                recording_object, sorting_object, wf_well_folder): unit for unit in units
            }

            # Collect results as they complete
            for future in as_completed(future_to_unit):
                unit, result = future.result()
                if result is not None:
                    spiking_metrics_by_unit[unit] = result

        # Store computed spiking metrics
        network_data['spiking_data']['spiking_metrics_by_unit'] = spiking_metrics_by_unit
        
        return network_data
    
def process_unit(unit, spike_times_by_unit, wf_extractor, sampling_rate, plot_wfs, recording_object, sorting_object, wf_well_folder):
    """Function to process a single unit."""
    try:
        print(f'Processing unit {unit}...')
        
        # Set path for unit plot
        #unit_wf_path = f"{wf_well_folder}/unit_{unit}_waveforms.png"
        
        # replace 'waveforms' w/ 'waveform_plots'
        unit_wf_path = wf_well_folder.replace('waveforms', 'waveform_plots') + f"/unit_{unit}_waveforms.png"
        
        # # create parent dir as needed
        # if not os.path.exists(os.path.dirname(unit_wf_path)):
        #     os.makedirs(os.path.dirname(unit_wf_path))        
        
        # Extract waveforms for each unit
        unit_wfs = wf_extractor.get_waveforms(unit)        
        
        # Compute mean waveform across all spikes
        avg_waveform = np.nanmean(unit_wfs, axis=0)  # Shape: (n_samples, n_channels)
        
        # Find the best channel (largest absolute amplitude)
        best_channel_idx = np.argmax(np.max(np.abs(avg_waveform), axis=0))  # Index of best channel
        
        # Collect waveforms from best channel
        best_channel_waveforms = unit_wfs[:, :, best_channel_idx]
        
        # Compute amplitude of the average waveform
        isi_diffs = np.diff(spike_times_by_unit[unit])
        
        return unit, {
            'num_spikes': len(spike_times_by_unit[unit]),
            'wf_metrics': compute_wf_metrics(best_channel_waveforms, sampling_rate, plot_wf=plot_wfs, save_fig=True, unit=unit, fig_name=unit_wf_path),
            'fr': get_unit_fr(recording_object, sorting_object, unit, sampling_rate),
            'isi': {
                'data': isi_diffs,
                'mean': np.nanmean(isi_diffs),
                'std': np.nanstd(isi_diffs),
                'median': np.nanmedian(isi_diffs),
                'cov': np.nanstd(isi_diffs) / np.nanmean(isi_diffs) if np.nanmean(isi_diffs) > 0 else np.nan,
                'max': np.nanmax(isi_diffs),
                'min': np.nanmin(isi_diffs),
            },
            'spike_times': spike_times_by_unit[unit],
        }

    except Exception as e:
        traceback.print_exc()
        print(f'Error processing unit {unit}: {e}')
        return unit, None  # Return None in case of an error

def extract_metrics_from_experimental_data_v3(spike_times, timeVector, spike_times_by_unit, **kwargs):
    """Parallelized function to extract spiking metrics from experimental data with 25 processes at a time."""
    
    # Extract objects from kwargs
    recording_object = kwargs['recording_object']
    sorting_object = kwargs['sorting_object']
    wf_extractor = kwargs['wf_extractor']
    sampling_rate = recording_object.get_sampling_frequency()    

    # Initialize network data
    network_data = {}
    network_data['source'] = 'experimental'
    network_data['timeVector'] = timeVector
    network_data['spiking_data'] = {
        'spike_times': spike_times,
        'spiking_times_by_unit': spike_times_by_unit
    }

    # Get unit IDs
    units = sorting_object.get_unit_ids()

    # Set directory for plot output
    plot_wfs = kwargs['plot_wfs']
    wf_posixpath = wf_extractor.folder
    wf_well_folder = str(wf_posixpath) + '_plots/unit_wf_plots'  

    spiking_metrics_by_unit = {}

    # Parallel processing with at most 25 processes at a time
    max_workers = 25  
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_unit = {
            executor.submit(process_unit, unit, spike_times_by_unit, wf_extractor, sampling_rate, plot_wfs, 
                            recording_object, sorting_object, wf_well_folder): unit for unit in units
        }

        # Collect results as they complete
        for future in as_completed(future_to_unit):
            unit, result = future.result()
            if result is not None:
                spiking_metrics_by_unit[unit] = result

    # Store computed spiking metrics
    network_data['spiking_data']['spiking_metrics_by_unit'] = spiking_metrics_by_unit

def compute_wf_metrics(best_channel_waveforms, sampling_rate, plot_wf=False, save_fig=False, fig_name="waveform_debug.png", unit=None):
    '''Get key waveform metrics from a single unit using a weighted average approach.'''
    
    unit_wfs = best_channel_waveforms
    num_samples = unit_wfs.shape[1]
    sampling_duration = num_samples / sampling_rate * 1000  # Convert to ms
    num_interp_samples = 1000  # Standardized interpolation length
    time_conversion_factor = sampling_duration / num_interp_samples  # Converts interp sample index → ms

    # Interpolate & smooth waveforms
    sigma = 9
    x = np.arange(num_samples)
    x_new = np.linspace(0, num_samples - 1, num_interp_samples)  # Standardize to 1000 samples
    interpolated_wfs = np.array([gaussian_filter1d(np.interp(x_new, x, wf), sigma) for wf in unit_wfs])
    
    # # remove wf if slope is too low
    # for wf in interpolated_wfs:
    #    # max_depolarization_slope = np.max(np.diff(interpolated_wfs[:, :np.argmin(interpolated_wfs[0])])) * (1000 / sampling_rate)  # uV/ms
    #     max_depolarization_slope = np.max(np.diff(wf[:np.argmin(wf)])) * (1000 / sampling_rate)  # uV/ms
    #     if max_depolarization_slope < 150:
    #         print("Warning: Low depolarization slope.")
    #         # remove from interpolated slope
    #         #interpolated_wfs = interpolated_wfs[:-1]
    #         # remove
    #         if wf in interpolated_wfs: interpolated_wfs = np.delete(interpolated_wfs, np.where(interpolated_wfs == wf), axis=0)
    
    # Compute mean waveform
    mean_wf = np.nanmean(interpolated_wfs, axis=0)

    # Compute Mahalanobis distances & weighted averaging
    if len(interpolated_wfs) < 2:
        print("Warning: Not enough waveforms to compute Mahalanobis distances.")
        return {
            'excluded': True,
            'reason': 'Not enough waveforms to compute Mahalanobis distances.'
        }
    try:
        covariance = np.cov(interpolated_wfs, rowvar=False)
        inv_covariance = np.linalg.pinv(covariance)  # Pseudo-inverse in case of singularity
        mean_vector = np.nanmean(interpolated_wfs, axis=0)
    except Exception as e:
        print(f"Error computing covariance matrix: {e}")
        pass

    distances = np.array([mahalanobis(wf, mean_vector, inv_covariance) for wf in interpolated_wfs])
    weights = np.exp(-distances / np.nanmedian(distances))  # Exponential weighting
    weighted_avg_wf = np.average(interpolated_wfs, axis=0, weights=weights)
    
    #measure difference between weighted and unweighted average
    difference = np.abs(weighted_avg_wf - mean_wf)
    
    # Compute weighted variability (biological variability indicators)
    weighted_variability = np.average(np.abs(interpolated_wfs - weighted_avg_wf), axis=0, weights=weights)
    weighted_coefficient_of_variation = weighted_variability / np.abs(weighted_avg_wf)
    
    # Compute weighted variability (biological variability indicators)
    weighted_variability = np.average(np.abs(interpolated_wfs - weighted_avg_wf), axis=0, weights=weights)

    # **Compute Vertical Smearing (Amplitude Variability)**
    vertical_std = np.std(interpolated_wfs, axis=0)  # Std across time points
    peak_amplitudes = np.max(interpolated_wfs, axis=1)  # Peak values per spike
    trough_amplitudes = np.min(interpolated_wfs, axis=1)  # Trough values per spike
    peak_to_trough = peak_amplitudes - trough_amplitudes  # Peak-to-trough amplitude per spike

    vertical_smearing = {
        "std_across_time": vertical_std,
        "cv_peak_amplitude": np.std(peak_amplitudes) / np.mean(peak_amplitudes),
        "cv_trough_amplitude": np.std(trough_amplitudes) / np.mean(trough_amplitudes),
        "peak_to_trough_variance": np.var(peak_to_trough),
    }

    # **Compute Horizontal Smearing (Timing Variability)**
    trough_times = np.argmin(interpolated_wfs, axis=1)  # Trough times per spike
    #peak_times = np.argmax(interpolated_wfs, axis=1)  # Peak times per spike
    # peak times must occur after trough times
    peak_times = np.array([np.argmax(wf[trough_idx:]) + trough_idx for wf, trough_idx in zip(interpolated_wfs, trough_times)])

    horizontal_smearing = {
        "trough_time_std_ms": np.std(trough_times) * time_conversion_factor,
        "peak_time_std_ms": np.std(peak_times) * time_conversion_factor,
        "trough_to_peak_jitter_ms": np.std((peak_times - trough_times)) * time_conversion_factor,
    }
    
    # Find key points (convert to ms)
    trough_idx = np.argmin(weighted_avg_wf)
    peak_idx = np.argmax(weighted_avg_wf[trough_idx:]) + trough_idx
    trough_time_ms = trough_idx * time_conversion_factor
    peak_time_ms = peak_idx * time_conversion_factor
    
    # **Zero-Crossings & Phase Markers**
    zero_crossings = np.where(np.diff(np.sign(weighted_avg_wf)))[0]

    #spike_start_idx = zero_crossings[0] if len(zero_crossings) > 0 else 0
    ap_start_idx = zero_crossings[zero_crossings < trough_idx][-1] if len(zero_crossings[zero_crossings < trough_idx]) > 0 else 0
    ap_end_idx = zero_crossings[zero_crossings > trough_idx][0] if len(zero_crossings[zero_crossings > trough_idx]) > 0 else len(weighted_avg_wf) - 1
    refractory_end_idx = zero_crossings[zero_crossings > ap_end_idx][0] if len(zero_crossings[zero_crossings > ap_end_idx]) > 0 else len(weighted_avg_wf) - 1
    #spike_end_idx = zero_crossings[-1] if len(zero_crossings) > 0 else num_interp_samples - 1
    
    flag_after_plot = False
    if ap_start_idx == 0 or refractory_end_idx == len(weighted_avg_wf) - 1:
        print('Warning: AP start or refractory end is at the beginning or end of the waveform.')
        flag_after_plot = True

    # Convert to ms
    #spike_start_ms = spike_start_idx * time_conversion_factor
    #spike_end_ms = spike_end_idx * time_conversion_factor
    ap_start_ms = ap_start_idx * time_conversion_factor
    ap_end_ms = ap_end_idx * time_conversion_factor
    refractory_end_ms = refractory_end_idx * time_conversion_factor

    # **Amplitude Metrics**
    peak_amplitude = weighted_avg_wf[peak_idx]
    trough_amplitude = weighted_avg_wf[trough_idx]
    peak_to_trough_amplitude = np.abs(peak_amplitude) + np.abs(trough_amplitude)

    # **Temporal Metrics**
    peak_to_trough_time_ms = (peak_idx - trough_idx) * time_conversion_factor
    ap_phase_duration_ms = (ap_end_idx - ap_start_idx) * time_conversion_factor
    refractory_phase_duration_ms = (refractory_end_idx - ap_end_idx) * time_conversion_factor
    #total_spike_duration_ms = (spike_end_idx - spike_start_idx) * time_conversion_factor

    # Find half-max points
    half_amplitude = (np.abs(peak_amplitude) + np.abs(trough_amplitude)) / 2
    half_max_value = trough_amplitude + half_amplitude
    half_max_idxs = np.where(weighted_avg_wf <= half_max_value)[0]

    if len(half_max_idxs) > 1:
        half_max_start_ms = half_max_idxs[0] * time_conversion_factor
        half_max_end_ms = half_max_idxs[-1] * time_conversion_factor
        half_max_width_ms = (half_max_end_ms - half_max_start_ms)
    else:
        half_max_start_ms, half_max_end_ms, half_max_width_ms = np.nan, np.nan, np.nan

    #total_spike_duration = (np.where(weighted_avg_wf > 0)[0][-1] - np.where(weighted_avg_wf < 0)[0][0]) / 1000 * sampling_duration

    # **Slope-Based Metrics**
    max_depolarization_slope = np.max(np.diff(weighted_avg_wf[:trough_idx])) * (1000 / sampling_duration)  # uV/ms
    max_repolarization_slope = np.min(np.diff(weighted_avg_wf[trough_idx:peak_idx])) * (1000 / sampling_duration)  # uV/ms
    slope_ratio = max_depolarization_slope / max_repolarization_slope

    # **Waveform Asymmetry**
    # trough_to_peak_ratio = peak_to_trough_time / total_spike_duration
    # waveform_asymmetry_index = (peak_to_trough_time - (total_spike_duration - peak_to_trough_time)) / total_spike_duration
    #trough_to_peak_ratio = peak_to_trough_time_ms / total_spike_duration_ms
    #waveform_asymmetry_index = (peak_to_trough_time_ms - (total_spike_duration_ms - peak_to_trough_time_ms)) / total_spike_duration_ms
    wf_duration = sampling_duration
    trough_to_peak_ratio = peak_to_trough_time_ms / wf_duration
    waveform_asymmetry_index = (peak_to_trough_time_ms - (wf_duration - peak_to_trough_time_ms)) / wf_duration

    # **Energy Metrics**
    ap_phase_power = np.mean(weighted_avg_wf[:trough_idx] ** 2)  # µV²
    refractory_phase_power = np.mean(weighted_avg_wf[trough_idx:] ** 2)  # µV²
    total_spike_power = ap_phase_power + refractory_phase_power  # µV²

    # if max_depolarization_slope < 150:
    #     print("Warning: Low depolarization slope.")
    #     # remove from interpolated slope
    #     interpolated_wfs = interpolated_wfs[:-1]
        
    
    # **Waveform Shape**
    waveform_skewness = np.mean((weighted_avg_wf - np.mean(weighted_avg_wf))**3) / np.std(weighted_avg_wf)**3
    waveform_kurtosis = np.mean((weighted_avg_wf - np.mean(weighted_avg_wf))**4) / np.std(weighted_avg_wf)**4

    # **Plot Waveform for Debugging**
    if plot_wf:
               
        plt.figure(figsize=(8, 6))

        # Convert x-axis to time in ms
        time_axis = np.linspace(0, sampling_duration, num_interp_samples)

        # **Plot Waveforms + Vertical Smearing Envelope**
        plt.subplot(2, 1, 1)
        #plt.fill_between(time_axis, weighted_avg_wf - vertical_std, weighted_avg_wf + vertical_std, color='red', alpha=0.3)

        #plt.plot(time_axis, mean_wf, color='black', linewidth=1, label="mean_wf")
        for wf in interpolated_wfs:
            plt.plot(time_axis, wf, color='gray', 
                     alpha=0.75, 
                     linewidth=0.25)
        plt.plot(time_axis, weighted_avg_wf, color='red', 
                 #linestyle='--', 
                 linewidth=1.0, label="weighted mean")
        #plt.xlabel("Time (ms)")
        #plt.ylabel("Amplitude (µV)")
        plt.fill_between(time_axis, weighted_avg_wf - vertical_std, weighted_avg_wf + vertical_std, color='red', alpha=0.5,
                         label="std", linestyle='--',  linewidth=1.0                         
                         )

        plt.ylabel("uV")
        #plt.title(f"Unit {unit}: Waveforms & Vertical Smearing")
        plt.title(f"Unit {unit}")
        plt.legend(fontsize='small')
        
        # plot zero line dotted
        plt.axhline(0, color="black", linestyle="--", alpha=0.5, linewidth=0.5)
        
        # **Adjust X-axis Limits so they are the same for both subplots**
        plt.xlim(0, max(time_axis))        
        #remove time ticks
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

        # **Plot Histogram of Trough & Peak Timing**
        plt.subplot(2, 1, 2)
        
        # format bins
        #bins = np.linspace(0, max(trough_times.max(), peak_times.max()) * time_conversion_factor, 50)
        bins = np.linspace(0, max(time_axis), 100)
        
        plt.hist(trough_times * time_conversion_factor, bins=bins, color="blue", alpha=0.6, label="trough timing")
        plt.hist(peak_times * time_conversion_factor, bins=bins, color="green", alpha=0.6, label="peak timing")
        plt.xlabel("time (ms)")
        plt.ylabel("spikes")
        #plt.title("Trough & Peak Timing Variability (Horizontal Smearing)")
        plt.legend(fontsize='small')
        
        # **Adjust X-axis Limits so they are the same for both subplots**
        plt.xlim(0, max(time_axis))

        plt.tight_layout()
        
        if save_fig:
            #plt.savefig(f'{str(unit)}_{fig_name}', dpi=300)
            fig_dir = os.path.dirname(fig_name)
            if not os.path.exists(fig_dir):
                os.makedirs(fig_dir)
            png_path = fig_name
            pdf_path = fig_name.replace('.png', '.pdf')
            plt.savefig(png_path, dpi=300)
            plt.savefig(pdf_path)
            print(f"Figure saved to {fig_name}")
            print(f"Figure saved to {pdf_path}")
            #plt.savefig(fig_name, dpi=300)
            #print(f"Figure saved to {fig_name}")
        else:
            plt.show()
            
        plt.close()

    # Store in dictionary
    wf_metrics = {
        'weighted_avg_wf': weighted_avg_wf,
        'bio_variability_metrics': {
            # NOTE: i have no idea if this is good practice or not - mostly only wanting single values to use for classification
            'weighted_variability': weighted_variability,
            'mean_variance': np.mean(weighted_variability),
            'std_variance': np.std(weighted_variability),
            'cov_variance': np.std(weighted_variability) / np.mean(weighted_variability),
            'weighted_coefficient_of_variation': weighted_coefficient_of_variation,
            'mean_cv': np.mean(weighted_coefficient_of_variation),
            'std_cv': np.std(weighted_coefficient_of_variation),
            'cov_cv': np.std(weighted_coefficient_of_variation) / np.mean(weighted_coefficient_of_variation),
            },
        'amplitude_metrics': {
            'peak_amplitude': peak_amplitude,
            'trough_amplitude': trough_amplitude,
            'peak_to_trough_amplitude': peak_to_trough_amplitude,
        },
        'temporal_metrics': {
            'trough_time_ms': trough_time_ms,
            'peak_time_ms': peak_time_ms,
            'peak_to_trough_time_ms': peak_to_trough_time_ms,
            'ap_phase_duration_ms': ap_phase_duration_ms,
            'ap_start_ms': ap_start_ms,
            'ap_end_ms': ap_end_ms,
            'refractory_phase_duration_ms': refractory_phase_duration_ms,
            'refractory_end_ms': refractory_end_ms,
            'spike_width_half_max_ms': half_max_width_ms,
            #'total_spike_duration_ms': total_spike_duration_ms,
        },
        'slope_metrics': {
            'max_depolarization_slope': max_depolarization_slope,
            'max_repolarization_slope': max_repolarization_slope,
            'slope_ratio': slope_ratio,
        },
        'waveform_asymmetry': {
            'trough_to_peak_ratio': trough_to_peak_ratio,
            'waveform_asymmetry_index': waveform_asymmetry_index,
        },
        'energy_metrics': {
            'ap_phase_power_uv2': ap_phase_power,
            'refractory_phase_power_uv2': refractory_phase_power,
            'total_spike_power_uv2': total_spike_power,
        },
        'waveform_shape': {
            'waveform_skewness': waveform_skewness,
            'waveform_kurtosis': waveform_kurtosis,
        },
        'vertical_smearing': vertical_smearing,
        'horizontal_smearing': horizontal_smearing,
    }
    return wf_metrics

def plot_single_wf_and_half_ap(wf_i, trough_idx, peak_idx, first_half_idx, second_half_idx, half_amplitude, ap_half_width, num_samples, num_interp_samples, sampling_duration, sigma):

    # to visualize what I'm doing, plot the waveform and the half-amplitude points
    plt.plot(wf_i, color='black')  # Average waveform
    plt.plot([first_half_idx, second_half_idx], [wf_i[first_half_idx], wf_i[second_half_idx]], 'ro')  # Half-amplitude points
    plt.plot([trough_idx, peak_idx], [wf_i[trough_idx], wf_i[peak_idx]], 'go')  # Peak and trough
    # plot horizontal line demonstrating ap_half_width
    plt.plot([first_half_idx, second_half_idx], [wf_i[trough_idx] + half_amplitude, wf_i[trough_idx] + half_amplitude], 'b-')
    # plot vertical lines on ends of horizontal line to connect to half-amplitude points
    plt.plot([first_half_idx, first_half_idx], [wf_i[trough_idx] + half_amplitude, wf_i[first_half_idx]], 'b--')
    plt.plot([second_half_idx, second_half_idx], [wf_i[trough_idx] + half_amplitude, wf_i[second_half_idx]], 'b--')
    
    # update x axis ticks following interpolation and sampling rate info
    xticks = np.linspace(0, len(wf_i)-1, 5)
    xticklabels = np.round(np.linspace(0, sampling_duration, 5), 2)
    plt.xticks(xticks, xticklabels)
    
    # label axes
    plt.xlabel('Time (ms)')
    plt.ylabel('Amplitude (uV)')
    
    # note interpolation and smoothing params on the plot - also include legend.
    plt.legend(['Waveform', 'Half-amplitude points', 'Peak and trough', 'Half-amplitude line'])
    #plt.text(0, 0, f'Interpolated from {num_samples} to {num_interp_samples} samples\nSmoothed with Gaussian filter (sigma = {sigma})')
    # move text just a above legend
    # legend_loc = plt.gca().get_legend().get_window_extent().get_points()
    # legend_y = legend_loc[1][1]
    # legend_x = legend_loc[1][0]
    # plt.text(legend_x, legend_y + 0.1, f'AP Half-width: {ap_half_width} ms')
    # plt.text(legend_x, legend_y + 0.2, f'Interpolated from {num_samples} to {num_interp_samples} samples\nSmoothed with Gaussian filter (sigma = {sigma})')
    
    # tight layout
    plt.tight_layout()            
    
    # make some space at the bottom of the plot for the text
    plt.subplots_adjust(bottom=0.2)
    
    # find min y point in plt
    #plt.get_current_fig_manager().
    # # aw 2025-02-12 00:22:03 #HACK: there's a smarter way to do this...
    
    # get y coord for text by putting it below the trough
    trough_y = wf_i[trough_idx]
    plt.text(0, trough_y - 30, f'AP Half-width: {ap_half_width} ms')
    plt.text(0, trough_y - 50, f'Interpolated from {num_samples} to {num_interp_samples} samples\nSmoothed with Gaussian filter (sigma = {sigma})')

    plt.savefig('waveform.png')
    plt.close()            

def compute_wf_metrics_too_detailed(best_channel_waveforms, sampling_rate):
    '''Get waveform metrics for a single unit.'''
    unit_wfs = best_channel_waveforms
    half_aps = []
    interpolated_wfs = []
    metrics = {}
    for i in range(unit_wfs.shape[0]):
        metrics[i] = {}
        try:
            # NOTE: We're analyzing extracellular action potentials (APs) here.
            '''
            Metric	        Intracellular AP	    Extracellular AP
            Depolarization	Positive Deflection	    Negative Deflection
            Repolarization	Negative Deflection	    Positive Deflection (sometimes weak)
            Signal Shape	Large (tens of mV)	    Small (~50-500 µV)
            Duration	    Longer (~1-2 ms)	    Shorter (~0.5-1 ms)
            '''
            
            # aw 2025-02-12 11:15:53 - to debug I'm going to plot components of this analysis as I do them - I will consolidate plotting later
            plot_individual_metrics = False
            
            # initialize metrics
            wf_i = unit_wfs[i, :]
            num_samples = len(wf_i)
            sampling_duration = num_samples / sampling_rate * 1000 #ms
            metrics[i]['wf'] = wf_i
            metrics[i]['num_samples'] = num_samples
            metrics[i]['sampling_duration_ms'] = sampling_duration        
            
            #interpolate and smooth the signal for better analysis/metrics
            sigma = 9
            x = np.arange(len(wf_i))
            x_new = np.linspace(0, len(wf_i)-1, 1000)
            wf_int = np.interp(x_new, x, wf_i)
            num_interp_samples = len(wf_int)
            wf_int = gaussian_filter1d(wf_int, sigma)
            metrics[i]['wf_int'] = wf_int
            metrics[i]['num_interp_samples'] = num_interp_samples
            metrics[i]['sigma'] = sigma   
            # interpolated_wfs.append(wf_int) # store interpolated waveforms
            
            # init plot for debugging
            #plt.plot(wf_int, color='black')  # current wf
            #xticks = np.linspace(0, len(wf_int)-1, 5)
            #xticklabels = np.round(np.linspace(0, sampling_duration, 5), 2)
            #plt.xticks(xticks, xticklabels)
            #plt.xlabel('Time (ms)')
            #plt.ylabel('Amplitude (uV)')
            #plt.tight_layout()
            #plt.savefig('waveform.png')
            
            # some quick and dirty wf exclusion criteria
            trough_idx = np.argmin(wf_int)
            if trough_idx == 0 or trough_idx == len(wf_int)-1:
                metrics[i] = {} # clear current metrics
                print('Warning: Trough is at the beginning or end of the waveform.')
                print('I think this means were actually looking at a spike that was cut off by the edge of the recording.')
                print('Plausibly, between spikes.')
                metrics[i]['exclude_reason'] = 'Trough cut off by edge of recording'
                metrics[i]['exclude'] = True
                continue
            
            # add to interp if not excluded
            interpolated_wfs.append(wf_int) # store interpolated waveforms
                        
            ## Amplitude metrics            
            # Peak-to-trough amplitude
            def get_amplitude(wf):
                # find min and max points on the waveform - max must follow min
                trough_idx = np.argmin(wf)
                peak_idx = np.argmax(wf[trough_idx:]) + trough_idx
                #amp = wf[peak_idx] - wf[trough_idx]
                amp = np.abs(wf[peak_idx]) + np.abs(wf[trough_idx]) # makes more sense to get absolute value here right...?
                
                if plot_individual_metrics:
                    # plot min and max points on the waveform and plot a vertical dotted line showing amplitude being measured
                    # also plot a horizontal line connecting trough to bottom of vertical line - for visual clarity
                    ##plt.plot(wf, color='black')  # current wf
                    plt.plot(peak_idx, wf[peak_idx], 'ro')  # Peak
                    plt.plot(trough_idx, wf[trough_idx], 'ro')  # Trough
                    #plt.plot([peak_idx, peak_idx], [wf[trough_idx], wf[peak_idx]], 'b--')  # Amplitude
                    #plt.plot([trough_idx, peak_idx], [wf[trough_idx], wf[trough_idx]], 'b--')  # Amplitude
                    #plt.savefig('waveform.png')
                    
                    
                    # plt.plot(wf, color='black')  # current wf
                    # plt.plot(np.argmax(wf), np.max(wf), 'ro')  # Peak
                    # plt.plot(np.argmin(wf), np.min(wf), 'go')  # Trough
                    # plt.plot([np.argmax(wf), np.argmax(wf)], [np.min(wf), np.max(wf)], 'b--')  # Amplitude
                    # plt.savefig('waveform.png')
                
                return amp
            amplitude = get_amplitude(wf_int)
            metrics[i]['amplitude'] = amplitude
            
            # Trough Depth
            def get_trough_depth(wf):
                return np.min(wf)
            trough_depth = get_trough_depth(wf_int)
            metrics[i]['trough_depth'] = trough_depth
            
            # Peak Height
            def get_peak_height(wf):
                trough_idx = np.argmin(wf)
                peak_idx = np.argmax(wf[trough_idx:]) + trough_idx
                return wf[peak_idx]
                #return np.max(wf)
            peak_height = get_peak_height(wf_int)
            metrics[i]['peak_height'] = peak_height
            
            # Total spike width, pre-hyper, AP phase duration, and refractory phase duration (post-hyper)
            def get_phases(wf):
                """
                Computes the total spike width for extracellular waveforms based on:
                - Identifying the **action potential phase** (from the first zero-crossing before the trough to the first zero-crossing after the trough).
                - Identifying the **refractory phase** (from the end of the AP phase to the next zero-crossing).
                - Summing both durations to get **total spike duration**.

                Parameters:
                - wf: numpy array, the waveform.

                Returns:
                - Total spike duration in samples.
                - AP phase duration in samples.
                - Refractory phase duration in samples.
                """
                # Find the trough (negative peak) - the AP event
                trough_idx = np.argmin(wf)

                # Find zero-crossings (where the waveform crosses zero or baseline)
                zero_crossings = np.where(np.diff(np.sign(wf)))[0]

                # Find the first zero-crossing **before** the trough (start of AP phase)
                ap_start_idx = zero_crossings[zero_crossings < trough_idx][-1] if len(zero_crossings[zero_crossings < trough_idx]) > 0 else 0
                if ap_start_idx == 0:
                    print('Warning: No zero-crossing before trough')

                # Find the first zero-crossing **after** the trough (end of AP phase)
                ap_end_idx = zero_crossings[zero_crossings > trough_idx][0] if len(zero_crossings[zero_crossings > trough_idx]) > 0 else len(wf) - 1

                # Find the first zero-crossing **after the AP phase** (end of refractory phase)
                refractory_end_idx = zero_crossings[zero_crossings > ap_end_idx][0] if len(zero_crossings[zero_crossings > ap_end_idx]) > 0 else len(wf) - 1
                
                # attempt to find pre-hyperpolarization phase start if it exists
                pre_hyper_start_idx = zero_crossings[zero_crossings < ap_start_idx][-1] if len(zero_crossings[zero_crossings < ap_start_idx]) > 0 else np.nan

                # Compute durations
                if pre_hyper_start_idx>0: pre_hyper_duration = ap_start_idx - pre_hyper_start_idx
                else: pre_hyper_duration = np.nan
                ap_phase_duration = ap_end_idx - ap_start_idx  # Samples
                refractory_phase_duration = refractory_end_idx - ap_end_idx  # Samples
                if pre_hyper_start_idx>0: total_spike_duration = ap_phase_duration + refractory_phase_duration + pre_hyper_duration  # Total duration
                else: total_spike_duration = ap_phase_duration + refractory_phase_duration  # Total duration
                
                #Biphasic or Triphasic
                if pre_hyper_start_idx>0: num_phases = 3
                else: num_phases = 2
                
                # spike_start and spike_end
                if pre_hyper_start_idx>0: spike_start_idx = pre_hyper_start_idx
                else: spike_start_idx = ap_start_idx
                spike_end_idx = refractory_end_idx
                
                # quality check
                try:
                    assert pre_hyper_start_idx <= ap_start_idx or np.isnan(pre_hyper_start_idx), f"Pre-hyper start index is after AP start index: {pre_hyper_start_idx} > {ap_start_idx}"
                    assert ap_start_idx <= ap_end_idx, f"AP start index is after AP end index: {ap_start_idx} > {ap_end_idx}"
                    assert ap_end_idx <= refractory_end_idx or np.isnan(refractory_end_idx), f"AP end index is after refractory end index: {ap_end_idx} > {refractory_end_idx}"
                except AssertionError as e:
                    print(e)
                    
                if trough_idx >= ap_end_idx:
                    print('Warning: Trough is at or after AP end - cannot compute max repolarization slope')
                    print("This warning shouldn't even be possible, but here we are...")
                    # plot this one because wtf
                    #plot spike start and end in green and red respectively
                    plt.axvline(spike_start_idx, color='green', linestyle='--', alpha=0.5, label="Spike Start")
                    plt.axvline(spike_end_idx, color='red', linestyle='--', alpha=0.5, label="Spike End")
                    
                    plt.axvline(ap_start_idx, color='gray', linestyle='--', alpha=0.5, label="AP Start")
                    plt.axvline(ap_end_idx, color='gray', linestyle='--', alpha=0.5, label="AP End")
                    plt.axvline(refractory_end_idx, color='gray', linestyle='--', alpha=0.5, label="Refractory End")
                    #plt.legend()
                    plt.plot(wf, color='black')  # current wf
                    plt.savefig('waveform.png')
                    plt.close()
                    return np.nan                

                # Plot key points
                if plot_individual_metrics:
                    #plt.plot(wf, label="Waveform")
                    # plot zero line
                    plt.plot([0, len(wf)], [0, 0], 'k--')  # Zero line
                    if pre_hyper_start_idx>0: plt.axvline(pre_hyper_start_idx, color='gray', linestyle='--', alpha=0.5, label="Pre-Hyper Start")
                    plt.axvline(ap_start_idx, color='gray', linestyle='--', alpha=0.5, label="AP Start")
                    plt.axvline(ap_end_idx, color='gray', linestyle='--', alpha=0.5, label="AP End")
                    plt.axvline(refractory_end_idx, color='gray', linestyle='--', alpha=0.5, label="Refractory End")
                    #plt.legend()
                    #plt.savefig('waveform.png')
                
                phases = {
                    'total_spike_duration_ms': total_spike_duration / num_interp_samples * sampling_duration,
                    'spike_start_ms': spike_start_idx / num_interp_samples * sampling_duration,
                    'spike_end_ms': spike_end_idx / num_interp_samples * sampling_duration,
                    
                    'pre_hyper_duration_ms': pre_hyper_duration / num_interp_samples * sampling_duration,
                    'pre_hyper_start_ms': pre_hyper_start_idx / num_interp_samples * sampling_duration,
                    
                    'ap_phase_duration_ms': ap_phase_duration / num_interp_samples * sampling_duration,
                    'ap_phase_start_ms': ap_start_idx / num_interp_samples * sampling_duration,
                    'ap_phase_end_ms': ap_end_idx / num_interp_samples * sampling_duration,
                    
                    'refractory_phase_duration_ms': refractory_phase_duration / num_interp_samples * sampling_duration,
                    'refractory_phase_end_ms': refractory_end_idx / num_interp_samples * sampling_duration,
                    
                    'num_phases': num_phases,
                    }

                return phases
            phases = get_phases(wf_int)
            metrics[i]['phase_metrics'] = phases

            # trough-to-peak duration
            def get_trough_to_peak_duration(wf):
                trough_idx = np.argmin(wf)
                peak_idx = np.argmax(wf[trough_idx:]) + trough_idx
                return peak_idx - trough_idx
                #return np.argmax(wf) - np.argmin(wf)
            trough_to_peak_duration = get_trough_to_peak_duration(wf_int) / num_interp_samples * sampling_duration #ms
            metrics[i]['trough_to_peak_duration_ms'] = trough_to_peak_duration
            
            # spike width at half-maximum (FWHM)
            def get_spike_width_at_half_max(wf):
                """
                Computes the full-width at half-maximum (FWHM) for extracellular spikes.

                Since extracellular spikes are negative-first, we measure the width 
                where the waveform crosses half of the peak-to-trough amplitude.

                Returns:
                - FWHM in samples
                """
                # trough = np.min(wf)
                # peak = np.max(wf)
                trough_idx = np.argmin(wf)
                peak_idx = np.argmax(wf[trough_idx:]) + trough_idx
                peak = wf[peak_idx]
                trough = wf[trough_idx]
                half_amplitude = abs(peak)+abs(trough) / 2  # Half of the peak-to-trough amplitude

                # For extracellular spikes, half-max is defined from the **trough upwards**
                half_max = trough + half_amplitude

                # Find where waveform recovers above half-max after the trough
                half_max_idx = np.where(wf <= half_max)[0]

                if len(half_max_idx) < 2:
                    return np.nan  # Not enough crossings
                
                if plot_individual_metrics:
                    # plt plot half max line and points where it crosses the waveform
                    plt.plot(wf, color='black')  # current wf
                    #plt.plot(half_max_idx, wf[half_max_idx], 'ro')  # Half-max points
                    plt.plot(half_max_idx[-1], half_max, 'go')  # Half-max points
                    plt.plot(half_max_idx[0], half_max, 'go')  # Half-max points
                    plt.plot([half_max_idx[0], half_max_idx[-1]], [half_max, half_max], 'g--')  # Half-max line
                    #plt.savefig('waveform.png')
                
                width_samples = half_max_idx[-1] - half_max_idx[0]  # Width in samples
                half_max_amp = half_max                

                return width_samples, half_max_amp
            #spike_width_at_half_max = get_spike_width_at_half_max(wf_int) / num_interp_samples * sampling_duration  # ms
            spike_width_at_half_max, half_max_amp = get_spike_width_at_half_max(wf_int)
            spike_width_at_half_max = spike_width_at_half_max / num_interp_samples * sampling_duration  # ms
            metrics[i]['spike_width_at_half_max_ms'] = spike_width_at_half_max
            metrics[i]['half_max_amp'] = half_max_amp

            # # Trough-to-half Recovery Time
            # def get_trough_to_half_recovery_time(wf):
            #     """
            #     Computes the time for the waveform to recover from trough 
            #     to half-amplitude level (extracellular waveform).

            #     Returns:
            #     - Recovery time in samples
            #     """
            #     trough_idx = np.argmin(wf)  # Find the trough
            #     half_max = get_trough_depth(wf) + get_half_amplitude(wf)

            #     # Find where the waveform recovers above half-max **after the trough**
            #     half_max_idx = np.where(wf[trough_idx:] >= half_max)[0]

            #     if len(half_max_idx) == 0:
            #         return np.nan  # No valid recovery

            #     return half_max_idx[0]  # Time in samples after the trough
            # trough_to_half_recovery_time = get_trough_to_half_recovery_time(wf_int) / num_interp_samples * sampling_duration  # ms
            # metrics[i]['trough_to_half_recovery_time_ms'] = trough_to_half_recovery_time
            
            # Trough-to-Peak Ratio: Ratio of the time from the trough to the peak relative to the total spike width.
            total_spike_width = phases['total_spike_duration_ms']
            def get_trough_to_peak_ratio(wf, total_spike_width):
                """
                Computes the ratio of trough-to-peak duration to total spike width.
                Handles cases where no peak is detected after the trough.

                Returns:
                - Ratio (float) or NaN if no peak exists.
                """
                trough_idx = np.argmin(wf)  # Find trough
                peak_idx = np.argmax(wf[trough_idx:]) + trough_idx  # Find peak after trough
                
                trough_ms = trough_idx / num_interp_samples * sampling_duration
                peak_ms = peak_idx / num_interp_samples * sampling_duration

                if peak_idx == trough_idx:  
                    return np.nan  # No peak detected, return NaN

                return (peak_ms - trough_ms) / total_spike_width  # Ratio ms/ms
            trough_to_peak_ratio = get_trough_to_peak_ratio(wf_int, total_spike_width)
            metrics[i]['trough_to_peak_ratio'] = trough_to_peak_ratio

            # Waveform Asymmetry Index
            def get_waveform_asymmetry_index(wf, trough_to_peak_duration_ms, total_spike_width_ms, sampling_rate):
                """
                Computes the Waveform Asymmetry Index (WAI) in **milliseconds**.

                Parameters:
                - wf: numpy array, the waveform data.
                - trough_to_peak_duration: int, time from trough to peak in samples.
                - total_spike_width: int, total spike width in samples.
                - sampling_rate: float, sampling rate in Hz.

                Returns:
                - Waveform Asymmetry Index (WAI) in ms.
                """
                if np.isnan(trough_to_peak_duration) or np.isnan(total_spike_width) or total_spike_width == 0:
                    return np.nan  # Avoid divide by zero errors

                # Convert from samples to milliseconds
                #trough_to_peak_duration_ms = (trough_to_peak_duration / sampling_rate) * 1000
                #total_spike_width_ms = (total_spike_width / sampling_rate) * 1000

                # Compute WAI in ms
                WAI = (trough_to_peak_duration_ms - (total_spike_width_ms - trough_to_peak_duration_ms)) / total_spike_width_ms

                return WAI
            waveform_asymmetry_index = get_waveform_asymmetry_index(wf_int, trough_to_peak_duration, total_spike_width, sampling_rate)
            metrics[i]['waveform_asymmetry_index'] = waveform_asymmetry_index

            ## Slope-based metrics
            # Maximum Depolarization Slope - Steepest slope on the way to trough
            def get_max_depolarization_slope(wf, phases):
                ap_start_idx = int(phases['ap_phase_start_ms'] * num_interp_samples / sampling_duration) # convert ms to samples
                
                trough_idx = np.argmin(wf)
                
                if trough_idx <= ap_start_idx:
                    print('Warning: Trough is at or before AP start - cannot compute max depolarization slope')
                    return np.nan
                
                
                slope = np.diff(wf[ap_start_idx:trough_idx]) # units: uV/sample
                #convert to uV/ms
                slope = slope * num_interp_samples / sampling_duration # uV/sample * samples/ms = uV/ms
                if len(slope) == 0:
                    return np.nan
                return np.max(slope)
            max_depolarization_slope = get_max_depolarization_slope(wf_int, phases)
            metrics[i]['max_depolarization_slope'] = max_depolarization_slope
            
            # Maximum Repolarization Slope - Steepest slope on the way to peak
            def get_max_repolarization_slope(wf, phases):
                ap_end_idx = int(phases['ap_phase_end_ms'] * num_interp_samples / sampling_duration) # convert ms to samples
                
                trough_idx = np.argmin(wf)
                
                if trough_idx >= ap_end_idx:
                    print('Warning: Trough is at or after AP end - cannot compute max repolarization slope')
                    print("This warning shouldn't even be possible, but here we are...")
                    # plot this one because wtf
                    plt.plot(wf, color='black')  # current wf
                    plt.savefig('waveform.png')
                    plt.close()
                    return np.nan
                
                slope = np.diff(wf[trough_idx:ap_end_idx]) # units: uV/sample
                #convert to uV/ms
                slope = slope * num_interp_samples / sampling_duration # uV/sample * samples/ms = uV/ms
                return np.min(slope)
            max_repolarization_slope = get_max_repolarization_slope(wf_int, phases)
            metrics[i]['max_repolarization_slope'] = max_repolarization_slope
            
            # Slope Ratio - Ratio of maximum depolarization slope to maximum repolarization slope
            def get_slope_ratio(max_depolarization_slope, max_repolarization_slope):
                return max_depolarization_slope / max_repolarization_slope
            slope_ratio = get_slope_ratio(max_depolarization_slope, max_repolarization_slope) # ms/ms i.e. unitless
            metrics[i]['slope_ratio'] = slope_ratio
            
            # Maximum voltage acceleration - useful for detecting fast spiking neurons
            def get_max_voltage_acceleration(wf):
                """
                Computes the maximum voltage acceleration (change in slope over time).

                Parameters:
                - wf: numpy array, the waveform data.
                - sampling_rate: float, sampling rate in Hz.

                Returns:
                - Maximum acceleration in μV/ms².
                """
                trough_idx = np.argmin(wf)  # Find trough
                
                if trough_idx == 0:
                    print('Warning: Trough is at the beginning of the waveform - cannot compute max voltage acceleration')
                    return np.nan
                
                acceleration = np.diff(np.diff(wf[:trough_idx]))  # Second derivative, units: μV/sample²
                # Convert from samples² to ms²
                acceleration = acceleration * (num_interp_samples ** 2) / (sampling_duration ** 2)  # dimensional analysis: μV/sample² * sample²/ms² = μV/ms²
                if len(acceleration) == 0:
                    return np.nan
                return np.max(acceleration)  # μV/ms²

                # Convert from samples² to ms²
                # time_scaling_factor = (sampling_rate ** 2) / (1000 ** 2)  # Converts from samples² to ms²
                # return np.max(acceleration) * time_scaling_factor  # μV/ms²
            max_voltage_acceleration = get_max_voltage_acceleration(wf_int)
            metrics[i]['max_voltage_acceleration'] = max_voltage_acceleration

            # Maximum voltage deceleration - useful for detecting fast spiking neurons
            def get_max_voltage_deceleration(wf):
                """
                Computes the maximum voltage deceleration (change in slope over time).

                Parameters:
                - wf: numpy array, the waveform data.
                - sampling_rate: float, sampling rate in Hz.

                Returns:
                - Maximum deceleration in μV/ms².
                """
                trough_idx = np.argmin(wf)
                
                # Find the peak after the trough
                peak_idx = np.argmax(wf[trough_idx:]) + trough_idx
                deceleration = np.diff(np.diff(wf[trough_idx:peak_idx]))  # Second derivative, units: μV/sample²
                # Convert from samples² to ms²
                deceleration = deceleration * (num_interp_samples ** 2) / (sampling_duration ** 2)  # dimensional analysis: μV/sample² * sample²/ms² = μV/ms²
                if len(deceleration) == 0:
                    return np.nan
                return np.min(deceleration)  # μV/ms²
                
                # Convert from samples² to ms²
                #time_scaling_factor = (sampling_rate ** 2) / (1000 ** 2)  # Converts from samples² to ms²                
                #return np.min(deceleration) * time_scaling_factor  # μV/ms²
            max_voltage_deceleration = get_max_voltage_deceleration(wf_int)
            metrics[i]['max_voltage_deceleration'] = max_voltage_deceleration
            
            ## Energe-Based and Spectral Metrics
            # calculate signal power as mean squared voltage for AP and refractory phases
            # NOTE: not technically power, but rather proportional to power. Need some metric of resistance to get power.
            def get_signal_power(wf, phases):
                """
                Computes the signal power (mean squared voltage, µV²) for the AP and Refractory phases.

                Steps:
                1. Identify the **AP phase** (first zero-crossing before the trough to first zero-crossing after trough).
                2. Identify the **Refractory phase** (from end of AP phase to next zero-crossing).
                3. Compute **power (µV²)** separately for each phase.

                Returns:
                - Pre-AP Phase Power (µV²) if exists
                - AP Phase Power (µV²)
                - Refractory Phase Power (µV²)
                - Total Spike Power (sum of both)
                """
                
                pre_ap_start_idx = int(phases['spike_start_ms'] * num_interp_samples / sampling_duration) # convert ms to samples
                ap_start_idx = int(phases['ap_phase_start_ms'] * num_interp_samples / sampling_duration) # convert ms to samples
                ap_end_idx = int(phases['ap_phase_end_ms'] * num_interp_samples / sampling_duration) # convert ms to samples
                refractory_end_idx = int(phases['refractory_phase_end_ms'] * num_interp_samples / sampling_duration) # convert ms to samples
                
                # compute power for each phase
                try:
                    pre_ap_phase_power = np.nanmean(wf[pre_ap_start_idx:ap_start_idx] ** 2)  # µV²
                    ap_phase_power = np.nanmean(wf[ap_start_idx:ap_end_idx] ** 2)  # µV²
                    refractory_phase_power = np.nanmean(wf[ap_end_idx:refractory_end_idx] ** 2)  # µV²
                    total_spike_power = np.nansum([pre_ap_phase_power, ap_phase_power, refractory_phase_power])  # µV²
                except:
                    traceback.print_exc()
                    print('Warning: Error computing power values')
                    # pre_ap_phase_power = np.nan
                    # ap_phase_power = np.nan
                    # refractory_phase_power = np.nan
                    # total_spike_power = np.nan
                
                # if any(np.isnan([pre_ap_phase_power, ap_phase_power, refractory_phase_power])):
                #     print('Warning: One or more power values are NaN')
                    
                if np.isnan(total_spike_power):
                    print('Critical: Total spike power is NaN')                

                return pre_ap_phase_power, ap_phase_power, refractory_phase_power, total_spike_power
            # ap_power, refractory_power, total_power = get_signal_power(wf_int, phases)
            pre_ap_power, ap_power, refractory_power, total_power = get_signal_power(wf_int, phases)
            metrics[i]['phase_metrics']['pre_ap_phase_power_uv2'] = pre_ap_power
            metrics[i]['phase_metrics']['ap_phase_power_uv2'] = ap_power
            metrics[i]['phase_metrics']['refractory_phase_power_uv2'] = refractory_power
            metrics[i]['phase_metrics']['total_spike_power_uv2'] = total_power
            
            # # Get the dominant frequency of the waveform
            # def get_dominant_frequency(wf, raw_sampling_rate, num_samples, num_interp_samples, phases):
            #     """
            #     Computes the dominant frequency of the waveform using power spectral density (PSD).
                
            #     Parameters:
            #     - wf: numpy array, interpolated waveform
            #     - raw_sampling_rate: float, original sampling rate of the raw data (Hz)
            #     - num_samples: int, number of samples in raw waveform
            #     - num_interp_samples: int, number of samples in interpolated waveform
                
            #     Returns:
            #     - Dominant frequency in Hz
            #     """
            #     #
            #     spike_start_idx = int(phases['spike_start_ms'] * num_interp_samples / sampling_duration) # convert ms to samples
            #     spike_end_idx = int(phases['spike_end_ms'] * num_interp_samples / sampling_duration) # convert ms to samples
            #     wf_crop = wf[spike_start_idx:spike_end_idx]
                
            #     # Compute the effective sampling rate after interpolation
            #     interp_factor = num_interp_samples / num_samples
            #     sampling_rate_interpolated = raw_sampling_rate * interp_factor

            #     # Compute the power spectral density (PSD)
            #     f, Pxx = signal.welch(wf_crop, fs=sampling_rate_interpolated, nperseg=len(wf_crop))

            #     return f[np.argmax(Pxx)]  # Frequency with max power
            # dominant_frequency = get_dominant_frequency(wf_int, sampling_rate, num_samples, num_interp_samples, phases)
            # metrics[i]['dominant_frequency'] = dominant_frequency
            
            # skewness of the waveform
            def get_waveform_skewness(wf, phases):
                spike_start_idx = int(phases['spike_start_ms'] * num_interp_samples / sampling_duration) # convert ms to samples
                spike_end_idx = int(phases['spike_end_ms'] * num_interp_samples / sampling_duration) # convert ms to samples
                
                wf_cropped = wf[spike_start_idx:spike_end_idx]
                return np.mean((wf_cropped - np.mean(wf_cropped))**3) / np.std(wf_cropped)**3
                #return np.mean((wf - np.mean(wf))**3) / np.std(wf)**3
            waveform_skewness = get_waveform_skewness(wf_int, phases)
            metrics[i]['waveform_skewness'] = waveform_skewness
            
            # kurtosis of the waveform
            def get_waveform_kurtosis(wf, phases):
                spike_start_idx = int(phases['spike_start_ms'] * num_interp_samples / sampling_duration)
                spike_end_idx = int(phases['spike_end_ms'] * num_interp_samples / sampling_duration)
                
                wf_cropped = wf[spike_start_idx:spike_end_idx]
                return np.mean((wf_cropped - np.mean(wf_cropped))**4) / np.std(wf_cropped)**4                
                #return np.mean((wf - np.mean(wf))**4) / np.std(wf)**4
            waveform_kurtosis = get_waveform_kurtosis(wf_int, phases)
            metrics[i]['waveform_kurtosis'] = waveform_kurtosis        
            
            # save all plotted metrics
            if plot_individual_metrics:
                # plot spike start and end in green and red respectively
                plt.axvline(phases['spike_start_ms'], color='green', linestyle='--', alpha=0.5, label="Spike Start")
                plt.axvline(phases['spike_end_ms'], color='red', linestyle='--', alpha=0.5, label="Spike End")                
                plt.plot(wf_int, color='black')  # current wf
                plt.tight_layout()
                plt.savefig('waveform.png')
            #print(f'Computed metrics for waveform {i}')
            
            # TODO: formaalize plots made here this, rn its just for debugging # aw 2025-02-12 09:40:36
            #plot_single_wf_and_half_ap(wf_i, trough_idx, peak_idx, first_half_idx, second_half_idx, half_amplitude, ap_half_width, num_samples, num_interp_samples, sampling_duration, sigma)
            
        except Exception as e:
            print(f'Error computing waveform metrics: {e}')
            traceback.print_exc()
            continue
            
    ## Collect data to compute stats
    errors = []
    
    ## General metrics
    amps = [] 
    troughs = []
    peaks = []
    
    ## Phase metrics
    total_spike_duration_ms = []
    spike_start_ms = []
    spike_end_ms = []
    pre_hyper_duration_ms = []
    pre_hyper_start_ms = []
    ap_phase_duration_ms = []
    ap_phase_start_ms = []
    ap_phase_end_ms = []
    refractory_phase_duration_ms = []
    refractory_phase_end_ms = []
    num_phases = []
    pre_ap_phase_power_uv2 = []
    ap_phase_power_uv2 = []
    refractory_phase_power_uv2 = []
    total_spike_power_uv2 = []
    
    ## advanced metrics
    trough_to_peak_durations = []
    spike_width_at_half_maxs = []
    half_max_amps = []
    trough_to_peak_ratios = []
    waveform_asymmetry_indices = []
    max_depolarization_slopes = []
    max_repolarization_slopes = []
    slope_ratios = []
    max_voltage_accelerations = []
    max_voltage_decelerations = []
    waveform_skewnesses = []
    waveform_kurtoses = []
    
    for wf_i, metric in metrics.items():
        
        ## General metrics
        try: amps.append(metric['amplitude'])
        except Exception as e: amps.append(np.nan); errors.append(e)
        try: troughs.append(metric['trough_depth'])
        except Exception as e: troughs.append(np.nan); errors.append(e)
        try: peaks.append(metric['peak_height'])
        except Exception as e: peaks.append(np.nan); errors.append(e)
        
        ## Phase metrics
        try: 
            phase_metric = metric['phase_metrics']
            try: total_spike_duration_ms.append(phase_metric['total_spike_duration_ms'])
            except Exception as e: total_spike_duration_ms.append(np.nan); errors.append(e)
            try: spike_start_ms.append(phase_metric['spike_start_ms'])
            except Exception as e: spike_start_ms.append(np.nan); errors.append(e)
            try: spike_end_ms.append(phase_metric['spike_end_ms'])
            except Exception as e: spike_end_ms.append(np.nan); errors.append(e)
            try: pre_hyper_duration_ms.append(phase_metric['pre_hyper_duration_ms'])
            except Exception as e: pre_hyper_duration_ms.append(np.nan); errors.append(e)
            try: pre_hyper_start_ms.append(phase_metric['pre_hyper_start_ms'])
            except Exception as e: pre_hyper_start_ms.append(np.nan); errors.append(e)
            try: ap_phase_duration_ms.append(phase_metric['ap_phase_duration_ms'])
            except Exception as e: ap_phase_duration_ms.append(np.nan); errors.append(e)
            try: ap_phase_start_ms.append(phase_metric['ap_phase_start_ms'])
            except Exception as e: ap_phase_start_ms.append(np.nan); errors.append(e)
            try: ap_phase_end_ms.append(phase_metric['ap_phase_end_ms'])
            except Exception as e: ap_phase_end_ms.append(np.nan); errors.append(e)
            try: refractory_phase_duration_ms.append(phase_metric['refractory_phase_duration_ms'])
            except Exception as e: refractory_phase_duration_ms.append(np.nan); errors.append(e)
            try: refractory_phase_end_ms.append(phase_metric['refractory_phase_end_ms'])
            except Exception as e: refractory_phase_end_ms.append(np.nan); errors.append(e)
            try: num_phases.append(phase_metric['num_phases'])
            except Exception as e: num_phases.append(np.nan); errors.append(e)
            try: pre_ap_phase_power_uv2.append(phase_metric['pre_ap_phase_power_uv2'])
            except Exception as e: pre_ap_phase_power_uv2.append(np.nan); errors.append(e)
            try: ap_phase_power_uv2.append(phase_metric['ap_phase_power_uv2'])
            except Exception as e: ap_phase_power_uv2.append(np.nan); errors.append(e)
            try: refractory_phase_power_uv2.append(phase_metric['refractory_phase_power_uv2'])
            except Exception as e: refractory_phase_power_uv2.append(np.nan); errors.append(e)
            try: total_spike_power_uv2.append(phase_metric['total_spike_power_uv2'])
            except Exception as e: total_spike_power_uv2.append(np.nan); errors.append(e)
        except Exception as e: print(e)
        
        ## advanced metrics
        try: trough_to_peak_durations.append(metric['trough_to_peak_duration_ms'])
        except Exception as e: trough_to_peak_durations.append(np.nan); errors.append(e)
        try: spike_width_at_half_maxs.append(metric['spike_width_at_half_max_ms'])
        except Exception as e: spike_width_at_half_maxs.append(np.nan); errors.append(e)
        try: half_max_amps.append(metric['half_max_amp'])
        except Exception as e: half_max_amps.append(np.nan); errors.append(e)
        try: trough_to_peak_ratios.append(metric['trough_to_peak_ratio'])
        except Exception as e: trough_to_peak_ratios.append(np.nan); errors.append(e)
        try: waveform_asymmetry_indices.append(metric['waveform_asymmetry_index'])
        except Exception as e: waveform_asymmetry_indices.append(np.nan); errors.append(e)
        try: max_depolarization_slopes.append(metric['max_depolarization_slope'])
        except Exception as e: max_depolarization_slopes.append(np.nan); errors.append(e)
        try: max_repolarization_slopes.append(metric['max_repolarization_slope'])
        except Exception as e: max_repolarization_slopes.append(np.nan); errors.append(e)
        try: slope_ratios.append(metric['slope_ratio'])
        except Exception as e: slope_ratios.append(np.nan); errors.append(e)
        try: max_voltage_accelerations.append(metric['max_voltage_acceleration'])
        except Exception as e: max_voltage_accelerations.append(np.nan); errors.append(e)
        try: max_voltage_decelerations.append(metric['max_voltage_deceleration'])
        except Exception as e: max_voltage_decelerations.append(np.nan); errors.append(e)
        try: waveform_skewnesses.append(metric['waveform_skewness'])
        except Exception as e: waveform_skewnesses.append(np.nan); errors.append(e)
        try: waveform_kurtoses.append(metric['waveform_kurtosis'])
        except Exception as e: waveform_kurtoses.append(np.nan); errors.append(e)        
        
        if len(errors) > 0:
            traceback.print_exc()
            print('An error occurred while collecting metrics')
            
    # Include/exclude
    included_wfs = [i for i, metric in metrics.items() if 'excluded' not in metric]
    excluded_wfs = [i for i, metric in metrics.items() if 'excluded' in metric]
    
    # Compute and store data and stats
    wf_metrics = {
        'unit_wfs': unit_wfs,
        'included_wfs': included_wfs,
        'excluded_wfs': excluded_wfs,
        'interpolated_wfs': interpolated_wfs,
        'num_wfs': len(interpolated_wfs),
        'note': (f'Interpolated from {num_samples} to {num_interp_samples} samples.'
                    f'Smoothed with Gaussian filter (sigma = {sigma}).'
                    f'Interpolated signal used to compute all metrics.'),
        'general_metrics': {
            'amps': {
                'data': amps,
                'mean': np.nanmean(amps),
                'std': np.nanstd(amps),
                'cov': np.nanstd(amps) / np.nanmean(amps),
                'median': np.nanmedian(amps),
                'min': np.nanmin(amps),
                'max': np.nanmax(amps),
            },
            'trough_depths': {
                'data': troughs,
                'mean': np.nanmean(troughs),
                'std': np.nanstd(troughs),
                'cov': np.nanstd(troughs) / np.nanmean(troughs),
                'median': np.nanmedian(troughs),
                'min': np.nanmin(troughs),
                'max': np.nanmax(troughs),
            },
            'peak_heights': {
                'data': peaks,
                'mean': np.nanmean(peaks),
                'std': np.nanstd(peaks),
                'cov': np.nanstd(peaks) / np.nanmean(peaks),
                'median': np.nanmedian(peaks),
                'min': np.nanmin(peaks),
                'max': np.nanmax(peaks),
            },
        },
        'phase_metrics': {
            'total_spike_duration_ms': {
                'data': total_spike_duration_ms,
                'mean': np.nanmean(total_spike_duration_ms),
                'std': np.nanstd(total_spike_duration_ms),
                'cov': np.nanstd(total_spike_duration_ms) / np.nanmean(total_spike_duration_ms),
                'median': np.nanmedian(total_spike_duration_ms),
                'min': np.nanmin(total_spike_duration_ms),
                'max': np.nanmax(total_spike_duration_ms),
            },
            'spike_start_ms': {
                'data': spike_start_ms,
                'mean': np.nanmean(spike_start_ms),
                'std': np.nanstd(spike_start_ms),
                'cov': np.nanstd(spike_start_ms) / np.nanmean(spike_start_ms),
                'median': np.nanmedian(spike_start_ms),
                'min': np.nanmin(spike_start_ms),
                'max': np.nanmax(spike_start_ms),
            },
            'spike_end_ms': {
                'data': spike_end_ms,
                'mean': np.nanmean(spike_end_ms),
                'std': np.nanstd(spike_end_ms),
                'cov': np.nanstd(spike_end_ms) / np.nanmean(spike_end_ms),
                'median': np.nanmedian(spike_end_ms),
                'min': np.nanmin(spike_end_ms),
                'max': np.nanmax(spike_end_ms),
            },
            'pre_hyper_duration_ms': {
                'data': pre_hyper_duration_ms,
                'mean': np.nanmean(pre_hyper_duration_ms),
                'std': np.nanstd(pre_hyper_duration_ms),
                'cov': np.nanstd(pre_hyper_duration_ms) / np.nanmean(pre_hyper_duration_ms),
                'median': np.nanmedian(pre_hyper_duration_ms),
                'min': np.nanmin(pre_hyper_duration_ms),
                'max': np.nanmax(pre_hyper_duration_ms),
            },
            'pre_hyper_start_ms': {
                'data': pre_hyper_start_ms,
                'mean': np.nanmean(pre_hyper_start_ms),
                'std': np.nanstd(pre_hyper_start_ms),
                'cov': np.nanstd(pre_hyper_start_ms) / np.nanmean(pre_hyper_start_ms),
                'median': np.nanmedian(pre_hyper_start_ms),
                'min': np.nanmin(pre_hyper_start_ms),
                'max': np.nanmax(pre_hyper_start_ms),
            },
            'ap_phase_duration_ms': {
                'data': ap_phase_duration_ms,
                'mean': np.nanmean(ap_phase_duration_ms),
                'std': np.nanstd(ap_phase_duration_ms),
                'cov': np.nanstd(ap_phase_duration_ms) / np.nanmean(ap_phase_duration_ms),
                'median': np.nanmedian(ap_phase_duration_ms),
                'min': np.nanmin(ap_phase_duration_ms),
                'max': np.nanmax(ap_phase_duration_ms),
            },
            'ap_phase_start_ms': {
                'data': ap_phase_start_ms,
                'mean': np.nanmean(ap_phase_start_ms),
                'std': np.nanstd(ap_phase_start_ms),
                'cov': np.nanstd(ap_phase_start_ms) / np.nanmean(ap_phase_start_ms),
                'median': np.nanmedian(ap_phase_start_ms),
                'min': np.nanmin(ap_phase_start_ms),
                'max': np.nanmax(ap_phase_start_ms),
            },
            'ap_phase_end_ms': {
                'data': ap_phase_end_ms,
                'mean': np.nanmean(ap_phase_end_ms),
                'std': np.nanstd(ap_phase_end_ms),
                'cov': np.nanstd(ap_phase_end_ms) / np.nanmean(ap_phase_end_ms),
                'median': np.nanmedian(ap_phase_end_ms),
                'min': np.nanmin(ap_phase_end_ms),
                'max': np.nanmax(ap_phase_end_ms),
            },
            'refractory_phase_duration_ms': {
                'data': refractory_phase_duration_ms,
                'mean': np.nanmean(refractory_phase_duration_ms),
                'std': np.nanstd(refractory_phase_duration_ms),
                'cov': np.nanstd(refractory_phase_duration_ms) / np.nanmean(refractory_phase_duration_ms),
                'median': np.nanmedian(refractory_phase_duration_ms),
                'min': np.nanmin(refractory_phase_duration_ms),
                'max': np.nanmax(refractory_phase_duration_ms),
            },
            'refractory_phase_end_ms': {
                'data': refractory_phase_end_ms,
                'mean': np.nanmean(refractory_phase_end_ms),
                'std': np.nanstd(refractory_phase_end_ms),
                'cov': np.nanstd(refractory_phase_end_ms) / np.nanmean(refractory_phase_end_ms),
                'median': np.nanmedian(refractory_phase_end_ms),
                'min': np.nanmin(refractory_phase_end_ms),
                'max': np.nanmax(refractory_phase_end_ms),
            },
            'num_phases': {
                'data': num_phases,
                'mean': np.nanmean(num_phases),
                'std': np.nanstd(num_phases),
                'cov': np.nanstd(num_phases) / np.nanmean(num_phases),
                'median': np.nanmedian(num_phases),
                'min': np.nanmin(num_phases),
                'max': np.nanmax(num_phases),
            },
            'pre_ap_phase_power_uv2': {
                'data': pre_ap_phase_power_uv2,
                'mean': np.nanmean(pre_ap_phase_power_uv2),
                'std': np.nanstd(pre_ap_phase_power_uv2),
                'cov': np.nanstd(pre_ap_phase_power_uv2) / np.nanmean(pre_ap_phase_power_uv2),
                'median': np.nanmedian(pre_ap_phase_power_uv2),
                'min': np.nanmin(pre_ap_phase_power_uv2),
                'max': np.nanmax(pre_ap_phase_power_uv2),
            },
            'ap_phase_power_uv2': {
                'data': ap_phase_power_uv2,
                'mean': np.nanmean(ap_phase_power_uv2),
                'std': np.nanstd(ap_phase_power_uv2),
                'cov': np.nanstd(ap_phase_power_uv2) / np.nanmean(ap_phase_power_uv2),
                'median': np.nanmedian(ap_phase_power_uv2),
                'min': np.nanmin(ap_phase_power_uv2),
                'max': np.nanmax(ap_phase_power_uv2),
            },
            'refractory_phase_power_uv2': {
                'data': refractory_phase_power_uv2,
                'mean': np.nanmean(refractory_phase_power_uv2),
                'std': np.nanstd(refractory_phase_power_uv2),
                'cov': np.nanstd(refractory_phase_power_uv2) / np.nanmean(refractory_phase_power_uv2),
                'median': np.nanmedian(refractory_phase_power_uv2),
                'min': np.nanmin(refractory_phase_power_uv2),
                'max': np.nanmax(refractory_phase_power_uv2),
            },
            'total_spike_power_uv2': {
                'data': total_spike_power_uv2,
                'mean': np.nanmean(total_spike_power_uv2),
                'std': np.nanstd(total_spike_power_uv2),
                'cov': np.nanstd(total_spike_power_uv2) / np.nanmean(total_spike_power_uv2),
                'median': np.nanmedian(total_spike_power_uv2),
                'min': np.nanmin(total_spike_power_uv2),
                'max': np.nanmax(total_spike_power_uv2),
            },
        },
        'advanced_metrics': {
            'trough_to_peak_durations': {
                'data': trough_to_peak_durations,
                'mean': np.nanmean(trough_to_peak_durations),
                'std': np.nanstd(trough_to_peak_durations),
                'cov': np.nanstd(trough_to_peak_durations) / np.nanmean(trough_to_peak_durations),
                'median': np.nanmedian(trough_to_peak_durations),
                'min': np.nanmin(trough_to_peak_durations),
                'max': np.nanmax(trough_to_peak_durations),
            },
            'spike_width_at_half_maxs': {
                'data': spike_width_at_half_maxs,
                'mean': np.nanmean(spike_width_at_half_maxs),
                'std': np.nanstd(spike_width_at_half_maxs),
                'cov': np.nanstd(spike_width_at_half_maxs) / np.nanmean(spike_width_at_half_maxs),
                'median': np.nanmedian(spike_width_at_half_maxs),
                'min': np.nanmin(spike_width_at_half_maxs),
                'max': np.nanmax(spike_width_at_half_maxs),
            },
            'half_max_amps': {
                'data': half_max_amps,
                'mean': np.nanmean(half_max_amps),
                'std': np.nanstd(half_max_amps),
                'cov': np.nanstd(half_max_amps) / np.nanmean(half_max_amps),
                'median': np.nanmedian(half_max_amps),
                'min': np.nanmin(half_max_amps),
                'max': np.nanmax(half_max_amps),
            },
            'trough_to_peak_ratios': {
                'data': trough_to_peak_ratios,
                'mean': np.nanmean(trough_to_peak_ratios),
                'std': np.nanstd(trough_to_peak_ratios),
                'cov': np.nanstd(trough_to_peak_ratios) / np.nanmean(trough_to_peak_ratios),
                'median': np.nanmedian(trough_to_peak_ratios),
                'min': np.nanmin(trough_to_peak_ratios),
                'max': np.nanmax(trough_to_peak_ratios),
            },
            'waveform_asymmetry_indices': {
                'data': waveform_asymmetry_indices,
                'mean': np.nanmean(waveform_asymmetry_indices),
                'std': np.nanstd(waveform_asymmetry_indices),
                'cov': np.nanstd(waveform_asymmetry_indices) / np.nanmean(waveform_asymmetry_indices),
                'median': np.nanmedian(waveform_asymmetry_indices),
            },
            'max_depolarization_slopes': {
                'data': max_depolarization_slopes,
                'mean': np.nanmean(max_depolarization_slopes),
                'std': np.nanstd(max_depolarization_slopes),
                'cov': np.nanstd(max_depolarization_slopes) / np.nanmean(max_depolarization_slopes),
                'median': np.nanmedian(max_depolarization_slopes),
            },
            'max_repolarization_slopes': {
                'data': max_repolarization_slopes,
                'mean': np.nanmean(max_repolarization_slopes),
                'std': np.nanstd(max_repolarization_slopes),
                'cov': np.nanstd(max_repolarization_slopes) / np.nanmean(max_repolarization_slopes),
                'median': np.nanmedian(max_repolarization_slopes),
            },
            'slope_ratios': {
                'data': slope_ratios,
                'mean': np.nanmean(slope_ratios),
                'std': np.nanstd(slope_ratios),
                'cov': np.nanstd(slope_ratios) / np.nanmean(slope_ratios),
                'median': np.nanmedian(slope_ratios),
            },
            'max_voltage_accelerations': {
                'data': max_voltage_accelerations,
                'mean': np.nanmean(max_voltage_accelerations),
                'std': np.nanstd(max_voltage_accelerations),
                'cov': np.nanstd(max_voltage_accelerations) / np.nanmean(max_voltage_accelerations),
                'median': np.nanmedian(max_voltage_accelerations),
            },
            'max_voltage_decelerations': {
                'data': max_voltage_decelerations,
                'mean': np.nanmean(max_voltage_decelerations),
                'std': np.nanstd(max_voltage_decelerations),
                'cov': np.nanstd(max_voltage_decelerations) / np.nanmean(max_voltage_decelerations),
                'median': np.nanmedian(max_voltage_decelerations),
            },
            'waveform_skewnesses': {
                'data': waveform_skewnesses,
                'mean': np.nanmean(waveform_skewnesses),
                'std': np.nanstd(waveform_skewnesses),
                'cov': np.nanstd(waveform_skewnesses) / np.nanmean(waveform_skewnesses),
                'median': np.nanmedian(waveform_skewnesses),
            },
            'waveform_kurtoses': {
                'data': waveform_kurtoses,
                'mean': np.nanmean(waveform_kurtoses),
                'std': np.nanstd(waveform_kurtoses),
                'cov': np.nanstd(waveform_kurtoses) / np.nanmean(waveform_kurtoses),
                'median': np.nanmedian(waveform_kurtoses),
            },
        },
        'individual_wf_metrics': metrics,
    }
    
    plot_avg_wf = True
    if plot_avg_wf:
        for i, wf in enumerate(interpolated_wfs):
            # plot individuals in grey and average in red
            plt.plot(wf, color='grey', alpha=0.5)
        plt.plot(np.nanmean(interpolated_wfs, axis=0), color='red')
        plt.tight_layout()
        plt.savefig('average_waveform.png')
        plt.close()  
    
    return wf_metrics

def get_experimental_network_activity_metrics_v2(sorting_object, recording_segment_object, wf_extractor, conv_params, mega_params, debug_mode = False, **kwargs):
    
    assert sorting_object is not None, 'No sorting object provided'
    
    #initialize network_data
    network_data = init_network_data_dict()
    
    #network data pathing info
    network_data['recording_path'] = recording_segment_object._kwargs.get('file_path', None)
    network_data['sorting_output'] = sorting_object._kwargs.get('folder_path', None)
    assert network_data['sorting_output'] is not None, 'No sorting object source found'
    
    # waveform info
    network_data['waveform_output'] = network_data['sorting_output'].replace('sorted', 'waveforms').replace('sorter_output', '')
    
    #get data
    sorting_object = sorting_object
    recording_object = recording_segment_object
    sampling_rate = recording_object.get_sampling_frequency()
    
    #convert time to seconds - get initially available data
    spike_times = get_spike_times(recording_object, sorting_object, sampling_rate=sampling_rate) #seconds
    timeVector = get_time_vector(recording_object, sampling_rate=sampling_rate) #seconds
    spike_times_by_unit = get_spike_times_by_unit(sorting_object, sampling_rate=sampling_rate) 
    
    #add conv_params to kwargs
    kwargs['conv_params'] = conv_params
    kwargs['mega_params'] = mega_params
    
    #add sorting_object and recording_object to kwargs if not already present
    kwargs['sorting_object'] = sorting_object
    kwargs['recording_object'] = recording_object
    kwargs['wf_extractor'] = wf_extractor
    
    #extract spiking metrics from simulated data
    # aw 2025-02-25 10:58:09 - plotting operations to be mindful of here:
    # 1. plot_wfs: whether to plot waveforms for individual units. Actual plotting occurrs in compute_wf_metrics
    try: 
        #extract_metrics_from_experimental_data_v2(spike_times, timeVector, spike_times_by_unit, **kwargs)
        #extract_metrics_from_experimental_data_v3(spike_times, timeVector, spike_times_by_unit, **kwargs) # aw 2025-02-21 01:12:14 - parallelized version
        extract_metrics_from_experimental_data_v3_threads(spike_times, timeVector, spike_times_by_unit, debug_mode = debug_mode, **kwargs) # aw 2025-02-21 01:12:14 - parallelized version w/ threads might be more efficient?
    except Exception as e:
        print(f'Error extracting metrics from simulated data: {e}')
        traceback.print_exc()
        pass
    
    #extract bursting metrics from simulated data (but this one works for both simulated and experimental data)
    try: 
        extract_bursting_activity_data(spike_times, spike_times_by_unit, debug_mode = debug_mode, **kwargs)
    except Exception as e:
        print(f'Error calculating bursting activity: {e}')
        traceback.print_exc()
        pass
    
    # classify neurons
    # from RBS_network_models.experimental import classify_neurons
    from MEA_Analysis.NeuronClassication.classify_neurons import classify_neurons_v2
    try:
        network_data = classify_neurons_v2(network_data, **kwargs)
    except Exception as e:
        print(f'Error classifying neurons: {e}')
        traceback.print_exc()
        pass
    
    # dynamic time warping (dtw) - burst analysis
    from MEA_Analysis.NetworkAnalysis.awNetworkAnalysis.dtw import dtw_burst_analysis
    try:
        network_data = dtw_burst_analysis(network_data, **kwargs)
    except Exception as e:
        print(f'Error in dtw burst analysis: {e}')
        traceback.print_exc()
        pass    
    
    return network_data

def extract_metrics_from_experimental_data_v2(spike_times, timeVector, spike_times_by_unit, **kwargs):
    
    #extract spiking metrics from experimental data
    recording_object = kwargs['recording_object']
    sorting_object = kwargs['sorting_object']
    wf_extractor = kwargs['wf_extractor']
    sampling_rate = recording_object.get_sampling_frequency()    
    
    #add immediately available data to network_data
    network_data['source'] = 'experimental'
    network_data['timeVector'] = timeVector
    network_data['spiking_data']['spike_times'] = spike_times
    network_data['spiking_data']['spiking_times_by_unit'] = spike_times_by_unit

    #spiking data by unit - just go through the simulated version of the data but de-identify pop
    #simulated_spiking_data_by_unit = network_data['simulated_data']['spiking_data_by_unit']
    #spiking_data_by_unit = {}
    spiking_metrics_by_unit = {}
    units = sorting_object.get_unit_ids()
    
    # set dir for plot output if needed
    #plot_wfs = kwargs.get('plot_wfs', False), # choose whether to plot waveforms
    plot_wfs = kwargs['plot_wfs']
    wf_extractor = kwargs['wf_extractor']
    wf_posixpath = wf_extractor.folder # this is posix path obj
    wf_well_folder = str(wf_posixpath)+'_plots/unit_wf_plots'  
    
    for unit in units:
        try:
            print(f'Processing unit {unit}...')
            
            # set path for unit plot
            unit_wf_path = wf_well_folder + f'/unit_{unit}_waveforms.png'
            
            #extract waveforms for each unit
            unit_wfs = wf_extractor.get_waveforms(unit)        
            
            # Compute mean waveform across all spikes
            avg_waveform = np.nanmean(unit_wfs, axis=0)  # Shape: (n_samples, n_channels)
            
            # Find the best channel (largest absolute amplitude)
            best_channel_idx = np.argmax(np.max(np.abs(avg_waveform), axis=0))  # Index of best channel
            
            # collect waveforms from best channel
            best_channel_waveforms = unit_wfs[:, :, best_channel_idx]
            
            # Compute amplitude of the average waveform
            
            #spiking_data_by_unit[unit] = {
            spiking_metrics_by_unit[unit] = {    
                # HACK: Hardcoded 6 second bin size
                # consistent with the 6 second bin size
                # user here: https://pmc.ncbi.nlm.nih.gov/articles/PMC3434456/
                
                # aw 2025-02-18 11:26:13 - as of now, I'm getting fano factor by getting spike count in distinct bursts - which I think is more flexible,
                # robust, and biologically informative than just binning spikes in 6 second bins.
                
                'num_spikes': len(spike_times_by_unit[unit]),
                'wf_metrics': compute_wf_metrics(best_channel_waveforms, sampling_rate, plot_wf = plot_wfs, save_fig=True, unit=unit, 
                                                 fig_name = unit_wf_path),
                'fr': get_unit_fr(recording_object, sorting_object, unit, sampling_rate),
                'isi': {
                    # 'data': spike_times_by_unit[unit],
                    # 'mean': np.nanmean(spike_times_by_unit[unit]),
                    # 'std': np.nanstd(spike_times_by_unit[unit]),
                    # 'cov': np.nanstd(spike_times_by_unit[unit]) / np.nanmean(spike_times_by_unit[unit]) if np.nanmean(spike_times_by_unit[unit]) > 0 else np.nan
                    'data': np.diff(spike_times_by_unit[unit]),
                    'mean': np.nanmean(np.diff(spike_times_by_unit[unit])),
                    'std': np.nanstd(np.diff(spike_times_by_unit[unit])),
                    'median': np.nanmedian(np.diff(spike_times_by_unit[unit])),
                    'cov': np.nanstd(np.diff(spike_times_by_unit[unit])) / np.nanmean(np.diff(spike_times_by_unit[unit])) if np.nanmean(np.diff(spike_times_by_unit[unit])) > 0 else np.nan,
                    'min': np.nanmin(np.diff(spike_times_by_unit[unit])),
                    'max': np.nanmax(np.diff(spike_times_by_unit[unit])),
                },
                'spike_times': spike_times_by_unit[unit],
            }
            
            # debug = True
            # if debug:
            #     print(f'Processed unit {unit}')
            #     if len(spiking_metrics_by_unit) > 2:
            #         break
            
        except Exception as e:            
            traceback.print_exc()
            print(f'Error processing unit {unit}: {e}')
            continue
    #network_data['spiking_data']['spiking_data_by_unit'] = spiking_data_by_unit
    network_data['spiking_data']['spiking_metrics_by_unit'] = spiking_metrics_by_unit

def process_burst(burst_id, burst_part):
    """Compute metrics for a given burst."""
    
    try:
        result = {}

        # Number of participating units
        result['num_units_participating'] = len(burst_part['participating_units'])
        
        # Spike counts by unit
        spike_counts = burst_part['spike_counts_by_unit']
        result['spike_counts_by_unit'] = {
            'data': spike_counts,
            'mean': np.nanmean(spike_counts),
            'std': np.nanstd(spike_counts),
            'cov': np.nanstd(spike_counts) / np.nanmean(spike_counts) if np.nanmean(spike_counts) > 0 else np.nan,
            'median': np.nanmedian(spike_counts),
            'min': np.nanmin(spike_counts),
            'max': np.nanmax(spike_counts),
        }
        
        # Burst duration
        sorted_spike_times = sorted(list(burst_part['spike_times']))
        duration = sorted_spike_times[-1] - sorted_spike_times[0]
        result['duration'] = duration
        
        # Spiking rate
        result['spike_rate'] = np.sum(spike_counts) / duration if duration > 0 else np.nan
        
        # ISI calculations
        isi_values = np.diff(sorted_spike_times)
        result['isi'] = {
            'data': isi_values,
            'mean': np.nanmean(isi_values),
            'std': np.nanstd(isi_values),
            'cov': np.nanstd(isi_values) / np.nanmean(isi_values) if np.nanmean(isi_values) > 0 else np.nan,
            'median': np.nanmedian(isi_values),
            'min': np.nanmin(isi_values),
            'max': np.nanmax(isi_values),
        }
        
        # Firing sequence
        sequence = []
        sequence_times = []
        for spike_time in sorted_spike_times:
            for unit, spike_times in zip(burst_part['participating_units'], burst_part['spike_times_by_unit']):
                if spike_time in spike_times:
                    sequence.append(unit)
                    sequence_times.append(spike_time)
                    #break  # Ensure each spike is only counted once
        #result['firing_sequence'] = sequence
        result['unit_seqeunce'] = sequence
        result['time_sequence'] = sequence_times
        result['relative_time_sequence'] = [time - sequence_times[0] for time in sequence_times]
        
        # #confirm sequence is at long as spike times
        # if len(sequence) != len(sorted_spike_times):
        #     print(f'Error: Sequence length does not match spike times')

        print(f'Burst {burst_id} sequenced')
        
        return burst_id, result
    except Exception as e:
        print(f'Error in processing burst {burst_id}: {e}')
        return burst_id, None

def compute_unit_burst_participation(unit_metrics, convolved_data, max_workers=4, debug_mode=False):
    '''For each burst, compute which units do and don't participate.'''
    indent_increase()
    left_base_times, right_base_times = convolved_data['left_base_times'], convolved_data['right_base_times']
    time_range = convolved_data['time_vector'][-1] - convolved_data['time_vector'][0]
    burst_starts = convolved_data['left_base_times']
    burst_ends = convolved_data['right_base_times']
    
    
    # Metrics ===============================================================
    participation_by_burst = {}
    for unit, metrics in unit_metrics.items():
        for burst_id, burst in metrics['bursts'].items():
            if burst_id not in participation_by_burst:
                participation_by_burst[burst_id] = {
                    #'duration': (burst[-1] - burst[0]),
                    'burst_start': burst_starts[burst_id],
                    'burst_end': burst_ends[burst_id],
                    'participating_units': [unit,],
                    'spike_count': len(burst),
                    'spike_times': set(burst),
                    'spike_counts_by_unit': [len(burst)],
                    'spike_times_by_unit': [burst,],
                    #'spike_times': {unit: burst,}                
                }
            else:
                participation_by_burst[burst_id]['participating_units'].append(unit)
                participation_by_burst[burst_id]['spike_count'] += len(burst)
                participation_by_burst[burst_id]['spike_times'].update(burst)
                participation_by_burst[burst_id]['spike_counts_by_unit'].append(len(burst))
                #participation_by_burst[burst_id]['spike_times'].append(set(burst))
                participation_by_burst[burst_id]['spike_times_by_unit'].append(burst)
    
    # sort participation by burst by key (burst_id)
    participation_by_burst = dict(sorted(participation_by_burst.items()))    
    
    # set max_workers
    if max_workers is None:
        max_workers = 1
    else:
        max_workers = max_workers
    
    # debug_mode
    if debug_mode:
        print(f'Debug mode: only processing first 5 burst parts')
        participation_by_burst = dict(list(participation_by_burst.items())[:5])
        max_workers = 5
    
    if max_workers is None:
        max_workers = 1
        
    if len(participation_by_burst) < max_workers:
        max_workers = len(participation_by_burst)
    
    if len(participation_by_burst) >= 1:    
        print(f'using {max_workers} cpus to process bursts')
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_burst_id = {executor.submit(process_burst, burst_id, burst_part): burst_id for burst_id, burst_part in participation_by_burst.items()}
            
            for future in as_completed(future_to_burst_id):
                burst_id, burst_result = future.result()
                participation_by_burst[burst_id].update(burst_result)
        
        # confirm that all bursts are accounted for and in order
        last_id = None
        for burst_id, burst_part in participation_by_burst.items():
            if last_id is None: 
                last_id = burst_id
                continue
            else:       
                if last_id > burst_id:
                    print(f'Error: Bursts are not in order')
                expected_val = last_id + 1
                if burst_id != expected_val:
                    print(f'Error: Missing burst {expected_val}')            
                last_id = burst_id                
        
        # return
        indent_decrease()
        return participation_by_burst
    else:
        print(f'No bursts to process')
        indent_decrease()
        return {}        

def compute_burst_metrics(unit_metrics, convolved_data, max_workers=4, debug_mode = False, **kwargs):
    ## unpack kwargs
    burst_sequencing = kwargs.get('burst_sequencing', False)
    
    # init
    warnings = None
    
    # get data
    peak_times = convolved_data['peak_times']
    
    # get unit participation metrics for each burst
    if burst_sequencing:
        burst_parts = compute_unit_burst_participation(unit_metrics, convolved_data, max_workers=max_workers, debug_mode = debug_mode)
        burst_parts = {**burst_parts} # HACK: Convert to dict replacing how its send to dict below...probably can do this better
    else:
        burst_parts = "Burst sequencing not enabled"
    
    #assemble dict
    burst_metrics = {
        'num_bursts': len(peak_times),
        'burst_rate': len(peak_times) / (convolved_data['time_vector'][-1] - convolved_data['time_vector'][0]) if len(peak_times) > 1 else np.nan,
        'burst_ids': list(range(len(peak_times))),
        'ibi': {
            'data': np.diff(peak_times),
            'mean': np.nanmean(np.diff(peak_times)),
            'std': np.nanstd(np.diff(peak_times)),
            'cov': np.nanstd(np.diff(peak_times)) / np.nanmean(np.diff(peak_times)) if np.nanmean(np.diff(peak_times)) > 0 else np.nan,
            'median': np.nanmedian(np.diff(peak_times)),
            'min': np.nanmin(np.diff(peak_times)),
            'max': np.nanmax(np.diff(peak_times)),
        },
        'burst_amp': {
            'data': convolved_data['peak_values'],
            'mean': np.nanmean(convolved_data['peak_values']),
            'std': np.nanstd(convolved_data['peak_values']),
            'cov': np.nanstd(convolved_data['peak_values']) / np.nanmean(convolved_data['peak_values']) if np.nanmean(convolved_data['peak_values']) > 0 else np.nan,
            'median': np.nanmedian(convolved_data['peak_values']),
            'min': np.nanmin(convolved_data['peak_values']),
            'max': np.nanmax(convolved_data['peak_values']),
        },
        'burst_duration': {
            'data': convolved_data['right_base_times'] - convolved_data['left_base_times'],
            'mean': np.nanmean(convolved_data['right_base_times'] - convolved_data['left_base_times']),
            'std': np.nanstd(convolved_data['right_base_times'] - convolved_data['left_base_times']),
            'cov': np.nanstd(convolved_data['right_base_times'] - convolved_data['left_base_times']) / np.nanmean(convolved_data['right_base_times'] - convolved_data['left_base_times']) if np.nanmean(convolved_data['right_base_times'] - convolved_data['left_base_times']) > 0 else np.nan,
            'median': np.nanmedian(convolved_data['right_base_times'] - convolved_data['left_base_times']),
            'min': np.nanmin(convolved_data['right_base_times'] - convolved_data['left_base_times']),
            'max': np.nanmax(convolved_data['right_base_times'] - convolved_data['left_base_times']),
        },
        #'burst_parts': {**burst_parts},
        'burst_parts': burst_parts,
        'num_units_per_burst': {
            'data': [burst_part['num_units_participating'] for burst_part in burst_parts.values()],
            'mean': np.nanmean([burst_part['num_units_participating'] for burst_part in burst_parts.values()]),
            'std': np.nanstd([burst_part['num_units_participating'] for burst_part in burst_parts.values()]),
            'cov': np.nanstd([burst_part['num_units_participating'] for burst_part in burst_parts.values()]) / np.nanmean([burst_part['num_units_participating'] for burst_part in burst_parts.values()]) if np.nanmean([burst_part['num_units_participating'] for burst_part in burst_parts.values()]) > 0 else np.nan,
            'median': np.nanmedian([burst_part['num_units_participating'] for burst_part in burst_parts.values()]),
            'min': np.nanmin([burst_part['num_units_participating'] for burst_part in burst_parts.values()]),
            'max': np.nanmax([burst_part['num_units_participating'] for burst_part in burst_parts.values()]),
        } if burst_parts != "Burst sequencing not enabled" else "Burst sequencing not enabled",
        'in_burst_fr': {
            'data': [burst_part['spike_rate'] for burst_part in burst_parts.values()],
            'mean': np.nanmean([burst_part['spike_rate'] for burst_part in burst_parts.values()]),
            'std': np.nanstd([burst_part['spike_rate'] for burst_part in burst_parts.values()]),
            'cov': np.nanstd([burst_part['spike_rate'] for burst_part in burst_parts.values()]) / np.nanmean([burst_part['spike_rate'] for burst_part in burst_parts.values()]) if np.nanmean([burst_part['spike_rate'] for burst_part in burst_parts.values()]) > 0 else np.nan,
            'median': np.nanmedian([burst_part['spike_rate'] for burst_part in burst_parts.values()]),
            'min': np.nanmin([burst_part['spike_rate'] for burst_part in burst_parts.values()]),
            'max': np.nanmax([burst_part['spike_rate'] for burst_part in burst_parts.values()]),
            } if burst_parts != "Burst sequencing not enabled" else "Burst sequencing not enabled",          
    }
    
    print(f'Burst Metrics Computed')
    return burst_metrics, warnings

def analyze_unit_activity(spike_times_by_unit, convolved_data):
    """Analyze bursts, quiet periods, firing rates, and ISIs for each unit."""
    total_number_of_bursts_in_convolved_data = len(convolved_data['peak_times'])
    left_base_times, right_base_times = convolved_data['left_base_times'], convolved_data['right_base_times']
    time_range = convolved_data['time_vector'][-1] - convolved_data['time_vector'][0]
    
    
    bursts_by_unit, non_bursts_by_unit, burst_durations, quiet_durations, warnings = {}, {}, {}, {}, []
    
    unit_data = {}
    for unit, spike_times in spike_times_by_unit.items():
        
        bursts, non_bursts = {}, {}
        burst_id = 0
        quiet_id = 0

        for left, right in zip(left_base_times, right_base_times):
            burst_spikes = spike_times[(spike_times >= left) & (spike_times <= right)]
            burst_durations[burst_id] = right - left
            if burst_spikes.size > 0:
                bursts[burst_id] = burst_spikes
                if np.any(np.diff(burst_spikes) < 0):
                    warnings.append(f'Negative ISI in unit {unit} (burst {burst_id})')
            burst_id += 1

        quiet_left = np.concatenate([[0], right_base_times[:-1]])
        quiet_right = np.concatenate([left_base_times[1:], [spike_times[-1]]])

        for left, right in zip(quiet_left, quiet_right):
            quiet_spikes = spike_times[(spike_times >= left) & (spike_times <= right)]
            quiet_durations[quiet_id] = right - left
            if quiet_spikes.size > 0:
                non_bursts[quiet_id] = quiet_spikes
                if np.any(np.diff(quiet_spikes) < 0):
                    warnings.append(f'Negative ISI in unit {unit} (non-burst {burst_id})')
            quiet_id += 1

        bursts_by_unit[unit] = bursts
        non_bursts_by_unit[unit] = non_bursts
        
        # 
        def isi_analytics():
            isi_data_in = {i: np.diff(burst) for i, burst in bursts.items()}
            isi_data_out = {i: np.diff(burst) for i, burst in non_bursts.items()}
            
            def compute_isi_metrics(isis):
                try:
                    all_isis = np.concatenate([list(isis) for isis in isis.values() if len(isis) > 0])
                    return {
                        'data': all_isis,
                        'mean': np.nanmean(all_isis),
                        'std': np.nanstd(all_isis),
                        'cov': np.nanstd(all_isis) / np.nanmean(all_isis) if np.nanmean(all_isis) > 0 else np.nan,
                        'median': np.nanmedian(all_isis),
                        'min': np.nanmin(all_isis),
                        'max': np.nanmax(all_isis),
                    }
                except Exception as e:
                    #print(f'Error in computing ISI metrics: {e}')
                    return {
                        'data': np.nan,
                        'mean': np.nan,
                        'std': np.nan,
                        'cov': np.nan,
                        'median': np.nan,
                        'min': np.nan,
                        'max': np.nan,
                    }
            return {
                'in_burst': {**compute_isi_metrics(isi_data_in)},
                'out_burst': {**compute_isi_metrics(isi_data_out)}
            }
            
        def fr_analytics():
            #use burst_durations to calculate firing rate
            fr_data_in = {i: len(burst) / burst_durations[i] if len(burst) > 1 else np.nan for i, burst in bursts.items()}
            fr_data_out = {i: len(burst) / quiet_durations[i] if len(burst) > 1 else np.nan for i, burst in non_bursts.items()}
            
            # remove inf values where burst duration is 0, replace with nan
            fr_data_in = {k: np.nan if v == np.inf else v for k, v in fr_data_in.items()}
            fr_data_out = {k: np.nan if v == np.inf else v for k, v in fr_data_out.items()}
            
            # assert that the only values in the dict are floats and np.nan values
            assert all(isinstance(v, (float, int)) or np.isnan(v) for v in fr_data_in.values())
            assert all(isinstance(v, (float, int)) or np.isnan(v) for v in fr_data_out.values())
            
            # Handle cases where burst or quiet periods are empty
            fr_data_in = {i: fr for i, fr in fr_data_in.items() if not (fr is None or np.isnan(fr))}
            fr_data_out = {i: fr for i, fr in fr_data_out.items() if not (fr is None or np.isnan(fr))}
            
            return {
                'in_burst': {
                    'data': fr_data_in,
                    'mean': np.nanmean(list(fr_data_in.values())) if fr_data_in else np.nan,
                    'std': np.nanstd(list(fr_data_in.values())) if fr_data_in else np.nan,
                    'cov': np.nanstd(list(fr_data_in.values())) / np.nanmean(list(fr_data_in.values())) if fr_data_in and np.nanmean(list(fr_data_in.values())) > 0 else np.nan,
                    'median': np.nanmedian(list(fr_data_in.values())) if fr_data_in else np.nan,
                    'min': np.nanmin(list(fr_data_in.values())) if fr_data_in else np.nan,
                    'max': np.nanmax(list(fr_data_in.values())) if fr_data_in else np.nan,
                    },
                'out_burst': {
                    'data': fr_data_out,
                    'mean': np.nanmean(list(fr_data_out.values())) if fr_data_out else np.nan,
                    'std': np.nanstd(list(fr_data_out.values())) if fr_data_out else np.nan,
                    'cov': np.nanstd(list(fr_data_out.values())) / np.nanmean(list(fr_data_out.values())) if fr_data_out and np.nanmean(list(fr_data_out.values())) > 0 else np.nan,
                    'median': np.nanmedian(list(fr_data_out.values())) if fr_data_out else np.nan,
                    'min': np.nanmin(list(fr_data_out.values())) if fr_data_out else np.nan,
                    'max': np.nanmax(list(fr_data_out.values())) if fr_data_out else np.nan,
                    },
                }

        def spike_count_analytics():
            spike_count_in_burst = {i: len(burst) for i, burst in bursts.items()} if bursts else {}
            spike_count_out_burst = {i: len(burst) for i, burst in non_bursts.items()} if non_bursts else {}

            def compute_stats(spike_counts):
                if not spike_counts or all(np.isnan(list(spike_counts.values()))):  # Check for empty or all NaN values
                    return {
                        'data': spike_counts,
                        'mean': np.nan,
                        'std': np.nan,
                        'cov': np.nan,
                        'median': np.nan,
                        'min': np.nan,
                        'max': np.nan
                    }
                
                values = np.array(list(spike_counts.values()), dtype=float)
                values = values[~np.isnan(values)]  # Remove NaNs if present
                
                if values.size == 0:  # If all values were NaN and got removed
                    return {
                        'data': spike_counts,
                        'mean': np.nan,
                        'std': np.nan,
                        'cov': np.nan,
                        'median': np.nan,
                        'min': np.nan,
                        'max': np.nan
                    }

                mean_val = np.nanmean(values)
                std_val = np.nanstd(values)
                return {
                    'data': spike_counts,
                    'mean': mean_val,
                    'std': std_val,
                    'cov': std_val / mean_val if mean_val > 0 else np.nan,
                    'median': np.nanmedian(values),
                    'min': np.nanmin(values),
                    'max': np.nanmax(values),
                }

            return {
                'in_burst': compute_stats(spike_count_in_burst),
                'out_burst': compute_stats(spike_count_out_burst)
            }

            
        def fano_factor_analytics():
            spike_count_in_burst = spike_count_analytics()['in_burst']['data']
            spike_count_out_burst = spike_count_analytics()['out_burst']['data']
            
            return {
                'in_burst': np.var(list(spike_count_in_burst.values())) / np.mean(list(spike_count_in_burst.values())) if np.mean(list(spike_count_in_burst.values())) > 0 else np.nan,
                'out_burst': np.var(list(spike_count_out_burst.values())) / np.mean(list(spike_count_out_burst.values())) if np.mean(list(spike_count_out_burst.values())) > 0 else np.nan,
                }
            
        fr_analytics_dict = fr_analytics()
        isi_analytics_dict = isi_analytics()
        spike_count_analytics_dict = spike_count_analytics()
        fano_factor_dict = fano_factor_analytics()
        
        try:
            # build unit data structure
            unit_data[unit] = {          
                'burst_id': [key for key in bursts.keys()],
                'quiet_id': [key for key in non_bursts.keys()],
                'bursts': bursts,
                'quiets': non_bursts,
                'burst_durations': burst_durations,
                'quiet_durations': quiet_durations,
                'burst_part_rate': len(bursts) / time_range,
                'quiet_part_rate': len(non_bursts) / time_range,
                'burst_part_perc': len(bursts) / total_number_of_bursts_in_convolved_data,
                'fr': fr_analytics_dict,
                'isi': isi_analytics_dict,
                'spike_counts': spike_count_analytics_dict,            
                'fano_factor': fano_factor_dict,
                'note': "The same spike may be represented in multipe bursts if they compound in overlapping epochs",
                }
        except Exception as e:
            print(f'Error in building unit data structure: {e}')
            return None
        
        # print(f'Unit {unit}: {len(bursts)} bursts, {len(non_bursts)} non-bursts')
    if len(warnings) == 0:
        warnings = None
    return unit_data, warnings
    # fr_by_unit = {unit: len(spikes) / time_range for unit, spikes in spike_times_by_unit.items()}
    # isis_by_unit = {unit: np.diff(spikes) for unit, spikes in spike_times_by_unit.items()}
    # return bursts_by_unit, non_bursts_by_unit, fr_by_unit, isis_by_unit, warnings

def compute_summary_metrics(unit_data):
    """Compute mean, standard deviation, and coefficient of variation for firing rates and ISIs."""
    def summarize(values):
        all_values = np.concatenate([list(fr.values()) for fr in values.values() if len(fr) > 0])
        return {
            'mean': np.nanmean(all_values),
            'std': np.nanstd(all_values),
            'cov': np.nanstd(all_values) / np.nanmean(all_values) if np.nanmean(all_values) > 0 else np.nan
        }
    
    return summarize(fr_by_unit), summarize(isis_by_unit)

def analyze_bursting_activity_v4(spike_times, spike_times_by_unit, min_peak_distance=1.0, binSize=0.001,
                                 gaussianSigma=0.01, thresholdBurst=0.5, prominence=1, title='Network Activity (v4)',
                                 debug_mode=False,
                                 **kwargs
                                 ):
    ''' 
    
    '''
    # init warnings
    unit_warnings = None
         
    #
    try:
        # 
        start = time.time()
        indent_increase()
        
        # Step 1: Generate convolved data from network activity plot
        try:
            fig, ax = plt.subplots()
            ax, convolved_data = plot_network_activity_aw(
                ax, spike_times_by_unit,
                binSize=binSize, gaussianSigma=gaussianSigma,
                thresholdBurst=thresholdBurst, prominence=prominence, title=title
            )
        except Exception as e:
            traceback.print_exc()
            print(f'Error in plotting network activity: {e}')
            #line_debug(ax)
            plt.close(fig)            
            return None

        # Step 2: Analyze individual units
        try:
            #max_workers = kwargs.get('max_workers', 4)
            unit_metrics, unit_warnings = analyze_unit_activity(spike_times_by_unit, convolved_data)
            if unit_warnings:
                print("\n".join(unit_warnings))
        except Exception as e:
            traceback.print_exc()
            print(f'Error in analyzing unit activity: {e}')
            #line_debug(ax)
            plt.close(fig)
            return {
                'ax' : ax,
                'convolved_data': convolved_data,
                'unit_metrics': None,
                'burst_metrics': None,
                'warnings': {
                    'unit': unit_warnings,
                    'burst': None
                }
            }
            
        # Step 3: Burst Summary Metrics
        try:
            max_workers = kwargs.get('max_workers', 4)
            burst_metrics, burst_warnings = compute_burst_metrics(unit_metrics, convolved_data, debug_mode=debug_mode, **kwargs)
            if burst_warnings:
                print("\n".join(burst_warnings))
        except Exception as e:
            traceback.print_exc()
            print(f'Error in computing burst metrics: {e}')
            indent_decrease()
            #debug - check if zero lines in ax
            #line_debug(ax)
            plt.close(fig) # NOTE: for some reason this is required for ax to persist later in the code -- I dont fully understand why
            return {
                'ax' : ax,
                'convolved_data': convolved_data,
                'unit_metrics': unit_metrics,
                'burst_metrics': None,
                'warnings': {
                    'unit': unit_warnings,
                    'burst': burst_warnings
                }
            }
        
        # return
        print(f'Elapsed time: {time.time() - start}')
        indent_decrease()
        
        #debug - check if zero lines in ax
        #line_debug(ax)
        plt.close(fig)
        return {
            'ax': ax,
            'convolved_data': convolved_data,
            'unit_metrics': unit_metrics,
            'burst_metrics': burst_metrics,
            'warnings': {
                'unit': unit_warnings,
                'burst': burst_warnings,
            }
        }
    except Exception as e:
        traceback.print_exc()
        print(f'Error in bursting activity analysis: {e}')
        indent_decrease()
        return None

def analyze_bursting_activity_v3(
    spike_times, 
    spike_times_by_unit, 
    min_peak_distance=1.0,
    binSize=0.001,
    gaussianSigma=0.01,
    thresholdBurst=0.5,
    prominence=1,
    title='Network Activity (v2)',
    ):
    try:
        # Step 1: Generate convolved data from network activity plot
        fig, ax = plt.subplots()
        ax, convolved_data = plot_network_activity_aw(
            ax, spike_times_by_unit,
            binSize=binSize, gaussianSigma=gaussianSigma,
            thresholdBurst=thresholdBurst, prominence=prominence, title=title
        )

        left_base_times, right_base_times = convolved_data['left_base_times'], convolved_data['right_base_times']
        time_range = convolved_data['time_vector'][-1] - convolved_data['time_vector'][0]

        ## --- Step 2: Parse Burst and Quiet Epochs --- ##
        def parse_epochs(spike_times_by_unit, left_base_times, right_base_times):
            bursts_by_unit, non_bursts_by_unit = {}, {}
            warnings = []

            for unit, spike_times in spike_times_by_unit.items():
                bursts, non_bursts = {}, {}
                burst_id = 0

                for left, right in zip(left_base_times, right_base_times):
                    burst_spikes = spike_times[(spike_times >= left) & (spike_times <= right)]
                    if burst_spikes.size > 0:
                        bursts[burst_id] = burst_spikes
                        if np.any(np.diff(burst_spikes) < 0):
                            warnings.append(f'Negative ISI in unit {unit} (burst {burst_id})')

                    burst_id += 1

                quiet_left = np.concatenate([[0], right_base_times[:-1]])
                quiet_right = np.concatenate([left_base_times[1:], [spike_times[-1]]])

                for left, right in zip(quiet_left, quiet_right):
                    quiet_spikes = spike_times[(spike_times >= left) & (spike_times <= right)]
                    if quiet_spikes.size > 0:
                        non_bursts[burst_id] = quiet_spikes
                        if np.any(np.diff(quiet_spikes) < 0):
                            warnings.append(f'Negative ISI in unit {unit} (non-burst {burst_id})')

                    burst_id += 1

                bursts_by_unit[unit] = bursts
                non_bursts_by_unit[unit] = non_bursts

            return bursts_by_unit, non_bursts_by_unit, warnings

        bursts_by_unit, non_bursts_by_unit, warnings = parse_epochs(spike_times_by_unit, left_base_times, right_base_times)
        if warnings:
            print("\n".join(warnings))

        ## --- Step 3: Compute Firing Rates --- ##
        fr_by_unit = {unit: len(spikes) / time_range for unit, spikes in spike_times_by_unit.items()}

        def compute_firing_rates(bursts_by_unit):
            fr_by_unit = {}
            for unit, bursts in bursts_by_unit.items():
                fr_by_unit[unit] = {
                    i: len(burst) / (burst[-1] - burst[0]) if len(burst) > 1 else np.nan
                    for i, burst in bursts.items()
                }
            return fr_by_unit

        fr_by_unit_in_burst = compute_firing_rates(bursts_by_unit)
        fr_by_unit_out_burst = compute_firing_rates(non_bursts_by_unit)

        def summarize_firing_rates(fr_by_unit):
            all_frs = np.concatenate([list(fr.values()) for fr in fr_by_unit.values()])
            return {
                'mean': np.nanmean(all_frs),
                'std': np.nanstd(all_frs),
                'cov': np.nanstd(all_frs) / np.nanmean(all_frs) if np.nanmean(all_frs) > 0 else np.nan
            }

        fr_summary_in_bursts = summarize_firing_rates(fr_by_unit_in_burst)
        fr_summary_out_bursts = summarize_firing_rates(fr_by_unit_out_burst)

        ## --- Step 4: Compute ISIs --- ##
        isis_by_unit = {unit: np.diff(spikes) for unit, spikes in spike_times_by_unit.items()}

        def compute_isi_metrics(isis_by_unit):
            all_isis = [isis for isis in isis_by_unit.values() if len(isis) > 0]  # Ensure only non-empty arrays
            
            # Repeat until all_isis is a flat list of floats
            try:
                while (any(isinstance(i, list) for i in all_isis) 
                        or any(isinstance(i, np.ndarray)) for i in all_isis
                        or any(isinstance(i, dict) for i in all_isis)
                        ):
                    if all(isinstance(i, (int,float)) for i in all_isis):
                        break
                    all_isis = [item for sublist in all_isis for item in sublist]
            except Exception as e:
                print(f'Error in flattening all_isis: {e}')

            if not all_isis:  # Check if the list is empty
                return {'mean': np.nan, 'std': np.nan, 'cov': np.nan, 'fano_factor': np.nan}
            
            #all_isis = np.concatenate(all_isis)  # Now safe to concatenate

            return {
                'mean': np.nanmean(all_isis),
                'std': np.nanstd(all_isis),
                'cov': np.nanstd(all_isis) / np.nanmean(all_isis) if np.nanmean(all_isis) > 0 else np.nan,
                'fano_factor': np.nanvar(all_isis) / np.nanmean(all_isis) if np.nanmean(all_isis) > 0 else np.nan
                }   

        isis_summary = compute_isi_metrics(isis_by_unit)

        def compute_isi_by_burst(bursts_by_unit):
            isis_by_unit = {unit: {i: np.diff(burst) for i, burst in bursts.items()} for unit, bursts in bursts_by_unit.items()}
            return isis_by_unit, compute_isi_metrics(isis_by_unit)

        isis_by_unit_in_burst, isis_summary_in_bursts = compute_isi_by_burst(bursts_by_unit)
        isis_by_unit_out_burst, isis_summary_out_bursts = compute_isi_by_burst(non_bursts_by_unit)

        ## --- Step 5: Compute Burst Statistics --- ##
        peak_times = convolved_data['peak_times']
        ibi = np.diff(peak_times)
        burst_stats = {
            'Number_Bursts': len(peak_times),
            'mean_IBI': np.nanmean(ibi),
            'std_IBI': np.nanstd(ibi),
            'cov_IBI': np.nanstd(ibi) / np.nanmean(ibi) if np.nanmean(ibi) > 0 else np.nan
        }

        ## --- Final Summary --- ##
        burst_summary_data = {
            #
            'MeanWithinBurstISI': isis_summary_in_bursts['mean'],
            'CoVWithinBurstISI': isis_summary_in_bursts['cov'],
            #
            'MeanOutsideBurstISI': isis_summary_out_bursts['mean'],
            'CoVOutsideBurstISI': isis_summary_out_bursts['cov'],
            #
            'MeanNetworkISI': isis_summary['mean'],
            'CoVNetworkISI': isis_summary['cov'],
            #
            'NumUnits': len(spike_times_by_unit),
            #
            'isi_fano_factor': isis_summary['fano_factor'],
            'isi_fano_factor_in_bursts': isis_summary_in_bursts['fano_factor'],
            'isi_fano_factor_out_bursts': isis_summary_out_bursts['fano_factor'],
        }

        burst_data_by_unit = {
            unit: {
                #
                'bursts': bursts_by_unit[unit],
                #
                'mean_isi_within': isis_summary_in_bursts['mean'],
                'cov_isi_within': isis_summary_in_bursts['cov'],
                #
                'mean_isi_outside': isis_summary_out_bursts['mean'],
                'cov_isi_outside': isis_summary_out_bursts['cov'],
                #
                'isis_within_bursts': isis_by_unit_in_burst[unit],
                'isis_outside_bursts': isis_by_unit_out_burst[unit],
                'isis_all': isis_by_unit[unit],
            } for unit in spike_times_by_unit.keys()
        }

        return {'burst_summary_data': burst_summary_data, 'burst_data_by_unit': burst_data_by_unit}

    except Exception as e:
        print(f'Error in bursting activity analysis: {e}')
        return None

def plot_network_activity_aw(ax,SpikeTimes, min_peak_distance=1.0, 
                             binSize=0.1, 
                             gaussianSigma=0.16, 
                             thresholdBurst=1.2, 
                             prominence=1, 
                             figSize=(10, 6),
                             title='Network Activity'
                             ):
    # TODO: stole this code from MEA_pipline ipn analysis. I'm going to make some changes for myself. 
            # I need to put them back into the MEA pipeline at somepoint
            
    # # aw 2025-01-16 13:26:25 - its been a while since I made teh note above. I'm pretty much ready to put changes
    # back into the MEA pipeline.
    
    # # aw 2025-02-10 11:23:30 recently realized that MP burst statistics function is
    # still based on threshold. I'm consolidating all analysis to use the prominence
    # idea.
    
    relativeSpikeTimes = []
    units = 0
    for unit_id, spike_times in SpikeTimes.items():
        temp_spike_times = spike_times
        if isinstance(temp_spike_times, np.ndarray): relativeSpikeTimes.extend(temp_spike_times) 
        elif isinstance(temp_spike_times, float): relativeSpikeTimes.append(temp_spike_times)
        elif isinstance(temp_spike_times, int): relativeSpikeTimes.append(temp_spike_times)
        else:
            Warning(f'unit {unit_id} has spike times that are not float or array')
            continue
        units += 1 # Set the first spike time to 0
    assert all(isinstance(x, (int, float)) for x in relativeSpikeTimes), 'All elements in relativeSpikeTimes must be ints or floats'
    
    relativeSpikeTimes = np.array(relativeSpikeTimes)
    relativeSpikeTimes.sort() # Convert to NumPy array

    # Step 1: Bin all spike times into small time windows
    timeVector = np.arange(min(relativeSpikeTimes), max(relativeSpikeTimes), binSize)  # Time vector for binning
    binnedTimes, _ = np.histogram(relativeSpikeTimes, bins=timeVector)  # Bin spike times
    binnedTimes = np.append(binnedTimes, 0)  # Append 0 to match MATLAB's binnedTimes length

    # Step 2: Smooth the binned spike times with a Gaussian kernel
    kernelRange = np.arange(-3*gaussianSigma, 3*gaussianSigma + binSize, binSize)  # Range for Gaussian kernel
    kernel = norm.pdf(kernelRange, 0, gaussianSigma)  # Gaussian kernel
    kernel *= binSize  # Normalize kernel by bin size
    firingRate = convolve(binnedTimes, kernel, mode='same') / binSize  # Convolve and normalize by bin size
    firingRate = firingRate / units  # Convert to Hz

    # Plot the smoothed network activity
    ax.plot(timeVector, firingRate, color='royalblue')
    ax.set_xlim([min(relativeSpikeTimes), max(relativeSpikeTimes)])
    ax.set_ylim([min(firingRate)*0.85, max(firingRate)*1.15])  # Set y-axis limits to min and max of firingRate
    ax.set_ylabel('Firing Rate [Hz]')
    ax.set_xlabel('Time [ms]')
    ax.set_title(title, fontsize=11)

    # Step 3: Peak detection on the smoothed firing rate curve
    #rmsFiringRate = np.sqrt(np.mean(firingRate**2))  # Calculate RMS of the firing rate

    peaks, properties = find_peaks(firingRate, prominence=prominence, distance=min_peak_distance)  # Find peaks above the threshold
    #print(properties.keys())
    burstPeakTimes = timeVector[peaks]  # Convert peak indices to times
    burstPeakValues = firingRate[peaks]  # Get the peak values

    # #Calculate the ISIs between spiketimes
    # # Calculate the intervals between consecutive peaks
    # intervals = np.diff(burstPeakTimes)

    # # Calculate the mean interval between peaks
    # mean_interburstinterval = np.mean(intervals)
    
    # covariance_interburstinterval = np.cov(intervals)
    # # Count the number of peaks
    # num_peaks = len(peaks)

    # # Calculate the mean peak height
    # mean_peak_height = np.mean(burstPeakValues)
    # cov_peak_height = np.cov(burstPeakValues)

    # #spikecounts= [len(x) for x in SpikeTimes.values()]
    # spikecounts = []
    # for x in SpikeTimes.values():
    #     try: count = len(x)
    #     except:
    #         if isinstance(x, float): count = 1
    #         else: count = 0
    #     spikecounts.append(count)
    
    # var_spikecounts = np.var(spikecounts)
    # mean_spikecounts = np.mean(spikecounts)
    # fanofact = var_spikecounts/mean_spikecounts
    
    # #mean burst rate
    # time_in_seconds = timeVector[-1] - timeVector[0]
    # mean_burst_rate = num_peaks/time_in_seconds

    # network_data = {
    #     "Number_Bursts": num_peaks,
    #     'mean_Burst_Rate': mean_burst_rate,
    #     "mean_IBI": mean_interburstinterval,
    #     "cov_IBI": covariance_interburstinterval,
    #     "mean_Burst_Peak": mean_peak_height, # TODO: add fitting func for this later
    #     "cov_Burst_Peak": cov_peak_height,
    #     "fano_factor": fanofact
    # }   


    # Plot the threshold line and burst peaks
    #ax.plot(np.arange(timeVector[-1]), thresholdBurst * rmsFiringRate * np.ones(np.ceil(timeVector[-1]).astype(int)), color='gray')
    ax.plot(burstPeakTimes, burstPeakValues, 'or')  # Plot burst peaks as red circles    
    #ax.plot(timeVector[peaks], firingRate[peaks], 'or')  # Plot burst peaks as red circles
    
    # save plot for debug
    #plt.savefig('network_activity_plot.png')
    #print(ax.getlines())
    #print(ax.get_lines()[0])
    
    convolved_data = {
        'convolved_FR': firingRate,
        'peak_idxs': peaks,
        'peak_times': burstPeakTimes,
        'peak_values': burstPeakValues,
        'prominences': properties['prominences'],
        'left_base_idxs': properties['left_bases'],
        'right_base_idxs': properties['right_bases'],
        'left_base_times': timeVector[properties['left_bases']],
        'right_base_times': timeVector[properties['right_bases']],
        'time_vector': timeVector,
        #'rms_FR': rmsFiringRate,
    }
    
    #return ax, network_data
    #return ax, signal, peaks
    return ax, convolved_data

def analyze_bursting_activity_v2_unfinished(
    spike_times, 
    spike_times_by_unit, 
    #isi_threshold,
    binSize=0.001,
    gaussianSigma=0.01,
    thresholdBurst=0.5,
    min_peak_distance=None,
    prominence=1,
    title = 'Network Activity (v2)',
    ):
    
    '''Slightly modified version of the code from mea_analysis_pipeline.py'''
    #spike_times = spike_times_by_unit
    
    # Retired part of Mandar Code
    # #burst_statistics = helper.detect_bursts_statistics(spike_times, isi_threshold=isi_threshold)
    # burst_statistics = detect_bursts_statistics(spike_times, isi_threshold=isi_threshold)
    # bursts_by_unit = [unit_stats['bursts'] for unit_stats in burst_statistics.values()]
    
    try: 
        # Step 1: Convolve the data? So that all downstream analsys is based
        # on the convolved data?
        fig, ax = plt.subplots()
        ax, convolved_data = plot_network_activity_aw(
            ax, 
            spike_times_by_unit,
            binSize=binSize, 
            gaussianSigma=gaussianSigma,
            thresholdBurst=thresholdBurst,
            prominence=prominence,
            title=title
            )   
        
        # spiketimes  ========================================
        # get spikes in and out of bursts parsed by bursting or non-bursting epochs
        left_base_times = convolved_data['left_base_times']
        right_base_times = convolved_data['right_base_times']
        def parse_burst_or_quiet_epochs(spike_times_by_unit, left_base_times, right_base_times):
            in_burst_mask_by_unit = {}
            out_burst_mask_by_unit = {}
            bursts_by_unit = {}
            non_bursts_by_unit = {}
            for unit, spike_times in spike_times_by_unit.items():
                in_burst_mask = np.zeros(len(spike_times), dtype=bool)
                bursts = {}
                non_bursts = {}
                burst_id = 0
                #current_burst = []
                #current_non_burst = []
                for left, right in zip(left_base_times, right_base_times):
                    in_burst_mask[(spike_times >= left) & (spike_times <= right)] = True
                    burst_spikes = spike_times[(spike_times >= left) & (spike_times <= right)]
                    #non_burst_spikes = spike_times[(spike_times < left) | (spike_times > right)]
                    
                    # if len(burst_spikes) > 0:
                    #     current_burst.extend(burst_spikes)
                    # else:
                    #     if len(current_burst) > 0:
                    #         bursts.append(np.array(current_burst))
                    #         current_burst = []
                    
                    current_burst = np.array(burst_spikes)
                    if len(current_burst) > 0:
                        #bursts.append(np.array(burst_spikes))
                        bursts[burst_id] = np.array(burst_spikes)
                            
                    # QUALITY
                    isis = np.diff(burst_spikes)
                    if any(isis < 0):
                        print(f'Warning: Negative ISI in unit {unit}')
                        break_switch = True
                        
                    isis = np.diff(current_burst)
                    if any(isis < 0):
                        print(f'Warning: Negative ISI in unit {unit}')
                        break_switch = True
                        
                    #update burst_id
                    burst_id += 1
                # invert left and right to get non-burst epochs
                quiet_left_base_times = np.concatenate([[0], right_base_times[:-1]])
                quiet_right_base_times = np.concatenate([left_base_times[1:], [spike_times[-1]]])
                quiet_zip = zip(quiet_left_base_times, quiet_right_base_times)       
                for left, right in quiet_zip:
                    #out_burst_mask[(spike_times >= left) & (spike_times <= right)] = True
                    #burst_spikes = spike_times[(spike_times >= left) & (spike_times <= right)]
                    quiet_spikes = spike_times[(spike_times >= left) & (spike_times <= right)]
                    #non_burst_spikes = spike_times[(spike_times < left) | (spike_times > right)]
                    # if len(quiet_spikes) > 0:
                    #     current_non_burst.extend(quiet_spikes)
                    # else:
                    #     if len(current_non_burst) > 0:
                    #         non_bursts.append(np.array(current_non_burst))
                    #         current_non_burst = []
                    
                    current_non_burst = np.array(quiet_spikes)
                    if len(current_non_burst) > 0:
                        #non_bursts.append(np.array(quiet_spikes))
                        non_bursts[burst_id] = np.array(quiet_spikes)
                    
                    # QUALITY
                    isis = np.diff(quiet_spikes)
                    if any(isis < 0):
                        print(f'Warning: Negative ISI in unit {unit}')
                        break_switch = True
                        
                    isis = np.diff(current_non_burst)
                    if any(isis < 0):
                        print(f'Warning: Negative ISI in unit {unit}')
                        break_switch = True
                    
                    #update burst_id
                    burst_id += 1 
                # if len(current_burst) > 0:
                #     bursts.append(np.array(current_burst))
                # if len(current_non_burst) > 0:
                #     non_bursts.append(np.array(current_non_burst))
                out_burst_mask = ~in_burst_mask
                in_burst_mask_by_unit[unit] = in_burst_mask
                out_burst_mask_by_unit[unit] = out_burst_mask
                bursts_by_unit[unit] = bursts
                non_bursts_by_unit[unit] = non_bursts
                
            burst_mask = in_burst_mask_by_unit
            return bursts_by_unit, non_bursts_by_unit, burst_mask
        bursts_by_unit, non_bursts_by_unit, burst_mask = parse_burst_or_quiet_epochs(spike_times_by_unit, left_base_times, right_base_times)
        spike_times_by_unit_in_burst = bursts_by_unit
        spike_times_by_unit_out_burst = non_bursts_by_unit
        
        ## firing rates ========================================
        fr_by_unit = {unit: len(x) / (convolved_data['time_vector'][-1] - convolved_data['time_vector'][0]) for unit, x in spike_times_by_unit.items()}
        
        #in_bursts
        fr_by_unit_in_burst = {}
        mean_fr_by_unit_in_bursts = {}
        cov_fr_by_unit_in_bursts = {}
        std_fr_by_unit_in_bursts = {}
        for unit, bursts in spike_times_by_unit_in_burst.items():
            fr_by_unit_in_burst[unit] = {}
            #for i, burst in enumerate(bursts):
            for i, burst in bursts.items():
                try:
                    epoch_start = burst[0] # in seconds
                    epoch_end = burst[-1]
                    assert epoch_start < epoch_end, f'epoch_start: {epoch_start}, epoch_end: {epoch_end}'
                    assert len(burst) > 1, f'len(burst): {len(burst)}'
                    fr_by_unit_in_burst[unit][i] = len(burst) / (epoch_end - epoch_start)
                except:
                    fr_by_unit_in_burst[unit][i] = np.nan
            
            try:
                mean_fr_by_unit_in_bursts[unit] = np.nanmean(list(fr_by_unit_in_burst[unit].values()))
                std_fr_by_unit_in_bursts[unit] = np.nanstd(list(fr_by_unit_in_burst[unit].values()))
                cov_fr_by_unit_in_bursts[unit] = (std_fr_by_unit_in_bursts[unit] / mean_fr_by_unit_in_bursts[unit]) if mean_fr_by_unit_in_bursts[unit] > 0 else np.nan
            except:
                mean_fr_by_unit_in_bursts[unit] = np.nan
                std_fr_by_unit_in_bursts[unit] = np.nan
                cov_fr_by_unit_in_bursts[unit] = np.nan
        all_fr_in_bursts = np.concatenate([list(fr_by_unit_in_burst[unit].values()) for unit in fr_by_unit_in_burst.keys()])
        mean_fr_in_bursts = np.nanmean(all_fr_in_bursts)
        std_fr_in_bursts = np.nanstd(all_fr_in_bursts)
        cov_fr_in_bursts = (std_fr_in_bursts / mean_fr_in_bursts) if mean_fr_in_bursts > 0 else np.nan
        
        #out_bursts
        fr_by_unit_out_burst = {}
        mean_fr_by_unit_out_bursts = {}
        cov_fr_by_unit_out_bursts = {}
        std_fr_by_unit_out_bursts = {}
        for unit, bursts in spike_times_by_unit_out_burst.items():
            fr_by_unit_out_burst[unit] = {}
            #for i, burst in enumerate(bursts):
            for i, burst in bursts.items():
                try:
                    epoch_start = burst[0] # in seconds
                    epoch_end = burst[-1]
                    assert epoch_start < epoch_end, f'epoch_start: {epoch_start}, epoch_end: {epoch_end}'
                    assert len(burst) > 1, f'len(burst): {len(burst)}'
                    fr_by_unit_out_burst[unit][i] = len(burst) / (epoch_end - epoch_start)
                except:
                    fr_by_unit_out_burst[unit][i] = np.nan
                    
            try:
                mean_fr_by_unit_out_bursts[unit] = np.nanmean(list(fr_by_unit_out_burst[unit].values()))
                std_fr_by_unit_out_bursts[unit] = np.nanstd(list(fr_by_unit_out_burst[unit].values()))
                cov_fr_by_unit_out_bursts[unit] = (std_fr_by_unit_out_bursts[unit] / mean_fr_by_unit_out_bursts[unit]) if mean_fr_by_unit_out_bursts[unit] > 0 else np.nan
            except:
                mean_fr_by_unit_out_bursts[unit] = np.nan
                std_fr_by_unit_out_bursts[unit] = np.nan
                cov_fr_by_unit_out_bursts[unit] = np.nan
                
        all_fr_out_bursts = np.concatenate([list(fr_by_unit_out_burst[unit].values()) for unit in fr_by_unit_out_burst.keys()])
        mean_fr_out_bursts = np.nanmean(all_fr_out_bursts)
        std_fr_out_bursts = np.nanstd(all_fr_out_bursts)
        cov_fr_out_bursts = (std_fr_out_bursts / mean_fr_out_bursts) if mean_fr_out_bursts > 0 else np.nan
        
        
        ## isis ========================================
        isis_by_unit = {unit: np.diff(x) for unit, x in spike_times_by_unit.items()}
        
        #in_bursts
        isis_by_unit_in_burst = {}
        mean_isi_by_unit_in_bursts = {}
        cov_isi_by_unit_in_bursts = {}
        std_isi_by_unit_in_bursts = {}
        mean_isi_by_unit_by_burst = {}
        std_isi_by_unit_by_burst = {}
        cov_isi_by_unit_by_burst = {}
        all_isis_by_unit = {}
        break_switch = False
        for unit, bursts in spike_times_by_unit_in_burst.items():
            isis_by_unit_in_burst[unit] = {}
            mean_isi_by_unit_by_burst[unit] = {}
            std_isi_by_unit_by_burst[unit] = {}
            cov_isi_by_unit_by_burst[unit] = {}
            #for i, burst in enumerate(bursts):
            for i, burst in bursts.items():
                try:
                    epoch_start = burst[0] # in seconds
                    epoch_end = burst[-1]
                    assert epoch_start < epoch_end, f'epoch_start: {epoch_start}, epoch_end: {epoch_end}'
                    assert len(burst) > 1, f'len(burst): {len(burst)}'
                    isis_by_unit_in_burst[unit][i] = np.diff(burst)
                    
                    #if any values in isis_by_unit_in_burst[unit][i] are negative, print warning
                    if any(isis_by_unit_in_burst[unit][i] < 0):
                        print(f'Warning: Negative ISI in unit {unit}, burst {i}')
                        break_switch = True
                    
                    #stats
                    mean_isi_by_unit_by_burst[unit][i] = np.nanmean(isis_by_unit_in_burst[unit][i])
                    std_isi_by_unit_by_burst[unit][i] = np.nanstd(isis_by_unit_in_burst[unit][i])
                    cov_isi_by_unit_by_burst[unit][i] = (std_isi_by_unit_by_burst[unit][i] / mean_isi_by_unit_by_burst[unit][i]) if mean_isi_by_unit_by_burst[unit][i] > 0 else np.nan
                except:
                    isis_by_unit_in_burst[unit][i] = np.nan
                    mean_isi_by_unit_by_burst[unit][i] = np.nan
                    std_isi_by_unit_by_burst[unit][i] = np.nan
                    cov_isi_by_unit_by_burst[unit][i] = np.nan
            # if break_switch: 
            #     print('BREAK')
            #     break                    

            #extend all isis into one list for this unit:
            #all_isis_by_unit = []
            all_isis_by_unit[unit] = []
            for burst in isis_by_unit_in_burst[unit].values():
                if isinstance(burst, float):
                    if np.isnan(burst): burst = [np.nan] #enable extend to work
                    else: burst = [burst]
                all_isis_by_unit[unit].extend(burst)
                
            #stats
            mean_isi_by_unit_in_bursts[unit] = np.nanmean(all_isis_by_unit[unit])
            std_isi_by_unit_in_bursts[unit] = np.nanstd(all_isis_by_unit[unit])
            cov_fr_by_unit_in_bursts[unit] = (std_isi_by_unit_in_bursts[unit] / mean_isi_by_unit_in_bursts[unit]) if mean_isi_by_unit_in_bursts[unit] > 0 else np.nan
            
        #concat all isis from all_isis_by_unit
        all_isis_in_bursts = []
        for unit in all_isis_by_unit.keys():
            all_isis_in_bursts.extend(all_isis_by_unit[unit])
        all_isis_in_bursts = np.array(all_isis_in_bursts)
        mean_isi_in_bursts = np.nanmean(all_isis_in_bursts)
        std_isi_in_bursts = np.nanstd(all_isis_in_bursts)
        cov_isi_in_bursts = (std_isi_in_bursts / mean_isi_in_bursts) if mean_isi_in_bursts > 0 else np.nan
        
        #out_bursts
        isis_by_unit_out_burst = {}
        mean_isi_by_unit_out_bursts = {}
        cov_isi_by_unit_out_bursts = {}
        std_isi_by_unit_out_bursts = {}
        mean_isi_by_unit_by_burst = {}
        std_isi_by_unit_by_burst = {}
        cov_isi_by_unit_by_burst = {}
        all_isis_by_unit = {}
        break_switch = False
        for unit, bursts in spike_times_by_unit_out_burst.items():
            isis_by_unit_out_burst[unit] = {}
            mean_isi_by_unit_by_burst[unit] = {}
            std_isi_by_unit_by_burst[unit] = {}
            cov_isi_by_unit_by_burst[unit] = {}
            #for i, burst in enumerate(bursts):
            for i, burst in bursts.items():
                try:
                    epoch_start = burst[0] # in seconds
                    epoch_end = burst[-1]
                    assert epoch_start < epoch_end, f'epoch_start: {epoch_start}, epoch_end: {epoch_end}'
                    assert len(burst) > 1, f'len(burst): {len(burst)}'
                    isis_by_unit_out_burst[unit][i] = np.diff(burst)
                    
                    #if any values in isis_by_unit_out_burst[unit][i] are negative, print warning
                    if any(isis_by_unit_out_burst[unit][i] < 0):
                        print(f'Warning: Negative ISI in unit {unit}, burst {i}')
                        break_switch = True
                    
                    #stats
                    mean_isi_by_unit_by_burst[unit][i] = np.nanmean(isis_by_unit_out_burst[unit][i])
                    std_isi_by_unit_by_burst[unit][i] = np.nanstd(isis_by_unit_out_burst[unit][i])
                    cov_isi_by_unit_by_burst[unit][i] = (std_isi_by_unit_by_burst[unit][i] / mean_isi_by_unit_by_burst[unit][i]) if mean_isi_by_unit_by_burst[unit][i] > 0 else np.nan
                except:
                    isis_by_unit_out_burst[unit][i] = np.nan
                    mean_isi_by_unit_by_burst[unit][i] = np.nan
                    std_isi_by_unit_by_burst[unit][i] = np.nan
                    cov_isi_by_unit_by_burst[unit][i] = np.nan
            # if break_switch: 
            #     print('BREAK')
            #     break                    

            #extend all isis into one list for this unit:
            #all_isis_by_unit = []
            all_isis_by_unit[unit] = []
            for burst in isis_by_unit_out_burst[unit].values():
                if isinstance(burst, float):
                    if np.isnan(burst): burst = [np.nan]
                    else: burst = [burst]
                all_isis_by_unit[unit].extend(burst)
                
            #stats
            mean_isi_by_unit_out_bursts[unit] = np.nanmean(all_isis_by_unit[unit])
            std_isi_by_unit_out_bursts[unit] = np.nanstd(all_isis_by_unit[unit])
            cov_fr_by_unit_out_bursts[unit] = (std_isi_by_unit_out_bursts[unit] / mean_isi_by_unit_out_bursts[unit]) if mean_isi_by_unit_out_bursts[unit] > 0 else np.nan
        
        #concat all isis from all_isis_by_unit
        all_isis_out_bursts = []
        for unit in all_isis_by_unit.keys():
            all_isis_out_bursts.extend(all_isis_by_unit[unit])
        all_isis_out_bursts = np.array(all_isis_out_bursts)
        mean_isi_out_bursts = np.nanmean(all_isis_out_bursts)
        std_isi_out_bursts = np.nanstd(all_isis_out_bursts)
        cov_isi_out_bursts = (std_isi_out_bursts / mean_isi_out_bursts) if mean_isi_out_bursts > 0 else np.nan
        
        
        #all
        all_isis = np.concatenate([isis for isis in isis_by_unit.values()])
        
        # fano factor ========================================
        isi_fano_factor = np.var(all_isis) / np.mean(all_isis) if np.mean(all_isis) > 0 else np.nan
        isi_fano_factor_in_bursts = np.var(all_isis_in_bursts) / mean_isi_in_bursts if mean_isi_in_bursts > 0 else np.nan
        isi_fano_factor_out_bursts = np.var(all_isis_out_bursts) / mean_isi_out_bursts if mean_isi_out_bursts > 0 else np.nan
        
        # Burst Stats ========================================
        #TODO: Complete this section tonight # aw 2025-02-10 16:47:44
        peak_times = convolved_data['peak_times']
        ibi = np.diff(peak_times)
        mean_ibi = np.mean(ibi)
        std_ibi = np.std(ibi)
        cov_ibi = std_ibi / mean_ibi if mean_ibi > 0 else np.nan
        
        burst_summary_data = {
            'MeanWithinBurstISI': mean_isi_in_bursts,
            'CoVWithinBurstISI': cov_isi_in_bursts,
            'MeanOutsideBurstISI': mean_isi_out_bursts,
            'CoVOutsideBurstISI': cov_isi_out_bursts,
            'MeanNetworkISI': mean_isi_by_unit_in_bursts,
            'CoVNetworkISI': cov_isi_by_unit_in_bursts,
            'NumUnits': len(spike_times_by_unit),
            'Number_Bursts': len(peak_times),
            'mean_IBI': mean_ibi,
            'cov_IBI': cov_ibi,
            'isi_fano_factor': isi_fano_factor,
            'isi_fano_factor_in_bursts': isi_fano_factor_in_bursts,
            'isi_fano_factor_out_bursts': isi_fano_factor_out_bursts,
        }
        
        burst_data_by_unit = {}
        for unit in spike_times_by_unit.keys():
            burst_data_by_unit[unit] = {
                'bursts': spike_times_by_unit_in_burst[unit],
                'mean_isi_within': mean_isi_by_unit_in_bursts[unit],
                'cov_isi_within': cov_isi_by_unit_in_bursts[unit],
                'mean_isi_outside': mean_isi_by_unit_out_bursts[unit],
                'cov_isi_outside': cov_isi_by_unit_out_bursts[unit],
                'mean_isi_all': mean_isi_by_unit_in_bursts[unit],
                'cov_isi_all': cov_isi_by_unit_in_bursts[unit],
                'isis_within_bursts': isis_by_unit_in_burst[unit],
                'isis_outside_bursts': isis_by_unit_out_burst[unit],
                'isis_all': isis_by_unit[unit],
            }
        
        # burst_data_by_unit = {
        #     "bursts": {unit: bursts for unit, bursts in spike_times_by_unit_in_burst.items()},
        #     "mean_isi_within": mean_isi_in_bursts,
        #     "cov_isi_within": 
        #     "mean_isi_outside":
        #     "cov_isi_outside":
        #     "mean_isi_all":
        #     "cov_isi_all":
        #     "isis_within_bursts":
        #     "isis_outside_bursts":
        #     "isis_all":
            
        #     #added # aw 2025-02-02 22:09:13
        #     "fano_factor_within": np.var(isis_within_bursts) / mean_isi_within if mean_isi_within > 0 else np.nan,
        #     "fano_factor_outside": np.var(isis_outside_bursts) / mean_isi_outside if mean_isi_outside > 0 else np.nan,
        #     "fano_factor_all": np.var(isis) / mean_isi_all if mean_isi_all > 0 else np.nan,
        # }
        
        
        
        # Old Code ========================================
        
        # #isis           
        # isis_by_unit = {unit: np.diff(x) for unit, x in spike_times_by_unit.items()}
        # isis_within_bursts_by_unit = {unit: np.diff(x) for unit, x in spike_times_by_unit_in_burst.items()}
        # isis_outside_bursts_by_unit = {unit: np.diff(x) for unit, x in spike_times_by_unit_out_burst.items()}       
        # mean_isi_within_combined = np.mean(all_isis_within_bursts) if all_isis_within_bursts.size > 0 else np.nan
        # cov_isi_within_combined = np.cov(all_isis_within_bursts) if all_isis_within_bursts.size > 0 else np.nan
        # #fano_factor_within_combined = np.concatenate([[stats['fano_factor_within']] for stats in burst_statistics.values() if stats['fano_factor_within'] is not None])

        # mean_isi_outside_combined = np.mean(all_isis_outside_bursts) if all_isis_outside_bursts.size > 0 else np.nan
        # cov_isi_outside_combined = np.cov(all_isis_outside_bursts) if all_isis_outside_bursts.size > 0 else np.nan
        # #fano_factor_outside_combined = np.concatenate([[stats['fano_factor_outside']] for stats in burst_statistics.values() if stats['fano_factor_outside'] is not None])

        # mean_isi_all_combined = np.mean(all_isis) if all_isis.size > 0 else np.nan
        # cov_isi_all_combined = np.cov(all_isis) if all_isis.size > 0 else np.nan    
        # #fano_factor_all_combined = np.concatenate([[stats['fano_factor_all']] for stats in burst_statistics.values() if stats['fano_factor'] is not None])

        # bursting_summary_data = {}
        # bursting_summary_data['MeanWithinBurstISI'] = mean_isi_within_combined
        # bursting_summary_data['CoVWithinBurstISI'] = cov_isi_within_combined
        # bursting_summary_data['MeanOutsideBurstISI'] = mean_isi_outside_combined   
        # bursting_summary_data['CoVOutsideBurstISI'] = cov_isi_outside_combined
        # bursting_summary_data['MeanNetworkISI'] = mean_isi_all_combined
        # bursting_summary_data['CoVNetworkISI'] = cov_isi_all_combined
        # #bursting_summary_data['fano_factor_within'] = fano_factor_within_combined
        # #bursting_summary_data['fano_factor_outside'] = fano_factor_outside_combined
        # #bursting_summary_data['fano_factor_all'] = fano_factor_all_combined
        # bursting_summary_data['NumUnits'] = len(spike_times)
        # bursting_summary_data['Number_Bursts'] = sum(len(unit_stats['bursts']) for unit_stats in burst_statistics.values())
        # bursting_summary_data['mean_IBI'] = np.mean(all_isis) if all_isis.size > 0 else np.nan
        # bursting_summary_data['cov_IBI'] = np.cov(all_isis) if all_isis.size > 0 else np.nan
        
        # bursting_data = {
        #     'bursting_summary_data': bursting_summary_data,
        #     'bursting_data_by_unit': burst_statistics,
        # }
    except Exception as e:
        print(f'Error in bursting activity analysis: {e}')
        bursting_data = None

    return bursting_data

'''Older Functions (Everything underthis predates # aw 2025-02-10 11:28:41)'''

'''derive optimization targets'''
    # def build_network_metric_targets_dict(network_metrics):
        
    #     #initialize network_metric_targets
    #     bursting_data = network_metrics['bursting_data']
    #     mega_bursting_data = network_metrics['mega_bursting_data']
    #     duration_seconds = network_metrics['timeVector'][-1]
        
    #     network_metric_targets = {
    #         #General Data
    #         'source': network_metrics['source'], # 'simulated' or 'experimental'
    #         #'timeVector': network_metrics['timeVector'],
            
    #         # Spiking Data
    #         'spiking_data': {
    #             #'spike_times': network_metrics['spiking_data']['spike_times'],
    #             #'spiking_times_by_unit': network_metrics['spiking_data']['spiking_times_by_unit'],
    #             #'spiking_data_by_unit': network_metrics['spiking_data']['spiking_data_by_unit'],
    #             'spiking_summary_data': {
    #                 'MeanFireRate': {
    #                     'target': network_metrics['spiking_data']['spiking_summary_data']['MeanFireRate'],
    #                     'min': get_min(network_metrics['spiking_data']['spiking_data_by_unit'], 'FireRate'),
    #                     'max': get_max(network_metrics['spiking_data']['spiking_data_by_unit'], 'FireRate'),
    #                     #'E_I_ratio_assumption': 0.7, #E_I_ratio = 5  # 1:5 ratio of E to I neurons #TODO: need scientific basis for this assumption
                        
    #                     ## TODO: Get scientific basis for these assumptions # aw 2025-01-26 21:01:23
    #                     #E
    #                     'max_E_assumption': 5, # Hz  #outside of bursts 
    #                     'min_E_assumption': 0.5, # Hz
    #                     'max_E_assumption_inBurst': 50, # Hz #inside of bursts
    #                     'min_E_assumption_inBurst': 10, # Hz
                        
    #                     #I
    #                     'max_I_assumption': 10, # Hz
    #                     'min_I_assumption': 1, # Hz
    #                     'max_I_assumption_inBurst': 100, # Hz
    #                     'min_I_assumption_inBurst': 20, # Hz                    
                                            
    #                     'weight': 1, # TODO: update these with Nfactors
    #                 },
    #                 'CoVFireRate': {
    #                     'target': network_metrics['spiking_data']['spiking_summary_data']['CoVFireRate'],
    #                     #'min': get_min(network_metrics['spiking_data']['spiking_data_by_unit'], 'fr_CoV'),
    #                     #'max': get_max(network_metrics['spiking_data']['spiking_data_by_unit'], 'fr_CoV'),
    #                     'weight': 1, # TODO: update these with Nfactors
    #                 },
    #                 'MeanISI': {
    #                     'target': network_metrics['spiking_data']['spiking_summary_data']['MeanISI'],
    #                     'min': get_min(network_metrics['spiking_data']['spiking_data_by_unit'], 'meanISI'),
    #                     'max': get_max(network_metrics['spiking_data']['spiking_data_by_unit'], 'meanISI'),
    #                     'weight': 1, # TODO: update these with Nfactors
    #                 },
    #                 'CoV_ISI': {
    #                     'target': network_metrics['spiking_data']['spiking_summary_data']['CoV_ISI'],
    #                     # 'min': get_min(network_metrics['spiking_data']['spiking_data_by_unit'], 'isi_CoV'),
    #                     # 'max': get_max(network_metrics['spiking_data']['spiking_data_by_unit'], 'isi_CoV'),
    #                     'weight': 1, # TODO: update these with Nfactors
    #                 },
    #             },
    #         },
            
            
    #         #Bursting Data
    #         'bursting_data': {
    #             'bursting_summary_data': {
    #                 'MeanBurstRate': {
    #                     'target': bursting_data['bursting_summary_data'].get('mean_Burst_Rate'),
    #                     # 'min': get_min_burst(bursting_data['bursting_data_by_unit'], 'bursts', duration_seconds), # NOTE: Calculated as individual unit burst participation rate.
    #                     # 'max': get_max_burst(bursting_data['bursting_data_by_unit'], 'bursts', duration_seconds), # NOTE: Calculated as individual unit burst participation rate.
    #                     'min': 0,
    #                     'max': None,
    #                     'weight': 1,
    #                 },                
    #                 'MeanWithinBurstISI': {
    #                     'target': bursting_data['bursting_summary_data'].get('MeanWithinBurstISI'),
    #                     'min': get_min(bursting_data['bursting_data_by_unit'], 'mean_isi_within'),
    #                     'max': get_max(bursting_data['bursting_data_by_unit'], 'mean_isi_within'),
    #                     'weight': 1,
    #                 },
    #                 'CovWithinBurstISI': {
    #                     'target': bursting_data['bursting_summary_data'].get('CoVWithinBurstISI'),
    #                     # 'min': get_min(bursting_data['bursting_data_by_unit'], 'cov_isi_within'),
    #                     # 'max': get_max(bursting_data['bursting_data_by_unit'], 'cov_isi_within'),
    #                     'weight': 1,
    #                 },
    #                 'MeanOutsideBurstISI': {
    #                     'target': bursting_data['bursting_summary_data'].get('MeanOutsideBurstISI'),
    #                     'min': get_min(bursting_data['bursting_data_by_unit'], 'mean_isi_outside'),
    #                     'max': get_max(bursting_data['bursting_data_by_unit'], 'mean_isi_outside'),
    #                     'weight': 1,
    #                 },
    #                 'CoVOutsideBurstISI': {
    #                     'target': bursting_data['bursting_summary_data'].get('CoVOutsideBurstISI'),
    #                     # 'min': get_min(bursting_data['bursting_data_by_unit'], 'cov_isi_outside'),
    #                     # 'max': get_max(bursting_data['bursting_data_by_unit'], 'cov_isi_outside'),
    #                     'weight': 1,
    #                 },
    #                 'MeanNetworkISI': {
    #                     'target': bursting_data['bursting_summary_data'].get('MeanNetworkISI'),
    #                     'min': get_min(bursting_data['bursting_data_by_unit'], 'mean_isi_all'),
    #                     'max': get_max(bursting_data['bursting_data_by_unit'], 'mean_isi_all'),
    #                     'weight': 1,
    #                 },
    #                 'CoVNetworkISI': {
    #                     'target': bursting_data['bursting_summary_data'].get('CoVNetworkISI'),
    #                     # 'min': get_min(bursting_data['bursting_data_by_unit'], 'cov_isi_all'),
    #                     # 'max': get_max(bursting_data['bursting_data_by_unit'], 'cov_isi_all'),
    #                     'weight': 1,
    #                 },
    #                 'NumUnits': {
    #                     'target': bursting_data['bursting_summary_data'].get('NumUnits'),
    #                     # 'min': 1,
    #                     # 'max': None,
    #                     # 'weight': 1,
    #                 },
    #                 'Number_Bursts': {
    #                     'target': bursting_data['bursting_summary_data'].get('Number_Bursts'),
    #                     # 'min': 1,
    #                     # 'max': None,
    #                     # 'weight': 1, #TODO: update these with Nfactors
    #                 },
    #                 'mean_IBI': {
    #                     'target': bursting_data['bursting_summary_data'].get('mean_IBI'),
    #                     'min': None,
    #                     'max': None,
    #                     'weight': 1, #TODO: update these with Nfactors
    #                 },
    #                 'cov_IBI': {
    #                     'target': bursting_data['bursting_summary_data'].get('cov_IBI'),
    #                     'min': None,
    #                     'max': None,
    #                     'weight': 1, #TODO: update these with Nfactors
    #                 },
    #                 'mean_Burst_Peak': {
    #                     'target': bursting_data['bursting_summary_data'].get('mean_Burst_Peak'),
    #                     'min': None,
    #                     'max': None,
    #                     'weight': 1, #TODO: update these with Nfactors
    #                 },
    #                 'cov_Burst_Peak': {
    #                     'target': bursting_data['bursting_summary_data'].get('cov_Burst_Peak'),
    #                     'min': None,
    #                     'max': None,
    #                     'weight': 1, #TODO: update these with Nfactors
    #                 },
    #                 'fano_factor': {
    #                     'target': bursting_data['bursting_summary_data'].get('fano_factor'),
    #                     'min': None,
    #                     'max': None,
    #                     'weight': 1, #TODO: update these with Nfactors
    #                 },
    #                 'baseline': {
    #                     'target': bursting_data['bursting_summary_data'].get('baseline'),
    #                     'min': None,
    #                     'max': None,
    #                     'weight': 1,
    #                 },
    #             },
    #         },
            
    #         #Mega Bursting Data
    #         'mega_bursting_data': {
    #             'bursting_summary_data': {
    #                 'MeanBurstRate': { # WHole network burst rate
    #                     'target': mega_bursting_data['bursting_summary_data'].get('mean_Burst_Rate'),
    #                     # 'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'bursts'),
    #                     # 'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'bursts'),
    #                     'min' : 0,
    #                     'max' : None,
    #                     'weight': 1,
                    

    #                     # 'min': get_min_burst(mega_bursting_data['bursting_data_by_unit'], 'bursts', duration_seconds), # NOTE: Calculated as individual unit burst participation rate.
    #                     # 'max': get_max_burst(mega_bursting_data['bursting_data_by_unit'], 'bursts', duration_seconds), # NOTE: Calculated as individual unit burst participation rate.
    #                     # 'weight': 1,
                        
    #                     # aw 2025-01-25 21:12:43 - #TODO: need to maybe rethink how I"m getting min and maxes for a lot of these data.
    #                     # like in this case, using min and max burst participation rate for individual units to set min and max for
    #                     # whole network burst rate doesnt make sense. Unless at least one unit participates in al bursts - max for individual
    #                     # units will always be less than burst rate for whole network.
    #                 },                  
    #                 'MeanWithinBurstISI': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('MeanWithinBurstISI'),
    #                     'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'mean_isi_within'),
    #                     'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'mean_isi_within'),
    #                     # 'min' : 0,
    #                     # 'max' : None,
    #                     'weight': 1,
    #                 },
    #                 'CovWithinBurstISI': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('CoVWithinBurstISI'),
    #                     # 'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'cov_isi_within'),
    #                     # 'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'cov_isi_within'),
    #                     'weight': 1,
    #                 },
    #                 'MeanOutsideBurstISI': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('MeanOutsideBurstISI'),
    #                     'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'mean_isi_outside'),
    #                     'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'mean_isi_outside'),
    #                     'weight': 1,
    #                 },
    #                 'CoVOutsideBurstISI': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('CoVOutsideBurstISI'),
    #                     # 'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'cov_isi_outside'),
    #                     # 'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'cov_isi_outside'),
    #                     'weight': 1,
    #                 },
    #                 'MeanNetworkISI': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('MeanNetworkISI'),
    #                     'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'mean_isi_all'),
    #                     'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'mean_isi_all'),
    #                     'weight': 1,
    #                 },
    #                 'CoVNetworkISI': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('CoVNetworkISI'),
    #                     # 'min': get_min(mega_bursting_data['bursting_data_by_unit'], 'cov_isi_all'),
    #                     # 'max': get_max(mega_bursting_data['bursting_data_by_unit'], 'cov_isi_all'),
    #                     'weight': 1,
    #                 },
    #                 'NumUnits': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('NumUnits'),
    #                     # 'min': 1,
    #                     # 'max': None,
    #                     # 'weight': 1,
    #                 },
    #                 'Number_Bursts': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('Number_Bursts'),
    #                     # 'min': 1,
    #                     # 'max': None,
    #                     # 'weight': 1, #TODO: update these with Nfactors
    #                 },
    #                 'mean_IBI': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('mean_IBI'),
    #                     'min': None,
    #                     'max': None,
    #                     'weight': 1, #TODO: update these with Nfactors
    #                 },
    #                 'cov_IBI': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('cov_IBI'),
    #                     'min': None,
    #                     'max': None,
    #                     'weight': 1, #TODO: update these with Nfactors
    #                 },
    #                 'mean_Burst_Peak': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('mean_Burst_Peak'),
    #                     'min': None,
    #                     'max': None,
    #                     'weight': 1, #TODO: update these with Nfactors
    #                 },
    #                 'cov_Burst_Peak': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('cov_Burst_Peak'),
    #                     'min': None,
    #                     'max': None,
    #                     'weight': 1, #TODO: update these with Nfactors
    #                 },
    #                 'fano_factor': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('fano_factor'),
    #                     'min': None,
    #                     'max': None,
    #                     'weight': 1, #TODO: update these with Nfactors
    #                 },
    #                 'baseline': {
    #                     'target': mega_bursting_data['bursting_summary_data'].get('baseline'),
    #                     'min': None,
    #                     'max': None,
    #                     'weight': 1,
    #                 },
    #             },        
    #         }
    #     }
    #     return network_metric_targets

def get_min(data, key):
    data_min = min(unit[key] for unit in data.values() if unit[key] is not None and unit[key] > 0 and not (
        unit[key] != unit[key] or unit[key] == float('inf') or unit[key] == float('-inf')
        ))
    
    #expand the code above into a more readable format
    # data_min = []
    # for unit in data.values():
    #     if unit[key] is not None and unit[key] > 0 and not (unit[key] != unit[key] or unit[key] == float('inf') or unit[key] == float('-inf')):
    #         data_min.append(unit[key])
    # data_min = min(data_min)
    
    
    
    #if data is an array with one value, return that value
    # if isinstance(data_min, list) and len(data_min) == 1:
    #     data_min = data_min[0]
    # try: 
    #     len(data_min)
    #     print('data_min is an array')
    # except TypeError: pass
    data_min = handle_numpy_float64(data_min)
    #print('data_min:', data_min)
    return data_min

def get_max(data, key):
    data_list = [unit[key] for unit in data.values() if not (unit[key] is None or unit[key] != unit[key] or unit[key] == float('inf') or unit[key] == float('-inf'))]
    data_max = max(unit[key] for unit in data.values() if not (unit[key] is None or unit[key] != unit[key] or unit[key] == float('inf') or unit[key] == float('-inf')))
    #if data is an array with one value, return that value
    # if type(data_max) == list:
    #     data_max = data_max[0]
    data_max = handle_numpy_float64(data_max)
    return data_max

def get_min_burst(data, key, duration_seconds):
    
    #start min_rate at positive infinity
    min_rate = float('inf')
    for unit in data.values():
        num_bursts = len(unit[key])
        burst_participation_rate = num_bursts / duration_seconds        
        try:
            if burst_participation_rate < min_rate:
                min_rate = burst_participation_rate
        except:
            continue
    
    return min_rate

def get_max_burst(data, key, duration_seconds):
    
    #start min_rate at positive infinity
    #min_rate = float('inf')
    max_rate = 0
    for unit in data.values():
        num_bursts = len(unit[key])
        burst_participation_rate = num_bursts / duration_seconds        
        try:
            if burst_participation_rate > max_rate:
                max_rate = burst_participation_rate
        except:
            continue
    
    return max_rate  

def save_network_metric_dict_with_timestamp(network_metrics, network_metrics_targets, output_dir):
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    date = datetime.now().strftime("%Y%m%d")
    #output_path = output_dir + '/fitness_args_' + timestamp + '.py'
    output_path = os.path.join(output_dir, f'{date}_features.py')
    
    output_path = output_path.replace('-', '_') #conver all '-' to '_' in the timestamp
    print('output_path:', output_path)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    def format_dict(d, indent=0):
        formatted_str = ''        
        if indent==0: formatted_str += 'import numpy as np\n' # add import statement to the top of the file
        for key, value in d.items():
            if isinstance(value, dict):
                formatted_str += ' ' * indent + f"'{key}': {{\n" + format_dict(value, indent + 4) + ' ' * indent + '},\n'
                
                # try: key = int(key)
                # except: pass
                
                # if isinstance(key, str):
                #     formatted_str += ' ' * indent + f"'{key}': {{\n" + format_dict(value, indent + 4) + ' ' * indent + '},\n'
                # elif isinstance(key, int):
                #     formatted_str += ' ' * indent + f"{key}: {{\n" + format_dict(value, indent + 4) + ' ' * indent + '},\n'
                # else:
                #     print('Warning: Unrecognized key type:', type(key), 'key:', key)
                #     formatted_str += ' ' * indent + f"'{key}': {{\n" + format_dict(value, indent + 4) + ' ' * indent + '},\n'
            elif isinstance(value, list):
                # convert list to string representation
                list_str = ', '.join([str(item) for item in value])
                formatted_str += ' ' * indent + f"'{key}': [{list_str}],\n"
                
                # formatted_str += ' ' * indent + f"'{key}': ["
                # for i, item in enumerate(value):
                #     formatted_str += ' ' + f"{item},"
                #     # if i < len(value) - 1:
                #     #     formatted_str += ','
                #     #formatted_str += '\n'
                # # remove last comma
                # formatted_str = formatted_str[:-1]
                # formatted_str += '],\n'
            elif isinstance(value, (int, float)):
                formatted_str += ' ' * indent + f"'{key}': {value},\n"
            elif isinstance(value, str):
                formatted_str += ' ' * indent + f"'{key}': '{value}',\n"
            elif isinstance(value, np.ndarray):
                #formatted_str += ' ' * indent + f"'{key}': np.array({value.tolist()}),\n"
                #convert to tuple
                formatted_str += ' ' * indent + f"'{key}': {tuple(value.tolist())},\n"
            elif isinstance(value, tuple):
                formatted_str += ' ' * indent + f"'{key}': {value},\n"
            elif isinstance(value, bool):
                formatted_str += ' ' * indent + f"'{key}': {value},\n"
            elif isinstance(value, None.__class__):
                formatted_str += ' ' * indent + f"'{key}': None,\n"
            else:
                print(type(value))
                print('Warning: Unrecognized data type for key:', key, 'value:', value)
                formatted_str += ' ' * indent + f"'{key}': {value},\n"
        return formatted_str

    # Convert network_metric_targets to a formatted string
    formatted_fitness_args = 'fitness_args = {\n' + format_dict(network_metrics_targets, 4) + '}'
    
    # append some additional info to the formatted string
    formatted_fitness_args += (
        f'\n\n#=== FITNESS FUNCTION ARGUMENTS ===\n'
        f'# fitness function\n'
        f'fitnessFuncArgs = {{}}\n'
        f'fitnessFuncArgs[\'targets\'] = fitness_args\n'
        f'number_of_units = fitness_args[\'bursting_data\'][\'bursting_summary_data\'][\'NumUnits\'][\'target\']\n'
        f'fitnessFuncArgs[\'features\'] = {{\n'
        f'  # tweaking the ratio a bit based on Roy\'s suggestion\n'
        f'  \'num_excite\': len(fitness_args[\'excit_units\']),\n'
        f'  \'num_inhib\': len(fitness_args[\'inhib_units\']),\n'
        f'}}\n'
        f'fitnessFuncArgs[\'maxFitness\'] = 1000\n'
        f'#fitnessFuncArgs[\'skip_existing\'] = False\n'
    )

    with open(output_path, 'w') as f:
        f.write(formatted_fitness_args)

    print('Updated fitness args saved to:', output_path)

def handle_numpy_float64(data):
    try: 
        if isinstance(data, np.ndarray):
            data = data.tolist()
    except: 
        pass
    return data

'''sim funcs'''
def get_simulated_network_activity_metrics(simData=None, popData=None, **kwargs):
    import time
    
    #candidate_label = kwargs.get('candidate_label', None)
    print('') #for formatting
    print('Calculating Network Activity Metrics for Simulated Data...')
    start_time = time.time()
    assert simData is not None, 'No simData provided'
    #this part should be useful for fitness during simulation
    # if simData is None:
    #     try: 
    #         simData = sim.simData
    #         print('Using simData from netpyne.sim.simData')
    #     except:
    #         print('No simData provided or found in netpyne.sim.simData')
    #         return None

    #initialize network_data
    #network_data = init_network_data_dict()
    # # HACK:
    # from RBS_network_models.network_analysis import init_network_data_dict, extract_bursting_activity_data
    # # HACK:
    network_data = init_network_data_dict()
    rasterData = simData.copy()
    
    #check if rasterData['spkt'] is empty
    if len(rasterData['spkt']) == 0:
        print('No spike times found in rasterData')
        return None
    
    #convert time to seconds - get initially available data
    spike_times = np.array(rasterData['spkt']) / 1000
    timeVector = np.array(rasterData['t']) / 1000
    spike_times_by_unit = {int(i): spike_times[rasterData['spkid'] == i] for i in np.unique(rasterData['spkid'])} #mea_analysis_pipeline.py expects spike_times as dictionary    
    
    #extract spiking metrics from simulated data
    try: 
        extract_metrics_from_simulated_data(spike_times, timeVector, spike_times_by_unit, rasterData, popData, **kwargs)
    except Exception as e:
        print(f'Error extracting metrics from simulated data: {e}')
        pass
    
    #extract bursting metrics from simulated data (but this one works for both simulated and experimental data)
    try: 
        extract_bursting_activity_data(spike_times, spike_times_by_unit, **kwargs)
    except Exception as e:
        print(f'Error calculating bursting activity: {e}')
        pass
    
    #return network_data
    
    def convert_single_element_arrays(data):
        if isinstance(data, dict):
            for key, value in data.items():
                data[key] = convert_single_element_arrays(value)
        elif isinstance(data, np.ndarray) and data.size == 1:
            return data.item()
        return data

    network_data = convert_single_element_arrays(network_data)    
    print('Network Activity Metrics Calculated and Extracted!')
    print(f'Elapsed time: {time.time() - start_time} seconds')
    print('') #for formatting    
    return network_data

def get_CoV_fr_simulated(spike_times, total_duration, window_size=1.0): 
    """
    Calculate the Coefficient of Variation of Firing Rate (CoV FR) over time windows.

    Parameters:
    - spike_times: List or array of spike times for a single unit (in seconds).
    - total_duration: Total duration of the simulation (in seconds).
    - window_size: Size of the time window for calculating firing rates (default: 1.0 second).

    Returns:
    - CoV_fr: Coefficient of Variation of Firing Rate (float).
    """
    #TODO import window from convolution_params?
    if len(spike_times) == 0:
        # If no spikes, CoV FR is undefined (return NaN)
        return np.nan

    # Divide the total duration into non-overlapping time windows
    #window_size = 0.1
    total_duration = total_duration[-1]
    num_windows = int(total_duration / window_size)
    window_edges = np.linspace(0, total_duration, num_windows + 1)

    # Count spikes in each window
    spike_counts = np.histogram(spike_times, bins=window_edges)[0]

    # Convert spike counts to firing rates (spikes per second)
    firing_rates = spike_counts / window_size

    # Compute CoV FR: standard deviation divided by mean
    mean_fr = np.mean(firing_rates)
    std_fr = np.std(firing_rates)

    # Avoid division by zero if mean_fr is 0
    if mean_fr == 0:
        return np.nan

    CoV_fr = std_fr / mean_fr
    return CoV_fr

def calculate_simulated_network_activity_metrics(spike_times_by_unit):
    # Initialize spiking data dictionary
    network_data['simulated_data']['spiking_data_by_unit'] = {}
    
    # Iterate over each unit to calculate individual metrics
    for unit, spike_times in spike_times_by_unit.items():
        if len(spike_times) < 2:
            meanISI = np.nan
            CoV_ISI = np.nan
        else:
            isi = np.diff(spike_times)
            meanISI = np.mean(isi)
            CoV_ISI = np.std(isi) / meanISI  # Correct CoV calculation for ISI
        
        fr = len(spike_times) / network_data['timeVector'][-1]  # Firing rate (spikes per second)
        
        # Calculate CoV FR using the get_CoV_fr function with a 1-second time window
        CoV_fr = get_CoV_fr_simulated(spike_times, network_data['timeVector'])

        # Store calculated metrics for each unit
        network_data['simulated_data']['spiking_data_by_unit'][unit] = {
            'FireRate': fr,
            'CoV_fr': CoV_fr,
            'meanISI': meanISI,
            'CoV_ISI': CoV_ISI,
            'spike_times': spike_times,
        }
        
        # Determine population type
        if unit in network_data['simulated_data']['E_Gids']:
            network_data['simulated_data']['spiking_data_by_unit'][unit]['pop'] = 'E'
        elif unit in network_data['simulated_data']['I_Gids']:
            network_data['simulated_data']['spiking_data_by_unit'][unit]['pop'] = 'I'
        else:
            raise ValueError(f'Unit {unit} not found in E_Gids or I_Gids')

    # Extract E and I gids
    E_Gids = network_data['simulated_data']['E_Gids']
    I_Gids = network_data['simulated_data']['I_Gids']

    # Calculate mean and CoV metrics for excitatory and inhibitory populations
    E_CoV_ISI = np.nanmean([network_data['simulated_data']['spiking_data_by_unit'][unit]['CoV_ISI'] 
                            for unit in E_Gids if unit in network_data['simulated_data']['spiking_data_by_unit']])
    I_CoV_ISI = np.nanmean([network_data['simulated_data']['spiking_data_by_unit'][unit]['CoV_ISI'] 
                            for unit in I_Gids if unit in network_data['simulated_data']['spiking_data_by_unit']])
    MeanISI_E = np.nanmean([network_data['simulated_data']['spiking_data_by_unit'][unit]['meanISI'] 
                            for unit in E_Gids if unit in network_data['simulated_data']['spiking_data_by_unit']])
    MeanISI_I = np.nanmean([network_data['simulated_data']['spiking_data_by_unit'][unit]['meanISI'] 
                            for unit in I_Gids if unit in network_data['simulated_data']['spiking_data_by_unit']])
    CoVFiringRate_E = np.nanmean([network_data['simulated_data']['spiking_data_by_unit'][unit]['CoV_fr'] 
                                  for unit in E_Gids if unit in network_data['simulated_data']['spiking_data_by_unit']])
    CoVFiringRate_I = np.nanmean([network_data['simulated_data']['spiking_data_by_unit'][unit]['CoV_fr'] 
                                  for unit in I_Gids if unit in network_data['simulated_data']['spiking_data_by_unit']])
    
    # Add population-level metrics to the network data
    network_data['simulated_data']['CoV_ISI_E'] = E_CoV_ISI
    network_data['simulated_data']['CoV_ISI_I'] = I_CoV_ISI
    network_data['simulated_data']['MeanISI_E'] = MeanISI_E
    network_data['simulated_data']['MeanISI_I'] = MeanISI_I
    network_data['simulated_data']['CoVFireRate_E'] = CoVFiringRate_E
    network_data['simulated_data']['CoVFireRate_I'] = CoVFiringRate_I
            
def extract_metrics_from_simulated_data(spike_times, timeVector, spike_times_by_unit, rasterData, popData, **kwargs):
    
    #add immediately available data to network_data
    network_data['source'] = 'simulated'
    network_data['timeVector'] = timeVector #seconds
    network_data['spiking_data']['spike_times'] = spike_times
    network_data['simulated_data']['soma_voltage'] = rasterData['soma_voltage']
    network_data['spiking_data']['spiking_times_by_unit'] = spike_times_by_unit

    #add uniquely simulated data to network_data
    network_data['simulated_data']['E_Gids'] = popData['E']['cellGids']
    network_data['simulated_data']['I_Gids'] = popData['I']['cellGids']
    network_data['simulated_data']['MeanFireRate_E'] = rasterData['popRates']['E']
    network_data['simulated_data']['MeanFireRate_I'] = rasterData['popRates']['I']
    
    #get simulated network activity metrics
    calculate_simulated_network_activity_metrics(spike_times_by_unit)
    
    #derive remaining spiking data from simulated data
    num_E = len(network_data['simulated_data']['E_Gids'])
    num_I = len(network_data['simulated_data']['I_Gids'])
    
    #overall mean firing rate
    E_popRates = rasterData['popRates']['E']
    I_popRates = rasterData['popRates']['I']
    network_data['spiking_data']['spiking_summary_data']['MeanFireRate'] = (
        (E_popRates * num_E + I_popRates * num_I) / (num_E + num_I)
        )
    
    #overall CoV firing rate
    E_CoV = network_data['simulated_data']['CoVFireRate_E']
    I_CoV = network_data['simulated_data']['CoVFireRate_I']
    network_data['spiking_data']['spiking_summary_data']['CoVFireRate'] = (
        (E_CoV * num_E + I_CoV * num_I) / (num_E + num_I)
        )
    
    #overall mean ISI
    E_meanISI = network_data['simulated_data']['MeanISI_E']
    I_meanISI = network_data['simulated_data']['MeanISI_I']
    network_data['spiking_data']['spiking_summary_data']['MeanISI'] = (
        (E_meanISI * num_E + I_meanISI * num_I) / (num_E + num_I)
        )
    
    #overall CoV ISI
    E_CoV_ISI = network_data['simulated_data']['CoV_ISI_E']
    I_CoV_ISI = network_data['simulated_data']['CoV_ISI_I']
    network_data['spiking_data']['spiking_summary_data']['CoV_ISI'] = (
        (E_CoV_ISI * num_E + I_CoV_ISI * num_I) / (num_E + num_I)
        )
    
    #spiking data by unit - just go through the simulated version of the data but de-identify pop
    simulated_spiking_data_by_unit = network_data['simulated_data']['spiking_data_by_unit']
    spiking_data_by_unit = {}
    for unit, unit_data in simulated_spiking_data_by_unit.items():
        spiking_data_by_unit[unit] = {
            'FireRate': unit_data['FireRate'],
            'CoV_fr': unit_data['CoV_fr'],
            'meanISI': unit_data['meanISI'],
            'CoV_ISI': unit_data['CoV_ISI'],
            'spike_times': unit_data['spike_times'],
        }
    network_data['spiking_data']['spiking_data_by_unit'] = spiking_data_by_unit

def get_simulated_network_activity_metrics(conv_params, mega_params, simData=None, popData=None, **kwargs):
    import time
    
    #initialize network_data
    network_data = init_network_data_dict()
    
    #add conv_params to kwargs
    kwargs['conv_params'] = conv_params
    kwargs['mega_params'] = mega_params
    
    #candidate_label = kwargs.get('candidate_label', None)
    print('') #for formatting
    print('Calculating Network Activity Metrics for Simulated Data...')
    start_time = time.time()
    assert simData is not None, 'No simData provided'
    
    # copy simData
    rasterData = simData.copy()
    
    #check if rasterData['spkt'] is empty
    if len(rasterData['spkt']) == 0:
        print('No spike times found in rasterData')
        return None
    
    #convert time to seconds - get initially available data
    spike_times = np.array(rasterData['spkt']) / 1000
    timeVector = np.array(rasterData['t']) / 1000
    spike_times_by_unit = {int(i): spike_times[rasterData['spkid'] == i] for i in np.unique(rasterData['spkid'])} #mea_analysis_pipeline.py expects spike_times as dictionary    
    
    #extract spiking metrics from simulated data
    try: 
        extract_metrics_from_simulated_data(spike_times, timeVector, spike_times_by_unit, rasterData, popData, **kwargs)
    except Exception as e:
        print(f'Error extracting metrics from simulated data: {e}')
        pass
    
    #extract bursting metrics from simulated data (but this one works for both simulated and experimental data)
    try: 
        extract_bursting_activity_data(spike_times, spike_times_by_unit, **kwargs)
    except Exception as e:
        print(f'Error calculating bursting activity: {e}')
        pass
    
    #return network_data
    
    def convert_single_element_arrays(data):
        if isinstance(data, dict):
            for key, value in data.items():
                data[key] = convert_single_element_arrays(value)
        elif isinstance(data, np.ndarray) and data.size == 1:
            return data.item()
        return data

    network_data = convert_single_element_arrays(network_data)    
    print('Network Activity Metrics Calculated and Extracted!')
    print(f'Elapsed time: {time.time() - start_time} seconds')
    print('') #for formatting    
    return network_data

'''Plotting Functions'''
def plot_network_summary(
    network_metrics,
    bursting_plot_path=None,
    bursting_fig_path=None,
    save_path=None,
    mode='2p',
    limit_seconds = None,
    plot_class = False,    
    ):
    # assertions
    assert save_path is not None, f"Error: save_path must be specified in plot_network_activity"
    
    #init paths
    if bursting_plot_path is None and bursting_fig_path is None:
        
        #HACK: dumb way to to allow for pdf or png save paths - fix later I guess.
        if '.pdf' in save_path: save_path = save_path.replace('.pdf', '.npy') #switch it back to .npy so that it gets handled correctly
        if '.png' in save_path: save_path = save_path.replace('.png', '.npy') 
        
        # summary plot path
        summary_plot_path = save_path.replace('.npy', '.pdf')
        summary_plot_path_png = save_path.replace('.npy', '.png')
        
 
    
    # import plotting functions
    from MEA_Analysis.NetworkAnalysis.awNetworkAnalysis.network_analysis import plot_raster_plot_experimental, plot_network_bursting_experimental
    
    # init metrics
    # spiking_data_by_unit = network_metrics['spiking_data']['spiking_data_by_unit']
    # bursting_ax = network_metrics['bursting_data']['bursting_summary_data']['ax']
    # mega_bursting_ax = network_metrics['mega_bursting_data']['bursting_summary_data']['ax']
    bursting_ax = network_metrics['bursting_data']['ax']
    mega_bursting_ax = network_metrics['mega_bursting_data']['ax'] 
    spiking_data_by_unit = network_metrics['spiking_data']['spiking_metrics_by_unit']    
    unit_types = network_metrics['unit_types']
    
    # choose mode:
    # mode = '2p' - 2 panels
    # mode = '3p' - 3 panels
    if mode == '2p':
        
        # init plot
        fig, ax = plt.subplots(2, 1, figsize=(16, 9))
        
        # plot raster plot
        print("Generating raster plot...")
        #spiking_data_by_unit = network_metrics['spiking_data']['spiking_data_by_unit']
        ax[0] = plot_raster_plot_experimental(ax[0], spiking_data_by_unit, unit_types=unit_types, plot_class=plot_class)

        # plot network bursting plots
        print("Generating network bursting plot...")
        ax[1] = plot_network_bursting_experimental(ax[1], bursting_ax, mega_ax=mega_bursting_ax)
        #mega_bursting_ax = mega_metrics['bursting_data']['bursting_summary_data']['ax']
        #ax[1] = plot_network_bursting_experimental(ax[1], mega_bursting_ax)
        
        # adjust axes
        print("Adjusting axes...")
        # ensure x-axes are the same for both plots
        # get the narrowes x-axis range shared by both plots
        x_min = min([ax[0].get_xlim()[0], ax[1].get_xlim()[0]])
        x_max = max([ax[0].get_xlim()[1], ax[1].get_xlim()[1]])
        ax[0].set_xlim(x_min, x_max)
        ax[1].set_xlim(x_min, x_max)
        
        ## override y axis to 17.5
        ## ==================
        ax[0].set_ylim(0, 17.5)
        ax[1].set_ylim(0, 17.5)
        #ax[2].set_ylim(0, 17.5)
        ##==================
        
        plt.tight_layout()
        
    elif mode == '3p':
        # init plot
        fig, ax = plt.subplots(3, 1, figsize=(16, 9))
        
        # plot raster plot
        print("Generating raster plot...")
        #spiking_data_by_unit = network_metrics['spiking_data']['spiking_data_by_unit']
        ax[0] = plot_raster_plot_experimental(ax[0], spiking_data_by_unit, unit_types=unit_types, plot_class=plot_class)
        
        # plot network bursting plots
        print("Generating network bursting plot...")
        ax[1] = plot_network_bursting_experimental(ax[1], bursting_ax)
        
        # plot network bursting plots
        print("Generating mega network bursting plot...")
        ax[2] = plot_network_bursting_experimental(ax[2], mega_bursting_ax, mode='mega')
        
        # adjust axes
        print("Adjusting axes...")
        # ensure x-axes are the same for both plots
        # get the narrowes x-axis range shared by both plots
        x_min = min([ax[0].get_xlim()[0], ax[1].get_xlim()[0], ax[2].get_xlim()[0]])
        x_max = max([ax[0].get_xlim()[1], ax[1].get_xlim()[1], ax[2].get_xlim()[1]])
        ax[0].set_xlim(x_min, x_max)
        ax[1].set_xlim(x_min, x_max)
        ax[2].set_xlim(x_min, x_max)
        
        ## override y axis to 17.5
        ## ==================
        ax[0].set_ylim(0, 17.5)
        ax[1].set_ylim(0, 17.5)
        ax[2].set_ylim(0, 17.5)
        ##==================
        
        plt.tight_layout()
        
    # # #for debugging, limit x-axis to 35 seconds (in ms)
    # try:
    #     ax[0].set_xlim(0, 35)
    #     ax[1].set_xlim(0, 35)
    #     ax[2].set_xlim(0, 35) # just using try except to catch the error if it's a 2 panel plot
    # except:
    #     pass
    
    # set limit on time as needed
    if limit_seconds is not None:
        try:
            ax[0].set_xlim(0, limit_seconds)
            ax[1].set_xlim(0, limit_seconds)
            ax[2].set_xlim(0, limit_seconds) # just using try except to catch the error if it's a 2 panel plot
        except:
            pass
            
    # make sure the plot path exists
    # if not os.path.exists(os.path.dirname(raster_plot_path)):
    #     os.makedirs(os.path.dirname(raster_plot_path), exist_ok=True)
    
    # # save plot as pdf
    # fig.savefig(bursting_plot_path) #save as pdf
    # print(f"Network summary plot saved to {bursting_plot_path}")
    
    # # save plot as png
    # fig.savefig(bursting_plot_path_png, dpi=600) #save as png
    # print(f"Network summary plot saved to {bursting_plot_path_png}")
    
    # save plot as pdf
    fig.savefig(summary_plot_path) #save as pdf
    print(f"Network summary plot saved to {summary_plot_path}")
    
    # save plot as png
    fig.savefig(summary_plot_path_png, dpi=600) #save as png
    print(f"Network summary plot saved to {summary_plot_path_png}")   

def plot_raster(ax, spiking_data_by_unit, unit_ids=None, E_gids=None, I_gids=None, data_type='simulated'):
        """
        Plot a raster plot of spiking activity.

        Parameters:
            ax (matplotlib.axes.Axes): The matplotlib axis to draw the plot on.
            spiking_data_by_unit (dict): A dictionary where keys are neuron IDs (GIDs) 
                                        and values are spike time data for each unit.
            E_gids (list): List of excitatory neuron IDs.
            I_gids (list): List of inhibitory neuron IDs.

        Returns:
            matplotlib.axes.Axes: The axis with the raster plot.
        """
        
        #main logic
        if data_type == 'simulated':
            list_of_gids = [E_gids, I_gids]
            list_of_colors = ['y.', 'b.']
            list_of_labels = ['Excitatory', 'Inhibitory']
            def plot_mixed(list_of_gids, list_of_colors, list_of_labels):
                assert E_gids is not None, "E_gids must be provided"
                assert I_gids is not None, "I_gids must be provided"

                # Create a dict of all gids and their respective colors and labels
                gid_color_label_dict = {gid: (color, label) for gids, color, label in zip(list_of_gids, list_of_colors, list_of_labels) for gid in gids}
                
                #sort all gids by fr
                sorted_spike_data = {k: v for k, v in sorted(spiking_data_by_unit.items(), key=lambda item: item[1]['FireRate'])}
                sorted_key = 0
                for gid, data in sorted_spike_data.items():
                    spike_times = data['spike_times']
                    #sometimes spike_times is a single value, not a list - make it a list
                    if isinstance(spike_times, (int, float)):
                        spike_times = [spike_times]
                    color, label = gid_color_label_dict[gid]
                    ax.plot(spike_times, [sorted_key] * len(spike_times), f'{color}', markersize=2, label=label)
                    sorted_key += 1
            plot_mixed(list_of_gids, list_of_colors, list_of_labels)
            ax.set_title('Simulated Raster Plot')
            # ax.set_xlabel('Time (s)')
            ax.set_ylabel('Neuron GID')
        elif data_type == 'experimental':            
            color = 'k.'
            label = 'Experimental'
            def plot_exprimental(unit_ids, color, label):           
                #sort
                sorted_spike_data = {k: v for k, v in sorted(spiking_data_by_unit.items(), key=lambda item: item[1]['FireRate'])}
                sorted_key = 0
                for gid, data in sorted_spike_data.items():
                    spike_times = data['spike_times']
                    ax.plot(spike_times, [sorted_key] * len(spike_times), f'{color}', markersize=1, label=label)
                    sorted_key += 1
                
                # #plot
                # for gid in unit_ids:
                #     if gid in spiking_data_by_unit:
                #         spike_times = spiking_data_by_unit[gid]['spike_times']
                #         if isinstance(spike_times, (int, float)):
                #             spike_times = [spike_times]
                #         #ax.plot(spike_times, [gid] * len(spike_times), f'{color}', markersize=1, label=label if gid == gids[0] else None)
                #         ax.plot(spike_times, [gid] * len(spike_times), f'{color}', markersize=1, label=label)
            plot_exprimental(unit_ids, color, label)
            #ax.set_xlabel('Time (s)')
            ax.set_ylabel('Unit ID')
            ax.set_title('Reference Raster Plot')
        else:
            raise ValueError(f"Invalid data_type: {data_type}")
        
        ax.set_title('Raster Plot')
        ax.legend()
        
        # Set the marker size in the legend
        legend = ax.legend()
        for handle in legend.legendHandles:
            #handle._legmarker.set_markersize(5)
            handle._markersize = 5

        # only show unique legend items
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
           
        plt.tight_layout()
        
        # HACK: there must be an easier way to do this...
        # fix x and y lims
        true_max_x = 0
        true_min_x = 0
        for line in ax.lines:
            try:
                x_data = line.get_xdata()
                max_x = max(x_data)
                min_x = min(x_data)
                if max_x > true_max_x:
                    true_max_x = max_x
                if min_x < true_min_x:
                    true_min_x = min_x
            except:
                continue
        ax.set_xlim(true_min_x, true_max_x)
        
        true_max_y = 0
        true_min_y = 0
        for line in ax.lines:
            try:
                y_data = line.get_ydata()
                max_y = max(y_data)
                min_y = min(y_data)
                if max_y > true_max_y:
                    true_max_y = max_y
                if min_y < true_min_y:
                    true_min_y = min_y
            except:
                continue
            ax.set_ylim(true_min_y, true_max_y)       
        
        return ax

def plot_network_summary_slide(sim_data_path, 
                                #convolution_params_path, 
                                #reference_data_npy, 
                                average_fitness,
                                trim_start=0,
                                ):
        start_network_metrics = time()
        #from CDKL5_DIV21.src.conv_params import conv_params
        def get_network_metrics():
            #convolution_params = import_module_from_path(CONVOLUTION_PARAMS)
            #from DIV21.src.conv_params import conv_params
            from RBS_network_models.CDKL5.DIV21.src.conv_params import conv_params
            #from RBS_network_models.developing.CDKL5.DIV21.src.conv_params import conv_params
            #get network metrics
            kwargs = {
                'simData': sim.allSimData,
                'simConfig': sim.cfg,
                #'conv_params': convolution_params.conv_params,
                'conv_params': conv_params,
                'popData': sim.net.allPops,
                'cellData': sim.net.allCells, #not actually used in the fitness calculation, but whatever
            }
            #from DIV21.utils.fitness_helper import calculate_network_metrics
            #from RBS_network_models.developing.utils.fitness_helper import calculate_network_metrics
            from RBS_network_models.utils.fitness_helper import calculate_network_metrics
            error, kwargs = calculate_network_metrics(kwargs)
            #save network metrics
            network_metrics_path = sim_data_path.replace('_data', '_network_metrics')
            if '.pkl' in network_metrics_path:
                network_metrics_path = network_metrics_path.replace('.pkl', '.npy')
            elif '.json' in network_metrics_path:
                network_metrics_path = network_metrics_path.replace('.json', '.npy')
            np.save(network_metrics_path, kwargs)
            # with open(network_metrics_path, 'w') as f:
            #     json.dump(kwargs, f, indent=4)
            return error, kwargs
        error, kwargs = get_network_metrics()
        if error: return error
        # privileged_print("\tNetwork metrics calculated - kwargs dict created.")
        # privileged_print(f'\tTime taken: {time()-start_network_metrics} seconds')
        print("\tNetwork metrics calculated - kwargs dict created.")
        print(f'\tTime taken: {time()-start_network_metrics} seconds')

        #plot raster
        start_raster = time()
        def plot_simulated_raster_wrapper(ax = None, subplot=False, trim_start = None, **kwargs):
            network_metrics = kwargs['network_metrics']
            spiking_data_by_unit = network_metrics['spiking_data']['spiking_data_by_unit']
            popData = sim.net.allPops
            E_gids = popData['E']['cellGids']
            I_gids = popData['I']['cellGids']
            
            if ax is None:
                fig, ax_raster = plt.subplots(1, 1, figsize=(16, 4.5))
            else:
                ax_raster = ax
                subplot = True
                
            ax_raster = plot_raster(ax_raster, spiking_data_by_unit, E_gids=E_gids, I_gids=I_gids, data_type='simulated')
            
            # if trim_start, trim first x seconds from the start of the simulation
            if trim_start is not None and trim_start > 0 and trim_start < ax_raster.get_xlim()[1]:
                ax_raster.set_xlim(trim_start, ax_raster.get_xlim()[1])
            elif trim_start is not None and trim_start > 0 and trim_start > ax_raster.get_xlim()[1]:
                modified_trim = ax_raster.get_xlim()[1]*0.1
                ax_raster.set_xlim(modified_trim, ax_raster.get_xlim()[1])
                print('boop')                
            
            #break if subplot
            if subplot:
                #plt.close()
                return ax_raster

            if DEBUG_MODE:
                #save local for debugging
                #dev_dir = '/pscratch/sd/a/adammwea/workspace/RBS_neuronal_network_models/optimization_projects/CDKL5_DIV21/scripts'
                dev_dir = os.path.dirname(os.path.realpath(__file__))
                plt.savefig(os.path.join(dev_dir, '_raster_plot.png'), dpi=300)
                #save local for debugging
                
            #save wherever data is saved
            #sim_data_path = SIMULATION_RUN_PATH
            raster_plot_path = sim_data_path.replace('_data', '_raster_plot')
            #remove file type and replace with png
            if '.json' in raster_plot_path:
                raster_plot_path = raster_plot_path.replace('.json', '.png')
            elif '.pkl' in raster_plot_path:
                raster_plot_path = raster_plot_path.replace('.pkl', '.png')
            #raster_plot_path = raster_plot_path.replace('.json', '.png')
            plt.savefig(raster_plot_path, dpi=300)
            print(f"Raster plot saved to {raster_plot_path}")
            
            #save as pdf
            raster_plot_path = sim_data_path.replace('_data', '_raster_plot')
            #raster_plot_path = raster_plot_path.replace('.json', '.pdf')
            if '.json' in raster_plot_path:
                raster_plot_path = raster_plot_path.replace('.json', '.pdf')
            elif '.pkl' in raster_plot_path:
                raster_plot_path = raster_plot_path.replace('.pkl', '.pdf')
            plt.savefig(raster_plot_path)
            print(f"Raster plot saved to {raster_plot_path}")
            plt.close()
            
            raster_plots_paths = [
                raster_plot_path,
                raster_plot_path.replace('.pdf', '.png'),
            ]
            
            return raster_plots_paths
        raster_plot_paths = plot_simulated_raster_wrapper(trim_start = trim_start, **kwargs)
        privileged_print("\tIndividual raster plots saved.")
        privileged_print(f'\tTime taken: {time()-start_raster} seconds')

        #plot bursting summary
        start_bursting = time()
        def plot_simulated_bursting_wrapper(ax = None, subplot=False, trim_start = None, **kwargs):
            if ax is None:
                fig, new_ax = plt.subplots(1, 1)
                fig.set_size_inches(16, 4.5)
            else:
                new_ax = ax
                subplot = True #if ax is passed in, then we are plotting on a subplot
            
            #
            conv_params = kwargs['conv_params']
            SpikeTimes = kwargs['network_metrics']['spiking_data']['spiking_times_by_unit']
            #from DIV21.utils.fitness_helper import plot_network_activity_aw
            #from RBS_network_models.developing.utils.analysis_helper import plot_network_activity_aw
            from MEA_Analysis.NetworkAnalysis.awNetworkAnalysis.network_analysis import plot_network_activity_aw
            #bursting_ax, _ = plot_network_activity_aw(new_ax, SpikeTimes, **conv_params) #TODO need to make sure this function agrees with mandar. would be best if we shared a function here.
            bursting_ax = kwargs['network_metrics']['bursting_data']['bursting_summary_data']['ax']
            mega_ax = kwargs['network_metrics']['mega_bursting_data']['bursting_summary_data']['ax']
            
            # # HACK
            # mega_conv_params = kwargs['conv_params'].copy()
            # mega_conv_params['binSize'] *= 5
            # mega_conv_params['gaussianSigma'] *= 15
            #mega_ax, _ = plot_network_activity_aw(new_ax, SpikeTimes, **mega_conv_params) #TODO need to make sure this function agrees with mandar. would be best if we shared a function here.
            
            from MEA_Analysis.NetworkAnalysis.awNetworkAnalysis.network_analysis import plot_network_bursting_experimental
            new_ax = plot_network_bursting_experimental(new_ax, bursting_ax, mega_ax=mega_ax)            
            
            # if trim_start, trim first x seconds from the start of the simulation
            if trim_start is not None and trim_start > 0 and trim_start < new_ax.get_xlim()[1]:
                new_ax.set_xlim(trim_start, new_ax.get_xlim()[1])
            elif trim_start is not None and trim_start > 0 and trim_start > new_ax.get_xlim()[1]:
                modified_trim = new_ax.get_xlim()[1]*0.1
                new_ax.set_xlim(modified_trim, new_ax.get_xlim()[1])    
            
            new_ax.set_title('Bursting summary')
            new_ax.set_xlabel('Time (s)')
            new_ax.set_ylabel('Fire rate (Hz)')
            
            #break if subplot
            if subplot:
                #plt.close()
                return new_ax
            
            plt.tight_layout()
            
            if DEBUG_MODE:
                # Save local for debugging
                dev_dir = os.path.dirname(os.path.realpath(__file__))
                plt.savefig(os.path.join(dev_dir, '_bursting_plot.png'), dpi=300)
                # Save local for debugging
            
            # save wherever data is saved
            #sim_data_path = SIMULATION_RUN_PATH
            bursting_plot_path = sim_data_path.replace('_data', '_bursting_plot')
            #remove file type and replace with png
            #bursting_plot_path = bursting_plot_path.replace('.json', '.png')
            if '.json' in bursting_plot_path:
                bursting_plot_path = bursting_plot_path.replace('.json', '.png')
            elif '.pkl' in bursting_plot_path:
                bursting_plot_path = bursting_plot_path.replace('.pkl', '.png')
            plt.savefig(bursting_plot_path, dpi=300)
            print(f"Bursting plot saved to {bursting_plot_path}")
            
            #save as pdf
            bursting_plot_path = sim_data_path.replace('_data', '_bursting_plot')
            #bursting_plot_path = bursting_plot_path.replace('.json', '.pdf')
            if '.json' in bursting_plot_path:
                bursting_plot_path = bursting_plot_path.replace('.json', '.pdf')
            elif '.pkl' in bursting_plot_path:
                bursting_plot_path = bursting_plot_path.replace('.pkl', '.pdf')
            plt.savefig(bursting_plot_path)
            print(f"Bursting plot saved to {bursting_plot_path}")
            plt.close()
            
            bursting_plot_paths = [
                bursting_plot_path,
                bursting_plot_path.replace('.pdf', '.png'),
            ]
            
            return bursting_plot_paths    
        bursting_plot_paths = plot_simulated_bursting_wrapper(trim_start = trim_start, **kwargs)
        privileged_print("\tIndividual bursting plots saved.")
        privileged_print(f'\tTime taken: {time()-start_bursting} seconds')

        # combine plots into a single summary plot
        start_summary = time()
        def plot_simulation_summary(trim_start = None, **kwargs):
            fig, ax = plt.subplots(2, 1, figsize=(16, 9))
            raster_plot_ax, bursting_plot_ax = ax
            subplot = True
            raster_plot_ax = plot_simulated_raster_wrapper(ax=raster_plot_ax, subplot=subplot, trim_start = trim_start, **kwargs)
            bursting_plot_ax = plot_simulated_bursting_wrapper(ax=bursting_plot_ax, subplot=subplot, trim_start = trim_start, **kwargs)
            
            #make both plots share the same x-axis
            raster_plot_ax.get_shared_x_axes().join(raster_plot_ax, bursting_plot_ax)    
            plt.tight_layout()
            
            if DEBUG_MODE:
                # Save local for debugging
                dev_dir = os.path.dirname(os.path.realpath(__file__))
                plt.savefig(os.path.join(dev_dir, '_summary_plot.png'), dpi=300)
                # Save local for debugging
                
            # save wherever data is saved
            #sim_data_path = SIMULATION_RUN_PATH
            summary_plot_path = sim_data_path.replace('_data', '_summary_plot')
            #remove file type and replace with png
            #summary_plot_path = summary_plot_path.replace('.json', '.png')
            if '.json' in summary_plot_path:
                summary_plot_path = summary_plot_path.replace('.json', '.png')
            elif '.pkl' in summary_plot_path:
                summary_plot_path = summary_plot_path.replace('.pkl', '.png')
            plt.savefig(summary_plot_path, dpi=300)
            print(f"Summary plot saved to {summary_plot_path}")
            
            #save as pdf
            summary_plot_path = sim_data_path.replace('_data', '_summary_plot')
            #summary_plot_path = summary_plot_path.replace('.json', '.pdf')
            if '.json' in summary_plot_path:
                summary_plot_path = summary_plot_path.replace('.json', '.pdf')
            elif '.pkl' in summary_plot_path:
                summary_plot_path = summary_plot_path.replace('.pkl', '.pdf')
            plt.savefig(summary_plot_path)  
            print(f"Summary plot saved to {summary_plot_path}")
            plt.close()
            
            summary_plot_paths = [
                summary_plot_path,
                summary_plot_path.replace('.pdf', '.png'),
            ]    
            return summary_plot_paths
        plot_simulation_summary(trim_start=trim_start, **kwargs)
        privileged_print("\tSimulation summary plot saved.")
        privileged_print(f'\tTime taken: {time()-start_summary} seconds')

        # comparison plot against reference data snippet
        start_comparison = time()
        activate_print()
        def plot_comparision_plot(ax_list=None, subplot=False, trim_start = None, **kwargs):
            #activate_print() #for debugging
            #fig, ax = plt.subplots(4, 1, figsize=(16, 9))
            if ax_list is None:
                fig, ax_list = plt.subplots(4, 1, figsize=(16, 9))
                sim_raster_ax, sim_bursting_ax, ref_raster_ax, ref_bursting_ax = ax_list
            else:
                #assert that there be 4 axes in the ax list
                assert len(ax_list) == 4, "There must be 4 axes in the ax_list."
                sim_raster_ax, sim_bursting_ax, ref_raster_ax, ref_bursting_ax = ax_list
                subplot = True
            
            #plot simulated raster
            sim_raster_ax = plot_simulated_raster_wrapper(ax=sim_raster_ax, subplot=True, trim_start=trim_start, **kwargs)
            
            #plot simulated bursting
            sim_bursting_ax = plot_simulated_bursting_wrapper(ax=sim_bursting_ax, subplot=True, trim_start=trim_start, **kwargs)
            
            # plot reference raster
            def plot_reference_raster_wrapper(ax=None, subplot=False, sim_data_length=None, **kwargs):
                #load npy ref data
                print(ax)
                print(subplot)
                #print(kwargs)
                
                #load npy ref data
                ref_data = np.load(
                    REFERENCE_DATA_NPY, 
                    allow_pickle=True
                    ).item()
                network_metrics = ref_data
                spiking_data_by_unit = network_metrics['spiking_data']['spiking_data_by_unit'].copy()
                unit_ids = [unit_id for unit_id in spiking_data_by_unit.keys()]
                
                # Trim reference data to match the length of simulated data
                if sim_data_length is not None:
                    for unit_id in unit_ids:
                        #spiking_data_by_unit[unit_id] = spiking_data_by_unit[unit_id][:sim_data_length]
                        #max_spike_time_index = np.argmax(spiking_data_by_unit[unit_id])
                        spike_times = spiking_data_by_unit[unit_id]['spike_times']
                        spike_times = spike_times[spike_times < sim_data_length]
                        spiking_data_by_unit[unit_id]['spike_times'] = spike_times
                
                if ax is None:
                    fig, ax_raster = plt.subplots(1, 1, figsize=(16, 4.5))
                else:
                    ax_raster = ax
                    subplot = True
                    
                ax_raster = plot_raster(ax_raster, spiking_data_by_unit, unit_ids=unit_ids, data_type='experimental')
                
                if subplot:
                    #plt.close()
                    print("subplot")
                    return ax_raster
                
                if DEBUG_MODE:
                    #save local for debugging
                    #dev_dir = '/pscratch/sd/a/adammwea/workspace/RBS_neuronal_network_models/optimization_projects/CDKL5_DIV21/scripts'
                    dev_dir = os.path.dirname(os.path.realpath(__file__))
                    plt.savefig(os.path.join(dev_dir, '_ref_raster_plot.png'), dpi=300)
                    #save as pdf
                    plt.savefig(os.path.join(dev_dir, '_ref_raster_plot.pdf'))
                    #save local for debugging
                
                
                #TODO: semi unfinished compared to similar function for simulated raster
            #sim_data_length = len(kwargs['network_metrics']['spiking_data']['spiking_data_by_unit'][0])
            sim_data_length = kwargs['simConfig']['duration']/1000 #in seconds
            ref_raster_ax = plot_reference_raster_wrapper(ax=ref_raster_ax, subplot=True, sim_data_length = sim_data_length, **kwargs)
            
            #plot reference bursting
            def plot_reference_bursting_wrapper(ax=None, subplot=False, sim_data_length=None, trim_start = None, **kwargs):
                #ref_data = np.load(REFERENCE_DATA_NPY, allow_pickle=True).item()
                ref_data = np.load(REFERENCE_DATA_NPY, allow_pickle=True).item()
                network_metrics = ref_data
                conv_params = kwargs['conv_params']
                SpikeTimes = network_metrics['spiking_data']['spiking_times_by_unit']
                if ax is None:
                    fig, new_ax = plt.subplots(1, 1)
                    fig.set_size_inches(16, 4.5)
                else:
                    new_ax = ax
                    subplot = True
                    
                if sim_data_length is not None:
                    for unit_id in SpikeTimes.keys():
                        spike_times = SpikeTimes[unit_id]
                        spike_times = spike_times[spike_times < sim_data_length]
                        SpikeTimes[unit_id] = spike_times
                
                #from DIV21.utils.fitness_helper import plot_network_activity_aw
                #from RBS_network_models.developing.utils.analysis_helper import plot_network_activity_aw
                from MEA_Analysis.NetworkAnalysis.awNetworkAnalysis.network_analysis import plot_network_activity_aw
                new_ax, _ = plot_network_activity_aw(new_ax, SpikeTimes, **conv_params) #TODO need to make sure this function agrees with mandar. would be best if we shared a function here.
                
                new_ax.set_title('Bursting summary')
                new_ax.set_xlabel('Time (s)')
                new_ax.set_ylabel('Fire rate (Hz)')
                
                #break if subplot
                if subplot:
                    #plt.close()
                    return new_ax
                
                plt.tight_layout()
                
                if DEBUG_MODE:
                    # Save local for debugging
                    dev_dir = os.path.dirname(os.path.realpath(__file__))
                    plt.savefig(os.path.join(dev_dir, '_ref_bursting_plot.png'), dpi=300)
                    # plot as pdf
                    plt.savefig(os.path.join(dev_dir, '_ref_bursting_plot.pdf'))
                    # Save local for debugging
                
                # TODO: semi unfinished compared to similar function for simulated bursting
            sim_data_length = kwargs['simConfig']['duration']/1000 #in seconds
            ref_bursting_ax = plot_reference_bursting_wrapper(ax=ref_bursting_ax, subplot=True, sim_data_length=sim_data_length, trim_start = trim_start, **kwargs)
            
            # # remove the first second of x-axis for all plots
            # def set_xlim_to_first_value_over_one(ax):
            #     x_data = ax.lines[0].get_xdata()
            #     first_value_over_one = next((x for x in x_data if x > 1), None)
            #     if first_value_over_one is not None:
            #         ax.set_xlim(first_value_over_one, ax.get_xlim()[1])
            
            # set_xlim_to_first_value_over_one(sim_raster_ax)
            # set_xlim_to_first_value_over_one(ref_raster_ax)
            # set_xlim_to_first_value_over_one(sim_bursting_ax)
            # set_xlim_to_first_value_over_one(ref_bursting_ax)
            
            # #print axes limits
            # print(sim_raster_ax.get_xlim())
            # print(ref_raster_ax.get_xlim())
            # print(sim_bursting_ax.get_xlim())
            # print(ref_bursting_ax.get_xlim())
            #ref_bursting_ax.append(ref_bursting_ax.get_xlim()[0] + 1)
            
            # # Ensure y-axis is the same for bursting plots
            # sim_bursting_ylim = sim_bursting_ax.get_ylim()
            # ref_bursting_ylim = ref_bursting_ax.get_ylim()
            
            # ensure xaxis of refernce plots matches simulated raster
            ref_raster_ax.set_xlim(sim_raster_ax.get_xlim())
            ref_bursting_ax.set_xlim(sim_raster_ax.get_xlim())
            sim_bursting_ax.set_xlim(ref_bursting_ax.get_xlim())
            
            # Ensure y-axis is the same for bursting plots
            def adjust_ylim_based_on_xlim(ax):
                x_data = ax.lines[0].get_xdata()
                y_data = ax.lines[0].get_ydata()
                xlim = ax.get_xlim()
                filtered_y_data = [y for x, y in zip(x_data, y_data) if xlim[0] <= x <= xlim[1]]
                if filtered_y_data:
                    ax.set_ylim(min(filtered_y_data) * 0.8, max(filtered_y_data) * 1.2)

            adjust_ylim_based_on_xlim(sim_bursting_ax)
            adjust_ylim_based_on_xlim(ref_bursting_ax)  
            
            
            # Calculate the combined y-axis limits
            sim_bursting_ylim = sim_bursting_ax.get_ylim()
            ref_bursting_ylim = ref_bursting_ax.get_ylim()
            # Calculate the combined y-axis limits
            combined_ylim = [
                min(sim_bursting_ylim[0], ref_bursting_ylim[0]) * 0.8,
                max(sim_bursting_ylim[1], ref_bursting_ylim[1]) * 1.2
            ]
            
            # Set the same y-axis limits for both axes
            sim_bursting_ax.set_ylim(combined_ylim)
            ref_bursting_ax.set_ylim(combined_ylim)
            
            # Make all four plots share the same x-axis as simulated raster
            sim_raster_ax.get_shared_x_axes().join(sim_raster_ax, sim_bursting_ax)
            sim_raster_ax.get_shared_x_axes().join(sim_raster_ax, ref_raster_ax)
            sim_raster_ax.get_shared_x_axes().join(sim_raster_ax, ref_bursting_ax)  

            #
            if subplot:
                #plt.tight_layout()
                #plt.close()
                return [sim_raster_ax, sim_bursting_ax, ref_raster_ax, ref_bursting_ax] 
            
            plt.tight_layout()
            
            if DEBUG_MODE:
                # Save local for debugging
                dev_dir = os.path.dirname(os.path.realpath(__file__))
                fig.savefig(os.path.join(dev_dir, '_comparison_plot.png'), dpi=300)
                # Save local for debugging
            
            # save wherever data is saved
            #sim_data_path = SIMULATION_RUN_PATH
            comparison_plot_path = sim_data_path.replace('_data', '_comparison_plot')
            #remove file type and replace with png
            #comparison_plot_path = comparison_plot_path.replace('.json', '.png')
            if '.json' in comparison_plot_path:
                comparison_plot_path = comparison_plot_path.replace('.json', '.png')
            elif '.pkl' in comparison_plot_path:
                comparison_plot_path = comparison_plot_path.replace('.pkl', '.png')
            fig.savefig(comparison_plot_path, dpi=300)
            print(f"Comparison plot saved to {comparison_plot_path}")
            
            #save as pdf
            comparison_plot_path = sim_data_path.replace('_data', '_comparison_plot')
            #comparison_plot_path = comparison_plot_path.replace('.json', '.pdf')
            if '.json' in comparison_plot_path:
                comparison_plot_path = comparison_plot_path.replace('.json', '.pdf')
            elif '.pkl' in comparison_plot_path:
                comparison_plot_path = comparison_plot_path.replace('.pkl', '.pdf')
            fig.savefig(comparison_plot_path)
            print(f"Comparison plot saved to {comparison_plot_path}")
            plt.close()
        plot_comparision_plot(trim_start=trim_start, **kwargs)
        privileged_print("\tComparison plot saved.")
        privileged_print(f'\tTime taken: {time()-start_comparison} seconds')

        # build comparison summary slide
        start_summary_slide = time()
        def build_comparision_summary_slide(sim_data_path, trim_start = None,):
            fig, ax_list = plt.subplots(4, 1, figsize=(16, 9))
            
            #plot comparison plot
            plot_comparision_plot(ax_list=ax_list, subplot=True, trim_start=trim_start, **kwargs)
            
            #remove x-axis labels from all plots
            for ax in ax_list:
                ax.set_xlabel('')
                
            #remove titles from all plots
            for ax in ax_list:
                ax.set_title('')
                
            #remove x-axis ticks from all plots except the bottom one
            for ax in ax_list[:-1]:
                ax.set_xticks([])
                
            # add 'simulated' and 'reference' labels above each raster plot
            ax_list[0].text(0.5, 1.1, 'Simulated', ha='center', va='center', transform=ax_list[0].transAxes, fontsize=12)
            ax_list[2].text(0.5, 1.1, 'Reference', ha='center', va='center', transform=ax_list[2].transAxes, fontsize=12)
            
            # add time(s) label to the bottom plot
            ax_list[-1].set_xlabel('Time (s)')
                    
            #create space to the right of plots for text
            fig.subplots_adjust(right=0.6)
            
            #create some space at the bottom for one line of text
            fig.subplots_adjust(bottom=0.1)
            
            #add text
            spiking_summary = kwargs['network_metrics']['spiking_data']['spiking_summary_data']
            bursting_summary = kwargs['network_metrics']['bursting_data']['bursting_summary_data']
            simulated_spiking_data = kwargs['network_metrics']['simulated_data']
            #simulated_bursting_data = kwargs['network_metrics']['simulated_data']
            
            #TODO: if meanBurstrate is not in the spiking summary, then the following will fail, add nan in that case for now
            if 'mean_Burst_Rate' not in bursting_summary:
                bursting_summary['mean_Burst_Rate'] = np.nan                
            
            text = (
                # f"Summary of comparison between simulated and reference data:\n"
                # f"Simulated raster plot is shown in the top left, reference raster plot is shown in the bottom left.\n"
                # f"Simulated bursting plot is shown in the top right, reference bursting plot is shown in the bottom right.\n"
                f"Simulated spiking summary:\n"
                f"  - Mean firing rate: {spiking_summary['MeanFireRate']} Hz\n"
                f"  - Coefficient of variation: {spiking_summary['CoVFireRate']}\n"
                f"  - Mean ISI: {spiking_summary['MeanISI']} s\n"
                f"  - Coefficient of variation of ISI: {spiking_summary['CoV_ISI']}\n"
                f"  - Mean firing rate E: {simulated_spiking_data['MeanFireRate_E']} Hz\n"
                f"  - Mean CoV FR E: {simulated_spiking_data['CoVFireRate_E']}\n"
                f"  - Mean firing rate I: {simulated_spiking_data['MeanFireRate_I']} Hz\n"
                f"  - Mean CoV FR I: {simulated_spiking_data['CoVFireRate_I']}\n"
                f"  - Mean ISI E: {simulated_spiking_data['MeanISI_E']} s\n"
                f"  - CoV ISI E: {simulated_spiking_data['CoV_ISI_E']}\n"
                f"  - Mean ISI I: {simulated_spiking_data['MeanISI_I']} s\n"
                f"  - CoV ISI I: {simulated_spiking_data['CoV_ISI_I']}\n"
                
                f"Simulated bursting summary:\n"        
                #mean burst rate
                #number of units
                f"  - Number of units: {bursting_summary['NumUnits']}\n"
                #baseline
                f"  - Baseline: {bursting_summary['baseline']} Hz\n"
                #fanofactor
                f"  - Fanofactor: {bursting_summary['fano_factor']}\n"
                #num bursts
                f"  - Number of bursts: {bursting_summary['Number_Bursts']}\n"
                #mean burst rate
                f"  - Mean burst rate: {bursting_summary['mean_Burst_Rate']} bursts/second\n"     
                #mean burst peak
                f"  - Mean burst peak: {bursting_summary['mean_Burst_Peak']} Hz\n"
                #burst peak CoV
                f"  - CoV burst peak: {bursting_summary['cov_Burst_Peak']}\n"
                #mean IBI
                f"  - Mean IBI: {bursting_summary['mean_IBI']} s\n"
                #CoV IBI
                f"  - CoV IBI: {bursting_summary['cov_IBI']}\n"
                #mean within burst isi
                f"  - Mean within burst ISI: {bursting_summary['MeanWithinBurstISI']} s\n"
                #CoV within burst isi
                f"  - CoV within burst ISI: {bursting_summary['CoVWithinBurstISI']}\n"
                #mean outside burst isi
                f"  - Mean outside burst ISI: {bursting_summary['MeanOutsideBurstISI']} s\n"
                #CoV outside burst isi
                f"  - CoV outside burst ISI: {bursting_summary['CoVOutsideBurstISI']}\n"
                #mean whole network isi
                f"  - Mean network ISI: {bursting_summary['MeanNetworkISI']} s\n"
                #CoV whole network isi
                f"  - CoV network ISI: {bursting_summary['CoVNetworkISI']}\n"
                )
            
            #append reference metrics to text
            ref_data = np.load(REFERENCE_DATA_NPY, allow_pickle=True).item()
            ref_spiking_summary = ref_data['spiking_data']['spiking_summary_data']
            ref_bursting_summary = ref_data['bursting_data']['bursting_summary_data']    
            text += (
                f"\n"
                f"\n"
                f"Reference spiking summary:\n"
                f"  - Mean firing rate: {ref_spiking_summary['MeanFireRate']} Hz\n"
                f"  - Coefficient of variation: {ref_spiking_summary['CoVFireRate']}\n"
                f"  - Mean ISI: {ref_spiking_summary['MeanISI']} s\n"
                f"  - Coefficient of variation of ISI: {ref_spiking_summary['CoV_ISI']}\n"
                
                f"Reference bursting summary:\n"
                #number of units
                f"  - Number of units: {ref_bursting_summary['NumUnits']}\n"
                #baseline
                f"  - Baseline: {ref_bursting_summary['baseline']} Hz\n"
                #fanofactor
                f"  - Fanofactor: {ref_bursting_summary['fano_factor']}\n"
                #num bursts
                f"  - Number of bursts: {ref_bursting_summary['Number_Bursts']}\n"
                #mean burst rate
                #f"  - Mean burst rate: {ref_bursting_summary['mean_Burst_Rate']} Hz\n" #TODO: this is not in the data        
                #mean burst peak
                f"  - Mean burst peak: {ref_bursting_summary['mean_Burst_Peak']} Hz\n"
                #burst peak CoV
                f"  - CoV burst peak: {ref_bursting_summary['cov_Burst_Peak']}\n"
                #mean IBI
                f"  - Mean IBI: {ref_bursting_summary['mean_IBI']} s\n"
                #CoV IBI
                f"  - CoV IBI: {ref_bursting_summary['cov_IBI']}\n"
                #mean within burst isi
                f"  - Mean within burst ISI: {ref_bursting_summary['MeanWithinBurstISI']} s\n"
                #CoV within burst isi
                f"  - CoV within burst ISI: {ref_bursting_summary['CoVWithinBurstISI']}\n"
                #mean outside burst isi
                f"  - Mean outside burst ISI: {ref_bursting_summary['MeanOutsideBurstISI']} s\n"
                #CoV outside burst isi
                f"  - CoV outside burst ISI: {ref_bursting_summary['CoVOutsideBurstISI']}\n"
                #mean whole network isi
                f"  - Mean network ISI: {ref_bursting_summary['MeanNetworkISI']} s\n"
                #CoV whole network isi
                f"  - CoV network ISI: {ref_bursting_summary['CoVNetworkISI']}\n"
                )
            
            #add average fitness to text
            text += (
                f"\n"
                f"\n"
                f"\nAverage fitness: {average_fitness}"
                ) 
            
            #add text to the right of the plots
            fig.text(0.65, 0.5, text, ha='left', va='center', fontsize=9)
            
            #add data path at the bottom of the slide
            #fig.text(0.5, 0.05, f"Data path: {SIMULATION_RUN_PATH}", ha='center', va='center', fontsize=10)
            #go a little lower and add data path
            fig.text(0.5, 0.02, f"Data path: {sim_data_path}", ha='center', va='center', fontsize=9) 
            
            if DEBUG_MODE:
                #save local for debugging
                dev_dir = os.path.dirname(os.path.realpath(__file__))
                fig.savefig(os.path.join(dev_dir, '_comparison_summary_slide.png'), dpi=300)
                #save as pdf
            
            # save wherever data is saved
            sim_data_path = sim_data_path
            comparison_summary_slide_path = sim_data_path.replace('_data', '_comparison_summary_slide')
            #remove file type and replace with png
            #comparison_summary_slide_path = comparison_summary_slide_path.replace('.json', '.png')
            if '.json' in comparison_summary_slide_path:
                comparison_summary_slide_path = comparison_summary_slide_path.replace('.json', '.png')
            elif '.pkl' in comparison_summary_slide_path:
                comparison_summary_slide_path = comparison_summary_slide_path.replace('.pkl', '.png')
            fig.savefig(comparison_summary_slide_path, dpi=300)
            privileged_print(f"Comparison summary slide saved to {comparison_summary_slide_path}")
            
            #save as pdf
            comparison_summary_slide_path = sim_data_path.replace('_data', '_comparison_summary_slide')
            #comparison_summary_slide_path = comparison_summary_slide_path.replace('.json', '.pdf')
            if '.json' in comparison_summary_slide_path:
                comparison_summary_slide_path = comparison_summary_slide_path.replace('.json', '.pdf')
            elif '.pkl' in comparison_summary_slide_path:
                comparison_summary_slide_path = comparison_summary_slide_path.replace('.pkl', '.pdf')
            fig.savefig(comparison_summary_slide_path)
            privileged_print(f"Comparison summary slide saved to {comparison_summary_slide_path}")
            plt.close()
        build_comparision_summary_slide(sim_data_path, trim_start = trim_start)
        privileged_print("\tComparison summary slide saved.")
        privileged_print(f'\tTime taken: {time()-start_summary_slide} seconds')

def plot_raster_plot_experimental(ax, spiking_data_by_unit, unit_types=None, plot_class=False):
    """Plot a raster plot for spiking data."""
    
    # Calculate the average firing rate for each unit
    firing_rates = {}
    for gid in spiking_data_by_unit:
        spike_times = spiking_data_by_unit[gid]['spike_times']
        spike_times = [spike_times] if isinstance(spike_times, (int, float)) else spike_times
        firing_rate = len(spike_times) / (max(spike_times) - min(spike_times)) if len(spike_times) > 1 else 0
        firing_rates[gid] = firing_rate
    
    # Sort the units based on their average firing rates
    sorted_units = sorted(firing_rates, key=firing_rates.get)
    
    # Create a mapping from original gid to new y-axis position
    gid_to_ypos = {gid: pos for pos, gid in enumerate(sorted_units)}
    #gid_to_ypos = {gid: pos for pos, gid in sorted_units.items()}
    
    # Plot the units in the sorted order
    #print(f'PLOT CLASS MODE: {plot_class}')
    if not plot_class: # typical raster plot
        for gid in sorted_units:
            spike_times = spiking_data_by_unit[gid]['spike_times']
            spike_times = [spike_times] if isinstance(spike_times, (int, float)) else spike_times
            ax.plot(spike_times, [gid_to_ypos[gid]] * len(spike_times), 'b.', markersize=2)
    else: # plot class
        #init legend for e/i units
        # e - blue, i - red
        ax.plot([], [], 'b.', markersize=2, label='E')
        ax.plot([], [], 'r.', markersize=2, label='I')
        ax.plot([], [], 'g.', markersize=2, label='U')
        ax.legend()
        #print(unit_types)
        for gid in sorted_units:
            spike_times = spiking_data_by_unit[gid]['spike_times']
            spike_times = [spike_times] if isinstance(spike_times, (int, float)) else spike_times
            
            # if gid not in unit_types, plot as green - unknown, avoid error
            if unit_types is None or gid not in unit_types:
                ax.plot(spike_times, [gid_to_ypos[gid]] * len(spike_times), 'g.', markersize=2)
                continue
            
            if unit_types[gid] == 'E':
                ax.plot(spike_times, [gid_to_ypos[gid]] * len(spike_times), 'b.', markersize=2)
            elif unit_types[gid] == 'I':
                ax.plot(spike_times, [gid_to_ypos[gid]] * len(spike_times), 'r.', markersize=2)
            #ax.plot(spike_times, [gid_to_ypos[gid]] * len(spike_times), 'b.', markersize=2)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Unit ID (sorted by firing rate)')
    ax.set_title('Raster Plot')
    plt.tight_layout()
    return ax

def plot_network_bursting_experimental(ax, bursting_ax, mega_ax=None, mode='normal'):
    """Plot network bursting activity with clear differentiation between lines."""
    
    if mode == 'normal':
        # Copy ax features to the new ax
        ax.set_xlim(bursting_ax.get_xlim())
        ax.set_ylim(bursting_ax.get_ylim())
        ax.set_ylabel('Firing Rate (Hz)')
        ax.set_xlabel('Time (s)')
        ax.set_title(bursting_ax.get_title())

        # Plot Line 1 (blue)
        ax.plot(
            bursting_ax.get_lines()[0].get_xdata(),
            bursting_ax.get_lines()[0].get_ydata(),
            color='blue',
            #label='Bursting'
        )
        
        # plot bursting peaks
        ax.plot(
            bursting_ax.get_lines()[1].get_xdata(),
            bursting_ax.get_lines()[1].get_ydata(),
            'o', color='red', 
            label='Bursts'
        )

        # # Mark peaks with vertical dotted lines
        # for x_peak in bursting_ax.get_lines()[1].get_xdata():
        #     ax.axvline(x=x_peak, linestyle='dotted', color='red', alpha=0.8)

        # If mega_ax is not None, plot filled area and peaks
        if mega_ax is not None:
            # Plot the filled area behind other lines
            ax.fill_between(
                mega_ax.get_lines()[0].get_xdata(),
                mega_ax.get_lines()[0].get_ydata(),
                color='red',
                alpha=0.2,  # Transparency for the filled area
                label='Mega Bursts'
            )
            
            #outline the filled aread, dotted line
            ax.plot(
                mega_ax.get_lines()[0].get_xdata(),
                mega_ax.get_lines()[0].get_ydata(),
                color='red',
                linestyle='dotted'
            )

            # Mark mega peaks with vertical dotted lines
            for x_peak in mega_ax.get_lines()[1].get_xdata():
                ax.axvline(x=x_peak, linestyle='dotted', 
                        color='red', 
                        #label='Mega Bursts'
                        #alpha=0.8
                        )
                
            # remove redundant legend items
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys())

        # Add a legend for clarity
        ax.legend()

        return ax
    elif mode == 'mega':
        # this is for the case where we're plotting mega bursts only
        # Copy ax features to the new ax
        ax.set_xlim(bursting_ax.get_xlim())
        ax.set_ylim(bursting_ax.get_ylim())
        ax.set_ylabel('Firing Rate (Hz)')
        ax.set_xlabel('Time (s)')
        ax.set_title(bursting_ax.get_title())
        
        # Plot Line 1 (red)
        ax.plot(
            bursting_ax.get_lines()[0].get_xdata(),
            bursting_ax.get_lines()[0].get_ydata(),
            color='red',
            #label='Mega Bursts'
        )
        
        # plot bursting peaks
        ax.plot(
            bursting_ax.get_lines()[1].get_xdata(),
            bursting_ax.get_lines()[1].get_ydata(),
            'o', 
            color='orange',
            label='Mega Bursts'
        )
        
        # Mark peaks with vertical dotted lines
        for x_peak in bursting_ax.get_lines()[1].get_xdata():
            ax.axvline(x=x_peak, linestyle='dotted', color='red', alpha=0.8)
            
        # Add a legend for clarity
        ax.legend()
        
        return ax

'''Experimental Data Functions'''
def analyze_bursting_activity(spike_times, spike_times_by_unit, isi_threshold):
    '''Slightly modified version of the code from mea_analysis_pipeline.py'''
    #spike_times = {int(i): rasterData['spkt'][rasterData['spkid'] == i] for i in np.unique(rasterData['spkid'])}
    
    # convert spike_times to dictionary
    #spike_times_dict = {int(i): spike_times[spike_times['spkid'] == i] for i in np.unique(spike_times['spkid'])}
    #spike_times = spike_times_dict
    spike_times = spike_times_by_unit
    
    
    #burst_statistics = helper.detect_bursts_statistics(spike_times, isi_threshold=isi_threshold)
    burst_statistics = detect_bursts_statistics(spike_times, isi_threshold=isi_threshold)
    bursts_by_unit = [unit_stats['bursts'] for unit_stats in burst_statistics.values()]
    
    all_isis_within_bursts = np.concatenate([stats['isis_within_bursts'] for stats in burst_statistics.values() if stats['isis_within_bursts'].size > 0])
    all_isis_outside_bursts = np.concatenate([stats['isis_outside_bursts'] for stats in burst_statistics.values() if stats['isis_outside_bursts'].size > 0])
    all_isis = np.concatenate([stats['isis_all'] for stats in burst_statistics.values() if stats['isis_all'].size > 0])

    mean_isi_within_combined = np.mean(all_isis_within_bursts) if all_isis_within_bursts.size > 0 else np.nan
    cov_isi_within_combined = np.cov(all_isis_within_bursts) if all_isis_within_bursts.size > 0 else np.nan
    #fano_factor_within_combined = np.concatenate([[stats['fano_factor_within']] for stats in burst_statistics.values() if stats['fano_factor_within'] is not None])

    mean_isi_outside_combined = np.mean(all_isis_outside_bursts) if all_isis_outside_bursts.size > 0 else np.nan
    cov_isi_outside_combined = np.cov(all_isis_outside_bursts) if all_isis_outside_bursts.size > 0 else np.nan
    #fano_factor_outside_combined = np.concatenate([[stats['fano_factor_outside']] for stats in burst_statistics.values() if stats['fano_factor_outside'] is not None])

    mean_isi_all_combined = np.mean(all_isis) if all_isis.size > 0 else np.nan
    cov_isi_all_combined = np.cov(all_isis) if all_isis.size > 0 else np.nan    
    #fano_factor_all_combined = np.concatenate([[stats['fano_factor_all']] for stats in burst_statistics.values() if stats['fano_factor'] is not None])

    bursting_summary_data = {}
    bursting_summary_data['MeanWithinBurstISI'] = mean_isi_within_combined
    bursting_summary_data['CoVWithinBurstISI'] = cov_isi_within_combined
    bursting_summary_data['MeanOutsideBurstISI'] = mean_isi_outside_combined   
    bursting_summary_data['CoVOutsideBurstISI'] = cov_isi_outside_combined
    bursting_summary_data['MeanNetworkISI'] = mean_isi_all_combined
    bursting_summary_data['CoVNetworkISI'] = cov_isi_all_combined
    #bursting_summary_data['fano_factor_within'] = fano_factor_within_combined
    #bursting_summary_data['fano_factor_outside'] = fano_factor_outside_combined
    #bursting_summary_data['fano_factor_all'] = fano_factor_all_combined
    bursting_summary_data['NumUnits'] = len(spike_times)
    bursting_summary_data['Number_Bursts'] = sum(len(unit_stats['bursts']) for unit_stats in burst_statistics.values())
    bursting_summary_data['mean_IBI'] = np.mean(all_isis) if all_isis.size > 0 else np.nan
    bursting_summary_data['cov_IBI'] = np.cov(all_isis) if all_isis.size > 0 else np.nan
    
    bursting_data = {
        'bursting_summary_data': bursting_summary_data,
        'bursting_data_by_unit': burst_statistics,
    }

    return bursting_data

def analyze_convolved_spiking_signal(spike_times, spike_times_by_unit, min_peak_distance, 
                                     binSize, gaussianSigma, 
                                     thresholdBurst, prominence, title = 'Network Activity'):
    fig, ax = plt.subplots()
    #spike_times = {i: rasterData['spkt'][rasterData['spkid'] == i] for i in np.unique(rasterData['spkid'])} #plot_network_activity expects spike_times as dictionary
    #spike_times_dict = {int(i): spike_times[spike_times['spkid'] == i] for i in np.unique(spike_times['spkid'])}
    try:
        #ax, convolved_signal_metrics = helper.plot_network_activity(ax, spike_times_by_unit, min_peak_distance=min_peak_distance, binSize=binSize, gaussianSigma=gaussianSigma, thresholdBurst=thresholdBurst)
        ax, convolved_signal_metrics = plot_network_activity_aw(ax, spike_times_by_unit, min_peak_distance=min_peak_distance, 
                                                                binSize=binSize, gaussianSigma=gaussianSigma, 
                                                                thresholdBurst=thresholdBurst,
                                                                prominence=prominence,
                                                                title=title,
                                                                )
        
        #pull y data out of ax
        y_data = ax.lines[0].get_ydata()
        signal = y_data
        
        #baseline = np.mean(signal)
        # convolved_data = {
        #     'mean_Burst_Rate': convolved_signal_metrics['mean_Burst_Rate'],
        #     'mean_Burst_Peak': convolved_signal_metrics['mean_Burst_Peak'],
        #     'cov_Burst_Peak': convolved_signal_metrics['cov_Burst_Peak'],
        #     'fano_factor': convolved_signal_metrics['fano_factor'],
        #     'mean_IBI': convolved_signal_metrics['mean_IBI'],
        #     'cov_IBI': convolved_signal_metrics['cov_IBI'],
        #     'Number_Bursts': convolved_signal_metrics['Number_Bursts'],
        #     #'baseline': convolve_signal_get_baseline(spike_times, binSize=binSize, gaussianSigma=gaussianSigma),
        #     'baseline': baseline,
        #     #'fig': fig,
        #     'ax': ax,
        # }
        # return convolved_data
        
        # aw 2025-02-10 11:37:18 - consolidating analysis to analyze_burst_activity_v2 - so just returning signal here.
        return ax, signal       
    except Exception as e:
        # #set all convolved data to nan
        # print(f'Error analyzing convolved spiking signal: {e}')
        # print(f'might be cause by all spike times being in the same bin')
        # convolved_data = {
        #     'mean_Burst_Rate': np.nan,
        #     'mean_Burst_Peak': np.nan,
        #     'cov_Burst_Peak': np.nan,
        #     'fano_factor': np.nan,
        #     'mean_IBI': np.nan,
        #     'cov_IBI': np.nan,
        #     'Number_Bursts': np.nan,
        #     'baseline': np.nan,
        # }
        # return convolved_data
        print(f'Error analyzing convolved spiking signal: {e}')
        return ax, None

def init_burst_analysis_params(isi_threshold=None):
    # try:
    #     #from simulate._config_files.burst_analysis_params import burst_analysis_params
    #     burst_analysis_params_temp = burst_analysis_params
    #     burst_analysis_params_temp['isi_threshold'] = isi_threshold if isi_threshold is not None else burst_analysis_params['isi_threshold']
    # except:
    #     burst_analysis_params_temp = {
    #         'isi_threshold': 0.1 if isi_threshold is None else isi_threshold,
    #     }
    
    burst_analysis_params_temp = {
        'isi_threshold': 0.1 if isi_threshold is None else isi_threshold,
    }    
    burst_analysis_params = burst_analysis_params_temp
    return burst_analysis_params

def extract_bursting_activity_data(spike_times, spike_times_by_unit, debug_mode = False, **kwargs):

    #conv_params = init_convolution_params()
    conv_params = kwargs.get('conv_params', None)
    try: conv_params = conv_params.conv_params if conv_params is not None else None
    except: conv_params = conv_params
    assert conv_params is not None, 'No convolution parameters provided'
    assert 'binSize' in conv_params, 'Convolution parameters must include "binSize"'
    assert 'gaussianSigma' in conv_params, 'Convolution parameters must include "gaussianSigma"'
    
    mega_conv_params = kwargs.get('mega_params', None)
    try: mega_conv_params = mega_conv_params.mega_conv_params if mega_conv_params is not None else None
    except: mega_conv_params = mega_conv_params
    assert mega_conv_params is not None, 'No mega convolution parameters provided'
    assert 'binSize' in mega_conv_params, 'Mega convolution parameters must include "binSize"'
    assert 'gaussianSigma' in mega_conv_params, 'Mega convolution parameters must include "gaussianSigma"'
    
    try:
        bursting_data = analyze_bursting_activity_v4(
                                spike_times, spike_times_by_unit,
                                binSize = conv_params['binSize'],
                                gaussianSigma = conv_params['gaussianSigma'],
                                thresholdBurst = conv_params['thresholdBurst'],
                                min_peak_distance = conv_params['min_peak_distance'],
                                prominence = conv_params['prominence'],
                                #title = 'Network Activity'
                                debug_mode = debug_mode,
                                **kwargs
                                )
    except Exception as e:
        print(f'Error analyzing bursting activity: {e}')
        bursting_data = None
    
    try:
        mega_bursting_data = analyze_bursting_activity_v4(
                                spike_times, spike_times_by_unit,
                                binSize = mega_conv_params['binSize'],
                                gaussianSigma = mega_conv_params['gaussianSigma'],
                                thresholdBurst = mega_conv_params['thresholdBurst'],
                                min_peak_distance = mega_conv_params['min_peak_distance'],
                                prominence = mega_conv_params['prominence'],
                                #title = 'Network Activity'
                                debug_mode = debug_mode,
                                **kwargs
                                )
    except Exception as e:
        print(f'Error analyzing mega bursting activity: {e}')
        mega_bursting_data = None
    
    # update network data
    network_data.update({
        'bursting_data': bursting_data,
        'mega_bursting_data': mega_bursting_data,
    })
    
    return network_data

def get_mean_fr(recording_object, sorting_object, sampling_rate=10000):
    total_fr = 0
    unit_count = 0
    units = sorting_object.get_unit_ids()
    for unit in units:
        spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
        assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
        if len(spike_train) > 1:
            duration = recording_object.get_total_duration()  # seconds
            fr = len(spike_train) / duration  # spikes per second
            if not np.isnan(fr) and not np.isinf(fr):
                total_fr += fr
                unit_count += 1
    return total_fr / unit_count if unit_count > 0 else None

def get_CoV_fr_experimental(recording_object, sorting_object, sampling_rate=10000, window_size=1.0):
    """
    Calculate the Coefficient of Variation of Firing Rate (CoV FR) for experimental data over time windows.

    Parameters:
    - recording_object: A recording extractor object containing the experiment's total duration.
    - sorting_object: A sorting extractor object with spike times for each unit.
    - sampling_rate: Sampling rate of the recording (in Hz, default: 10000).
    - window_size: Size of the time window for calculating firing rates (default: 1.0 second).

    Returns:
    - CoV_fr: Coefficient of Variation of Firing Rate (float).
    """
    # Get total duration of the recording (in seconds)
    total_duration = recording_object.get_total_duration()

    # Initialize list to hold firing rates across all windows for all units
    firing_rates = []

    # Get unit IDs from the sorting object
    units = sorting_object.get_unit_ids()

    # Process each unit's spike train
    for unit in units:
        spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate  # Convert spike times to seconds
        assert all(spike_time >= 0 for spike_time in spike_train), f"Spike times for unit {unit} contain negative values"

        if len(spike_train) > 1:
            # Divide the total duration into non-overlapping time windows
            num_windows = int(np.ceil(total_duration / window_size))
            window_edges = np.linspace(0, total_duration, num_windows + 1)

            # Count spikes in each window for the current unit
            spike_counts = np.histogram(spike_train, bins=window_edges)[0]

            # Convert spike counts to firing rates (spikes per second) for this unit
            unit_firing_rates = spike_counts / window_size

            # Append non-NaN and non-inf firing rates to the global list
            firing_rates.extend([fr for fr in unit_firing_rates if not np.isnan(fr) and not np.isinf(fr)])

    # Calculate CoV FR if there are valid firing rates
    if len(firing_rates) > 1:
        mean_fr = np.mean(firing_rates)
        std_fr = np.std(firing_rates)
        CoV_fr = std_fr / mean_fr
        return CoV_fr
    else:
        # Return None if there are no valid firing rates
        return None

def get_mean_isi(recording_object, sorting_object, sampling_rate=10000):
    all_means = []
    units = sorting_object.get_unit_ids()
    for unit in units:
        spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
        assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
        if len(spike_train) > 1:
            isi = np.diff(spike_train)
            mean_isi = np.mean(isi)
            if not np.isnan(mean_isi) and not np.isinf(mean_isi):
                all_means.append(mean_isi)
    
    overall_mean = np.mean(all_means) if all_means else None
    return overall_mean if overall_mean is not None else None  # already in seconds

def get_CoV_isi(recording_object, sorting_object, sampling_rate=10000):
    all_CoVs = []
    units = sorting_object.get_unit_ids()
    for unit in units:
        spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
        assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
        if len(spike_train) > 1:
            isi = np.diff(spike_train)  # already in seconds
            if len(isi) > 1:
                CoV = np.std(isi) / np.mean(isi)
                if not np.isnan(CoV) and not np.isinf(CoV):
                    all_CoVs.append(CoV)
    
    overall_mean_CoV = np.mean(all_CoVs) if all_CoVs else None
    return overall_mean_CoV

def get_unit_fr(recording_object, sorting_object, unit, sampling_rate=10000):
    spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
    assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
    if len(spike_train) > 1:
        duration = recording_object.get_total_duration()  # seconds
        fr = len(spike_train) / duration  # spikes per second
        return fr
    else:
        return 0

def get_unit_fr_CoV(recording_object, sorting_object, unit, sampling_rate=10000):
    spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
    assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
    if len(spike_train) > 1:
        isi = np.diff(spike_train)  # inter-spike intervals
        if len(isi) > 1:
            CoV = np.std(isi) / np.mean(isi)
            return CoV
        else:
            return None
    else:
        return None
    
def get_unit_mean_isi(recording_object, sorting_object, unit, sampling_rate=10000):
    spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
    assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
    if len(spike_train) > 1:
        isi = np.diff(spike_train)
        mean_isi = np.mean(isi)
        return mean_isi  # already in seconds
    else:
        return None
    
def get_unit_isi_CoV(recording_object, sorting_object, unit, sampling_rate=10000):
    spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
    assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
    if len(spike_train) > 2:
        isi = np.diff(spike_train)  # already in seconds
        if len(isi) > 1:
            CoV = np.std(isi) / np.mean(isi)
            return CoV
        else:
            return None
    else:
        return None

def get_unit_fano_factor(recording_object, sorting_object, unit, sampling_rate=10000, bin_size=6.0):
    """
    Computes the Fano Factor (FF) for a given unit, which is the variance-to-mean ratio
    of spike counts in a 6-second time window.

    Parameters:
    - recording_object: SpikeInterface recording object.
    - sorting_object: SpikeInterface sorting object.
    - unit: The unit ID to analyze.
    - sampling_rate (int): Sampling rate in Hz (default: 10000 Hz).
    - bin_size (float): Window size in seconds for spike count bins (default: 6s).

    Returns:
    - Fano Factor (float) or None if not enough spikes.
    """

    # Get spike train and convert to seconds
    spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate  
    assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'

    # Ensure sufficient spikes
    if len(spike_train) < 2:
        return None  # Not enough data

    # Get total duration of recording (in seconds)
    duration = recording_object.get_total_duration()
    
    # Define bin edges for 6-second bins
    num_bins = max(1, int(duration / bin_size))  # Ensure at least 1 bin
    bin_edges = np.linspace(0, duration, num_bins + 1)
    
    # Compute spike counts in bins
    spike_counts, _ = np.histogram(spike_train, bins=bin_edges)
    
    # Compute mean and variance of spike counts
    mean_count = np.mean(spike_counts)
    var_count = np.var(spike_counts)

    # Compute Fano Factor (variance-to-mean ratio)
    fano_factor = var_count / mean_count if mean_count > 0 else None

    return fano_factor

def extract_metrics_from_experimental_data(spike_times, timeVector, spike_times_by_unit, **kwargs):
    
    #extract spiking metrics from experimental data
    recording_object = kwargs['recording_object']
    sorting_object = kwargs['sorting_object']
    sampling_rate = recording_object.get_sampling_frequency()    
    
    #add immediately available data to network_data
    network_data['source'] = 'experimental'
    network_data['timeVector'] = timeVector
    network_data['spiking_data']['spike_times'] = spike_times
    network_data['spiking_data']['spiking_times_by_unit'] = spike_times_by_unit
    
    #overall mean firing rate
    network_data['spiking_data']['spiking_summary_data']['MeanFireRate'] = get_mean_fr(
        recording_object, sorting_object, sampling_rate=sampling_rate
        )
    
    #overall CoV firing rate
    network_data['spiking_data']['spiking_summary_data']['CoVFireRate'] = get_CoV_fr_experimental(
        recording_object, sorting_object, sampling_rate=sampling_rate
        )
    
    #overall mean ISI
    network_data['spiking_data']['spiking_summary_data']['MeanISI'] = get_mean_isi(
        recording_object, sorting_object, sampling_rate=sampling_rate
    )
    
    #overall CoV ISI
    network_data['spiking_data']['spiking_summary_data']['CoV_ISI'] = get_CoV_isi(
        recording_object, sorting_object, sampling_rate=sampling_rate
    )
    
    #spiking data by unit - just go through the simulated version of the data but de-identify pop
    #simulated_spiking_data_by_unit = network_data['simulated_data']['spiking_data_by_unit']
    spiking_data_by_unit = {}
    units = sorting_object.get_unit_ids()
    for unit in units:
        spiking_data_by_unit[unit] = {
            #'unitProperty': sorting_object.get_unit_property(unit),
            #'SpikeTrain': sorting_object.get_unit_spike_train(unit),
            'FireRate': get_unit_fr(recording_object, sorting_object, unit),
            'fr_CoV': get_unit_fr_CoV(recording_object, sorting_object, unit),
            'meanISI': get_unit_mean_isi(recording_object, sorting_object, unit),
            'isi_CoV': get_unit_isi_CoV(recording_object, sorting_object, unit),
            'fano_factor': get_unit_fano_factor(recording_object, sorting_object, unit), # HACK: Hardcoded 6 second bin size
                                                                                        # consistent with the 6 second bin size
                                                                                        # user here: https://pmc.ncbi.nlm.nih.gov/articles/PMC3434456/
            'spike_times': spike_times_by_unit[unit],
        }
    network_data['spiking_data']['spiking_data_by_unit'] = spiking_data_by_unit

def get_spike_times(recording_object, sorting_object, sampling_rate=10000):
    spike_times = []
    units = sorting_object.get_unit_ids()
    total_duration = recording_object.get_total_duration() #seconds
    for unit in units:
        spike_train = sorting_object.get_unit_spike_train(unit) #in samples
        spike_times.extend(spike_train)
    spike_times.sort()
    spike_times = np.array(spike_times) / sampling_rate #convert to seconds
    
    # quality control
    # spike_times_max need to be rounded to the nearest thousandths place to avoid floating point errors
    spike_times_max = np.round(max(spike_times), 3)
    
    print(f'Max Spike Time: {max(spike_times)}')
    #assert max(spike_times) <= total_duration, 'Spike times are not in seconds'
    assert spike_times_max <= total_duration, 'Spike times are not in seconds'
    assert all(spike_time >= 0 for spike_time in spike_times), 'Spike times contain negative values'
    
    #return
    return spike_times

def get_time_vector(recording_object, sampling_rate=10000):
    duration = recording_object.get_total_duration() #seconds
    time_vector = np.linspace(0, duration, int(duration * sampling_rate))
    assert len(time_vector) == int(duration * sampling_rate), 'Time vector length mismatch'
    return time_vector

def get_spike_times_by_unit(sorting_object, sampling_rate=10000):
    spike_times_by_unit = {}
    units = sorting_object.get_unit_ids()
    for unit in units:
        spike_train = sorting_object.get_unit_spike_train(unit) / sampling_rate # convert to seconds
        assert all(spike_time >= 0 for spike_time in spike_train), f'Spike times for unit {unit} contain negative values'
        spike_times_by_unit[unit] = spike_train
    return spike_times_by_unit

def init_network_data_dict():
    # global network_data
    # from RBS_network_models.network_data_struct import network_data_dict
    # empty_network_data = network_data_dict
    # network_data = empty_network_data
    # return empty_network_data
    
    empty_network_data = {
        #General Data
        'source': None, # 'simulated' or 'experimental'
        'timeVector': None,
        
        #Simulated Data
        'simulated_data': {
            'soma_voltage': None,
            'E_Gids': None,
            'I_Gids': None,
            'MeanFireRate_E': None,
            'CoVFireRate_E': None,
            'MeanFireRate_I': None,
            'CoVFireRate_I': None,
            'MeanISI_E': None,
            'MeanISI_I': None,
            'CoV_ISI_E': None,
            'CoV_ISI_I': None,
            'spiking_data_by_unit': None, 
        },
        
        #Spiking Data
        'spiking_data': {
            'spike_times': None,
            'spiking_summary_data': {
                #'spike_times': None,
                'MeanFireRate': None,
                'CoVFireRate': None,
                'MeanISI': None,
                'CoV_ISI': None,         
            },
            'spiking_times_by_unit': None,
            'spiking_data_by_unit': None,
        },
        
        #Bursting Data
        'bursting_data': {
            'bursting_summary_data': {
                'baseline': None,
                'MeanWithinBurstISI': None,
                'CovWithinBurstISI': None,
                'MeanOutsideBurstISI': None,
                'CoVOutsideBurstISI': None,
                'MeanNetworkISI': None,
                'CoVNetworkISI': None,
                'NumUnits': None,
                'Number_Bursts': None,
                'mean_IBI': None,
                'cov_IBI': None,
                'mean_Burst_Peak': None,
                'cov_Burst_Peak': None,
                'fano_factor': None,
            },
            'bursting_data_by_unit': None,
        },
        
        # Neuron Classification Data
        'classification_data': {
            'waveform_data': None,
            'classification_data': None,
            'excitatory_neurons': None,
            'inhibitory_neurons': None,            
        },
    }
    global network_data
    network_data = empty_network_data
    return empty_network_data
    
    ## Keep for reference
    
            # return {
        #     # 'burstPeakValues': None,
        #     # 'IBIs': None,
        #     # 'baseline': None,
        #     # 'peakFreq': None,
        #     # 'firingRate': None,
        #     # 'burstPeakTimes': None,
        #     # 'timeVector': None,
        #     # 'threshold': None,
            
        #     'Number_Bursts': None,
        #     'mean_IBI': None,
        #     'cov_IBI': None,
        #     'mean_Burst_Peak': None,
        #     'cov_Burst_Peak': None,
        #     'fano_factor': None,
        #     'MeanWithinBurstISI': None,
        #     'CoVWithinBurstISI': None,
        #     'MeanOutsideBurstISI': None,
        #     'CoVOutsideBurstISI': None,
        #     'MeanNetworkISI': None,
        #     'CoVNetworkISI': None,
        #     'NumUnits': None,
        #     #'fileName': None
        # }

def get_experimental_network_activity_metrics(sorting_object, recording_segment_object, conv_params, mega_params, **kwargs):
    
    assert sorting_object is not None, 'No sorting object provided'
    
    #initialize network_data
    network_data = init_network_data_dict()
    
    #network data pathing info
    network_data['recording_path'] = recording_segment_object._kwargs.get('file_path', None)
    network_data['sorting_output'] = sorting_object._kwargs.get('folder_path', None)
    assert network_data['sorting_output'] is not None, 'No sorting object source found'
    
    # waveform info
    network_data['waveform_output'] = network_data['sorting_output'].replace('sorted', 'waveforms').replace('sorter_output', '')
    
    #get data
    sorting_object = sorting_object
    recording_object = recording_segment_object
    sampling_rate = recording_object.get_sampling_frequency()
    
    # analyze each neuron. classify if possible. Excitatory (RS) or Inhibitory (FS)
    #from RBS_network_models.experimental import classify_neurons
    #network_data['waveform_output'] = network_data['sorting_output'].replace('sorted', 'neurons').replace('sorter_output', '')
    #network_data['waveform_output'] = network_data['waveform_output'].replace('sorter', 'waveform')
    #classifications = classify_neurons(sorting_object, recording_object, sampling_rate=sampling_rate, output_path=network_data['waveform_output']) 
    
    #convert time to seconds - get initially available data
    spike_times = get_spike_times(recording_object, sorting_object, sampling_rate=sampling_rate) #seconds
    timeVector = get_time_vector(recording_object, sampling_rate=sampling_rate) #seconds
    spike_times_by_unit = get_spike_times_by_unit(sorting_object, sampling_rate=sampling_rate) 
    
    #add conv_params to kwargs
    kwargs['conv_params'] = conv_params
    kwargs['mega_params'] = mega_params
    
    #add sorting_object and recording_object to kwargs if not already present
    kwargs['sorting_object'] = sorting_object
    kwargs['recording_object'] = recording_object
    
    #extract spiking metrics from simulated data
    try: 
        extract_metrics_from_experimental_data_v2(spike_times, timeVector, spike_times_by_unit, **kwargs)
    except Exception as e:
        print(f'Error extracting metrics from simulated data: {e}')
        pass
    
    #extract bursting metrics from simulated data (but this one works for both simulated and experimental data)
    try: 
        extract_bursting_activity_data(spike_times, spike_times_by_unit, **kwargs)
    except Exception as e:
        print(f'Error calculating bursting activity: {e}')
        pass
    
    #classify neurons
    from RBS_network_models.experimental import classify_neurons
    try:
        network_data = classify_neurons(network_data, **kwargs)
    except Exception as e:
        print(f'Error classifying neurons: {e}')
        pass
    
    return network_data

def get_simulated_network_activity_metrics_nvm(conv_params, mega_params, **kwargs):
    
    #initialize network_data
    network_data = init_network_data_dict()
    
    #add conv_params to kwargs
    kwargs['conv_params'] = conv_params
    kwargs['mega_params'] = mega_params
    
    #extract spiking metrics from simulated data
    try: 
        extract_metrics_from_experimental_data(spike_times, timeVector, spike_times_by_unit, **kwargs)
    except Exception as e:
        print(f'Error extracting metrics from simulated data: {e}')
        pass
    
    #extract bursting metrics from simulated data (but this one works for both simulated and experimental data)
    try: 
        extract_bursting_activity_data(spike_times, spike_times_by_unit, **kwargs)
    except Exception as e:
        print(f'Error calculating bursting activity: {e}')
        pass
    
    #classify neurons
    from RBS_network_models.experimental import classify_neurons
    try:
        network_data = classify_neurons(network_data, **kwargs)
    except Exception as e:
        print(f'Error classifying neurons: {e}')
        pass
    
    return network_data
