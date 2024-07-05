# MEA_Analysis
MEA_Analysis is a python pipeline which utilizes the [spikeinterface](https://github.com/SpikeInterface) and the kilorsort to get the desired spikes from the mice brain.Our research focused on how the healthy neuron spikes differs from a diseased neuron.
---
# Working:

We take the readings we have 3 different types _'Network Today'_, _'Network last'_, _'Network Best'_. <br />
1.**Network Today:** This refers to the recording of the neuron activity we took on that particular day. <br />
2.**Network last:** This refers to the number of days since the last recording. <br />
3.**Network Best:** This refers to the best available neuron spike recording that exists. Usually after 20 days we get the best recording.<br />
---
# Usage
### Data Preprocessing
1. **Data Preprocessing**: Here we make sure we get the correct the sequence of data as in the data and the runids will match so we can recongnize on what date the neuron activity was recorded.
```python
def preprocess(recording):  ## some hardcoded stuff.
    """
    Does the bandpass filtering and the common average referencing of the signal
    
    """
    recording_bp = spre.bandpass_filter(recording, freq_min=300, freq_max=6000)
    

    recording_cmr = spre.common_reference(recording_bp, reference='global', operator='median')

    recording_cmr.annotate(is_filtered=True)

    return recording_cmr
```
---
## Process Block 
this is our main code where the preprocessed data is being fed into for our processing and extracting the sorter output and to extract the waveforms. 
```python
dir_name = sorting_folder
        #os.mkdir(dir_name,0o777,)
        os.chdir(dir_name)
        logging.debug(f"currentdirectory: {current_directory}") #changes made here 
        #rohan made changes here 
        kilosort_output_folder = f"{current_directory}/../AnalyzedData/{desired_pattern}/kilosort2_{rec_name}"
        print("ks folder:"+ kilosort_output_folder) #changes made here 
        start = timer()
        sortingKS3 = run_kilosort(recording_chunk,output_folder='./Kilosort_tmp')
        logging.debug("Sorting complete")
        sortingKS3 = sortingKS3.remove_empty_units()
        sortingKS3 = spikeinterface.curation.remove_excess_spikes(sortingKS3,recording_chunk)
```
For our purposed we have set the folder to match our recording sequence for easier analysis 
---
In the kilosort function we have utilized spikeinterface's kilosort2, we are using kilosort2 specifically because our data is coming from a immobile mice brains. We have utilized the functionality of docker images to create a custom docker image over the pre-existing kilosort2 docker image to meet our techincal requiremnt. 
```python
def run_kilosort(recording,output_folder):
    logging.debug("run_kilosort_output folder:"+output_folder) #rohan made changes here
    default_KS2_params = ss.get_default_sorter_params('kilosort2')
    default_KS2_params['keep_good_only'] = True
    # default_KS2_params['detect_threshold'] = 12
    # default_KS2_params['projection_threshold']=[18, 10]
    # default_KS2_params['preclust_threshold'] = 14
    #run_sorter=run_sorter_local(sorter_name="kilosort3",recording=recording, output_folder=output_folder, delete_output_folder=False,verbose=True,with_output=True,**default_KS2_paramsdocker_image= "rohanmalige/rohan_si-98:v8")
    #sorting=run_sorter(sorter_name="kilosort2",recording=recording,output_folder=output_folder,remove_existing_folder=True, delete_output_folder=False,verbose=True,docker_image="rohanmalige/rohan_si-98:v8",with_output=True, **default_KS2_params)
    run_sorter = ss.run_sorter('kilosort2',recording=recording, output_folder=output_folder,docker_image= "si-98-ks2-maxwell",verbose=True, **default_KS2_params)
    #run_sorter = ss.run_kilosort2(recording, output_folder=output_folder, docker_image= "si-98-ks2-maxwell",verbose=True, **default_KS2_params) #depreciation warning 
    #sorting_KS3 = ss.Kilosort3Sorter._get_result_from_folder(output_folder+'/sorter_output/')
    return run_sorter
```
---
The next step after the kilosort function has run is to extract the waveforms we send the recording with the kilosort output folder 
```python
def extract_waveforms(recording,sorting_KS3,folder):
   
    folder = Path(folder)
    logging.debug("waveforms_folder"+folder) #rohan made changes here 
    global_job_kwargs = dict(n_jobs=24) 
    si.set_global_job_kwargs(**global_job_kwargs)
    waveforms = si.extract_waveforms(recording=recording,sorting=sorting_KS3,sparse=False,folder=folder,max_spikes_per_unit=500,overwrite=True)
    #waveforms = si.extract_waveforms(recording,sorting_KS3,folder=folder,overwrite=True, sparse = True, ms_before=1., ms_after=2.,allow_unfiltered=True,**job_kwargs)
    return wavefor
```
---
After extracting the waveforms we use the Quality metrics to filter out the waveforms, the quality metrics depends on the firing rate, SNR, RV.
1. **Firing rate**: Firing rate is the number of times the neuron has been active, a few of the neurons might have low firing rate and hence cannot be used to analyze the data.
2. **Signal to Noise Ratio**: When we record we will have noise that corrupts the data so we want to take the data where the signal strenght is strong and the noise is low, hence we want a ratio of the signal to noise ratio high enough.
3. **RV**: we know that a neuron should become active again after certain period of time and shouldn't be active prematuraly before the set time for our data we have kept the threshold as 2ms. So if a neuron becomes active between 2ms and 3ms. We consider this particular recording of that neuron to be corrputed by some other unit.
```python
def compute_quality_metrics(
    waveform_extractor,
    load_if_exists=False,
    metric_names=None,
    qm_params=None,
    peak_sign=None,
    seed=None,
    sparsity=None,
    skip_pc_metrics=False,
    verbose=False,
    **job_kwargs,
):
    """Compute quality metrics on waveform extractor.

    Parameters
    ----------
    waveform_extractor: WaveformExtractor
        The waveform extractor to compute metrics on.
    load_if_exists : bool, default: False
        Whether to load precomputed quality metrics, if they already exist.
    metric_names : list or None
        List of quality metrics to compute.
    qm_params : dict or None
        Dictionary with parameters for quality metrics calculation.
        Default parameters can be obtained with: `si.qualitymetrics.get_default_qm_params()`
    sparsity : dict or None
        If given, the sparse channel_ids for each unit in PCA metrics computation.
        This is used also to identify neighbor units and speed up computations.
        If None (default) all channels and all units are used for each unit.
    skip_pc_metrics : bool
        If True, PC metrics computation is skipped.
    n_jobs : int
        Number of jobs (used for PCA metrics)
    verbose : bool
        If True, output is verbose.
    progress_bar : bool
        If True, progress bar is shown.

    Returns
    -------
    metrics: pandas.DataFrame
        Data frame with the computed metrics
    """
    if load_if_exists and waveform_extractor.is_extension(QualityMetricCalculator.extension_name):
        qmc = waveform_extractor.load_extension(QualityMetricCalculator.extension_name)
    else:
        qmc = QualityMetricCalculator(waveform_extractor)
        qmc.set_params(
            metric_names=metric_names,
            qm_params=qm_params,
            peak_sign=peak_sign,
            seed=seed,
            sparsity=sparsity,
            skip_pc_metrics=skip_pc_metrics,
        )
        qmc.run(verbose=verbose, **job_kwargs)

    metrics = qmc.get_data()

    return metrics
```
---
Lastly we save the waveforms that are good enough for analysis and curation inside a folder called 'waveforms_good'. We also save the qualitymetrics and waveforms in excel for human-readibility and data analysis.
```python
 waveform_good = waveforms.select_units(non_violated_units_new,new_folder=f"{current_directory}/../AnalyzedData/{desired_pattern}/waveforms_good")
        template_metrics = sp.compute_template_metrics(waveform_good)
        #template_metrics = template_metrics.loc[update_qual_metrics.index.values]
        qual_metrics = qm.compute_quality_metrics(waveform_good ,metric_names=['num_spikes','firing_rate', 'presence_ratio', 'snr',
                                                       'isi_violation', 'amplitude_cutoff','amplitude_median'])  ## to do : have to deal with NAN values
        #rohan made change here
        template_metrics.to_excel(f"{current_directory}/../AnalyzedData/{desired_pattern}/template_metrics.xlsx")
        locations = sp.compute_unit_locations(waveform_good)
        qual_metrics['location_X'] = locations[:,0]
        qual_metrics['location_Y'] = locations[:,1]
        #rohan made change here 
        qual_metrics.to_excel(f"{current_directory}/../AnalyzedData/{desired_pattern}/quality_metrics.xlsx")
        ax = plt.subplot(111)
```
----
# License 
 pending 

