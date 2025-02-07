import sys
import os

from modules import mea_processing_library as MPL
import modules.lib_helper_functions as helper
from modules.axon_reconstructor import AxonReconstructor

def run_pipeline(h5_parent_dirs, mode='normal', **kwargs):
    '''
    Main function to run the pipeline.
    Supports two modes: 'normal' and 'lean'.
    '''
    h5_files = helper.get_list_of_h5_files(h5_parent_dirs, **kwargs)
    
    if mode == 'normal':
        '''
        Run the pipeline normally
        '''
        reconstructor = AxonReconstructor(h5_files, **kwargs)
        reconstructor.run_pipeline(**kwargs)
    
    elif mode == 'lean':
        '''
        Hacky way to run the pipeline on max two plates, one well at a time to minimize size of temp data
        '''
        for h5_file in h5_files:
            max_two_wells = 6  # 6 wells total
            max_two_wells_analyzed = 0
            
            while max_two_wells_analyzed < max_two_wells:
                # Get reconstructor ID
                h5_details = MPL.extract_recording_details([h5_file])
                date = h5_details[0]['date']
                chipID = h5_details[0]['chipID']
                runID = h5_details[0]['runID']
                scanType = h5_details[0]['scanType']
                
                if 'stream_select' not in kwargs:
                    #kwargs['stream_select'] = max_two_wells_analyzed  # should go 0 through 5
                    stream_select = max_two_wells_analyzed
                else:
                    stream_select = kwargs['stream_select']
                    if max_two_wells_analyzed != kwargs['stream_select']:
                        max_two_wells_analyzed += 1
                        continue
                
                projectName = h5_details[0]['projectName']
                reconstructorID = f'{date}_{chipID}_{runID}_well00{stream_select}'
                
                # Set log files to be unique for each reconstructor
                
                #old naming scheme
                # kwargs['log_file'] = f'{kwargs["output_dir"]}/{projectName}/{reconstructorID}_axon_reconstruction.log'
                # kwargs['error_log_file'] = f'{kwargs["output_dir"]}/{projectName}/{reconstructorID}_axon_reconstruction_error.log'
                # kwargs['recordings_dir'] = os.path.join(kwargs['output_dir'], projectName, 'recordings')
                # kwargs['sortings_dir'] = os.path.join(kwargs['output_dir'], projectName, 'sortings')
                # kwargs['waveforms_dir'] = os.path.join(kwargs['output_dir'], projectName, 'waveforms')
                # kwargs['templates_dir'] = os.path.join(kwargs['output_dir'], projectName, 'templates')
                # kwargs['recon_dir'] = os.path.join(kwargs['output_dir'], projectName, 'reconstructions')
                # kwargs['reconstructor_dir'] = os.path.join(kwargs['output_dir'], projectName, 'reconstructors')
                
                #new naming scheme
                #kwargs['log_file'] = f'{kwargs["output_dir"]}/{projectName}/{date}/{chipID}/{scanType}/{runID}/{wellID}/{reconstructorID}_axon_reconstruction.log'
                wellID = f'well00{stream_select}'
                well_path = f'{kwargs["output_dir"]}/{projectName}/{date}/{chipID}/{scanType}/{runID}/{wellID}'
                kwargs['log_file'] = f'{well_path}/{wellID}_axon_reconstruction.log'
                kwargs['error_log_file'] = f'{well_path}/{wellID}_axon_reconstruction_error.log'
                kwargs['recordings_dir'] = os.path.join(well_path, 'recordings')
                kwargs['sortings_dir'] = os.path.join(well_path, 'sortings')
                kwargs['waveforms_dir'] = os.path.join(well_path, 'waveforms')
                kwargs['templates_dir'] = os.path.join(well_path, 'templates')
                kwargs['recon_dir'] = os.path.join(well_path, 'reconstructions')
                kwargs['reconstructor_dir'] = well_path
                
                print(f'Starting reconstruction for {reconstructorID}...')
                print(f'Log file: {kwargs["log_file"]}')
                print(f'Error log file: {kwargs["error_log_file"]}')
                print(f'Recordings dir: {kwargs["recordings_dir"]}')
                print(f'Sortings dir: {kwargs["sortings_dir"]}')
                print(f'Waveforms dir: {kwargs["waveforms_dir"]}')
                print(f'Templates dir: {kwargs["templates_dir"]}')
                print(f'Reconstructions dir: {kwargs["recon_dir"]}')
                print(f'Reconstructor dir: {kwargs["reconstructor_dir"]}')
                
                # Run the pipeline
                kwargs['stream_select'] = stream_select
                reconstructor = AxonReconstructor([h5_file], **kwargs)
                try: reconstructor.run_pipeline(**kwargs)
                except Exception as e:
                    print(f'Error running pipeline for {reconstructorID}: {e}')
                    print(f'Continuing to next well...')
                kwargs.pop('stream_select', None)
                
                max_two_wells_analyzed += 1