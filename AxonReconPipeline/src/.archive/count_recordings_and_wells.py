import spikeinterface.extractors as se
import spikeinterface.extractors as se

def count_recordings_and_wells(file_path, scan_type, verbose=False):
    recording_count = 0
    well_num = 0
    rec_num = 0
    well_names = []
    recordings_per_well = {}
    wells_per_chip = 0

    while True:
        well_name = 'well' + str(well_num).zfill(3)
        well_recording_count = 0
        well_found = False

        try:
            while True:
                rec_name = 'rec' + str(well_recording_count).zfill(4)
                recording = se.read_maxwell(file_path,rec_name=rec_name, stream_id=well_name)
                recording_count += 1
                rec_num += 1
                well_recording_count += 1
                well_found = True
                if scan_type == 'Network':
                    break
        except:
            pass

        if well_found:
            well_num += 1
            wells_per_chip += 1
            well_names.append(well_name)
            recordings_per_well[well_name] = well_recording_count
        else:
            try:
                rec_name = 'rec' + str(rec_num).zfill(4)
                recording = se.read_maxwell(file_path,rec_name=rec_name, stream_id=None)
                recording_count += 1
                rec_num += 1
            except:
                break

    return wells_per_chip, recordings_per_well, recording_count 

# Path: src/axon_trace_classes.py
