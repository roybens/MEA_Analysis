import maxlab
import maxlab.util
import maxlab.chip
import maxlab.characterize

import random
import time

def stimulation_pulse(dac_code):
    seq = maxlab.Sequence()
    for amplitude in range(10, 26, 5):
        for rep in range(1,10):
            seq.append( maxlab.chip.DAC(0, dac_code+amplitude) )
            seq.append( maxlab.system.DelaySamples(4) )
            seq.append( maxlab.chip.DAC(0, dac_code-amplitude) )
            seq.append( maxlab.system.DelaySamples(4) )
            seq.append( maxlab.chip.DAC(0, dac_code) )
    
        # Wait 100 ms between two pulses (2000 samples)
        seq.append( maxlab.system.DelaySamples(2000) )
    # Wait 1 second between two sequences (20000 samples
    seq.append( maxlab.system.DelaySamples(20000) )
    
    return seq

maxlab.util.initialize()

# Create the class to characterize the electrodes.
characterizer = maxlab.characterize.StimulationUnitCharacterizer()

# Select 1024 electrodes at random
electrodes = random.sample(range(0, 26400), 1024)
stimulation_electrodes = random.sample(electrodes, 2)

# Smaller test case. Uncomment to test these electrodes
#electrodes = [1000,1001,1002,1003]
#stimulation_electrodes = [1220,1221,1222,1223] #random.sample(electrodes, 32)

array = maxlab.chip.Array('stimulation')
array.reset()
array.clear_selected_electrodes( )
array.select_electrodes( electrodes )
array.select_stimulation_electrodes( stimulation_electrodes )
array.route()
array.download()
time.sleep(2)

# Run the offset compensation
maxlab.util.offset()
time.sleep(5)

for stim_el in stimulation_electrodes:
    # Connect the electorde to stimulation unit
    array.connect_electrode_to_stimulation( stim_el )
    stim = array.query_stimulation_at_electrode( stim_el )
    if stim:
        stimulation_units.append( stim )
    else:
        print("No stimulation channel can connect to electrode: " + str(stim_el))
        continue
        
    # Get the 0V DAC Code
    dac_code = characterizer.charactierize(stim)
    
    # It's now important to correctly set the gain and download the array
    # Set the gain to 512
    amplifier = maxlab.chip.Amplifier().set_gain(512)
    maxlab.send(amplifier)
    
    array.download()
    time.sleep(5)
    
    # Send the stimulation pulse
    print("Power up stimulation unit " + str(stimulation_unit))
    stim = maxlab.chip.StimulationUnit(stimulation_unit)
    stim.power_up(True).connect(True).set_current_mode().dac_source(0)
    maxlab.send(stim)
    print("Send pulse")
    stimulation_pulse(dac_code).send()
    # Sleep for some time to let the stimulation settle
    time.sleep(2)
    print("Power down stimulation unit " + str(stimulation_unit))
    stim = maxlab.chip.StimulationUnit(stimulation_unit).power_up(False) 
    maxlab.send(stim)
    time.sleep(2)
    

