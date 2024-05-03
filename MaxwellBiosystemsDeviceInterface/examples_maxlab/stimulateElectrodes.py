#!/usr/bin/python

import maxlab
import maxlab.system
import maxlab.chip
import maxlab.util

import random
import time

## Stimulation script to stimulate 32 electrodes with biphasic pulses
#
# 0. Initialize system into a defined state
# 
#    It is best to initialize the system into a defined state before
#    Starting any script. This way, one can be sure that the system
#    is always in the same state while running the script, regardless
#    of what has been done before with it.
#
#    In this step, we also issue the following command:
# 
#       maxlab.send(maxlab.chip.Core().enable_stimulation_power(True))
# 
#    This is needed, because by default the stimulation units are
#    powered off globally...
# 
# 
# 1. Select 1024 electrodes at random.
# 
#    For illustration purposes we select 1024 electrodes at random. Other
#    schemes are possible. For example large consecutive blocks of
#    electrodes or specific electrodes based on prior data analysis
#    done e.g. in Matlab.
# 
# 2. Connect stimulation units to 32 electrodes.
# 
#    Once an array configuration has been obtained, either through routing
#    or through loading a previous configuration, the stimulation units
#    can be connected to the desired electrodes.
#    With this step, one needs to be careful, in rare cases it can happen
#    that an electrode cannot be stimulated. For example, the electrode
#    could not be routed (due to routing constraints).
#    
# 3. Power up and configure the stimulation units
# 
#    Once the electrodes are connected to the stimulation units (this was
#    done in the step above), the stimulation units need to be configured
#    and powererd up.
#    
# 4. Apply the same pulses to all 32 units at the same time
# 
#    The stimulation units can be controlled through three independent
#    sources, what we call DAC channels (for digital analog converter).
#    By programming a DAC channel with digital values, we can control
#    the output the stimulation units. DAC inputs are in the range between
#    0 to 1023 bits, whereas 512 corresponds to zero volt and one bit
#    corresponds to 2.9 mV. Thus, to give a pulse of 100mV, the DAC
#    channel temporarily would need to be set to 512 + 34 (100mV/2.9)
#    and back again to 512.
#    In this example, all 32 units are controlled through the same DAC
#    channel ( dac_source(0) ), thus by programming a buphasic pulse
#    on DAC channel 0, all the stimulation units exhibit the biphasic
#    pulse.
# 
# 5. Apply different pulses to the units sequentially
#
#    In this example, we go through the stimulation units one by one
#    and program a series of pulses on each one sequentially
#


# 0. Initialize system into a defined state
maxlab.util.initialize()
maxlab.send(maxlab.chip.Core().enable_stimulation_power(True))



######################################################################
# 1. Select 1024 electrodes at random.
######################################################################

electrodes = random.sample(range(0, 26400), 1024)
stimulation_electrodes = random.sample(electrodes, 32)

# Smaller test case. Uncomment to test these electrodes
#electrodes = [1000,1001,1002,1003]
#stimulation_electrodes = [1220,1221,1222,1223] #random.sample(electrodes, 32)

array = maxlab.chip.Array('stimulation')
array.reset()
array.clear_selected_electrodes( )
array.select_electrodes( electrodes )
array.select_stimulation_electrodes( stimulation_electrodes )
array.route()


######################################################################
# 2. Connect stimulation units to the stimulation_electrodes
######################################################################

stimulation_units = []

for stim_el in stimulation_electrodes:
    array.connect_electrode_to_stimulation( stim_el )
    stim = array.query_stimulation_at_electrode( stim_el )
    if stim:
        stimulation_units.append( stim )
    else:
        print("No stimulation channel can connect to electrode: " + str(stim_el))


array.download()

maxlab.util.offset()


######################################################################
# 3. Power up and configure the stimulation units
######################################################################

print("Setup stim units")
stimulation_unit_commands = []

for stimulation_unit in stimulation_units:
    # Stimulation Unit
    stim = maxlab.chip.StimulationUnit(stimulation_unit)
    stim.power_up(True)
    stim.connect(True)
    stim.set_voltage_mode()
    stim.dac_source(0)
    stimulation_unit_commands.append(stim)
    maxlab.send(stim)


######################################################################
# 4. Apply the same pulses to all 32 units at the same time
######################################################################

# Create a stimulation pulse
# When the stimulation buffers are set to voltage mode, they act like
# an inverting amplifier. Thus here we need to program a negative first
# biphasic pulse, to get a positive first pulse on the chip.
def append_stimulation_pulse(seq, amplitude):
    seq.append( maxlab.chip.DAC(0, 512-amplitude) )
    seq.append( maxlab.system.DelaySamples(4) )
    seq.append( maxlab.chip.DAC(0, 512+amplitude) )
    seq.append( maxlab.system.DelaySamples(4) )
    seq.append( maxlab.chip.DAC(0, 512) )
    return seq


#for pulse in range(1,10):
#    seq = maxlab.Sequence()
#    for rep in range(1,30):
#        append_stimulation_pulse(seq, 200)
#        seq.append( maxlab.system.DelaySamples(200) )
#    print("Send pulse")
#    seq.send()
#    time.sleep(2)


######################################################################
# 5. Apply different pulses to the units sequentially
######################################################################

# First, poweroff all the stimulation units
for stimulation_unit in range(0,32):
    # Stimulation Unit
    stim = maxlab.chip.StimulationUnit(stimulation_unit)
    stim.power_up(False)
    stim.connect(False)
    maxlab.send(stim)

# Prepare sequence of stimulation pulses with increasing amplitudes
# Below we'll repeatedly deliver this sequence of pulses for each
# electrode/stimulation_unit pair
seq = maxlab.Sequence()
for amplitude in range(10, 26, 5):
    for rep in range(1,30):
        append_stimulation_pulse(seq, amplitude)
        # Wait 100 ms between two pulses (2000 samples)
        seq.append( maxlab.system.DelaySamples(2000) )
    # Wait 1 second between two sequences (20000 samples
    seq.append( maxlab.system.DelaySamples(20000) )


# Iterate through all units:
# Power up one unit
# Deliver a series of stimulation pulses
# Power down the unit
for stimulation_unit in stimulation_units:
    print("Power up stimulation unit " + str(stimulation_unit))
    stim = maxlab.chip.StimulationUnit(stimulation_unit)
    stim.power_up(True).connect(True).set_voltage_mode().dac_source(0)
    maxlab.send(stim)
    print("Send pulse")
    seq.send()
    print("Power down stimulation unit " + str(stimulation_unit))
    stim = maxlab.chip.StimulationUnit(stimulation_unit).power_up(False)
    maxlab.send(stim)
    time.sleep(2)

