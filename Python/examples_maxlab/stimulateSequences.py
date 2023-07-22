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
# 1. Load a previously created configuration
# 
#    For this example, we assume that a previous configuration has
#    been generated. Either through the Python API or through the Scope GUI.
#    We need to know two electrodes from this configuration which are
#    routed, so that we can attempt to stimulate them.
# 
# 2. Connect two electrodes to stimulation units and power up stimulation units
# 
#    In rare cases it can happen that the selected electrode cannot be
#    stimulated. So, always check whether the electrode got properly
#    connected. As it is done in this example.
#    If the electrode is not properly connected, the stimulation has no effect.
#    
# 3. Prepare two different sequences of pulse trains
# 
#    An almost infinite amount (only limited by the computer memory) of
#    independent stimulation pulse trains can be prepared without actually
#    deliver them yet.
#    
# 4. Deliver the two pulse trains randomly
# 
#    The previously prepared pulse trains can be delivered whenever
#    seems reasonable, or following a specific stimulation schedudle.
# 


name_of_configuration = "/home/user/my_config.cfg"
electrode1 = 13083
electrode2 = 13095

######################################################################
# 0. Initialize system into a defined state
######################################################################

maxlab.util.initialize()
maxlab.send(maxlab.chip.Core().enable_stimulation_power(True))


######################################################################
# 1. Load a previously created configuration
######################################################################

array = maxlab.chip.Array('stimulation')
array.load_config( name_of_configuration )

######################################################################
# 2. Connect two electrodes to stimulation units and power up stimulation units
######################################################################

array.connect_electrode_to_stimulation( electrode1 )
array.connect_electrode_to_stimulation( electrode2 )

stimulation1 = array.query_stimulation_at_electrode( electrode1 )
stimulation2 = array.query_stimulation_at_electrode( electrode2 )

if not stimulation1:
    print("Error: electrode: " + str(electrode1) + " cannot be stimulated")

if not stimulation2:
    print("Error: electrode: " + str(electrode2) + " cannot be stimulated")

# Download the prepared array configuration to the chip
array.download()

# Prepare commands to power up and power down the two stimulation units
cmd_power_stim1 = maxlab.chip.StimulationUnit(stimulation1).power_up(True).connect(True).set_voltage_mode().dac_source(0)
cmd_power_down_stim1 = maxlab.chip.StimulationUnit(stimulation1).power_up(False)
cmd_power_stim2 = maxlab.chip.StimulationUnit(stimulation2).power_up(True).connect(True).set_voltage_mode().dac_source(0)
cmd_power_down_stim2 = maxlab.chip.StimulationUnit(stimulation2).power_up(False)


######################################################################
# 3. Prepare two different sequences of pulse trains
######################################################################

def append_stimulation_pulse(seq, amplitude):
    seq.append( maxlab.chip.DAC(0, 512-amplitude) )
    seq.append( maxlab.system.DelaySamples(4) )
    seq.append( maxlab.chip.DAC(0, 512+amplitude) )
    seq.append( maxlab.system.DelaySamples(4) )
    seq.append( maxlab.chip.DAC(0, 512) )
    return seq


# Create a sequence with 4 bursts with increasing amplitude.
# 5 pulses per burst
sequence1 = maxlab.Sequence()
sequence1.append( cmd_power_stim1 )
for amplitude in range(10, 26, 5):
    for rep in range(1,5):
        append_stimulation_pulse(sequence1, amplitude) 
        # Wait 10 ms between two pulses (200 samples)
        sequence1.append( maxlab.system.DelaySamples(200) )
    # Wait 100 ms between two sequences (2000 samples
    sequence1.append( maxlab.system.DelaySamples(2000) )
sequence1.append( cmd_power_down_stim1 )


# Create a sequence with 4 bursts with increasing time between the pulses.
# 5 pulses per burst
sequence2 = maxlab.Sequence()
sequence2.append( cmd_power_stim2 )
for delta_time in range(200, 501, 100):
    for rep in range(1,5):
        append_stimulation_pulse(sequence2, 15) 
        # Wait 10 ms between two pulses (200 samples)
        sequence2.append( maxlab.system.DelaySamples(delta_time) )
    # Wait 100 ms between two sequences (2000 samples
    sequence2.append( maxlab.system.DelaySamples(2000) )
sequence2.append( cmd_power_down_stim2 )


######################################################################
# 4. Deliver the two pulse trains randomly
######################################################################


sequence1.send()
time.sleep(4)

sequence2.send()
time.sleep(3)

sequence1.send()
time.sleep(5)

sequence2.send()
time.sleep(3)


