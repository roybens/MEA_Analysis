import maxlab
import maxlab.chip
import maxlab.system
import maxlab.pycompat as pycompat

import time
import struct

def power_down_all_readout_channels_except(channel):
        block = int(channel / 128)
        value = (2**8-1) - 2**block
        maxlab.send(maxlab.chip.Bias(0,0,value))

def get_mean():
    ret = maxlab.send_raw("system_mean ")
    size_str = ret.split(',')[0]
    size = int(size_str)
#    buf = ret[len(size_str)+1:]
    buf = pycompat.encode(ret[len(size_str)+1:])
    fmt = '=%df' % size
    return list(struct.unpack(fmt, buf))
    
def binary_search_code(readout_channel, target, min_code=0, max_code=1024, sleeptime=1.2, debug=False):
        span = int((max_code - min_code)/2)
        last_code = int((max_code - min_code)/2)
        current_code = last_code
        
        while True:
            # Set the DAC code
            maxlab.send(maxlab.chip.DAC(0, current_code))
            
            # Sleep for the given amount
            time.sleep(sleeptime)
            
            # Get the deviation from target
            mean = get_mean()
            
            diff = (mean[readout_channel] - target)           
            span = span/2.0
            current_code = int(current_code + span) if diff > 0.5 else int(current_code - span)
            
            if debug:
                print(span , current_code, mean[readout_channel], target, diff)
                
            if current_code == last_code:
                break
                
            last_code = current_code
            
        return last_code


class StimulationUnitCharacterizer(object):
    
    def __init__(self):
        core = maxlab.chip.Core()
        core.enable_stimulation_power(True)
        core.use_external_port(True)
        maxlab.send(core)
        self.array = maxlab.chip.Array('current')
        
    def connect_stimulation_unit_to_external_port(self, unit_no):
        stim_unit = maxlab.chip.StimulationUnit(unit_no)
        readout_channel = stim_unit.get_readout_channel()
        
        # Power down all NOT needed blocks of readout channels, to avoid kickback noise
        # power_down_all_readout_channels_except(readout_channel)
        
        # Power down all stimulation buffers
        maxlab.chip.power_down_all_stimulation_buffers()
        
        # Connect the corresponding readout channel to the stimBuffer
        self.array.connect_amplifier_to_stimulation(readout_channel)
        
        # Connect the corresponding readout channel to the VRef ring of the array
        self.array.connect_amplifier_to_ringnode(readout_channel)
        
        self.array.download()
        
        # Power up the correct stimulation unit
        stim_unit.power_up(True).connect(True).set_current_mode().set_large_current_range().dac_source(0)
        maxlab.send(stim_unit)
        return stim_unit
        
    def disconnect_stimulation_unit_from_external_port(self, unit_no):
        stim_unit = maxlab.chip.StimulationUnit(unit_no)
        readout_channel = stim_unit.get_readout_channel()
        
        # Power down all stimulation buffers
        maxlab.chip.power_down_all_stimulation_buffers()
        
        # Disconnect the corresponding readout channel to the stimBuffer
        self.array.disconnect_amplifier_from_stimulation(readout_channel)
        
        # Disconnect the corresponding readout channel to the VRef ring of the array
        self.array.disconnect_amplifier_from_ringnode(readout_channel)
        
    def characterize(self, unit_no):
        stim_unit = maxlab.chip.StimulationUnit(unit_no)
        readout_channel = stim_unit.get_readout_channel()
        
        # Connect the 1M Ohm resistor to ProbeIn
        sw = maxlab.system.Switches(1, 0, 0, 0, 1, 0, 0, 1)
        maxlab.send(sw)
        
        # Connect stimulation unit to external port
        stim_unit = self.connect_stimulation_unit_to_external_port(unit_no)
        
        # Disconnect the stimulation unit to read tearget mean
        stim_unit.connect(False)
        maxlab.send(stim_unit)
        
        # Get the target mean value to achieve
        time.sleep(2)
        target = get_mean()[readout_channel]
        time.sleep(2)
        
        # Power up the stimulation unit
        stim_unit.connect(True)
        maxlab.send(stim_unit)
        
        # Binary search for the correct code
        code = binary_search_code(readout_channel, target)
        
        # Disconnect stimulation unit from external port
        self.disconnect_stimulation_unit_from_external_port(unit_no)
        
        # Load the original switch configuration
        sw = maxlab.system.Switches(1, 0, 0, 0, 1, 0, 0, 0)
        maxlab.send(sw)
        
        return code
        
        
        
    
    
