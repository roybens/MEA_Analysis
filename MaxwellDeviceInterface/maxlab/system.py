#!/usr/bin/python
# -*- coding: utf-8


class DelaySamples:
    """Insert a delay between two consecutive commands

    Args:
        samples (int): How many samples to delay between two consecutive
            commands. Samples here correspond to the recording
            samples. Time between two recording samples is 50 µs.
    """
    def __init__(self, samples=0):
        self.samples=samples
    def set(self):
        return "system_delay_samples " + str(self.samples)

class MidSupply:
    def __init__(self, reference_voltage=128):
        self.reference_voltage=reference_voltage
    def set(self):
        return "system_set_mid_supply " + str(self.reference_voltage)

class ReferenceStimulationHigh:
    """Sets the high reference voltage for the stimulation buffers.

    Args:
        reference_voltage (int): Value in bits to set the reference
            voltage. Allowed values are between 0 and 4095, whereas
            4095 corresponds to 3.3 volts. When doing voltage controlled
            stimulation, the recommended default value for this reference
            voltage is 2.25 volts, i.e. 2792 bits.

    Refer to the system manual for more details about this reference
    voltage.
    """
    def __init__(self, reference_voltage=2820):
        self.reference_voltage=reference_voltage
    def set(self):
        return "system_set_reference_stimulation_high " + str(self.reference_voltage)

class ReferenceStimulationMiddle:
    def __init__(self, reference_voltage=2048):
        self.reference_voltage=reference_voltage
    def set(self):
        return "system_set_reference_stimulation_middle " + str(self.reference_voltage)

class ReferenceStimulationLow:
    """Sets the low reference voltage for the stimulation buffers.

    Args:
        reference_voltage (int): Value in bits to set the reference
            voltage. Allowed values are between 0 and 4095, whereas
            4095 corresponds to 3.3 volts. When doing voltage controlled
            stimulation, a typically recommended value for this reference
            voltage is 1.05 volts, i.e. 1304.

    Refer to the system manual for more details about this reference
    voltage.
    """
    def __init__(self, reference_voltage=1304):
        self.reference_voltage=reference_voltage
    def set(self):
        return "system_set_reference_stimulation_low " + str(self.reference_voltage)

class ReferenceADCStart:
    def __init__(self, reference_voltage=6):
        self.reference_voltage=reference_voltage
    def set(self):
        return "system_set_reference_adc_start " + str(self.reference_voltage)

class ReferenceADCStop:
    def __init__(self, reference_voltage=45):
        self.reference_voltage=reference_voltage
    def set(self):
        return "system_set_reference_adc_stop " + str(self.reference_voltage)

class ReferenceRampGen:
    def __init__(self, reference_voltage=2048):
        self.reference_voltage=reference_voltage
    def set(self):
        return "system_set_reference_rampgen " + str(self.reference_voltage)

class ReferenceMOSResistor:
    def __init__(self, reference_voltage=4095):
        self.reference_voltage=reference_voltage
    def set(self):
        return "system_set_reference_mos_resistor " + str(self.reference_voltage)

class ReferenceVoltage:
    def __init__(self, reference_voltage=2048):
        self.reference_voltage=reference_voltage
    def set(self):
        return "system_set_reference_voltages " + str(self.reference_voltage)

class VariableReference:
    def __init__(self, reference_voltage=2048):
        self.reference_voltage=reference_voltage
    def set(self):
        return "system_set_variable_reference " + str(self.reference_voltage)

class Switches:
    def __init__(self, sw_0=1, sw_1=1, sw_2=0, sw_3=0, sw_4=0, sw_5=0, sw_6=0, sw_7=0):
        self.sw = [sw_0, sw_1, sw_2, sw_3, sw_4, sw_5, sw_6, sw_7]
    def get(self):
        return "system_get_switches"
    def set(self):
        return "system_set_switches " + " ".join([str(sw) for sw in self.sw])

class GPIODirection:
    """Set the direction of the GPIOs.

    Args:
        direction (int):
            An 8-bit value, where each bit of `direction` controls the
            state of one GPIO, whether one GPIO is either in input or
            in output mode.

    * bit = '0': setting the corresponding GPIO channel into input mode.
    * bit = '1': setting the corresponding GPIO channel into output mode.

    For example:

    >>> # This sets GPIO3 to output mode, GPIO7 - GPIO4 and GPIO2 - GPIO0 to input mode
    >>> maxlab.system.GPIODirection(0b1000)

    When the GPIO is set to input mode, the value applied to the GPIO
    pin is sampled together with each sampling frame of the MaxOne chip
    and the values are recorded into the recording file. Refer to the
    systems manual how to access these values from the recording file.
    """
    def __init__(self, direction=0):
        self.direction=direction
    def set(self):
        return "system_set_gpio_direction " + str(self.direction)

class GPIOOutput:
    """Control the output of the GPIOs.

    Args:
        output (int):
            An 8-bit value specifying the output state of each GPIO channel.

    When the corresponding bit is set to output mode with
    :func:`maxlab.system.GPIODirection`, the output of that channel can
    be controlled through this function.

    The following example generates a pulse on the GPIO channel 4,
    by toggling bit 4:

    >>> maxlab.send( maxlab.system.GPIOOutput(0b1000))
    >>> maxlab.send( maxlab.system.DelaySamples(100))
    >>> maxlab.send( maxlab.system.GPIOOutput(0b0000))
    """
    def __init__(self, output=0):
        self.output=output
    def set(self):
        return "system_set_gpio_output " + str(self.output)

class StatusLED:
    """Control the color of the status LED.

    Args:
        color (3-bit integer):
            Value to control the color of the LED on the Recording Unit.

    Values for the color parameter can be:
     * 0: light blue.
     * 1: cyan.
     * 2: pink
     * 3: dark blue
     * 4: light green
     * 5: dark green
     * 6: red
     * 7: off

    When the experiment is very light sensitive, for example during a retina experiment,
    this LED can be switched off by setting the `color` parameter to 0.
    """
    def __init__(self, color=0):
        self.color=color
    def set(self):
        return "system_set_status_led " + str(self.color)


class StatusOut:
    """Send status bits

    """
    def __init__(self, status=0):
        self.status=status
    def set(self):
        return "system_status_out " + str(self.status)


class Event:
    """ Triggers the sending of event which is then converted to status_out command
    """

    def __init__(self, well_id, event_type, user_id, properties):
        self.well_id = well_id
        self.event_type = event_type
        self.user_id = user_id
        self.properties = properties

    def set(self):
        return "system_trigger_event " + str(self.well_id) + " " + str(self.event_type) + " " + str(self.user_id) + " " + self.properties


