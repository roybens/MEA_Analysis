#!/usr/bin/python
# -*- coding: utf-8

import maxlab.apicomm as comm
import maxlab.config

def power_down_all_stimulation_buffers():
    return StimulationUnit(-1).power_up(False)


class Amplifier:
    """Program gain of the MaxOne on-chip amplifiers.

    Args:
        gain (int): The gain to be set.

    Default gain for this command is 512. Use :func:`maxlab.chip.Amplifier.set_gain` to set a different gain.
    """
    def __init__(self, gain=512):
        self = self.set_gain(gain)

    def params(self,  stage1_bypass=0,
                       stage1_gain=1,
                       stage1_reset_mode=0,
                       stage2_bypass=0,
                       stage2_gain=5,                   # 5 == x16
                       stage3_bypass=0,
                       stage3_gain=0):
        self.settings = [ stage1_bypass,
                          stage1_gain,
                          stage1_reset_mode,
                          stage2_bypass,
                          stage2_gain,                   # 5 == x16
                          stage3_bypass,
                          stage3_gain]
        return self

    def set(self):
        return "mea_set_amplifier " + " ".join([str(settings) for settings in self.settings])

    def set_gain(self,gain):
        """Set different gain values for the amplifier.

        Args:
            gain (int): The gain to be set.

        Raises:
            ValueError: If `gain` is not a valid gain setting.

        Possible gain values are:
         * 1
         * 7
         * 112
         * 512
         * 1024
         * 1025
         * 2048

        Note:
            Other values are not valid and will raise an exception. For example trying to set the gain to 1000:

        >>> try:
        >>>     maxlab.chip.Amplifier().set_gain(1000)
        >>> except ValueError as error:
        >>>     print(error.args[0])
        Not a valid gain parameter: 1000
        """
        if gain==1:
            return self.params(1,0,0,1,0,1,0)
        if gain==7:
            return self.params(0,0,0,1,0,1,0)
        if gain==112:
            return self.params(0,0,0,0,4,0,0)
        if gain==512:
            return self.params(0,1,0,0,5,0,0)
        if gain==1024:
            return self.params(0,1,0,0,5,0,1)
        if gain==1025:
            return self.params(0,1,0,0,6,0,0)
        if gain==2048:
            return self.params(0,1,0,0,6,0,1)
        raise ValueError("Not a valid gain parameter: " + str(gain))


class Core:
    """Control core settings of the MaxOne chip.

    """
    def __init__(self):
        self.params(enable_external_port=False,
                    stimulation_power=True,
                    controller_multiplication=0,
                    output_enable=True,
                    tx_mode=0,
                    reset_mode=0,
                    reset_speed=7
                    )

    def params(self, enable_external_port=False,
                     stimulation_power=True,
                     controller_multiplication=0,
                     output_enable=True,
                     tx_mode=0,
                     reset_mode=0,
                     reset_speed=7):
         self.enable_external_port=enable_external_port
         self.stimulation_power=stimulation_power
         self.controller_multiplication=controller_multiplication
         self.output_enable=output_enable
         self.tx_mode=tx_mode
         self.reset_mode=reset_mode
         self.reset_speed=reset_speed
         return self

    def set(self):
        return "mea_set_core  " + str(int(self.enable_external_port)) + \
                            " " + str(int(self.stimulation_power)) + \
                            " " + str(int(self.controller_multiplication)) + \
                            " " + str(int(self.output_enable)) + \
                            " " + str(int(self.tx_mode)) + \
                            " " + str(int(self.reset_mode)) + \
                            " " + str(int(self.reset_speed))

    def use_external_port(self, enable):
        """Enables the external port access to the array."""
        self.enable_external_port=enable
        return self

    def enable_stimulation_power(self, enable):
        self.stimulation_power=enable
        return self

    def set_controller_multiplication(self, multiplication):
        """Sets the multiplication stage in the Controller.

        Arg:
            `multiplication`: possible values are: 0 (1x), 1 (2x), 2 (4x)
        """
        self.controller_multiplication=multiplication
        return self

    def enable_digital_output(self, enable):
        self.output_enable=enable
        return self

    def set_tx_mode(self, mode):
        self.tx_mode=mode
        return self

    def set_reset_mode(self, mode):
        self.reset_mode=mode
        return self

    def set_reset_speed(self, speed):
        self.reset_speed=speed
        return self


class RampGen:
    """Control RampGen settings of the MaxOne chip.

    """
    def __init__(self, enable=True, amplifier=4):
        self.enable = enable
        self.amplifier = amplifier

    def set(self):
        return "mea_set_rampgen  " + str(int(self.enable)) + \
                               " " + str(int(self.amplifier))


class Controller:
    """Control Controller (SwitchCap) of the MaxOne chip.

    """
    def __init__(self, enable=True, ampMult=0, ampFreq=3, clkMon1=0, clkMon2=0):
        self.enable = enable
        self.ampMult = ampMult
        self.ampFreq = ampFreq
        self.clkMon1 = clkMon1
        self.clkMon2 = clkMon2

    def set(self):
        return "mea_set_controller  " + str(int(self.enable)) + \
                                  " " + str(int(self.ampMult)) + \
                                  " " + str(int(self.ampFreq)) + \
                                  " " + str(int(self.clkMon1)) + \
                                  " " + str(int(self.clkMon2))


class Bias:
    def __init__(self, value=0, address=0, power_down=0):
        self.value=value
        self.address=address
        self.power_down=power_down

    def set(self):
        return "mea_set_bias " +str(self.value) + " " + str(self.address) + " " + str(self.power_down)

class ResetDisconnect:
    def __init__(self, enable=0, broadcast=0, channel=0):
        self.enable=enable
        self.broadcast=broadcast
        self.channel=channel

    def set(self):
        return "mea_set_reset_disconnect " + str(self.enable) + " " + \
                                             str(self.broadcast) + " " + \
                                             str(self.channel)


class DAC:
    """Program the three on-chip DAC (digital-analog converter) channels.

    Args:
        dac_no          (int): Which DAC channel to control.
        dac_code        (int): DAC code in bits.
        dac_code_second (int): this parameter is only relevant when `dac_no=3`.

    Depending on which DAC channel is specified in the command, the
    parameter `dac_code` has a different meaning.

    The meaning of `dac_code` depending on `dac_no` (DAC channel) is:

     * 0: `dac_code` controls DAC0
     * 1: `dac_code` controls DAC1
     * 2: `dac_code` controls DAC2
     * 3: `dac_code` controls DAC0 and `dac_code_second` controls DAC1 at the very same time.

    The last mode `dac_code=3` is useful for applications where the
    outputs of two DAC channels need to be controlled/changed at the
    exact same time.

    Each DAC channel can be programmed with 10 bits, i.e. possible
    values range between 0 and 1023.  To convert between DAC bits and
    stimulation amplitude, one needs to multiply the DAC bits with the
    value from :func:`maxlab.system.query_DAC_lsb()`.

    When the stimulation buffers are in voltage mode, the conversion
    between DAC code and stimulation voltage goes as follows:

    >>> stimVoltage = (dac_code - 512) * maxlab.query_DAC_lsb()

    Note:
        This conversion does not hold when doing current
        stimulation. Refer the the system manual for more details about
        this process.
    """

    def __init__(self, dac_no=0, dac_code=512, dac_code_second=512):
        self.dac_no = dac_no
        self.dac_code = dac_code
        self.dac_code_second = dac_code_second

    def set(self):
        return "mea_set_dac " + str(self.dac_no) + " " + str(self.dac_code) + " " + str(self.dac_code_second)



class StimulationUnit:
    """Program the 32 on-chip stimulation units.

    Args:
        unit_no (int): The stimulation unit to configure.

    Warning:
        When doing voltage stimulation, the amplifier in the stimulation
        unit acts as an inverting amplifier.  A DAC code of 512
        corresponds to mid-supply, i.e. 0-volt stimulation. 512+100
        results in a negative voltage.  512-100 results in a positive
        voltage.

    To convert between bits and stimulation voltage, one needs to know
    the DAC lsb:

    >>> lsb = maxlab.query_DAC_lsb_mV()
    >>> amplitude = lsb * bits
    """
    def __init__(self, unit_no):
        self.unit_no=unit_no
        self.power=0
        self.do_connect=0
        self.current_mode=False
        self.current_range=False
        self.dac=0
        self.ref_source=False
        self.mapping = \
              [ ( 0, 796), ( 1, 640), ( 2, 556), ( 3, 384), ( 4, 364), ( 5, 128),
                ( 6,  84), ( 7,  16), ( 8,  50), ( 9, 130), (10, 366), (11, 386),
                (12, 558), (13, 642), (14, 798), (15, 898), (16, 896), (17, 893),
                (18, 773), (19, 641), (20, 533), (21, 385), (22, 129), (23,  85),
                (24,   1), (25,  51), (26, 131), (27, 387), (28, 535), (29, 643),
                (30, 775), (31, 895)]

    def set(self):
        return "mea_set_stimulation_unit " + str(int(self.unit_no)) +\
                                         " " + str(int(self.power)) +\
                                         " " + str(int(self.do_connect)) +\
                                         " " + str(int(self.current_mode)) +\
                                         " " + str(int(self.current_range)) +\
                                         " " + str(int(self.dac)) +\
                                         " " + str(int(self.ref_source))

    def power_up(self, power_up):
        """Power-up this stimulation unit.

        Args:
            power_up (bool): Enable on-chip power for this particular
                stimulation unit.
        """
        self.power=power_up
        return self

    def connect(self, connect): # i.e. autozeroing
        """Connect or disconnect the stimulation unit.

        Args:
            connect (bool): connect or disconnect the output of the unit.
        """
        self.do_connect=connect
        return self

    def set_current_mode(self):
        """Set the stimulation unit into current mode."""
        self.current_mode=True
        return self

    def set_voltage_mode(self):
        """Set the stimBuffer unit into voltage mode."""
        self.current_mode=False
        return self

    def set_large_current_range(self): #'large' / 'small'):
        """Set current range to `large` currents when in current mode.

        Refer to the system manual to what large means.
        """
        self.current_range=True
        return self

    def set_small_current_range(self): #'large' / 'small'):
        """Set current range to `small` currents when in current mode.

        Refer to the system manual to what small means.
        """
        self.current_range=False
        return self

    def dac_source(self, dac):
        """Choose the DAC channel for the stimulation unit.

        Args:
            dac (int) : The DAC channel to choose for the unit.

        Possible values for DAC channels are:
         * 0=DAC0
         * 1=DAC1
         * 2=DAC2
        """
        self.dac=dac
        return self

    def external_reference(self, ref_source=True):
        """Use an external reference voltage for this stimulation unit
        instead of a DAC channel.

        Args:
            ref_source (bool): If `true`, use the external reference
                source instead of a DAC channel.
        """
        self.ref_source=ref_source
        return self

    def get_readout_channel(self):
        """Get the readout channel for this stimulation unit
        """
        return next(sCh for sCh in self.mapping if sCh[0]==int(self.unit_no))[1]



class Offset:
    def __init__(self):
        pass


class Array:
    """Control the electrode array of the MaxOne chip.

    Args:
        token (str):
            Token id to identify this instance on the server. Default is
            'online'
        persistent (bool):
            Should the instance of this array get deleted on the server
            when 'close' gets called

    Similar as for :class:`maxlab.Sequence`, for each *Array* object,
    there is a counter part on the server, which can be uniquely
    identified through the token id parameter. However, different than
    for :class:`maxlab.Sequence`, creating such Array objects on the
    server is computationally expensive. It takes multiple seconds and
    consumes significant amounts of RAM.

    This is why we can control whether we want to use one of the already
    available array objects, or whether we want to create a new array
    object on the server. The names of the available objects are
    'online' and 'offline'. The online array object is the same as
    which is controlled through MaxLab Live. I.e. when downloading a
    configuration from the scope GUI, the 'online' array will have the
    same configuration set. If it is required to manipulate the array
    without affecting the online version, it is recommended to use the
    'offline' version. By specifying a different name, arbitrary amounts
    of arrays can be controlled at the same time.

    All operations on the array, such as loading a configuration
    or connecting amplifiers to stimulation units are all done in
    software and no changes are applied to the hardware. Only once the
    :func:`Array.download` function is executed, the state of the array
    gets downloaded to the actual electrode array on the chip.

    When creating configurations through the :func:`Array.route` method,
    all other configurations, such as connecting stimulation channels etc.
    will be lost. So the best procedure is to:

     1. Select electrodes
     2. Route
     3. Connect stimulation channels
     4. Save configuration

    >>> array = maxlab.chip.Array('online')
    >>> array.reset()
    >>> array.select_electrodes( [11110, 11111, 11112] )
    >>> array.select_stimulation_electrodes( [11330, 11331, 11332] )
    >>> array.route()
    >>> array.connect_electrode_to_stimulation( 11330 )
    >>> array.connect_electrode_to_stimulation( 11331 )
    >>> array.connect_electrode_to_stimulation( 11332 )
    >>> array.save_config("myConfig.cfg")
    >>>
    >>> stim1 = array.query_stimulation_at_electrode( 11330 )
    >>> stim2 = array.query_stimulation_at_electrode( 11331 )
    >>> stim3 = array.query_stimulation_at_electrode( 11332 )
    """
    def __init__(self, token='online', persistent=False):
        self.token = token
        self.persistent = persistent
        with comm.api_context() as api:
            api.send("mea_array_new " + self.token)

    def close(self):
        """Close this `Array` object.

        If `persistent` is not set, delete the array with the same token
        from the server.
        """
        if not self.persistent:
            with comm.api_context() as api:
                api.send("mea_array_delete " + self.token)

    def shutdown(self):
        self.close()

    def __send(self, command):
        with comm.api_context() as api:
            return api.send("mea_array_command " + self.token + "\n" + command)

    def reset(self):
        """Resets the array into a defined state.

        This function disconnects all electrodes and all amplifiers.
        """
        return self.__send("mea_array_reset")

    # Electrode selection & routing functions
    def select_stimulation_electrodes(self, electrodes):
        """Select the given electrodes for stimulation.

        Args:
            electrodes (list of int): List of stimulation electrodes.

        The selected electrodes get automatically a high priority when
        routed.

        Note:
            Make sure not to select more than a few hundreds of these
            electrodes.  Otherwise, routing will not work (converge) well.
        """
        return self.select_electrodes(electrodes, 1000)

    def select_electrodes(self, electrodes, weight=1):
        """Select the given electrodes for routing.

        Args:
            electrodes (list of int): List of recording electrodes.

        By passing a weight parameter, the routing priority for the
        electrodes can be adjusted.  The higher the weight, the higher
        the routing priority during routing.

        Only one weight can be set for the electrode ids in the function
        argument. Usually, different priorities need to be assigned for
        this to make sense at all. To achieve this, call this function
        for each set of priorities. See below for an example:

        >>> array = maxlab.chip.Array()
        >>> array.select_electrodes([1,2,3,4], 10) # electrodes with priority of '10'
        >>> array.select_electrodes([5,6,7,8], 15) # other electrodes with a higher priority
        """
        return self.__send("mea_array_select_electrodes " + " ".join([str(el)+"/"+str(weight) for el in electrodes]) )

    def clear_selected_electrodes(self):
        """Clear all selected electrodes.
        """
        return self.__send("mea_array_clear_selected_electrodes")

    def route(self):
        """Route electrode configuration.

        Note:
            Be aware that any manual switch settings, such as stimulation
            connections are lost after the routing. This is because
            routing starts from an 'empty' array.
        """
        return self.__send("mea_array_route")

    def load_config_data(self, config_file_data):
        return self.__send("mea_array_set_config " + config_file_data)

    def load_config(self, config_file_name):
        """Loads an electrode configuration from disk and sends it to the API.

        Args:
            config_file_name (str): File name with the config. Usually
                ends with ".cfg".
        """
        with open(config_file_name, 'r') as config_file:
            config=config_file.read()
            return self.__send("mea_array_set_config " + config)
        return -1

    def save_config(self, config_file_name):
        """Save current array configuration to disk.

        Args:
            config_file_name (str): File name where to save the
                configuration to. Should end with ".cfg"
        """
        with open(config_file_name, 'w') as config_file:
            config = self.__send("mea_array_get_config")
            config_file.write(config)
        return 0

    def get_config(self):
        return maxlab.config.Config(self.__send("mea_array_get_config"))

    def download(self, wells=None):
        """Download the electrode configuration to the chip."""
        if wells is None:
            self.__send("mea_array_download")
        else:
            self.__send("mea_array_download_wells " + ','.join([str(w) for w in wells]))

    def set(self):
        return "mea_array_download"

    # Connect functions
    def connect_amplifier_to_stimulation(self, amplifier_channel): # chip.readoutToStim( self.getReadoutCh(stimBuffer) )
        """Connect the amplifier to a stimulation unit.

        Args:
            amplifier_channel (int): Amplifier to be connected to a
                stimulation unit.
        """
        return self.__send("mea_array_connect_amplifier_to_stimulation " + str(amplifier_channel))

    def connect_amplifier_to_ringnode(self, amplifier_channel):    # chip.readoutToVRef( self.getReadoutCh(stimBuffer) )
        """Connect amplifier to the common ring node.

        Args:
            amplifier_channel (int): Amplifier to be connected to the
                ring node.

        Note:
            By enabling :func:`Core.enable_external_port`, the ring node
            can be connected to the external port.
        """
        return self.__send("mea_array_connect_amplifier_to_ringnode " + str(amplifier_channel))

    def connect_electrode_to_stimulation(self, electrode_no):
        """Connect electrode to stimulation unit.

        Args:
            electrode_no (int): Electrode ID to be connected to a
                stimulation unit.

        For this function to work, the selected electrode ID needs
        already be routed to an amplifier.
        """
        return self.__send("mea_array_connect_electrode_to_stimulation " + str(electrode_no))

    # Connect electrode to direct stimulation. Do not use this
    def connect_electrode_to_direct_stimulation(self, electrode_no):
        return self.__send("mea_array_connect_electrode_to_direct_stimulation " + str(electrode_no))

    def connect_all_floating_amplifiers(self):
        """Connect all floating amplifiers to a common node.

        Note:
            By enabling :func:`Core.enable_external_port`, the floating
                amplifiers can be connected to the external port
        """
        return self.__send("mea_array_connect_all_floating_amplifiers")

    # Query functions
    def query_amplifier_at_stimulation(self, stimulation_channel):
        """Query amplifier at stimulation unt.

        Args:
            stimulation_channel (int): Which stimulation unit to query.
        """
        return self.__send("mea_array_query_amplifier_at_stimulation " + str(stimulation_channel))

    def query_stimulation_at_amplifier(self, amplifier_channel):
        """Query stimulation unit at amplifier.

        Args:
            amplifier_channel (int): Which amplifier channel to query.
        """
        return self.__send("mea_array_query_stimulation_at_amplifier " + str(amplifier_channel))

    def query_amplifier_at_electrode(self, electrode_no):
        """Query amplifier at the electrode.

        Args:
            electrode_no (int): Which electrode ID to query.
        """
        return self.__send("mea_array_query_amplifier_at_electrode " + str(electrode_no))

    def query_amplifier_at_ringnode(self):
        """Query which amplifiers are connected to the ring node."""
        return self.__send("mea_array_query_amplifier_at_ringnode")

    def query_stimulation_at_electrode(self, electrode_no):
        """Query stimulation unit at the electrode.

        Args:
            electrode_no (int): Which electrode ID to query.
        """
        return self.__send("mea_array_query_stimulation_at_electrode " + str(electrode_no))

    # Disconnect functions
    def disconnect_amplifier_from_stimulation(self, amplifier_channel):
        """Disconnect amplifier from stimulation.

        Args:
            amplifier_channel (int): Which amplifier channel to disconnect
                from a stimulation unit.
        """
        return self.__send("mea_array_disconnect_amplifier_from_stimulation " + str(amplifier_channel))

    def disconnect_electrode_from_stimulation(self, electrode):
        """Disconnect electrode from stimulation.

        Args:
            electrode (int): Which electrode ID to disconnect from a
                stimulation unit.
        """
        return self.__send("mea_array_disconnect_electrode_from_stimulation " + str(electrode))

    def disconnect_amplifier_from_ringnode(self, amplifier_channel):
        """Disconnect amplifier from ring node.

        Args:
            amplifier_channel (int): Which amplifier channel to disconnect
                from the ring node.
        """
        return self.__send("mea_array_disconnect_amplifier_from_ringnode " + str(amplifier_channel))

    # Aux functions
    def connect_all(self):
        return self.__send("mea_array_connect_all")

    def connect_electrode(self, electrode_no):
        return self.__send("mea_array_connect_electrode " + str(electrode_no))

    def disconnect_electrode(self, electrode_no):
        return self.__send("mea_array_disconnect_electrode " + str(electrode_no))

