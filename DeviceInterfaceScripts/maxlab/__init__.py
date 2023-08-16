#!/usr/bin/python
# -*- coding: utf-8
"""
Python module for programmatic control of MaxLab Live.

.. moduleauthor:: Jan Mueller <jan.mueller@mxwbio.com>

"""

import maxlab.system
import maxlab.chip
import maxlab.apicomm as comm

#from sequence import Sequence
#from Loop import Loop
#query_DAC_lsb

def send(obj):
    """Utility function to send command objects to the server API.

    Args:
        obj (command object): a command object from either
            :class:`maxlab.system`
            or :class:`maxlab.chip`. The command object needs to have a
            method called `send()`, which contains all the information
            needed to program the API.

    """
    with comm.api_context() as api:
        return api.send(obj.set())

def send_raw(msg):
    with comm.api_context() as api:
        return api.send(msg)
    
def error():
    return send_raw("get_errors")

def interrupt_routing():
    """Utility function to interrupt long running routing jobs in
    the server.

    When a routing job or a sequence of routing jobs are ongoing, the
    user can choose to prematurely terminate the routing.
    """
    return send_raw("system_interrupt_routing")


def query_DAC_lsb_mV():
    """Query LSB value of the DAC channels.

    Note:
        This value is only valid when doing voltage stimulation.
    """
    return send_raw("system_query_dac_lsb_mv")


class Sequence:
    """Holds a sequence of commands.

    Args:
        initial_delay (samples) : how many samples to delay the first
            command.

    The commands contained in a sequence are executed immediately after
    each other.

    Thus, the timing between consecutive commands can be precisely
    controlled, down to a resolution of 50 µs. The jitter between
    consecutive commands is much less than than.  Particularly for
    sequences implementing electrical stimulation pulses, this is
    relevant. For a more detailed discussion about that, please see the
    introduction of this document.

    It can happen, that the first command of the sequence gets executed
    before all remaining commands were successfully downloaded. In that
    case the timing would not be precise anymore. To prevent this from
    happening, it is good practise to insert a delay command right at
    the beginning of the sequence. This is done through  the parameter
    `initial_delay`. Delaying for 100 samples (i.e. 5 ms) should be
    enough in virtually all cases.

    The sequences are managed on the server. As many sequences as needed
    can be generated. They will be maintained on the server. Whenever
    a sequence object in the python code gets deleted, either through
    explicitely calling del or if the code goes out of scope, the
    respective sequence in the server will also be removed:

    >>> s1 = maxlab.Sequence()
    >>> s1.append( maxlab.system.DelaySamples(4) )
    >>> del s1 # All information about Sequence s1 is removed from the server

    >>> def func(samples):
    >>>     s = maxlab.Sequence()
    >>>     s.append(maxlab.system.DelaySamples(samples))
    >>>     return 0 # As soon as the function finishes, the sequence s is removed from the server


    """
    def __init__(self, token=None, initial_delay=100, persistent=False):
        """Create a Sequence object.

        Args:
            initial_delay (int): Delay for this samples. This allows
                this sequence to buffer in memory before the first command
                is executed.
        """
        self.persistent = persistent
        
        with comm.api_context() as api:
            # request a new sequence and save its token
            if token is None:
                self.token = api.send("sequence_new")
            else:
                self.token = api.send("sequence_new " + token)

        if initial_delay>0:
            self.append(maxlab.system.DelaySamples(initial_delay))
            
    def shutdown(self):
        with comm.api_context() as api:
            if not self.persistent:
                api.send("sequence_delete " + self.token)

    def __del__(self):
        self.shutdown()

    def reset(self):
        """Resets the sequence.

        Calling this function removes all previously appended commands
        from the sequence.
        """
        with comm.api_context() as api:
            api.send("sequence_reset " + str(self.token) )
            return self

    def append(self, obj):
        """Append the command object `obj` to the sequence.

        Args:
            obj (command object): The command object to be appended to
                the sequence.
        """
        with comm.api_context() as api:
            api.send("sequence_append " + self.token + "\n" + obj.set())
            return self

    def send(self):
        """Sends the sequence to the system.

        This function downloads the sequence to the MaxHub, from where
        it gets executed, command by command.
        """
        with comm.api_context() as api:
            api.send("sequence_send " + self.token)
            return self


class Loop:
    @staticmethod
    def prepare():
        with comm.api_context() as api:
            api.send("system_loop_prepare")

    @staticmethod
    def finish():
        with comm.api_context() as api:
            api.send("system_loop_finish")

    @staticmethod
    def download():
        with comm.api_context() as api:
            api.send("system_loop_download")

    @staticmethod
    def start():
        with comm.api_context() as api:
            api.send("system_loop_start")

    @staticmethod
    def stop():
        with comm.api_context() as api:
            api.send("system_loop_stop")

    @staticmethod
    def run_once():
        with comm.api_context() as api:
            api.send("system_loop_run_once")

    @staticmethod
    def append_delay(samples):
        with comm.api_context() as api:
            api.send("system_loop_append_delay " + str(samples))

    @staticmethod
    def append_dac(code, value):
        with comm.api_context() as api:
            api.send("system_loop_append_dac " + str(code) + " " + str(value))

    @staticmethod
    def append_event(well_id, event_type, user_id, properties):
        with comm.api_context() as api:
            api.send("system_loop_event " + str(well_id) + " " + str(event_type) + " " + str(user_id) + " " + properties)
