#!/usr/bin/python
# -*- coding: utf-8

import struct
import maxlab.apicomm as comm
import maxlab.pycompat as pycompat

class Stream(object):

    @staticmethod
    def _deserialize(buffer):
        values = {}
        for line in buffer.split('\n##\n'):
            if ':' in line:
                well_str = line.split(':')[0]
                well = int(well_str)
                ret = line[len(well_str)+1:]
                size_str = ret.split(',')[0]
                size = int(size_str)
                buf = pycompat.encode(ret[len(size_str)+1:])
                fmt = '=%df' % size
                values[well] = list(struct.unpack(fmt, buf))
        return values

    @staticmethod
    def start_demodulate(freq):
        with comm.api_context() as api:
            api.send("stream_set_demodulate " + str(freq))

    @staticmethod
    def stop_demodulate():
        with comm.api_context() as api:
            api.send("stream_set_demodulate 0")

    @staticmethod
    def get_amplitudes(wells):
        with comm.api_context() as api:
            return Stream._deserialize(api.send("stream_get_amplitudes " + ','.join([str(w) for w in wells])))

    @staticmethod
    def get_mean():
        with comm.api_context() as api:
            return Stream._deserialize(api.send("system_mean "))
