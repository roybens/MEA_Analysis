#!/usr/bin/python
# -*- coding: utf-8

import math
import time
from collections import defaultdict

import maxlab.apicomm as comm
import maxlab.stream as stream

def initialize(wells=None):
    with comm.api_context() as api:
        if wells is None:
            api.send("system_initialize")
        else:
            api.send("system_initialize_wells " + ','.join([str(w) for w in wells]))

def activate(wells):
    with comm.api_context() as api:
        api.send("system_set_activated_wells " + ','.join([str(w) for w in wells]))

def offset():
    with comm.api_context() as api:
        api.send("system_offset")
        time.sleep(15)

def hpf(cutoff):
    if (cutoff=="1Hz"):
        cutoff_value=4095
    elif (cutoff=="300Hz"):
        cutoff_value=1100
    else:
        raise ValueError("Not a valid high-pass parameter: " + str(cutoff) + "\n" + \
                         "Allowed values are: '1Hz' and '300Hz'" )
    with comm.api_context() as api:
        api.send("system_hpf " + str(cutoff))

def set_gain(gain):
    with comm.api_context() as api:
        api.send("system_gain " + str(gain))

def get_mean():
    return stream.Stream.get_mean()

def percentile(values, percent):
    """Find the percentile of a list of values.

    @parameter values - is a list of values. Note values MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.

    @return - the percentile of the values
    """
    if not values:
        return None
    k = (len(values)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return values[int(k)]
    d0 = values[int(f)] * (c-k)
    d1 = values[int(c)] * (k-f)
    return d0+d1

def group_wells_per_bank(wells):
    """
    Group the given list of wells by banks.
    @return - A list of lists of wells
    """
    with comm.api_context() as api:
        mapping = api.send("wellplate_query_well_to_bank_mapping").strip(',').split(',', )
        if len(mapping) % 2 != 0:
            return []

        result = defaultdict(list)
        for bank, well in zip(mapping[1::2], mapping[0::2]):
            if int(well) in wells:
                result[int(bank)].append(int(well))
            else:
                continue
        return list(result.values())

def get_no_of_banks(wells):
    """
    Get the number of banks occupied by the given list of wells.
    """
    with comm.api_context() as api:
        mapping = api.send("wellplate_query_well_to_bank_mapping").strip(',').split(',', )
        if len(mapping) % 2 != 0:
            return -1

        result = defaultdict(list)
        for bank, well in zip(mapping[1::2], mapping[0::2]):
            if int(well) in wells:
                result[int(bank)].append(int(well))
            else:
                continue
        return len(result.keys())

def set_primary_well(well):
    with comm.api_context() as api:
        api.send("system_set_primary_well " + str(well))

def init_wells_for_mx2(wells, assay, array, total_time):
    """
    Initialize the given wells.
    """

    def _progress(seconds):
        return int(seconds * 100.0 / total_time)

    bank_grouping = group_wells_per_bank(wells)
    no_of_banks = len(bank_grouping)
    with comm.api_context() as api:
        api.send("stream_blanking_on " + str((Timing.waitInit + Timing.waitAfterOffset) * no_of_banks + 1))
    for grouped_wells in bank_grouping:
        initialize(grouped_wells)
        assay.sleep(Timing.waitInit, _progress(Timing.waitInit))

        array.connect_all()
        array.download()
        assay.sleep(1)
        offset()
        array.reset()
        assay.progress = assay.progress + _progress(Timing.waitInOffset)
        assay.sleep(Timing.waitAfterOffset, _progress(Timing.waitAfterOffset))
    with comm.api_context() as api:
        api.send("stream_blanking_off")

class Timing:
    waitInit = 2
    waitAfterDownload = 5
    waitAfterOffset = 5
    waitInOffset = 15
    waitAfterBankSwitch = 2
    waitAfterRecording = 2