#!/usr/bin/python
# -*- coding: utf-8

import maxlab.apicomm as comm


class Saving:
    # COMPAT
    def start(self, file_name=""):
        with comm.api_context() as api:
            api.send("saving_start " + file_name)

    # COMPAT
    def start_spikes_only(self, file_name=""):
        with comm.api_context() as api:
            api.send("saving_spikes_only " + file_name)

    # COMPAT
    def stop(self):
        with comm.api_context() as api:
            api.send("saving_stop")

    def set_legacy_format(self, use=True):
        with comm.api_context() as api:
            api.send("saving_use_old_file_format " + str(int(use)))

    def start_file(self, file_name):
        with comm.api_context() as api:
            api.send("saving_start_file " + file_name)

    def stop_file(self):
        with comm.api_context() as api:
            api.send("saving_stop_file")

    def start_recording(self, wells=None):
        if wells is None:
            wells = [0]

        with comm.api_context() as api:
            api.send("saving_set_recording_wells " + ','.join([str(w) for w in wells]))
            api.send("saving_start_recording")

    def stop_recording(self):
        with comm.api_context() as api:
            api.send("saving_stop_recording")

    def open_directory(self, directory_name):
        with comm.api_context() as api:
            api.send("saving_open_dir " + directory_name)

    def record_wells(self, wells):
        """Set the wells that should be recorded in the files
        """
        with comm.api_context() as api:
            api.send("saving_set_recording_wells " + ','.join([str(w) for w in wells]))

    def group_define(self, well, name, channels=[]):
        channels = [c for c in channels if 0 <= c < 1024]
        with comm.api_context() as api:
            return api.send("saving_group_define %d %s %s" % (well, name, ','.join([str(c) for c in channels])))

    def group_set_trigger(self, well, name, channels, pre=100, post=100, min_amp=-10000, max_amp=10000):
        channels = [c for c in channels if 0 <= c < 1024]
        with comm.api_context() as api:
            return api.send("saving_group_trigger %d %s %s %d %d %f %f" % (well, name, ','.join([str(c) for c in channels]), pre, post, min_amp, max_amp))

    def group_clear_trigger(self, well, name):
        with comm.api_context() as api:
            return api.send("saving_group_trigger %d %s clear" % (well, name))

    def group_delete(self, well, name):
        with comm.api_context() as api:
            return api.send("saving_group_delete %d %s" % (well, name))

    def group_delete_well(self, well):
        with comm.api_context() as api:
            return api.send("saving_group_clear %d" % well)

    def group_delete_all(self):
        with comm.api_context() as api:
            return api.send("saving_group_clear")

    def group_info(self, well):
        with comm.api_context() as api:
            return api.send("saving_group_info %d" % well)

    def write_assay_property(self, key, value):
        with comm.api_context() as api:
            return api.send("saving_set_assay_property " + key + " " + value)

    def write_assay_input(self, key, value):
        with comm.api_context() as api:
            return api.send("saving_set_assay_input " + key + " " + value)
