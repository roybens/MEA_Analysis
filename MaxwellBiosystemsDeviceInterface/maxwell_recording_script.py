#!/usr/bin/python

import maxlab
import maxlab.system
import maxlab.chip
import maxlab.util
import maxlab.saving

import time

## How to load configurations and run recordings
#
#
# 0. Initialize system into a defined state
#
#    It is best to initialize the system into a defined state beforefile:///usr/share/doc/HTML/index.html
#    Starting any script. This way, one can be sure that the system
#    is always in the same state while running the script, regardless
#    of what has been done before with it.
#
#
# 1. Load your electrode selection
#
#    We use an existing config file to route the electrodes of interest.
#
# 2. Start a recording
#
#    Once we are ready to start an experiment e.g. a stimulation protocol
#    we can start a recording.
#
# 3. Stop and store the recording
#
#    When our experiment is done, we can stop the recording and store the data
#    at a location of our choice.
#
#

######################################################################
# User Input
######################################################################

# Input the path to the config file you want to use
name_of_configuration = ''
# Input the path of the directory where the recording and sequence log should be stored
data_path = ''
# Input the name the recording file should have
recording_file_name = ''
# Input your choice of data format (True for legacy format)
use_legacy_write = True
# Input wheter you want to record spikes only, or signals as well
record_only_spikes = False
# Input how many wells you want to record from (range(1) for MaxOne, range(6) for MaxTwo)
wells = range(6)


# 0. Initialize system into a defined state
maxlab.util.initialize()

# If we want to stimulate during our experiment, we also issue the following command:
#
#    maxlab.send(maxlab.chip.Core().enable_stimulation_power(True))
#
# This is needed, because by default the stimulation units are
# powered off globally...


######################################################################
# 1. Load your electrode selection
######################################################################

array = maxlab.chip.Array('online')
array.load_config(name_of_configuration)

# After loading the configuration, we could also add more recording or stimulation
# electrodes. Once everything is setup according to our needs, we need to download
# the prepared array configuration to the chip
array.download()
maxlab.util.offset()


######################################################################
# 2. Start a recording
######################################################################
# Now that the chip is ready, we start the recording and perform our experiment.
s = maxlab.saving.Saving()
s.open_directory(data_path)
s.set_legacy_format(use_legacy_write)
# If we use the new data format and also store the signal traces, we must declare which electrodes
# we want to store data from. This can be set through the group_define function which has the following form:
# s.group_define(well_nr, "any_name", list_of_channels_to_record_from)
# e.g s.group_define(0, "all_channels", range(1024)) will store data from all 1024 channels
# for the first well and store in the group 'all_channels'.
# In case you need, you can also define multiple groups for same recording, with a set of channels.
# These groups can also have overlapping channels. For example if you stimulate on channel 200:
# s.group_define(0, "stim_channel", [200])
# s.group_define(0, "all", range(1024))
s.group_delete_all()
if not record_only_spikes:
	for well in wells:
		s.group_define(well, "routed")

s.start_file(recording_file_name)
s.start_recording(wells)
# Sleep for the amount of recording time + buffer time during which the experiment will run
time.sleep(2)


######################################################################
# 3. Stop and store the recording
######################################################################
# After our experiment is over, we just need to stop the recording and store everything.
# The recording should now be available  at data_path/recording_file_name and can be loaded
# into MaxLab Live.
s.stop_recording()
s.stop_file()
s.group_delete_all()
