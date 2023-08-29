import maxlab
import sys
import maxlab.util
import maxlab.chip
import maxlab.saving
import maxlab.characterize
import maxlab.config
import time
from datetime import datetime
import maxlab.apicomm
import os


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

###############################################################################
# User Input 
###############################################################################

CHIP_ID = 16848
#Path of the directory where the recording and sequence log should be stored

data_path = '~/recording/'

# Input your choice of data format (True for legacy format)
use_legacy_write = True
# Input wheter you want to record spikes only, or signals as well
record_only_spikes = False
# Input how many wells you want to record from (range(1) for MaxOne, range(6) for MaxTwo)
wells = range(1)


def main():
	try:
		print("Script begin")
		#unit = 16
		#if len(sys.argv) > 1:
		    #unit = int(sys.argv[1])

		# Initalize
		maxlab.util.initialize()
		print("Init complete")

		# Set the gain to 1

		#maxlab.util.set_gain(512)
		#maxlab.util.hpf("1Hz")

		#select and record block by block.
		#going with 35 blocks of size to (44 x 20) to cover the whole chip
		XMULTIPLIER = 44
		YMULTIPLIER = 20

		block_number = 0
		for y in range(6):  #TODO : change the hardcoding
			for x in range(5):
				block_number += 1
				print(f"For the block {block_number}")
				#print("(",XMULTIPLIER*x,YMULTIPLIER*y,XMULTIPLIER*x+XMULTIPLIER-1,YMULTIPLIER*y+YMULTIPLIER-1,")")
				selected_electrodes = maxlab.config.electrode_rectangle_indices(XMULTIPLIER*x,YMULTIPLIER*y,XMULTIPLIER*x+XMULTIPLIER-1,YMULTIPLIER*y+YMULTIPLIER-1)
				array = maxlab.chip.Array('online')
				array.reset()
				array.clear_selected_electrodes( )
				array.select_electrodes(selected_electrodes)
				array.route()
				array.download()
				time.sleep(2)
				# Run the offset compensation
				maxlab.util.offset()
				time.sleep(5)

				#creating a recording object
				now = datetime.now()
				strfmt = now.strftime("%Y%m%d_%H_%M_%S")
				filename = "recording_chip_id_"+str(CHIP_ID)+"_block_id_"+str(block_number)+"_"+strfmt
				saving_object = maxlab.saving.Saving()
				try:
					
					saving_object.group_delete_all()
					saving_object.open_directory("~/recordings/")
					#saving_object.start_recording()
					#saving_object.start(filename)
					#
					saving_object.start_file(filename)
					saving_object.start_recording()
					
					print("recording start")

					#sleep for recording time
					time.sleep(60)
					
					saving_object.stop_file()
					#saving_object.stop_recording()
					#saving_object.stop_recording()
					#saving_object.stop()
					#
					print("recording stop")
				except Exception as e:
					print(e)
					sys.exit(0)
		print("Script end")

	except Exception as e:
		print(e)
            

if __name__ == "__main__":

    main()
