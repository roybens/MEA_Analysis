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

CHIP_ID = 16848

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
		#amplifier = maxlab.chip.Amplifier().set_gain(512)
		#maxlab.send(amplifier)
		#maxlab.util.set_gain(512)
		#maxlab.util.hpf("1Hz")

	#	characterizer = maxlab.characterize.StimulationUnitCharacterizer()

		#print("\n************ Characterizing unit: %d *************\n" % unit) 
		#dac_code = characterizer.characterize(unit)
		#print("\n0V DAC CODE = %d\n" % dac_code)

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
					
					#
					#saving_object.open_directory("~/recordings/")
					#saving_object.start_recording()
					#saving_object.start(filename)
					saving_object.start_recording()
					saving_object.start_file(filename)
					
					print("recording start")

					#sleep for recording time
					time.sleep(60)
					
					saving_object.stop_file()
					saving_object.stop_recording()
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
