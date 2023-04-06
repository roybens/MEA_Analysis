import maxlab
import sys
import maxlab.util
import maxlab.chip
import maxlab.saving
import maxlab.characterize
import time
from datetime import datetime

CHIP_ID = 16848

def main():

    unit = 16
    if len(sys.argv) > 1:
        unit = int(sys.argv[1])

    # Initalize
    maxlab.util.initialize()

    # Set the gain to 1
    amplifier = maxlab.chip.Amplifier().set_gain(512)
    maxlab.send(amplifier)

    characterizer = maxlab.characterize.StimulationUnitCharacterizer()

    print("\n************ Characterizing unit: %d *************\n" % unit) 
    dac_code = characterizer.characterize(unit)
    print("\n0V DAC CODE = %d\n" % dac_code)

    #select and record block by block.
    #going with 35 blocks of size to (44 x 20) to cover the whole chip
    XMULTIPLIER = 44
    YMULTIPLIER = 20

    block_number = 0
    for y in range(6):  #TODO : change the hardcoding
        for x in range(5):
            block_number += 1
            #print("(",XMULTIPLIER*x,YMULTIPLIER*y,XMULTIPLIER*x+XMULTIPLIER-1,YMULTIPLIER*y+YMULTIPLIER-1,")")
            selected_electrodes = maxlab.chip.electrode_rectangle_indices(XMULTIPLIER*x,YMULTIPLIER*y,XMULTIPLIER*x+XMULTIPLIER-1,YMULTIPLIER*y+YMULTIPLIER-1)
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
            now = datetime.datetime.now()
            strfmt = now.strftime("%Y%m%d_%H_%M_%S")
            filename = "recording_chip_id_"+CHIP_ID+"_block_id_"+block_number+"_"+strfmt+".raw.h5"
            saving_object = maxlab.saving.Saving()
            try :
                saving_object.start(filename)
                print("recording start")

                #sleep for recording time
                time.sleep(60)
                saving_object.stop()
            except Exception as e:
                print(e)
                sys.exit(0)



            

if __name__ == "__main__":

    main()



