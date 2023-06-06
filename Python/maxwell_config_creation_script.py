import maxlab
import maxlab.chip
import sys
import maxlab.util
import maxlab.characterize
import json
import time

CONFIG_FILENAME = "electrodes4min.cfg"


def main():

    #read the location of the electrodes found out by 
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

    #read the electrode indices file
    with open('electrode_locations_4min.json',) as fileptr:
        electrode_data = json.load(fileptr)
    
    #electode data is in terms of indices convert it to electrodes

    selected_electrodes = [220 * locations[1] + locations[0] for locations in  electrode_data.values()]  # x is location[0] and y the other

    array = maxlab.chip.Array('online')
    array.reset()
    array.clear_selected_electrodes( )
    array.select_electrodes(selected_electrodes)
    array.route()
    array.save_config(CONFIG_FILENAME)




if __name__ == "__main__":

    main()