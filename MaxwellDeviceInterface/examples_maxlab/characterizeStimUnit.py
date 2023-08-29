import maxlab
import maxlab.chip
import maxlab.util
import maxlab.characterize

import sys

unit = 16
if len(sys.argv) > 1:
    unit = int(sys.argv[1])

# Initalize
maxlab.util.initialize()

# Set the gain to 1
amplifier = maxlab.chip.Amplifier().set_gain(1)
maxlab.send(amplifier)

characterizer = maxlab.characterize.StimulationUnitCharacterizer()

print("\n************ Characterizing unit: %d *************\n" % unit) 
dac_code = characterizer.characterize(unit)
print("\n0V DAC CODE = %d\n" % dac_code)
