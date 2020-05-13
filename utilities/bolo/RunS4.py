import numpy as np
import sys
import optical_calcs

import DetectorDefs
import optical_calcs
import tes_calcs
import print_optics
import print_totals

''' RunS4.py 
Usage:   RunS4 system_name
where system name must be defined in DetectorDefs.py
'''


system_name = sys.argv[1]
print(system_name)

#

#detector, optics = DetectorDefs.S4_SAT_30()
detector, optics = getattr(DetectorDefs,system_name)()
detector, optics = optical_calcs.optical_calcs(detector, optics)
detector['DO_FREQ_RESPONSE'] = False
detector = tes_calcs.tes_calcs(detector)

print_optics.print_optics(optics)
print_totals.print_totals(detector)

