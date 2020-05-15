import numpy as np
import optical_calcs

import DetectorDefs
import optical_calcs
import tes_calcs
import print_optics
import print_totals


detector, optics = DetectorDefs.generic150()

detector = optical_calcs.optical_calcs(detector, optics)
detector['DO_FREQ_RESPONSE'] = False
detector = tes_calcs.tes_calcs(detector)

print_optics.print_optics(optics)
print_totals.print_totals(detector)

