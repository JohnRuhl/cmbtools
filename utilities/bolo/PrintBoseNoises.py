import matplotlib.pyplot as plt
import numpy as np

import optical_calcs
import DetectorDefs
import tes_calcs


# 
detector, optics = DetectorDefs.S4_SAT_155()
detector['BoseFactor'] = 1.0
detector, optics = optical_calcs.optical_calcs(detector, optics)



