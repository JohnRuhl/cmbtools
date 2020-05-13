# ndf.py
# A script for plotting ndf parameters, 
# using eccosorb MF-110, like CR-110 from Halpern and gush materials artile.

# CR110 materials properties from
# Applied Optics, Vol. 25, Issue 4, pp. 565-570 (1986)
# http://dx.doi.org/10.1364/AO.25.000565

import numpy as np
import pylab

k = 1.38e-23
dT = 300-77.

a = 0.3   # for eccosorb CR110
b = 1.2   # for eccosorb CR110
d_in = 1.0
d = d_in*2.54  # cm

bands = dict()
bands['90band'] = np.arange(75,105,1.0)
bands['150band'] = np.arange(120,160,1.0)
bands['220band'] = np.arange(200,250,1.0)


for obsband in ['90band','150band','220band']:
    icm = bands[obsband]/30.
    alpha = a*icm**b
    avg_transmission = np.mean(np.exp(-alpha*d))
    bandwidth = bands[obsband][-1] - bands[obsband][0]
    expected_power = k*dT*bandwidth*1e9*avg_transmission
    print obsband, bandwidth, avg_transmission, expected_power

