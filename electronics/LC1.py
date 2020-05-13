# Plot LC comb

import numpy as np
import pylab

pylab.ion()


fcenter = 6e6  # 3MHz
f = np.arange(fcenter-40000., fcenter+2000, 10.)
foffset = f-fcenter
w = 2*np.pi*f
w0 = fcenter*2*np.pi

L = 60e-6
C = 1/(L*w0**2)

pylab.figure(1)
pylab.clf()

for R in [0.0, 0.1, 0.2, 0.3, 0.4] :
    Z = R + 1j*w*L + 1./(1j*w*C)
    Zinv =  1/Z
    pylab.subplot(211)
    pylab.plot(foffset,np.abs(Z))
    pylab.subplot(212)
    Rstr = 'Rp = ' + str(R)
    pylab.plot(foffset,(180./np.pi)*np.angle(Z),label=Rstr)

pylab.subplot(211)
pylab.grid()
pylab.ylabel('abs(Z)')
pylab.subplot(212)
pylab.legend(loc=2)
pylab.ylabel('Phase (deg)')
pylab.grid()
pylab.xlabel('Frequency offset')


pylab.show()




