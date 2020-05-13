# Plot amp and phase for varying L.

import numpy as np
import pylab

pylab.ion()


L = 60e-6
C = 1/(L*w0**2)
w0 = 1/np.sqrt(L*C)
f0 = w0/(2*np.pi)

R = 0   # parasitic resistance




f = np.arange(fcenter-100., fcenter+100, 10.)
foffset = f-fcenter
w = 2*np.pi*f
w0 = fcenter*2*np.pi


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




