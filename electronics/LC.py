# Plot LC comb

import numpy as np
import pylab

pylab.ion()

f = np.arange(2e6,2.3e6,101.)
w = 2*np.pi*f

fcenters = np.arange(2e6,2.3e6,16.e3)
w0 = fcenters*2*np.pi

R =  3*0.025
L = 3e-6
Carr = 1/(L*w0**2)


Zinv = np.zeros(len(w))
for ii in range(len(Carr)):
    Cx = Carr[ii]
    Z = R + 1j*w*L + 1./(1j*w*Cx)
    Zinv = Zinv + 1/Z

Ztot = 1/Zinv
I = 1/Ztot

pylab.figure(1)
pylab.clf()

pylab.plot(f,I)
pylab.grid()


pylab.show()




