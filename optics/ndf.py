# ndf.py
# A script for plotting ndf parameters, 
# using eccosorb MF-110, like CR-110 from Halpern and gush materials artile.

# CR110 materials properties from
# Applied Optics, Vol. 25, Issue 4, pp. 565-570 (1986)
# http://dx.doi.org/10.1364/AO.25.000565

import numpy as np
import pylab
pylab.ion()


a = 0.3   # for eccosorb CR110
b = 1.2   # for eccosorb CR110
d_in = np.arange(0.125,1.0,0.001)
d = d_in*2.54  # cm

transmission = {}
for nu in [90,150,220,280]:
    icm = nu/30 #
    alpha = a*icm**b   
    freqstr = str(nu)
    transmission[freqstr] = np.exp(-alpha*d) 


pylab.figure(1)
pylab.clf()
pylab.semilogy(d_in,transmission['90'],label='90 GHz')
pylab.semilogy(d_in,transmission['150'],label='150 GHz')
pylab.semilogy(d_in,transmission['220'],label='220 GHz')
pylab.semilogy(d_in,transmission['280'],label='280 GHz')
pylab.semilogy(d_in,0*d_in + 0.05,'-k')
pylab.legend()
pylab.grid()
pylab.xlabel('CR110 thickness in inches')
pylab.ylabel('Transmission')
pylab.title('Transmission at 4K, CR110 (Halpern,Gush,Wishnow and De Cosmo)')

pylab.show()



