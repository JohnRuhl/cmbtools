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
nu = np.arange(50.,250.,1.0)
icm = nu/30. #
alpha = a*icm**b   
#d_in = np.arange(0.125,1.0,0.001)
#d = d_in*2.54  # cm

transmission = {}
for d_in in [0.25,0.375,0.5,.75,1.0]:
    d = d_in*2.54  # 3/8" thick NDF
    dstr = str(d_in)
    transmission[dstr] = np.exp(-alpha*d) 

pylab.figure(1)
pylab.clf()
pylab.semilogy(nu,transmission['0.25'],label='0.25inch')
pylab.semilogy(nu,transmission['0.375'],label='0.375inch')
pylab.semilogy(nu,transmission['0.5'],label='0.5inch')
pylab.semilogy(nu,transmission['0.75'],label='0.75inch')
pylab.semilogy(nu,transmission['1.0'],label='1.0inch')
pylab.semilogy(nu,0*nu + 0.05,'-k')
pylab.legend()
pylab.grid()
pylab.xlabel('Frequency (GHz)')
pylab.ylabel('Transmission')
pylab.title('Transmission at 4K, CR110 (Halpern,Gush,Wishnow and De Cosmo)')

pylab.figure(2)
pylab.clf()
pylab.plot(nu,transmission['0.25'],label='0.25inch')
pylab.plot(nu,transmission['0.375'],label='0.375inch')
pylab.plot(nu,transmission['0.5'],label='0.5inch')
pylab.plot(nu,transmission['0.75'],label='0.75inch')
pylab.plot(nu,transmission['1.0'],label='1.0inch')
pylab.plot(nu,0*nu + 0.05,'-k')
pylab.legend()
pylab.grid()
pylab.ylim([0,0.5])
pylab.xlabel('Frequency (GHz)')
pylab.ylabel('Transmission')
pylab.title('Transmission at 4K, CR110 (Halpern,Gush,Wishnow and De Cosmo)')

pylab.show()



