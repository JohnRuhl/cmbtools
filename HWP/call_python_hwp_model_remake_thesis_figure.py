import python_hwp_model
import numpy as nm
import pylab

# this code remakes Figure 8.1 in Sean Bryan's thesis about half-wave plates
# http://arxiv.org/pdf/1402.2591v1.pdf

# calculate HWP mueller matrix elements for nominal thicknesses of everything
freq,T,rho,c,s,R,tau,h,q = python_hwp_model.calculate_hwp_mueller_matrix(use_preset_thickness_and_index='150GHz', \
                                                                         cooled=True, \
                                                                         freq=nm.arange(0.0,300.0,0.1))

# calculate for an ideal retarder designed for 150 GHz
T_ideal = nm.ones_like(freq)
rho_ideal = nm.zeros_like(freq)
c_ideal = nm.cos(nm.pi*(freq/150.0))
s_ideal = nm.sin(nm.pi*(freq/150.0))

# plot the mueller matrix elements
pylab.figure(1,figsize=(9,6))
pylab.clf()
pylab.subplot(2,2,1)
pylab.plot(freq,T,'b',label='Sapphire Slab + AR')
pylab.plot(freq,T_ideal,'--r',label='Ideal Retarder')
pylab.ylim((-1.05,1.05))
pylab.xlabel('Frequency [GHz]')
pylab.title('T')
pylab.fill([130.0,160.0,160.0,130.0], [-1.05,-1.05,1.05,1.05], 'k', alpha=0.2, edgecolor='k',label='Spider 150 GHz Band')
pylab.legend(loc='lower right')

pylab.subplot(2,2,2)
pylab.plot(freq,rho,'b',label='Sapphire Slab + AR')
pylab.plot(freq,rho_ideal,'--r',label='Ideal Retarder')
pylab.ylim((-1.05,1.05))
pylab.xlabel('Frequency [GHz]')
pylab.title(r'$\rho$')
pylab.fill([130.0,160.0,160.0,130.0], [-1.05,-1.05,1.05,1.05], 'k', alpha=0.2, edgecolor='k',label='Spider 150 GHz Band')

pylab.subplot(2,2,3)
pylab.plot(freq,c,'b',label='Sapphire Slab + AR')
pylab.plot(freq,c_ideal,'--r',label='Ideal Retarder')
pylab.ylim((-1.05,1.05))
pylab.xlabel('Frequency [GHz]')
pylab.title('c')
pylab.fill([130.0,160.0,160.0,130.0], [-1.05,-1.05,1.05,1.05], 'k', alpha=0.2, edgecolor='k',label='Spider 150 GHz Band')

pylab.subplot(2,2,4)
pylab.plot(freq,s,'b',label='Sapphire Slab + AR')
pylab.plot(freq,s_ideal,'--r',label='Ideal Retarder')
pylab.ylim((-1.05,1.05))
pylab.xlabel('Frequency [GHz]')
pylab.title('s')
pylab.fill([130.0,160.0,160.0,130.0], [-1.05,-1.05,1.05,1.05], 'k', alpha=0.2, edgecolor='k',label='Spider 150 GHz Band')

pylab.tight_layout()

pylab.ion()
pylab.show()
