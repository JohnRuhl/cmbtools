import numpy as np
import pylab

pylab.ion()

#AM Stevies
f1,Trj1,B1 = np.loadtxt('Spole_output.txt',unpack=True)
#Trj1 = Trj1/2

#AM mine
#f2,Trj2,B2 = np.loadtxt('SPole_winter_256u.txt',unpack=True)

pylab.figure(1)
#pylab.plot(f1,Trj1,f2,Trj2)
pylab.plot(f1,Trj1)
pylab.ylim([0,300])
pylab.xlim([0,400])
pylab.grid()
pylab.xlabel('Frequency (GHz)')
pylab.ylabel('T_RJ')

#pylab.errorbar(SPT3gfs,SPT3gTrj,xerr=SPT3gdfs,elinewidth=4,fmt='None')
#pylab.legend(['ATM_250u','AM_250u','AM_500u','Quad','Acbar','Bicep2','SPT3G'],loc='upper left')
pylab.legend(['SB','JR'],loc='upper left')
pylab.title('Spole, winter')

