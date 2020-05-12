import pylab 
import numpy as np
import bolo_module

global c, h, k, T_cmb, numpy

pylab.ion()

#physical constants in SI units
c = float(299792456)
h = 6.626068e-34
k = 1.3806503e-23
T_cmb = 2.725 #Kelvin

#set up data structure
name = 'bolo_dataSPT150.txt'
import imp
f = open(name)
global data
experiment = imp.load_source('experiment', 'r', f)
f.close()
data = experiment.data
datasrc = experiment.datasrc
print(name)

# calling optical_calcs
print('Calling optical_calcs')
bolo_module.optical_calcs(data, datasrc)
data['W'] = 3.5*data['Qtot'] #total power to put bolo at operating point

R_bolo = np.arange(0.5,5.0,0.1)

Rplot = np.empty(0)
N_tot = np.empty(0)
N_johnson = np.empty(0)
N_phonon = np.empty(0)
N_photon = np.empty(0)
N_readout = np.empty(0)
for R in R_bolo:
    data['R_bolo'] = R
    Rplot = np.append(Rplot,R)
    # call bolo_calcs
    bolo_module.bolo_calcs(data)   # modifies data

    NoiseType = 'NET_CMB'
    N_tot = np.append(N_tot,data['NoiseTotal'][NoiseType])
    N_phonon = np.append(N_phonon,data['Phonon'][NoiseType])
    N_photon = np.append(N_photon,data['Photon'][NoiseType])
    N_readout = np.append(N_readout,data['Readout'][NoiseType])
    N_johnson = np.append(N_johnson,data['Johnson'][NoiseType])

print(data['W'])
print(data['Qtot'])

#ynorm = 1
ynorm = 1e6/np.sqrt(2)  # convert K/rtHz to uKrtsec
toplabel = name[9:-4] + ',  Ptot={:5.2e},  Popt={:5.2e}'.format(data['W'],data['Qtot'])
pylab.figure(1)
pylab.clf()
pylab.plot(Rplot, N_tot*ynorm, label='Total')
pylab.plot(Rplot, N_phonon*ynorm,label='Phonon')
pylab.plot(Rplot, N_photon*ynorm,label='Photon')
pylab.plot(Rplot, N_readout*ynorm,label='Readout')
pylab.plot(Rplot, N_johnson*ynorm,label='Johnson')
pylab.xlabel('R_bolo')
pylab.ylabel('NEP_tot')
#pylab.ylim([0,16e-17])
pylab.title(toplabel)
pylab.grid()
pylab.legend()
pylab.show()





