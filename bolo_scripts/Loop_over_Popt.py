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
name = 'Generic150.txt'
import imp
f = open(name)
global data
experiment = imp.load_source('experiment', 'r', f)
f.close()
data = experiment.data
datasrc = experiment.datasrc
print(name)

eps_vec = np.empty(0)
Qtot_vec = np.empty(0)
W_vec = np.empty(0)
N_tot = np.empty(0)
N_johnson = np.empty(0)
N_phonon = np.empty(0)
N_photon = np.empty(0)
N_readout = np.empty(0)
#
for eps_lens in np.arange(0.0,0.2,0.01):  # loss in lens, dumped to same temp as lens

    # change all three lens emissivities
    datasrc[2]['eps'] = eps_lens
    datasrc[5]['eps'] = eps_lens
    datasrc[6]['eps'] = eps_lens

    # calling optical_calcs
    print('Calling optical_calcs'+', eps_lens =' + str(eps_lens))
    bolo_module.optical_calcs(data, datasrc)
    data['W'] = data['PsatFactor']*data['Qtot'] #total power to put bolo at operating point

    # call bolo_calcs and save vector elements to plot
    bolo_module.bolo_calcs(data)   # modifies data
    #
    NoiseType = 'NEP' #'NET_CMB' # string, NET_CMB, NET_RJ, NEP, NEI
    eps_vec = np.append(eps_vec, eps_lens)
    Qtot_vec = np.append(Qtot_vec, data['Qtot'])
    W_vec = np.append(W_vec, data['W'])
    N_tot = np.append(N_tot,data['NoiseTotal'][NoiseType])
    N_phonon = np.append(N_phonon,data['Phonon'][NoiseType])
    N_photon = np.append(N_photon,data['Photon'][NoiseType])
    N_readout = np.append(N_readout,data['Readout'][NoiseType])
    N_johnson = np.append(N_johnson,data['Johnson'][NoiseType])

#toplabel = name[9:-4] + ',  Ptot={:5.2e},  Popt={:5.2e}'.format(data['W'],data['Qtot'])
toplabel = 'NET_CMB vs eps_lens'
xvec = eps_vec
#ynorm = 1e6/np.sqrt(2);  # convert from K/rtHz to uKrtsec
ynorm = 1. # for NEP, remain in W/sqrt(Hz)
pylab.figure(1)
pylab.clf()
pylab.plot(xvec, N_tot*ynorm, label='Total')
pylab.plot(xvec, N_phonon*ynorm,label='Phonon')
pylab.plot(xvec, N_photon*ynorm,label='Photon')
pylab.plot(xvec, N_readout*ynorm,label='Readout')
pylab.plot(xvec, N_johnson*ynorm,label='Johnson')
pylab.xlabel('eps_lens')
pylab.ylabel('NET_CMB (uKrtsec)')
#pylab.ylim([0,16e-17])
pylab.title(toplabel)
pylab.grid()
pylab.legend()
pylab.show()

pylab.figure(2)
pylab.clf()
Rel_MappingSpeed = (N_tot[0]/N_tot)**2
pylab.plot(xvec,Rel_MappingSpeed)
pylab.ylabel('Relative Mapping Speed')
pylab.grid()




