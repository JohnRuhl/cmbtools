# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 10:40:17 2014

@author: sam
"""
global c, h, k, T_cmb, numpy
import numpy, scipy, math, bolo_module
import sys
#physical constants in SI units
c = float(299792456)
h = 6.626068e-34
k = 1.3806503e-23
T_cmb = 2.725 #Kelvin

#optical properties:
#sources (eg CMB, atmos, 100K): eps(nu_vector), T
#datasrc1[name]
#datasrc1[eps]
#datasrc1[T]
#datasrc1[tau] = tau
    
#tau (perhaps for each source)
#data[eta]
#data[Nmodes] (or A*Omega as f(nu))
#data[Npol]
#data[nu]
#data[nuGHZ]
#data[band]
    
#set up data structure
band_width = 25.8  # GHz
band = 94.6  #GHz
data = {'nuGHZ':numpy.arange(band - 0.5*band_width, band + 0.5*band_width,band_width/1000., dtype = float)}
data['nu'] = data['nuGHZ']*float(1e9)

#for i in range(len(data['nu'])):
#    data['nu'][i] = Decimal(data['nu'][i])
#nu = []
#for i in range(len(data['nuGHZ'])):    
#    nu.append(Decimal(data['nuGHZ'][i]*1e9))
#data['nu'] = nu
data['band'] = numpy.ones(len(data['nu']),float)
data['Npol'] = 1.0
data['Nmodes'] = 1.0
data['eta'] = .176
data['tau'] = 1.0 #.018/data['eta']
data['L'] = 0.6 #scattering
data['R_bolo'] = 0.9 #ohms (operating point)
data['R_load'] = 0.03 #ohms (shunt)
data['alpha'] = 10.0 #d(logR) / d(logT) at operating point
data['beta'] = 0.0 #d(logR) / d(logT) at fixed T
data['tau_0'] = 0.010 #seconds, = C/G
data['tau_electronics'] = 0.0001 #wild guess - need to figure this out
data['NEI_squid'] = float(5e-9) #wild guess
 #Bolometer parameters
data['n'] = 3.0 #G index - 3 for insulators
data['T_bolo'] = 0.5
data['T_base'] = 0.28 #note: we use W, n, and these T's to calculate dynamic G

# These need to be in order, from detector to CMB
datasrc1 = {'name':'antenna', 'eps':0., 'T':0.25, 'tau':0.948}
datasrc2 = {'name':'lenslet', 'eps':0., 'T':0.25, 'tau':0.95}
datasrc3 = {'name':'collimating lens', 'eps':0.052, 'T':5.0, 'tau':0.928}
datasrc4 = {'name':'stop', 'eps':0.684, 'T':5.0, 'tau':0.316}
datasrc5 = {'name':'lens filter 1', 'eps':0.05, 'T':5.0, 'tau':0.95}
datasrc6 = {'name':'lens filter 2', 'eps':0.05, 'T':5.0, 'tau':0.95}
datasrc7 = {'name':'aperture lens', 'eps':0.052, 'T':5.0, 'tau':0.948}
datasrc8 = {'name':'nylon filter', 'eps':0.02, 'T':5.0, 'tau':0.98}
datasrc9 = {'name':'field lens', 'eps':0.052, 'T':6.0, 'tau':0.948}
datasrc10 = {'name':'IR shader', 'eps':0.02, 'T':50.0, 'tau':0.98}
datasrc11 = {'name':'window alumina', 'eps':0.013, 'T':60.0, 'tau':0.987}
datasrc12 = {'name':'flat mirror', 'eps':0.01, 'T':280., 'tau':0.989}
datasrc13 = {'name':'secondary', 'eps':0.01, 'T':280., 'tau':0.984}
datasrc14 = {'name':'primary', 'eps':0.01, 'T':220., 'tau':0.99}
datasrc15 = {'name':'atm', 'eps':0.097, 'T':230.0, 'tau':0.903} 
datasrc16 = {'name':'CMB', 'eps':1., 'T':2.725, 'tau':1.0}
datasrc = [datasrc1, datasrc2, datasrc3, datasrc4, datasrc5, datasrc6, datasrc7, datasrc8,\
datasrc9, datasrc10, datasrc11, datasrc12, datasrc13, datasrc14, datasrc15, datasrc16]
    
bolo_module.optical_calcs(data, datasrc)


data['W'] = 2*data['Qtot'] #total power to put bolo at operating point
#data['Qtot'] (calculated by optical_calcs, or explicitly defined if you don't run that)
    #note: P_elec = data['W'] - data['Qtot']
bolo_module.bolo_calcs(data)

print '               Source',
print '  eps ',
print '   T_src',
print '   eta_to_bolo',
print '   Q',
print '       T_RJ',
print '  NEP_photon  NET_photon_RJ NET_photon_cmb'
for i in range(len(datasrc)):
    print "%20s" %datasrc[i]['name'],
    print "%7.3f" %datasrc[i]['eps'],
    print "%7.1f" %datasrc[i]['T'],
    print "%7.3f" %datasrc[i]['eta_to_bolo'],
    print "%10.3e" %datasrc[i]['Q'],
    print "%7.1f" %datasrc[i]['T_RJ'],
    print "%10.3e" %datasrc[i]['NEP_photon'],
    print "%12.3e" %datasrc[i]['NET_photon_RJ'],
    print "%12.3e" %datasrc[i]['NET_photon_cmb']
    
print 'Total: ----------------------------------------------' 
print '                  Qtot = ' "%10.3e" %data['Qtot']
print '              T_RJ_tot = ' "%10.1f" %data['T_RJ_tot']
print '               dPdT_RJ = ' "%10.3e" %data['dPdT_RJ']
print '              dPdT_cmb = ' "%10.3e" %data['dPdT_cmb']
print '      NEP_photon_total = ' "%10.3e" %data['NEP_photon_total']
print '   NET_photon_total_RJ = ' "%10.3e" %data['NET_photon_total_RJ']
print '  NET_photon_total_cmb = ' "%10.3e" %data['NET_photon_total_cmb']
print '-----------------------------------------------------' 
    
