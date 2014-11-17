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
band_width = 25.8 # GHz
band = 94.6 #GHz
data = {'nuGHZ':numpy.arange(band - 0.5*band_width, band + 0.5*band_width,band_width/1000., dtype = float)}
data['nu'] = data['nuGHZ']*float(1e9)
data['band'] = numpy.ones(len(data['nu']),float)
# Set up bolo parameters
data['Npol'] = 1.0
data['Nmodes'] = 1.0
data['eta'] = .352 # optical efficiency from sky to bolo
#data['eta'] = .173 # optical efficiency from sky to bolo
data['tau'] = 1.0 #.018/data['eta']
data['L'] = 0.6 #scattering?
data['R_bolo'] = 0.9 #ohms (operating point)
data['R_load'] = 0.03 #ohms (shunt)
data['alpha'] = 30.0 #d(logR) / d(logT) at operating point
data['beta'] = 0.0 #d(logR) / d(logT) at fixed T
data['tau_0'] = 0.020 #seconds, = C/G
data['tau_electronics'] = 0.001 #wild guess - need to figure this out
data['NEI_squid'] = float(5e-15) #wild guess
#Bolometer parameters
data['n'] = 2.7 #G index - 3 for insulators
data['T_bolo'] = 0.5
data['T_base'] = 0.28 #note: we use W, n, and these T's to calculate dynamic G
# These need to be in order, from detector to CMB
datasrc1 = {'name':'antenna', 'eps':0., 'T':0.25, 'tau':0.948}
datasrc2 = {'name':'lenslet', 'eps':0., 'T':0.25, 'tau':0.95}
datasrc3 = {'name':'collimating lens', 'eps':0.052, 'T':5.0, 'tau':0.928}
datasrc4 = {'name':'stop', 'eps':.37, 'T':5.0, 'tau':0.63}
#datasrc4 = {'name':'stop', 'eps':.69, 'T':5.0, 'tau':0.31}
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
# Print out a table of values for each source
print ' Source',
print ' eps ',
print ' T_src',
print 'eta_to_bolo',
print ' Q',
print ' T_RJ',
print ' NEP_photon NET_photon_RJ NET_photon_cmb'
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
NEP_total = numpy.sqrt(data['NEP_photon_total']**2 + data['Phonon']['NEP']**2 + data['Johnson']['NEP']**2 + data['Readout']['NEP']**2)
NET_cmb_total = NEP_total/data['dPdT_cmb']/numpy.sqrt(2) # uKrtsec rather than uK/rtHz
# print out the totals
print 'Total: ----------------------------------------------'
print ' Qtot = ' "%10.3e" %data['Qtot']
print ' T_RJ_tot = ' "%10.1f" %data['T_RJ_tot']
print ' dPdT_RJ = ' "%10.3e" %data['dPdT_RJ']
print ' dPdT_cmb = ' "%10.3e" %data['dPdT_cmb']
print ' NEP_photon_total = ' "%10.3e" %data['NEP_photon_total']
print ' NEP_phonon_total = ' "%10.3e" %data['Phonon']['NEP']
print ' NEP_Johnson_total = ' "%10.3e" %data['Johnson']['NEP'][0]
print ' NEP_Readout_total = ' "%10.3e" %data['Readout']['NEP'][0]
print ' NET_photon_total_RJ = ' "%10.3e" %data['NET_photon_total_RJ']
print ' NET_photon_total_cmb = ' "%10.3e" %data['NET_photon_total_cmb']
print ' NEP_total = ' "%10.3e" %NEP_total[0]
print ' NET_cmb_total = ' "%10.3e" %NET_cmb_total[0],
print 'uKrtsec'
print '-----------------------------------------------------' 
