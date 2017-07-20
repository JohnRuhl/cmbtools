# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 10:40:17 2014

@author: sam
Modified by JR
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
band_width = 84.2  # GHz
band = 288.6  #GHz
#band = 300  #GHz
data = {'nuGHZ':numpy.arange(band - 0.5*band_width, band + 0.5*band_width,band_width/1000., dtype = float)}
data['nu'] = data['nuGHZ']*float(1e9)
data['band'] = numpy.ones(len(data['nu']),float)

# Set up bolo parameters
data['Npol'] = 1.0
data['Nmodes'] = 1.0
data['eta'] = .3  # optical efficiency from sky to bolo
data['tau'] = 1.0 #.018/data['eta']
data['L'] = 16e-6 #inductance for electrical bandwidth
data['R_bolo'] = 0.9 #ohms (operating point)
data['R_load'] = 0.03 #ohms (shunt)
data['alpha'] = 100.0 #d(logR) / d(logT) at operating point
data['beta'] = 0.0 #d(logR) / d(logT) at fixed T
data['tau_0'] = 0.020 #seconds, = C/G
#data['tau_electronics'] = 0.001 #wild guess - need to figure this out
data['tau_electronics'] = 1*data['L']/data['R_bolo']
data['NEI_squid'] = float(7e-12) #wild guess
 #Bolometer parameters
data['n'] = 2.7 #G index - 3 for insulators
data['T_bolo'] = 0.5
data['T_base'] = 0.28 #note: we use W, n, and these T's to calculate dynamic G

# atmos 0.5, 4.0 pW
# stop 0.0238pW
# 4K 0.0637 pW
# cmb 0.0876 pW
# VCS1 0.268pW 
# VCS2 0.518pW
# window 0.465 pW (273, 0.5% emissive)

# These need to be in order, from detector to CMB
datasrc1 = {'name':'antenna', 'eps':.5, 'T':0.25, 'tau':0.5}
datasrc2 = {'name':'lens 1', 'eps':0.05, 'T':5.0, 'tau':0.95}
datasrc3 = {'name':'lens 2', 'eps':0.05, 'T':5.0, 'tau':0.95}
datasrc4 = {'name':'stop', 'eps':0.2, 'T':3.0, 'tau':0.8}
datasrc5 = {'name':'nylon filter', 'eps':0.02, 'T':5.0, 'tau':0.98}
datasrc6 = {'name':'snout filter', 'eps':0.03, 'T':5.0, 'tau':0.97}
datasrc7 = {'name':'hwp', 'eps':0.05, 'T':20., 'tau':0.95}
datasrc8 = {'name':'Ade filter VCS1', 'eps':0.02, 'T':45., 'tau':0.97}
datasrc9 = {'name':'IR shaders VCS1 (3)', 'eps':0.003, 'T':45.0, 'tau':0.997}
datasrc10 = {'name':'IR shaders VCS2 (4)', 'eps':0.004, 'T':155.0, 'tau':0.996}
datasrc11 = {'name':'IR shaders 300K (2)', 'eps':0.002, 'T':270.0, 'tau':0.998}
datasrc12 = {'name':'window', 'eps':0.005, 'T':250.0, 'tau':0.995}
#datasrc13 = {'name':'atm', 'eps':0.006, 'T':240.0, 'tau':0.977}  # low
datasrc13 = {'name':'atm', 'eps':0.047, 'T':240.0, 'tau':0.977}  # high
datasrc14 = {'name':'CMB', 'eps':1., 'T':2.725, 'tau':1.0}
datasrc = [datasrc1, datasrc2, datasrc3, datasrc4, datasrc5, datasrc6, datasrc7, datasrc8,\
datasrc9, datasrc10, datasrc11, datasrc12, datasrc13, datasrc14]
    
bolo_module.optical_calcs(data, datasrc)


data['W'] = 2*data['Qtot'] #total power to put bolo at operating point
#data['Qtot'] (calculated by optical_calcs, or explicitly defined if you don't run that)
    #note: P_elec = data['W'] - data['Qtot']
bolo_module.bolo_calcs(data)


# Print out a table of values for each source
print '               Source',
print '  eps ',
print '   T_src',
print 'eta_to_bolo',
print ' Q',
print '     T_RJ',
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

NEP_total = numpy.sqrt(data['NEP_photon_total']**2 + data['Phonon']['NEP']**2 + data['Johnson']['NEP']**2 + data['Readout']['NEP']**2)

NET_cmb_total = NEP_total/data['dPdT_cmb']/numpy.sqrt(2)  # uKrtsec rather than uK/rtHz
    
# print out the totals
print 'Total: ----------------------------------------------' 
print '                  Qtot = ' "%10.3e" %data['Qtot']
print '              T_RJ_tot = ' "%10.1f" %data['T_RJ_tot']
print '               dPdT_RJ = ' "%10.3e" %data['dPdT_RJ']
print '              dPdT_cmb = ' "%10.3e" %data['dPdT_cmb']
print '      NEP_photon_total = ' "%10.3e" %data['NEP_photon_total']
print '      NEP_phonon_total = ' "%10.3e" %data['Phonon']['NEP']
print '     NEP_Johnson_total = ' "%10.3e" %data['Johnson']['NEP'][0]
print '     NEP_Readout_total = ' "%10.3e" %data['Readout']['NEP'][0]
print '   NET_photon_total_RJ = ' "%10.3e" %data['NET_photon_total_RJ']
print '  NET_photon_total_cmb = ' "%10.3e" %data['NET_photon_total_cmb']
print '             NEP_total = ' "%10.3e" %NEP_total[0]
print '         NET_cmb_total = ' "%10.3e" %NET_cmb_total[0],
print 'Krtsec'
print '-----------------------------------------------------' 

    
