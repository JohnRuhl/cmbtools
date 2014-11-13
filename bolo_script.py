# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 10:40:17 2014

@author: sam
"""
global c, h, k, T_cmb, numpy
import numpy, scipy, math, bolo_module
import sys
#physical constants in SI units

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
import imp
f = open('bolo_data.py')
global data
experiment = imp.load_source('experiment', 'r', f)
f.close()
data = experiment.data
datasrc = experiment.datasrc
    
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