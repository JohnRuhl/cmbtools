# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 10:40:17 2014

@author: sam
"""
global c, h, k, T_cmb, numpy
import numpy, scipy, math, bolo_module
#physical constants in SI units
c = float(299792456)
h = 6.626068e-34
k = 1.3806503e-23
T_cmb = 2.725 #Kelvin
#optical properties:
#sources (eg CMB, atmos, 100K): eps(nu_vector), T
# note: to investigate band edge placement, need eps(nu) for atmosphere, ie line effects
#data.
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
band_width = 25.8
band = 94.6
data = {'nuGHZ':numpy.arange(band - 0.5*band_width, band + 0.5*band_width,0.1, dtype = float)}
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
data['T_base'] = 0.28
    #note: we use W, n, and these T's to calculate dynamic G
#source parameters for CMB
datasrc1 = {'name':'CMB', 'eps':1., 'T':2.725, 'tau':data['tau']}
datasrc2 = {'name':'atm', 'eps':0.097, 'T':230.0, 'tau':data['tau']} 
datasrc3 = {'name':'antenna', 'eps':0., 'T':0.25, 'tau':data['tau']}
datasrc4 = {'name':'lenslet', 'eps':0., 'T':0.25, 'tau':data['tau']}
datasrc5 = {'name':'collimating lens', 'eps':0.052, 'T':5.0, 'tau':data['tau']}
datasrc6 = {'name':'stop', 'eps':0.684, 'T':5.0, 'tau':data['tau']}
datasrc7 = {'name':'lens filter 1', 'eps':0.05, 'T':5.0, 'tau':data['tau']}
datasrc8 = {'name':'lens filter 2', 'eps':0.05, 'T':5.0, 'tau':data['tau']}
datasrc9 = {'name':'aperture lens', 'eps':0.052, 'T':5.0, 'tau':data['tau']}
datasrc10 = {'name':'nylon filter', 'eps':0.02, 'T':5.0, 'tau':data['tau']}
datasrc11 = {'name':'field lens', 'eps':0.052, 'T':6.0, 'tau':data['tau']}
datasrc12 = {'name':'IR shader', 'eps':0.02, 'T':50.0, 'tau':data['tau']}
datasrc13 = {'name':'window alumina', 'eps':0.013, 'T':60.0, 'tau':data['tau']}
datasrc14 = {'name':'flat mirror@', 'eps':0.01, 'T':280., 'tau':data['tau']}
datasrc15 = {'name':'secondary@', 'eps':0.01, 'T':280., 'tau':data['tau']}
datasrc16 = {'name':'primary', 'eps':0.01, 'T':220., 'tau':data['tau']}
datasrc = [datasrc1, datasrc2, datasrc3, datasrc4, datasrc5, datasrc6, datasrc7, datasrc8,\
datasrc9, datasrc10, datasrc11, datasrc12, datasrc13, datasrc14, datasrc15, datasrc16]
    
bolo_module.optical_calcs(data, datasrc)
data['W'] = 2*data['Qtot'] #total power to put bolo at operating point
#data['Qtot'] (calculated by optical_calcs, or explicitly defined if you don't run that)
    #note: P_elec = data['W'] - data['Qtot']
bolo_module.bolo_calcs(data)

print 'Qtot=' "%.3e" %data['Qtot']
print 'T_RJ_tot = ' "%.3e" %data['T_RJ_tot']
print 'dPdT_RJ = ' "%.3e" %data['dPdT_RJ']
print 'dPdT_cmb = ' "%.3e" %data['dPdT_cmb']
print 'NEP_photon_total = ' "%.3e" %data['NEP_photon_total']
print 'NET_photon_total_RJ = ' "%.3e" %data['NET_photon_total_RJ']
print 'NET_photon_total_cmb = ' "%.3e" %data['NET_photon_total_cmb']
    
for i in range(len(datasrc)):
    print 'Source:' + datasrc[i]['name']
    print 'Q: '"%.3e" %datasrc[i]['Q']
    print 'T_RJ: ' "%.3e" %datasrc[i]['T_RJ']
    print 'NEP_photon = ' "%.3e" %datasrc[i]['NEP_photon']
    print 'NET_photon_RJ = ' "%.3e" %datasrc[i]['NET_photon_RJ']
    print 'NET_photon_cmb = ' "%.3e" %datasrc[i]['NET_photon_cmb']