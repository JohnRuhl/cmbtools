# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 10:51:20 2014

@author: sam
"""

def plot_temp(band, band_width, eta):
    import pylab as plot
    import numpy, scipy, bolo_module
    global c, h, k, T_cmb, numpy
    c = float(299792456)
    h = 6.626068e-34
    k = 1.3806503e-23
    T_cmb = 2.725
    data = {'nuGHZ':numpy.arange(band - 0.5*band_width, band + 0.5*band_width,0.1, dtype = float)}
    data['T_bolos'] = numpy.arange(0.4,0.6,0.001)
    data['T_base'] = 0.28
    data['nu'] = data['nuGHZ']*float(1e9)
    data['band'] = numpy.ones(len(data['nu']),float)
    data['Npol'] = 1.0
    data['Nmodes'] = 1.0
    data['eta'] = 0.088
    data['tau'] = .018/eta
    data['L'] = 0.6
    datasrc1 = {'name':'CMB', 'eps':1, 'T':2.725, 'tau':data['tau']}
    datasrc2 = {'name':'atm', 'eps':1, 'T':230.0, 'tau':data['tau']}
    data['n'] = 3.0
    data['R_bolo'] = 0.9 #ohms (operating point)
    data['R_load'] = 0.03 #ohms (shunt)
    data['alpha'] = 10.0 #d(logR) / d(logT) at operating point
    data['beta'] = 0.0 #d(logR) / d(logT) at fixed T
    data['tau_0'] = 0.010 #seconds, = C/G
    data['tau_electronics'] = 0.0001 #wild guess - need to figure this out
    data['NEI_squid'] = float(5e-9) #wild guess
    datasrc = [datasrc1, datasrc2]
    empty = numpy.empty(len(data['T_bolos']), dtype = float)
    data['NET_cmb_phonon'] = empty
    data['NET_cmb_Johnson'] = empty
    data['NET_cmb_readout'] = empty
    data['NET_rj_phonon'] = empty
    data['NET_rj_Johnson'] = empty
    data['NET_rj_readout'] = empty
    bolo_module.optical_calcs(data, datasrc)
    data['W'] = 2*data['Qtot']
    for i in range(len(data['T_bolos'])):
        data['T_bolo'] = data['T_bolos'][i]
    bolo_calcs(data)
    data['NET_cmb_phonon'][i] = data['Phonon']['NET_CMB']
    data['NET_cmb_Johnson'][i] = sum(data['Johnson']['NET_CMB'])/len(data['Johnson']['NET_CMB'])
    data['NET_cmb_readout'][i] = sum(data['Readout']['NET_CMB'])/len(data['Readout']['NET_CMB'])
    data['NET_rj_phonon'][i] = data['Phonon']['NET_RJ']
    data['NET_rj_Johnson'][i] = sum(data['Johnson']['NET_RJ'])
    data['NET_rj_readout'][i] = sum(data['Readout']['NET_RJ'])
    plot.gca().set_color_cycle(['red', 'green', 'blue'])
    plot.plot(data['T_bolos'], data['NET_cmb_phonon'])
    print data['NET_cmb_phonon']
    plot.plot(data['T_bolos'], data['NET_cmb_Johnson'])
    print data['NET_cmb_Johnson']
    plot.plot(data['T_bolos'], data['NET_cmb_readout'])
    print data['NET_cmb_readout']
    plot.xlabel('Bolo temperature in K')
    plot.ylabel('NET in uK/sqrt(Hz)')
    plot.title('NET vs. T_bolo over band '+str(band)+' GHz')
    plot.show()
    #print data['NET_photon_totals']
