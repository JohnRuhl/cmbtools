# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 11:31:39 2014

@author: sam
"""
def bolo_play(band, band_width, eta):
    import numpy, scipy, math, bolo_functions
    global c, h, k, T_cmb, data, datasrc, numpy, scipy, math, decimal
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
    data['eta'] = eta
    data['tau'] = .018/data['eta']
    data['L'] = 0.6 #scattering
    
    #source parameters for CMB
    datasrc1 = {'name':'CMB', 'eps':1, 'T':2.725, 'tau':data['tau']}
    datasrc2 = {'name':'atm', 'eps':0.1, 'T':230.0, 'tau':data['tau']}
    datasrc = [datasrc1, datasrc2]
    
    bolo_functions.optical_calcs(data, datasrc)
    
    #Bolometer parameters
    data['W'] = 2*data['Qtot'] #total power to put bolo at operating point
    #data['Qtot'] (calculated by optical_calcs, or explicitly defined if you don't run that)
        #note: P_elec = data['W'] - data['Qtot']
    data['n'] = 3.0 #G index - 3 for insulators
    data['T_bolo'] = 0.5
    data['T_base'] = 0.28
    #note: we use W, n, and these T's to calculate dynamic G
    
    data['R_bolo'] = 0.9 #ohms (operating point)
    data['R_load'] = 0.03 #ohms (shunt)
    data['alpha'] = 10.0 #d(logR) / d(logT) at operating point
    data['beta'] = 0.0 #d(logR) / d(logT) at fixed T
    data['tau_0'] = 0.010 #seconds, = C/G
    data['tau_electronics'] = 0.0001 #wild guess - need to figure this out
    data['NEI_squid'] = float(5e-9) #wild guess
    
    #do calculations...
    
    bolo_functions.bolo_calcs(data)
    
    #optical report:
    #report interesting scalars:
    print 'Qtot='+str(data['Qtot'])
    print 'T_RJ_tot = '+str(data['T_RJ_tot'])
    print 'dPdT_RJ = '+str(data['dPdT_RJ'])
    print 'dPdT_cmb = '+str(data['dPdT_cmb'])
    print 'NEP_photon_total = '+str(data['NEP_photon_total'])
    print 'NET_photon_total_RJ = '+str(data['NET_photon_total_RJ'])
    print 'NET_photon_total_cmb = '+str(data['NET_photon_total_cmb'])
    
    for i in range(len(datasrc)):
        print 'Source: '+str(datasrc[i]['name'])
        print 'Q: '+str(datasrc[i]['Q'])
        print 'T_RJ: '+str(datasrc[i]['T_RJ'])
        print 'NEP_photon = '+str(datasrc[i]['NEP_photon'])
        print 'NET_photon_RJ = '+str(datasrc[i]['NET_photon_RJ'])
        print 'NET_photon_cmb = '+str(datasrc[i]['NET_photon_cmb'])
        
def plot_power():
    import pylab as plot
    bands = [94.6, 147.9, 223.0]
    band_widths = [25.8, 38.6, 52.3]
    for i in range(len(bands)):
        bolo_play(bands[i], band_widths[i])
        plot.plot(data['nuGHZ'], n)
    plot.show()
    
def plot_temp(band, band_width, eta):
    import pylab as plot
    import numpy, scipy
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
    optical_calcs(data, datasrc)
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
