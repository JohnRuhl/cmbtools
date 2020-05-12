# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 10:40:17 2014
@author: sam
"""
def run(name):
    global c, h, k, T_cmb, numpy
    import numpy, scipy, math, bolo_module
    import sys
    #physical constants in SI units
    c = float(299792456)
    h = 6.626068e-34
    k = 1.3806503e-23
    T_cmb = 2.725 #Kelvin

    #set up data structure
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
    data['W'] = data['PsatFactor']*data['Qtot'] #total power to put bolo at operating point
    
    # call bolo_calcs
    print('Calling bolo_calcs')
    bolo_module.bolo_calcs(data)
    NEP_total = 0

    # Print out a table of values for each source
    print ' \t   Source\t',
    print ' eps \t',
    print ' T_src',
    print 'eta_to_bolo',
    print ' Q\t'
    #print '    T_RJ',
    #print ' NEP_photon NET_photon_RJ NET_photon_cmb'
    for i in range(len(datasrc)):
        print "%20s" %datasrc[i]['name'],
        print "%7.3f" %datasrc[i]['eps'],
        print "%8.1f" %datasrc[i]['T'],
        print "%9.3f" %datasrc[i]['eta_to_bolo'],
        print "%12.3e" %datasrc[i]['Q']
        #print "%7.1f" %datasrc[i]['T_RJ'],
        #print "%10.3e" %datasrc[i]['NEP_photon'],
        #print "%12.3e" %datasrc[i]['NET_photon_RJ'],
        #print "%12.3e" %datasrc[i]['NET_photon_cmb']

    NEP_total = numpy.sqrt(data['NEP_photon_total']**2 + data['Johnson']['NEP']**2 + data['Phonon']['NEP']**2 + data['Readout']['NEP']**2)
    NEI_total = NEP_total*abs(data['s_w'])
    NET_cmb_total = NEP_total/data['dPdT_cmb']/numpy.sqrt(2) # uKrtsec rather than uK/rtHz
    # print out the totals
    print 'Total: ----------------------------------------------'
    print ' P_optical = ' "%10.3e" %data['Qtot']
    print ' P_total = ' "%10.3e" %data['W']
    print ' T_RJ_tot = ' "%10.1f" %data['T_RJ_tot']
    print ' dPdT_RJ = ' "%10.3e" %data['dPdT_RJ']
    print ' dPdT_cmb = ' "%10.3e" %data['dPdT_cmb']
    print('')
    #
    print(' I_bolo = {0:3.2e} Amps'.format(data['I_0']))
    print(' R_bolo = {0:3.2e} Ohms'.format(data['R_bolo']))
    print(' V_bolo = {0:3.2e} Volts'.format(data['V_0']))
    print(' S_DC = {0:3.2e} '.format(data['Sdc']))
    print('')
    #
    print ' NEP_photon_total = ' "%10.3e" %data['Photon']['NEP']
    print ' NEP_phonon_total = ' "%10.3e" %data['Phonon']['NEP']
    print ' NEP_Johnson_total = ' "%10.3e" %data['Johnson']['NEP']
    print ' NEP_Readout_total = ' "%10.3e" %data['Readout']['NEP']
    print ' NEP_total = ' "%10.3e" %NEP_total
    print(' NEP_total = {0:3.2e}'.format(data['NoiseTotal']['NEP']))
    print('')
    #
    print ' NEI_photon_total = ' "%10.3e" %data['Photon']['NEI']
    print ' NEI_phonon_total = ' "%10.3e" %data['Phonon']['NEI']
    print ' NEI_Johnson_total = ' "%10.3e" %data['Johnson']['NEI']
    print ' NEI_Readout_total = ' "%10.3e" %data['Readout']['NEI']
    print ' NEI_total = ' "%10.3e" %NEI_total
    print('')
    #
    print ' NET_photon_total_RJ = ' "%10.3e" %data['NET_photon_total_RJ']
    print ' NET_photon_total_cmb = ' "%10.3e" %data['NET_photon_total_cmb']
    print ' NET_cmb_total = ' "%10.3e" %NET_cmb_total#[0],
    print 'uKrtsec'
    print('')
    print '-----------------------------------------------------' 
    print data['Sdc']

    return data
