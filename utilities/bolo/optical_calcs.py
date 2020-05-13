# -*- coding: utf-8 -*-
"""
optical_calcs(detector,optics)
detector and optics are typically from an experiment file.

In this code, NEP is in units of Watts/sqrt(Hz), while NEP2 is in units of Watts**2/Hz.
"""

import numpy as np

c = 299792456.0
h = 6.626068e-34
k = 1.3806503e-23
T_cmb = 2.725

def optical_calcs(detector, optics):

    # Calculate power load and photon NEP, by summing from the detector outward.
    # initialize eta_tot and the things we'll be summing into
    NEP2_photon_total = 0.0
    detector['Qtot'] = 0.0
    eta_tot = 1.0 # this is the optical efficiency from the bolo to the element in current loop.
    power_prefactor = h*detector['Nmodes']*detector['Npol']

    # Loop over optical sources, from detector to sky.
    Ptot = 0.
    for i in range(len(optics)):

        # Note that detector['nu'] is typically a vector.
        x = h*detector['nu']/(k*optics[i]['T'])

        #occupation number; this is eps*eta_tot*band/(e^x -1).
        n = eta_tot*optics[i]['eps']*detector['bandshape']\
          /np.longdouble((np.exp(x)-float(1)))     

        # photon power on bolometer at each freq
        P = power_prefactor*detector['nu']*n #power per mode per Hz, times Nmodes*Npol

        Ptot = Ptot + P

        #total power integrated across band
        optics[i]['Q'] = np.trapz(P, detector['nu'])

        # update eta_tot for the next loop through here
        # eta_tot is the optical effiency between the relevant element and the bolo.
        optics[i]['eta_to_bolo'] = eta_tot
        if optics[i]['name'] != 'CMB':
            eta_tot = eta_tot*(1-optics[i]['eps'])

        #sum powers and NEP^2 to get total
        detector['Qtot'] = detector['Qtot'] + optics[i]['Q']
        ####### end for loop over sources

    #photon noise integral - single-mode, single-polarization.
    # Note this does have Nmodes and Npol in the prefactor.
    noise_prefactor = 2 #*detector['Nmodes']*detector['Npol']) # have to handle other modes in quadrature!
    integrand = noise_prefactor*(h*detector['nu']*Ptot + detector['BoseFactor']*Ptot**2)
    NEP2_photon_total = np.trapz(integrand, detector['nu'])


    # Reset the total power on the detector based on the calculated optical power.
    detector['W'] = detector['PsatFactor']*detector['Qtot']

    # Now that we have eta_tot, we can do conversions to T_RJ, T_cmb units.
    #
    # dP_bolo/dT_RJ, referenced to celstial_sphere (ie what dT_RJ do you need there to create dP_bolo)
    prefactor = eta_tot*detector['Npol']*detector['Nmodes']
    RJintegrand = k*prefactor*detector['bandshape']   
    detector['dPdT_RJ'] = np.trapz(RJintegrand, detector['nu'])

    # dP/dT_cmb, again referenced to T at celestial sphere
    x = h*detector['nu']/(k*T_cmb)
    prefactor =eta_tot*(detector['Npol']*h**2)/(k*c**2*T_cmb**2) 
    AOmega = detector['Nmodes']*c**2/detector['nu']**2
    integrand2 = AOmega*(detector['nu']**4)*np.exp(x)/((np.exp(x)-1.0)**2)
    detector['dPdT_cmb'] = prefactor*np.trapz(integrand2, detector['nu'])

    # Record all the converted total quantities
    detector['eta_tot'] = eta_tot
    detector['NEP_photon_total'] = np.sqrt(NEP2_photon_total)  
    detector['Photon'] = {'NEP': detector['NEP_photon_total']}  # New arrangement
    detector['T_RJ_tot'] = detector['Qtot']/detector['dPdT_RJ']
    detector['NET_photon_total_RJ'] = detector['NEP_photon_total']/detector['dPdT_RJ']
    detector['NET_photon_total_cmb'] = detector['NEP_photon_total']/detector['dPdT_cmb']
    return detector, optics

