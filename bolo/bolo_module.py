# -*- coding: utf-8 -*-
"""
Created on Tue Aug 5 09:03:22 2014
Translated, debugged, and completed by Samuel Musilli from code originally in Matlab
by James T. Sayre as adapted from John Ruhl's code.
Most comments are copied from JT's original code and are not my own.
However, they may be adjusted if they refer to Matlab-specific grammar.

In this code, NEP is in units of Watts/sqrt(Hz), while NEP2 is in units of Watts**2/Hz.

"""
def optical_calcs(data, datasrc):
    c = 299792456.0
    h = 6.626068e-34
    k = 1.3806503e-23
    T_cmb = 2.725
    import numpy
    #from operator import add
    #input optical properties are in datasrc dictionaries.
    #One source per dictionary index, eg datasrc[i]['field']eps
    # sources (eg CMB or atmosphere or optical elements...): eps(nu_vector), T
    #outputs go to both
    # (a) the datasrc dictionaries (to store something associated with that source) and
    # (b) the data dictionary.
    #
    #data[dPdT_cmb] (band average)
    #data[dPdT_RJ] (band average)
    #datasrcX[Q] (optical power)
    #datasrcX[NEP_photon]
    #datasrcX[net_photon_cmb]
    #datasrcX[NET_photon_RJ]
    #data[Qtot]
    #data[]

    # First we calculate dP/dT_RJ and dP/dT_cmb, where these are defined as being
    # the change in optical power absorbed by the bolometer for a small temperature
    # change in a source in the sky.
    #
    

    # Calculate power load and photon NEP, by summing from the detector outward.
    # initialize eta_tot and the things we'll be summing into
    NEP2_photon_total = 0.0
    data['Qtot'] = 0.0
    eta_tot = 1.0 # this is the optical efficiency from the bolo to the element in current loop.
    noise_prefactor = 2*((h**2)*(data['nu']**2))*(2*data['Nmodes']*data['Npol'])
    for i in range(len(datasrc)):
        x = h*data['nu']/(k*datasrc[i]['T'])
        #occupation number; this is eps*eta_tot*band/(e^x -1).
        # band can be freq-dependent
        n = eta_tot*datasrc[i]['eps']*data['bandshape']\
          /numpy.longdouble((numpy.exp(x)-float(1)))     

        # photon power on bolometer
        P = h*data['nu']*n #power per mode per Hz
        #total power integrated across band
        datasrc[i]['Q'] = numpy.trapz(P, data['nu'])
        #datasrc[i]['T_RJ'] = datasrc[i]['Q']/data['dPdT_RJ']

        #photon noise integral
        integrand = noise_prefactor*(n + data['BoseFactor']*n**2) # if you want that term... )
        NEP2_photon = numpy.trapz(integrand, data['nu'])
        datasrc[i]['NEP_photon'] = numpy.sqrt(NEP2_photon)
        #datasrc[i]['NET_photon_cmb'] = datasrc[i]['NEP_photon'] / data['dPdT_cmb']
        #datasrc[i]['NET_photon_RJ'] = datasrc[i]['NEP_photon'] / data['dPdT_RJ']

        # update eta_tot for the next loop through here
        # eta_tot is the optical effiency between the relevant element and the bolo.
        datasrc[i]['eta_to_bolo'] = eta_tot
        if datasrc[i]['name'] != 'CMB':
            eta_tot = eta_tot*(1-datasrc[i]['eps'])

        #sum powers and NEP^2 to get total
        data['Qtot'] = data['Qtot'] + datasrc[i]['Q']
        NEP2_photon_total = NEP2_photon_total + NEP2_photon

        # end for loop over sources

    # Now that we have eta_tot, we can do conversions to T_RJ, T_cmb units.
    #
    # dP_bolo/dT_RJ, referenced to celstial_sphere (ie what dT_RJ do you need there to create dP_bolo)
    prefactor = eta_tot*data['Npol']*data['Nmodes']
    RJintegrand = k*prefactor*data['bandshape']   
    data['dPdT_RJ'] = numpy.trapz(RJintegrand, data['nu'])

    # dP/dT_cmb, again referenced to T at celestial sphere
    x = h*data['nu']/(k*T_cmb)
    AOmega = data['Nmodes'] * c**2 / data['nu']**2
    integrand2 = AOmega*data['bandshape']*(prefactor*h**2*(numpy.power(data['nu'],4)))*numpy.exp(x)/\
      (k*(c**2)*(T_cmb**2)*(numpy.exp(x)-1.0)**2)
    data['dPdT_cmb'] = numpy.trapz(integrand2, data['nu'])

    # Record all the converted total quantities
    data['eta_tot'] = eta_tot
    data['NEP_photon_total'] = numpy.sqrt(NEP2_photon_total)
    data['Photon'] = {'NEP': data['NEP_photon_total']}
    data['T_RJ_tot'] = data['Qtot']/data['dPdT_RJ']
    data['NET_photon_total_RJ'] = data['NEP_photon_total']/data['dPdT_RJ']
    data['NET_photon_total_cmb'] = data['NEP_photon_total']/data['dPdT_cmb']
    return data


def bolo_calcs(data):
    #does calculations relevant to bolometer responsivity, NEP and NEI for thermal
    # fluctuation (G), Johnson, and readout noise.
    #Inputs and outputs all live in 'data' structure
    import math, numpy
    #physical constants in SI units
    c = 299792456.0
    h = 6.626068e-34
    k = 1.3806503e-23
    #input bolo properties (from experiment file)
    #data['W'] (total power to put bolo at T_bolo)
    #data['Q'] (optical power - could be from optical_calcs or separate estimate, P_elec = W-Q)
    #data['n'] (G index - 1 for metals, 3 for crystalline dielectrics/superconductors)
    #data['T_bolo'] (bolometer operating temperature (>T_base!))
    #Tbase (bolometer bath temperature)
    #data['R_bolo'] (operating point resistance of TES)
    #data['R_load'] (shunt/load Thevenin equivalent series resistance)
    #data['alpha'] (T/R d_R/d_T, const I -- Irwin & Hilton 2.3)
    #data['beta'] (I/R d_R/d_I, const T --Irwin & Hilton 2.3)
    #data['tau_0'] (C/G. bolometer natural time constant)
    #data['tau_el'] (bolometer electrical time constant)
    #data['NEI_readout'] (readout/squid noise)
    #data['f'] (audio frequency bandwidth over which to calculate noise)
    #outputs:
    #data['s_w'] (current responsivity S(omega) = (dI/dP) -- I&H 2.3)
    #data['tau_eff'] (effective time constant (sped up by TES ETF))
    #data['Lg'] (loopgain due to ETF)
    #data['Gbar'] (average thermal conductivity)
    #data['G'] (dynamic thermal conductivity)
    #data['P_el'] (electrical power through the TES)
    #data['I_O'] (current through TES at P_el and R_0)
    #data[Johnson[NEP]] (Johnson noise, referred as NEP on bolo)
    #data[Johnson[NEI]] (Johnson noise)
    #data[Johnson[NET]] (NET of Johnson noise (vs. data['f']))
    #data[phonon[NEP]] (thermal fluctuation/G/phonon noise)
    #data[phonon[NEI]] (phonon nosie, referred as NEI at readout)
    #data[phonon[NET]] (NET of phonon noise (vs. data['f']))
    #relations between W, G, and T
    # Gbar = W/(Tbolo-Tbase)
    # kappa = W/(Tbolo**n - Tbase**n)
    # G = n*kappa*Tbolo^n-1
    #unpack some of the structure variables to make them more compact and
    #easier to follow in the calculations

    # Bolometer and electronics properties
    R0 = data['R_bolo']
    RL = data['R_load']
    tau_0 = data['tau_0']
    tau_el = data['tau_electronics']
    L = data['L'] #inductance. Not sure what that is for fmux, but I think it's the full series inductor.
    Tbolo = data['T_bolo']
    Tbase = data['T_base']
    beta = data['beta']
    alpha = data['alpha']
    n = data['n']
    #ratio of equivalent load resistance to Rbolo
    R_rat = RL/R0

    #set up T vector for phonon noise integral
    T = numpy.arange(Tbase, Tbolo, (Tbolo-Tbase)/1000.)

    # set up (audio) frequency vector, eg relevant for time constants/etc
    # Note, 8/9/2017, JR simplified to DC.
    data['f'] = 0.0 #numpy.arange(1.0,10.0,0.1)   
    w = 2*math.pi*data['f']   # angular frequency vector

    #do some calculatin'
    # note that data['W'] is Psat (total, ie optical + electrical)
    Gbar = data['W']/(Tbolo-Tbase)
    kappa = data['W']/(Tbolo**n - Tbase**n)
    G = n*kappa*Tbolo**(n-1)   # dynamic G
    P_el = data['W'] - data['Qtot'] # Electrical power
    Lg = (P_el*alpha/(G*Tbolo)) # loop gain
    I_0 = numpy.sqrt(P_el/R0)
    #effective time consant (for zero bias inductance ?????)
    tau_eff = tau_0*(1+beta+R_rat)/(1+beta+R_rat+(1-R_rat)*Lg)
    #time constant for constant-current fluctuations
    tau_i = tau_0/(1-Lg)
    #calculate tau+/- from Irwin & Hilton section 2.3, table 1
    #### XXX CHECK THESE OUT AND WHAT L SHOULD BE
    A = 1/(2*tau_el)+1/(2*tau_i)
    B = .5*numpy.sqrt((1/tau_el -1/tau_i)**2 - 4*R0/L*(Lg*(2+beta)/tau_0))
    tau_p = 1/(A+B)
    tau_m = 1/(A-B)
    #calculate responsivity S(omega) from Irwin and Hilton sec 2.3
    T1 = -1/(I_0*R0)
    T2 = 1/(2+beta)
    T3 = (1-tau_p/tau_i)/(1+1j*w*tau_p)
    T4 = (1-tau_m/tau_i)/(1+1j*w*tau_m)
    s_w = T1 #T1*T2*T3*T4
    #### This was causing problems, T4 appears to be way too low...
    #s_w = T1*numpy.ones(len(w))
    data['I_0'] = I_0
    data['V_0'] = I_0*R0
    data['Sdc'] = T1
    data['T2'] = T2
    data['T3'] = T3
    data['T4'] = T4
    data['s_w'] = T1

    #Phonon (Thermal Fluctuation) noise
    #k integral = int{[(t*k(t))/(T*k(T))]^2dt}/int{[k(t)/k(T)]}
    #between Tbase and T (Tbolo)
    #note that k(t) = kappa*(T_bolo^n), but the kappas cancel, giving k(t) = (T_bolo^n)
    #heat link integral (Mather 1982, eq 33)
    b = n-1
    F_Tbolo_Tbase = sum((T*(T**b)/(Tbolo*(Tbolo**b)))**2)/sum((T**b)/(Tbolo**b))
    data['Phonon'] = {'NEP':numpy.sqrt(4*k*Tbolo**2*G*F_Tbolo_Tbase)}
    data['Phonon']['NEI'] = data['Phonon']['NEP']*abs(s_w)
    
    # Convert photon noise to NEI
    data['Photon']['NEI'] = data['Photon']['NEP']*abs(s_w)
    
    #TES Johnson noise
    Si_tes = 4*k*Tbolo*R0*I_0**2*(1+w**2*tau_0**2)*abs(s_w)**2/Lg**2
    #load resistor noise (assume single resistor at 4K)
    T_L = 4.2
    Si_L = 4*k*T_L*I_0**2*RL*(Lg-1)**2*(1+w**2*tau_i**2)*abs(s_w)**2/Lg**2
    data['Johnson'] = {'NEI':numpy.sqrt(Si_tes+Si_L)}
    data['Johnson']['NEP'] = data['Johnson']['NEI']/abs(s_w)

    #readout noise
    data['Readout'] = {'NEI':data['NEI_squid']} #originally NEI_readout, assumed to be NEI_squid
    data['Readout']['NEP'] = data['Readout']['NEI']/abs(s_w) #Si_tes/abs(s_w)**2 + Si_L/abs(s_w)**2

    # total NEP
    nep2tot = 0.
    for source in ['Johnson','Photon','Phonon','Readout']:
        nep2tot = nep2tot + data[source]['NEP']**2
    data['NoiseTotal']= {'NEP': numpy.sqrt(nep2tot)}

    #If we've done the optical calculations and have dPdT values, we can get
    #NET on the sky. Otherwise, return NET as 0.
    if 'dPdT_RJ' in data:
        data['Phonon']['NET_RJ'] = data['Phonon']['NEP']/data['dPdT_RJ']
        data['Photon']['NET_RJ'] = data['Photon']['NEP']/data['dPdT_RJ']
        data['Johnson']['NET_RJ'] = data['Johnson']['NEP']/data['dPdT_RJ']
        data['Readout']['NET_RJ'] = data['Readout']['NEP']/data['dPdT_RJ']
        data['NoiseTotal']['NET_RJ'] = data['NoiseTotal']['NEP']/data['dPdT_RJ']
    else:
        data['Phonon']['NET_RJ'] = 0
        data['Photon']['NET_RJ'] = 0
        data['Johnson']['NET_RJ'] = 0
        data['Readout']['NET_RJ'] = 0
        data['NoiseTotal']['NET_RJ'] = 0
    if 'dPdT_cmb' in data:
        data['Phonon']['NET_CMB'] = data['Phonon']['NEP']/data['dPdT_cmb']
        data['Photon']['NET_CMB'] = data['Photon']['NEP']/data['dPdT_cmb']
        data['Johnson']['NET_CMB'] = data['Johnson']['NEP']/data['dPdT_cmb']
        data['Readout']['NET_CMB'] = data['Readout']['NEP']/data['dPdT_cmb']
        data['NoiseTotal']['NET_CMB'] = data['NoiseTotal']['NEP']/data['dPdT_cmb']
    else:
        data['Phonon']['NET_CMB'] = 0
        data['Photon']['NET_CMB'] = 0
        data['Johnson']['NET_CMB'] = 0
        data['Readout']['NET_CMB'] = 0
        data['NoiseTotal']['NET_CMB'] = 0

    #pack up all the variables to be returned into the data structure
    data['G'] = G
    data['Lg'] = Lg
    data['P_el'] = P_el
    data['I_0'] = I_0
    data['tau_eff'] = tau_eff
    data['Gbar'] = Gbar
    data['kappa'] = kappa
    data['tau_i'] = tau_i
    data['s_w'] = s_w
    #return data
    #print data['Phonon']['NEP'], data['Johnson']['NEP'], data['NEP_photon_total']
