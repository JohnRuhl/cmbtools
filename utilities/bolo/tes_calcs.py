# -*- coding: utf-8 -*-
"""
tes_calcs(data)
(used to be called bolo_calcs, but this is for tes's, so renamed.  8/18/2018.)

inputs and outputs all live in the "data" dictionary, which is typically loaded from an experiment configuation file.
You may want to run optical_calcs first, if you want NETs and or photon noise.
In this code, NEP is in units of Watts/sqrt(Hz), while NEP2 is in units of Watts**2/Hz.
Physical constants and all calculations are done in SI units.
"""
import numpy

def tes_calcs(data):
    c = 299792456.0
    h = 6.626068e-34
    k = 1.3806503e-23

    # Bolometer and electronics properties
    R0 = data['R_bolo']
    RL = data['R_load']
    tau_0 = data['tau_0']
    tau_el = data['tau_electronics']
    L = data['L'] #inductance. Not sure what that is for fmux, but I think it's the full series inductor.
    Tbolo = data['T_bolo']
    Tbase = data['T_base']
    beta = data['beta']
    alpha = data['alpha']   # dln(R)/dln(T)
    n = data['n']  # thermal leg index
    #ratio of equivalent load resistance to Rbolo
    R_rat = RL/R0

    #set up T vector for phonon noise integral
    T = numpy.arange(Tbase, Tbolo, (Tbolo-Tbase)/1000.)


    # Calculations
    # Note that data['W'] is Psat (total, ie optical + electrical)
    Gbar = data['W']/(Tbolo-Tbase)
    kappa = data['W']/(Tbolo**n - Tbase**n)
    G = n*kappa*Tbolo**(n-1)   # dynamic G
    P_el = data['W'] - data['Qtot'] # Electrical power
    Lg = (P_el*alpha/(G*Tbolo)) # loop gain
    I_0 = numpy.sqrt(P_el/R0)

    T1 = -1/(I_0*R0)
    T2 = 1/(2+beta)
    #print(data['DO_FREQ_RESPONSE'])
    if data['DO_FREQ_RESPONSE']==True :
        # set up (audio) frequency vector, eg relevant for time constants/etc
        # Note, 8/9/2017, JR simplified to DC.  Need to revise this if we ever go back to doing this as a function of audio frequency.
        data['f'] = 0.0 #numpy.arange(1.0,10.0,0.1)   
        w = 2*numpy.pi*data['f']   # angular frequency vector

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
        T3 = (1-tau_p/tau_i)/(1+1j*w*tau_p)
        T4 = (1-tau_m/tau_i)/(1+1j*w*tau_m)
        s_w = T1 #T1*T2*T3*T4
        #### This was causing problems, T4 appears to be way too low...
        #s_w = T1*numpy.ones(len(w))
        data['tau_eff'] = tau_eff
        data['tau_i'] = tau_i
    else:
        T3 = 1.
        T4 = 1.
        s_w = T1
        data['s_w'] = 1.

    data['Sdc'] = T1
    data['I_0'] = I_0
    data['V_0'] = I_0*R0
    data['T2'] = T2
    data['T3'] = T3
    data['T4'] = T4

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
    #Si_tes = 4*k*Tbolo*R0*(s_w**2)*I_0**2/Lg**2  # Irwin and Hilton, but with chsi set to 1.
    Si_tes = 4*k*Tbolo/R0/Lg**2   # easy version.
    #*(1+w**2*tau_0**2)*abs(s_w)**2/Lg**2
    #load resistor noise (assume single resistor at 4K)
    T_L = 4.2
    Si_L = 0# 4*k*T_L*I_0**2*RL*(Lg-1)**2/Lg**2
    #*(1+w**2*tau_i**2)*abs(s_w)**2/Lg**2
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
        for ntype in ['Phonon','Photon','Johnson','Readout','NoiseTotal']:
            data[ntype]['NET_RJ'] = data[ntype]['NEP']/data['dPdT_RJ']
    else:
        for ntype in ['Phonon','Photon','Johnson','Readout','NoiseTotal']:
            data[ntype]['NET_RJ'] = 0.

    if 'dPdT_cmb' in data:
        for ntype in ['Phonon','Photon','Johnson','Readout','NoiseTotal']:
            data[ntype]['NET_CMB'] = data[ntype]['NEP']/data['dPdT_cmb']
    else:
        for ntype in ['Phonon','Photon','Johnson','Readout','NoiseTotal']:
            data[ntype]['NET_CMB'] = 0.

    #pack up all the variables to be returned into the data structure
    data['G'] = G
    data['Lg'] = Lg
    data['P_el'] = P_el
    data['I_0'] = I_0
    data['Gbar'] = Gbar
    data['kappa'] = kappa
    data['s_w'] = s_w
    return data
