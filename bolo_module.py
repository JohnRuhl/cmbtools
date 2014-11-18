# -*- coding: utf-8 -*-
"""
Created on Tue Aug 5 09:03:22 2014
Translated, debugged, and completed by Samuel Musilli from code originally in Matlab
by James T. Sayre as adapted from John Ruhl's code.
Most comments are copied from JT's original code and are not my own.
However, they may be adjusted if they refer to Matlab-specific grammar.
"""
def optical_calcs(data, datasrc):
    c = 299792456.0
    h = 6.626068e-34
    k = 1.3806503e-23
    T_cmb = 2.725
    import numpy
    from operator import add
    #input optical properties are in datasrc dictionaries.
    #One source per dictionary index, eg datasrc[i]['field']
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
    # AOmega = n_modes * Lambda**2
    #AOmega = data['Nmodes'] * c**2 / numpy.power(data['nu'],2)
    AOmega = data['Nmodes'] * c**2 / data['nu']**2
    # dP/dT_RJ
    # Note that eta*tau is the optical efficiency from the sky to the bolo.
    # data['eta'] is the "sky to bolo" total optical efficiency.
    # data['tau'] is vestigial, set to one.
    prefactor = float(data['tau']*data['eta']*data['Npol']*data['Nmodes'])
    RJintegrand = k*prefactor*numpy.ones(len(data['nu']))
    data['dPdT_RJ'] = numpy.trapz(RJintegrand, data['nu'])
    # dP/dT_cmb
    x = h*data['nu']/(k*T_cmb)
    integrand2 = AOmega*(prefactor*h**2*(numpy.power(data['nu'],4)))*numpy.exp(x)/\
    (k*(c**2)*(T_cmb**2)*(numpy.exp(x)-1.0)**2)
    data['dPdT_cmb'] = numpy.trapz(integrand2, data['nu'])
    # calculate NEP's etc for each source.
    NEP2_photon_total = 0.0
    data['Qtot'] = 0.0
    eta_tot = 1.0
    for i in range(len(datasrc)):
    x = h*data['nu']/(k*datasrc[i]['T'])
    #occupation number; this is eps*eta_tot*band/(e^x -1).
    # band*tau is freq-dependent, and n the case of tau element-dependent, optical efficiency
    n = datasrc[i]['eps']*eta_tot*data['band']\
    /numpy.longdouble((numpy.exp(x)-float(1)))
    #power per mode per Hz
    P = h*data['nu']*n
    #total power integrated across band
    datasrc[i]['Q'] = numpy.trapz(P, data['nu'])
    datasrc[i]['T_RJ'] = datasrc[i]['Q']/data['dPdT_RJ']
    #phonon noise integral
    integrand = 2*((h**2)*(data['nu']**2))*(data['Nmodes']*data['Npol'])*(n+n**2)
    #print n, n**2, map(add,n,n**2)
    NEP2_photon = numpy.trapz(integrand, data['nu'])
    datasrc[i]['NEP_photon'] = numpy.sqrt(NEP2_photon)
    datasrc[i]['NET_photon_cmb'] = datasrc[i]['NEP_photon'] / data['dPdT_cmb']
    datasrc[i]['NET_photon_RJ'] = datasrc[i]['NEP_photon'] / data['dPdT_RJ']
    datasrc[i]['eta_to_bolo'] = eta_tot
    # update eta_tot for the next loop through here
    # eta_tot is the optical effiency between the relevant element and the bolo.
    eta_tot = eta_tot*datasrc[i]['tau']
    #sum to get total
    data['Qtot'] = data['Qtot'] + datasrc[i]['Q']
    NEP2_photon_total = add(NEP2_photon_total, NEP2_photon)
    data['T_RJ_tot'] = data['Qtot']/data['dPdT_RJ']
    data['NEP_photon_total'] = numpy.sqrt(NEP2_photon_total)
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
    #input bolo properties:
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
    #global data
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
    data['f'] = numpy.arange(1.0,10.0,0.01)
    #convert frequency vector to radians for responsivity calculation
    w = 2*math.pi*data['f']
    #do some calculatin'
    Gbar = data['W']/(Tbolo-Tbase)
    kappa = data['W']/(Tbolo**n - Tbase**n)
    G = n*kappa*Tbolo**(n-1)
    P_el = data['W'] - data['Qtot']
    Lg = (P_el*alpha/(G*Tbolo))
    I_0 = numpy.sqrt(P_el/R0)
    #effective time consant (for zero bias inductance ?????)
    tau_eff = tau_0*(1+beta+R_rat)/(1+beta+R_rat+(1-R_rat)*Lg)
    #time constant for constant-current fluctuations
    tau_i = tau_0/(1-Lg)
    #calculate tau+/- from Irwin & Hilton section 2.3, table 1
    A = 1/(2*tau_el)+1/(2*tau_i)
    B = .5*numpy.sqrt((1/tau_el -1/tau_i)**2 - 4*R0/L*(Lg*(2+beta)/tau_0))
    tau_p = 1/(A+B)
    tau_m = 1/(A-B)
    #calculate responsivity S(omega) from Irwin and Hilton sec 2.3
    T1 = -1/(I_0*R0)
    T2 = 1/(2+beta)
    T3 = (1-tau_p/tau_i)/(1+1j*w*tau_p)
    T4 = (1-tau_m/tau_i)/(1+1j*w*tau_m)
    s_w = T1*T2*T3*T4
    #### This was causing problems, T4 appears to be way too low...
    #s_w = T1*numpy.ones(len(w))
    data['Sdc'] = T1
    data['T2'] = T2
    data['T3'] = T3
    data['T4'] = T4
    data['s_w'] = s_w
    #Phonon (Thermal Fluctuation) noise
    #k integral = int{[(t*k(t))/(T*k(T))]^2dt}/int{[k(t)/k(T)]}
    #between Tbase and T (Tbolo)
    #note that k(t) = kappa*(T_bolo^n), but the kappas cancel, giving k(t) = (T_bolo^n)
    #heat link integral (Mather 1982, eq 33)
    b = n-1
    F_Tbolo_Tbase = sum((T*(T**b)/(Tbolo*(Tbolo**b)))**2)/sum((T**b)/(Tbolo**b))
    data['Phonon'] = {'NEP':numpy.sqrt(4*k*Tbolo**2*G*F_Tbolo_Tbase)}
    data['Phonon']['NEI'] = data['Phonon']['NEP']*abs(s_w)
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
    #If we've done the optical calculations and have dPdT values, we can get
    #NET on the sky. Otherwise, return NET as 0.
    if 'dPdT_RJ' in data:
    data['Phonon']['NET_RJ'] = data['Phonon']['NEP']/data['dPdT_RJ']
    data['Johnson']['NET_RJ'] = data['Johnson']['NEP']/data['dPdT_RJ']
    data['Readout']['NET_RJ'] = data['Readout']['NEP']/data['dPdT_RJ']
    else:
    data['Phonon']['NET_RJ'] = 0
    data['Johnson']['NET_RJ'] = 0
    data['Readout']['NET_RJ'] = 0
    if 'dPdT_cmb' in data:
    data['Phonon']['NET_CMB'] = data['Phonon']['NEP']/data['dPdT_cmb']
    data['Johnson']['NET_CMB'] = data['Johnson']['NEP']/data['dPdT_cmb']
    data['Readout']['NET_CMB'] = data['Readout']['NEP']/data['dPdT_cmb']
    else:
    data['Phonon']['NET_CMB'] = 0
    data['Johnson']['NET_CMB'] = 0
    data['Readout']['NET_CMB'] = 0
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
    return data
def bolo_play(band, band_width, eta):
    import numpy, scipy, math
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
    # data['nu'][i] = Decimal(data['nu'][i])
    #nu = []
    #for i in range(len(data['nuGHZ'])):
    # nu.append(Decimal(data['nuGHZ'][i]*1e9))
    #data['nu'] = nu
    data['band'] = numpy.ones(len(data['nu']),float)
    data['Npol'] = 1.0
    data['Nmodes'] = 1.0
    data['eta'] = eta
    data['tau'] = 1.0 #.018/data['eta']
    data['L'] = 0.6 #scattering
    #source parameters for CMB
    datasrc1 = {'name':'CMB', 'eps':1, 'T':2.725, 'tau':data['tau']}
    datasrc2 = {'name':'atm', 'eps':0.097, 'T':230.0, 'tau':data['tau']}
    datasrc = [datasrc1, datasrc2]
    optical_calcs(data, datasrc)
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
    data['NEI_squid'] = float(5e-12) #wild guess
    #do calculations...
    bolo_calcs(data)
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
