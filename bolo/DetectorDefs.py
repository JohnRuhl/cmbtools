"""
Experiment configuration for bolo_calcs and tes_calcs.

defines "detector" dictionary and "optics" list of dictionaries.

Use SI units.

Generic 150GHz instrument, 100mK base, Tc = 200mK
Rbolo = 1.5 Ohm (fmux)
Bosefactor = 0.5
PsatFactor = 2.0
"""

import numpy


def S4_SAT_155():

    T_coldstage = 0.1  #K
    T_stop = 1.0 #K

    # Define the optical elements, for calculating optical efficiency, optical loading, and photon noise,
    # using optical_calcs().
    # These need to be in order, from detector to CMB
    optics = []
    optics.append({'name':'Cold stage feed', 'eps':0.60, 'T':T_coldstage})
    optics.append({'name':'5K lenses',    'eps':0.15, 'T':5.0 })
    optics.append({'name':'30K filter',   'eps':0.02, 'T':30.0 })
    optics.append({'name':'70K filter',   'eps':0.01, 'T':70.0 })
    optics.append({'name':'150K filter',  'eps':0.01, 'T':150.0 })
    optics.append({'name':'window',       'eps':0.03, 'T':280.0 })  
    optics.append({'name':'atm',          'eps':0.03, 'T':250.0 })  
    optics.append({'name':'CMB',          'eps':1.,   'T':2.73 })


    # Set up detector and other parameters.
    detector = {}

    # Set up optical properties;  divide into bins of width band_resolution
    # Generic:  30% banwidth
    #
    band_center = 155.0  #GHz
    band_width = 0.22*band_center
    band_resolution = band_width/100.
    detector['nuGHZ']=numpy.linspace(band_center-0.5*band_width, band_center+0.5*band_width, num=100, endpoint=True)
    detector['nu'] = detector['nuGHZ']*1.e9
    detector['bandshape'] = numpy.ones(len(detector['nu']),float)  # this can be a bandshape
    detector['Npol'] = 1.0
    detector['Nmodes'] = 1.0
    detector['BoseFactor'] = 1.0

    # bolometer parameters;  
    # these are fmux appropriate
    #
    # Thermal properties
    detector['T_bolo'] = 1.8*T_coldstage #note: we use W, n, and these T's to calculate dynamic G
    detector['T_base'] = T_coldstage
    detector['PsatFactor'] = 2.5
    # This should be overwritten by optical_calcs
    detector['Qtot'] = 6e-12   # optical power
    detector['W'] = detector['PsatFactor']*detector['Qtot']   # total power
    detector['n'] = 2.8 #G index - 3 for insulators
    #
    # Electrical properties
    detector['DO_FREQ_RESPONSE'] = False
    detector['R_bolo'] = 1.0 #ohms (operating point)
    detector['R_load'] = 0.03 #ohms (shunt)
    detector['alpha'] = 20.0 #d(logR) / d(logT) at operating point
    detector['beta'] = 0.0 #d(logR) / d(logT) at fixed T
    detector['tau_0'] = 0.010 #seconds, = C/G
    detector['tau_electronics'] = 0.0001 #wild guess - need to figure this out
    detector['NEI_squid'] = 5.e-12 # pA/sqrtHz 
    detector['L'] = 60.e-6 # inductance, set's L/R time constant and thus bandwidth of electronics

    return detector, optics

def S4_SAT_145():
    detector, optics = S4_SAT_155()
    band_center = 145.0  #GHz
    band_width = 0.22*band_center
    band_resolution = band_width/100.
    detector['nuGHZ']=numpy.linspace(band_center-0.5*band_width, band_center+0.5*band_width, num=100, endpoint=True)
    detector['nu'] = detector['nuGHZ']*1.e9
    #
    return detector, optics

def S4_SAT_85():
    detector, optics = S4_SAT_155()
    band_center = 85.0  #GHz
    band_width = 0.24*band_center
    band_resolution = band_width/100.
    detector['nuGHZ']=numpy.linspace(band_center-0.5*band_width, band_center+0.5*band_width, num=100, endpoint=True)
    detector['nu'] = detector['nuGHZ']*1.e9
    #
    optics = []
    optics.append({'name':'Cold stage feed', 'eps':0.60, 'T':T_coldstage})
    optics.append({'name':'5K lenses',    'eps':0.15, 'T':5.0 })
    optics.append({'name':'30K filter',   'eps':0.02, 'T':30.0 })
    optics.append({'name':'70K filter',   'eps':0.01, 'T':70.0 })
    optics.append({'name':'150K filter',  'eps':0.01, 'T':150.0 })
    optics.append({'name':'window',       'eps':0.03, 'T':280.0 })
    optics.append({'name':'atm',          'eps':0.03, 'T':250.0 })
    optics.append({'name':'CMB',          'eps':1.,   'T':2.73 })

    optics['window']['eps'] = 0.02
    optics['atm']['eps'] = 0.05
    optics['5K lenses']['eps']=
    #
    return detector, optics

def S4_SAT_95():
    detector, optics = S4_SAT_85()
    band_center = 95.0  #GHz
    band_width = 0.24*band_center
    band_resolution = band_width/100.
    detector['nuGHZ']=numpy.linspace(band_center-0.5*band_width, band_center+0.5*band_width, num=100, endpoint=True)
    detector['nu'] = detector['nuGHZ']*1.e9
    #
    return detector, optics

def S4_SAT_30():
    detector, optics = S4_SAT_85()
    band_center = 30.0  #GHz
    band_width = 0.30*band_center
    band_resolution = band_width/100.
    detector['nuGHZ']=numpy.linspace(band_center-0.5*band_width, band_center+0.5*band_width, num=100, endpoint=True)
    detector['nu'] = detector['nuGHZ']*1.e9
    #
    for x in optics:
        if x['name']=='30K filter':
            x['eps'] = 0.01
        if x['name']=='window':
            x['eps'] = 0.01
    #
    return detector, optics

def S4_SAT_40():
    detector, optics = S4_SAT_30()
    #
    band_center = 40.0  #GHz
    band_width = 0.30*band_center
    band_resolution = band_width/100.
    detector['nuGHZ']=numpy.linspace(band_center-0.5*band_width, band_center+0.5*band_width, num=100, endpoint=True)
    detector['nu'] = detector['nuGHZ']*1.e9
    return detector, optics


def S4_90():
    detector, optics = S4_150()
    #
    for item in optics: 
        if item['name'] == 'atm':
            item['eps'] = 0.05
        if item['name'] == 'window':
            item['eps'] = 0.01
    band_center = 90.0  #GHz
    band_width = 0.3*band_center
    band_resolution = band_width/100.
    detector['nuGHZ']=numpy.linspace(band_center-0.5*band_width, band_center+0.5*band_width, num=100, endpoint=True)
    detector['nu'] = detector['nuGHZ']*1.e9
    detector['bandshape'] = numpy.ones(len(detector['nu']),float)  # this can be a bandshapekkkkkk
    return detector, optics


