"""
Experiment configuration for bolo_calcs and tes_calcs.

defines "detector" dictionary and "optics" list of dictionaries.

Use SI units.

"""

import numpy

def S4_150():
    ''' Ground-based 150GHz fmux bolos at South Pole (approximately).
    '''

    T_coldstage = 0.1  #K
    T_stop = 1.0 #K

    # Define the optical elements, for calculating optical efficiency, optical loading, and photon noise,
    # using optical_calcs().
    # These need to be in order, from detector to CMB
    optics = []
    optics.append({'name':'300mK feed',   'eps':0.05, 'T':T_coldstage})
    optics.append({'name':'300mK filter', 'eps':0.05, 'T':T_coldstage})
    optics.append({'name':'4K lenses',    'eps':0.05, 'T':4.0 })
    optics.append({'name':'Cold stop',    'eps':0.10, 'T':T_stop })   # changes in loop
    optics.append({'name':'5K lens',      'eps':0.10, 'T':5.0 })
    optics.append({'name':'5K filter',    'eps':0.05, 'T':5.0 })
    optics.append({'name':'50K filter',   'eps':0.05, 'T':50.0 })
    optics.append({'name':'window',       'eps':0.02, 'T':280.0 })  
    optics.append({'name':'atm',          'eps':0.07, 'T':230.0 })
    optics.append({'name':'CMB',          'eps':1.,   'T':2.73 })

    #optics.append({'name':'260K mirrors', 'eps':0.01, 'T':260.0 })

    # Set up detector and other parameters.
    detector = {}

    # Set up optical properties;  divide into bins of width band_resolution
    # Generic:  30% banwidth
    #
    band_center = 150.0  #GHz
    band_width = 0.3*band_center
    band_resolution = band_width/1000.   
    detector['nuGHZ']=numpy.arange(band_center-0.5*band_width, band_center+0.5*band_width, band_resolution)
    detector['nu'] = detector['nuGHZ']*1.e9
    detector['bandshape'] = numpy.ones(len(detector['nu']),float)  # this can be a bandshape
    detector['Npol'] = 1.0
    detector['Nmodes'] = 1.0
    detector['BoseFactor'] = 0.5

    # bolometer parameters;  
    # these are fmux appropriate
    #
    # Thermal properties
    detector['T_bolo'] = 2*T_coldstage #note: we use W, n, and these T's to calculate dynamic G
    detector['T_base'] = T_coldstage
    detector['PsatFactor'] = 2.0
    # This should be overwritten by optical_calcs
    detector['Qtot'] = 6e-12   # optical power
    detector['W'] = detector['PsatFactor']*detector['Qtot']   # total power
    detector['n'] = 2.8 #G index - 3 for insulators
    #
    # Electrical properties
    detector['DO_FREQ_RESPONSE'] = False
    detector['R_bolo'] = 1.5 #ohms (operating point)
    detector['R_load'] = 0.03 #ohms (shunt)
    detector['alpha'] = 20.0 #d(logR) / d(logT) at operating point
    detector['beta'] = 0.0 #d(logR) / d(logT) at fixed T
    detector['tau_0'] = 0.010 #seconds, = C/G
    detector['tau_electronics'] = 0.0001 #wild guess - need to figure this out
    detector['NEI_squid'] = 5.e-12 # pA/sqrtHz 
    detector['L'] = 60.e-6 # inductance, set's L/R time constant and thus bandwidth of electronics

    return detector, optics


def S4_90():
    ''' Similar to S4_90, but 90GHz
    '''

    detector, optics = S4_150()
    #
    for item in optics: 
        if item['name'] == 'atm':
            item['eps'] = 0.05
        if item['name'] == 'window':
            item['eps'] = 0.01
    band_center = 90.0  #GHz
    band_width = 0.3*band_center
    band_resolution = band_width/1000.
    detector['nuGHZ']=numpy.arange(band_center-0.5*band_width, band_center+0.5*band_width, band_resolution)
    detector['nu'] = detector['nuGHZ']*1.e9
    detector['bandshape'] = numpy.ones(len(detector['nu']),float)  # this can be a bandshapekkkkkk
    return detector, optics


