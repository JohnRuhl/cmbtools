# -*- coding: utf-8 -*-
# ModelEquations

'''
All 'master' ModelEquations are designated with Signal and must only accept 2 inputs
This file holds all the functions for making the theoretical CMB Models
'''

import numpy as np
from Functions import parameterSplit


##########################################################################################
                                    # Total Signal

def totalSignal((lList, const, frequency, BMode_data), parameters):
    ''' totalSignal((lList, const, frequency, BMode_data), parameters)
        Creates the total theoretical signal
        Calls the BModeSignal and dustSignal functions
        Adds them
    '''
    # creating parameter inputs
    params_BM, params_d = parameterSplit(parameters)
    # getting functions
    BMode = BModeSignal(BMode_data, params_BM)
    Dust = dustSignal((lList, const, frequency), params_d)
    # adding into one signal
    signal = np.array(BMode+Dust)
    return signal



##########################################################################################
                                    # BMode

def BModeSignal(BMode, parameters):
    ''' BModeSignal(BMode, parameters)
        Fits BMode with equation a*BMode
        parameters:
            - integer
            - amplitude = params[:]
    '''
    # getting parameters
    R = parameters[0]
    # making BMode model
    BModeSignal = (R)*BMode
    return np.array(BModeSignal)



##########################################################################################
                                    # Dust

def dustSignal((lList, const, frequency), parameters):
    ''' dustSignal((lList, const, frequency), parameters)
        Creates the dust signal
        const holds the constants
        parameters is the list of fit parameters.
            parameters = [amplitude, exponent, dust temperature]
    '''
    # getting parameter
    ampl, exp = parameters[:]  # , dustT
    # getting prefactor
    prefactor = dustPreFactor(const, frequency, exp)  # , dustT)
    # making dust model
    dustSignal = [(ampl)*(5300.*1.15*10.**(-20.))*prefactor*dust(l) for l in lList]  # *** multiplying by the constant in temporary ***
    return np.array(dustSignal)


def dust(lList):  # *** change this name ***
    ''' dust(lList)
        Exponential dust function
    '''
    lList = np.float64(lList)  # precision and datatype correct
    # making model
    dust = (lList/80.)**(-0.42)
    return dust


def dustPreFactor(const, frequency, exp, dustT=19.6):
    ''' dustPreFactor(const, frequency, exp, dustT):
        Dust frequency powerlaw and blackbody curve
    '''
    # getting out the constants
    nu0, List, Tcmb = const['nu0'], const['List'], const['Tcmb']
    # finding the dust frequency powerlaw
    powerlaw = dustFreqPowLaw(frequency, nu0, exp)
    # finding the blackbody curve
    blackbody = blackbody_nu(frequency, List, dustT)
    # getting the conversion for the blackbody
    blackbody_conversion = blackbody_convert(frequency, List, Tcmb)
    # makeing the dust pre-factor
    dustPreFactor = powerlaw*blackbody/blackbody_conversion
    return dustPreFactor


def dustFreqPowLaw(frequency, nu0, exp):
    ''' dustFreqPowLaw(frequency, nu0, exp)
        the exponent is added to the theoretical value of 1.59
    '''
    frequency, nu0 = np.float64(frequency), np.float64(nu0)
    powerlaw = (frequency/nu0)**(exp)  # theory exp is 1.59
    return powerlaw


def blackbody_nu(nu, List, dustT=19.6):
    ''' blackbody_nu(nu, List, dustT)
        Another Planck's Law equation
        nu, h=6.62606957*(10**-34), c=299792458, k=1.3806488*(10**-23), T=2.7
        the T starts at the mean dust value of 19.6 K
    '''
    h, c, k = List[0:3]
    blackbody_Nu = (2.*(h*(nu**3.))/(c**2.))*(1./(np.exp((h*(nu))/(k*(dustT))) - 1.))
    return blackbody_Nu


def blackbody_convert(nu, List, T):
    ''' blackbody_convert(nu, List, T)
        A Plancks Law frequency function to convert the Black Body equation to the right units
        That's why it uses Tcmb instead of TDust, it is converting to the same metric
        nu,  const = [h, c, k, Tcmb, TDust]
        TDust=19.6, h=6.62606957*(10**-34), c=299792458, k=1.3806488*(10**-23), Tcmb = 2.7
        It is the derivative of blackbody_nu, converting to CMB
    '''
    h, c, k, Tcmb, TDust = List[:]
    numerator = (2.*(h**2.)*(nu**4.)*np.exp((h*nu)/(k*Tcmb)))
    denominator = (k*(Tcmb**2.)*(c**2.)*((np.exp((h*nu)/(k*Tcmb)) - 1.)**2.))
    return numerator/denominator
