# -*- coding: utf-8 -*-
# Equations

'''
THERE IS SOMETHING VERY WRONG WITH THESE EQUATIONS
Make sure the two parameters exp and dustT are going to the right places
     dustT could be actually teh parameter for dust**-0.42
** All 'master' equations are designated with Signal and must only accept 2 inputs
This file holds all the functions for making the theoretical CMB dust data
'''

import numpy as np
from Functions import parameterSplit



##########################################################################################
                                    # BMode

def BModeSignal(BMode, parameters):
    ''' Fits BMode with equation a*BMode
        parameters:
            - must be a list of length 1
            - amplitude = paramss[:]

       add +c to the parameters?
    '''
    ampl = parameters
    BModeSignal = (ampl)*BMode
    return np.array(BModeSignal)



##########################################################################################
                                    # Dust

def dustSignal((lList, const, frequency), parameters):
    ''' Creates the dust signal
        format: ((lList, const), parameters)
        const holds the constants
        parameters is the list of fit parameters.
            parameters = [amplitude, exponent, dust temperature]
            They are *added* to their respetive theoretical values
            ex: (1+amplitude)*somefunction()

        add +c to the parameters?
    '''
    ampl, exp, dustT = parameters[:]
    prefactor = dustPreFactor(const, frequency, exp, dustT)
    dustSignal = [(ampl)*(1*10**-20)*prefactor*dust(l) for l in lList]  # *** multiplying by the constant in temporary ***
    return np.array(dustSignal)


def dust(lList):  # CHANGE THIS NAME
    ''' Exponential dust function '''
    lList = np.float64(lList)
    dust = (lList/80)**(-0.42)
    return dust


def dustPreFactor(const, frequency, exp, dustT):
    ''' Dust frequency powerlaw and blackbody curve
    '''
    # getting out the constants
    nu0, List, TVac = const['nu0'], const['List'], const['TVac']
    # finding the dust frequency powerlaw
    powerlaw = dustFreqPowLaw(frequency, nu0, exp)
    # finding the blackbody curve
    blackbody = blackbody_nu(frequency, List, dustT)
    blackbody_conversion = blackbody_convert(frequency, List, TVac)
    # makeing the dust pre-factor
    dustPreFactor = powerlaw*blackbody/blackbody_conversion
    return dustPreFactor


def dustFreqPowLaw(frequency, nu0, exp):
    '''
        the exponent is added to the theoretical value of 1.59
    '''
    frequency, nu0 = np.float64(frequency), np.float64(nu0)
    powerlaw = (frequency/nu0)**(exp)  # theory exp is 1.59
    return powerlaw


def blackbody_nu(nu, List, dustT):
    ''' Another Planck's Law equation
        nu, h=6.62606957*(10**-34), c=299792458, k=1.3806488*(10**-23), T=2.7
        the T is added to the mean dust value of 19.6 K
    '''
    h, c, k = List[0:3]
    blackbody_Nu = (2*(h*(nu**3))/(c**2))*(1/(np.exp((h*(nu))/(k*(dustT))) - 1))
    return blackbody_Nu


def blackbody_convert(nu, List, T):
    ''' A Plancks Law frequency function to convert the Black Body equation to the right units
        nu,  const = [h, c, k, TVac, TDust]
        TDust=19.6, h=6.62606957*(10**-34), c=299792458, k=1.3806488*(10**-23), TVac = 2.7
        It is the derivative of blackbody_nu
    '''
    h, c, k, TVac, TDust = List[:]
    numerator = (2*(h**2)*(nu**4)*np.exp((h*nu)/(k*TVac)))
    denominator = (k*(TVac**2)*(c**2)*((np.exp((h*nu)/(k*TVac)) - 1)**2))
    return numerator/denominator



##########################################################################################
                                    # Total Signal

def totalSignal((lList, const, frequency, BMode_data), parameters):  # Change name of inputs
    ''' Creates the total theoretical signal
    '''
    # creating parameter inputs
    params_BM, params_d = parameterSplit(parameters)
    # getting functions
    BMode = BModeSignal(BMode_data, params_BM)
    Dust = dustSignal((lList, const, frequency), params_d)
    # adding into one signal
    signal = np.array(BMode+Dust)
    return signal
