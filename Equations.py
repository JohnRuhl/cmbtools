# -*- coding: utf-8 -*-
# Equations

'''
THERE IS SOMETHING VERY WRONG WITH THESE EQUATIONS
This file holds all the functions for making the theoretical CMB data
'''

import numpy as np
import matplotlib.pyplot as plt


def dustFreqPowLaw(nu, nu0):
    ''' The power law used in some dust equations (mostly a filler function) '''
    return (nu/nu0)**1.59


def blackbodyConvertofNu(nu, List):
    ''' A Plancks Law frequency function to convert the Black Body equation to the right units
        nu,  const = [h, c, k, TVac, TDust]
        TDust=19.6, h=6.62606957*(10**-34), c=299792458, k=1.3806488*(10**-23), TVac = 2.7
    '''
    h, c, k, TVac, TDust = List[:]
    return 2*(h**2)*(nu**4)*np.exp((h*nu)/(k*TVac))/(k*(TVac**2)*(c**2)*((np.exp((h*nu)/(k*TVac)) - 1)**2))


def blackbody(nu, List):                                                     # Why is blackbody function never called?
    ''' Another Plankcs Law equation
        nu, h=6.62606957*(10**-34), c=299792458, k=1.3806488*(10**-23), T=2.7
    '''
    h, c, k, TVac, TDust = List[:]
    return (2*(h*(nu**3))/(c**2))*(1/(np.exp((h*(nu))/(k*TDust)) - 1))


def dust(l):
    ''' Exponential dust function '''
    return (l/80.0)**(-0.42)


def dustRatio(const):
    ''' Provides the ratio to multiply the dust by when given two frequencies of blackbody
        This is one of the two "outermost" functions that is called
        const stands for a List of constants given in main.py
    '''
    # extracting values from list of constants
    nu1, nu2, nu0, List = const['nu1'], const['nu2'], const['nu0'], const['List']
    # finding the dust ratio
    dustRatio = (dustFreqPowLaw(nu1, nu0)*blackbody(nu2, List) * blackbodyConvertofNu(nu2, List)) \
              / (dustFreqPowLaw(nu2, nu0)*blackbody(nu1, List) * blackbodyConvertofNu(nu1, List))
    return dustRatio


# # Seeing some outputs:
# List = [6.62606957*(10**-34), 299792458, 1.3806488*(10**-23), 2.7, 19.6]
# nuList = np.arange(1.0*(10**9), 601*(10**9), 1*(10**9))
# test = [blackbodyConvertofNu(nu, List) for nu in nuList]
# fig = plt.figure()
# plt.plot(np.log10(nuList), np.log10(test))
# plt.savefig('blackbody.png')
# plt.close(fig)