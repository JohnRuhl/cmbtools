# -*- coding: utf-8 -*-
# Equations
import numpy as np

# Equation Functions

# The power law used in some dust equations (mostly a filler function)
def dustFreqPowLaw(nu, nu0 = 1):
    return (nu/nu0)**1.59

# A Plancks Law equation
def blackbodyConvertofNu(nu, h = 6.62606957*(10**-34), c = 299792458, k = 1.3806488*(10**-23), T = 2.7):
    return 2*(h**2)*(nu**4)*np.exp((h*nu)/(k*T))/(k*(T**2)*(c**2)*((np.exp((h*nu)/(k*T)) - 1)**2))

# Another Plankcs Law equation
def blackbody(nu, h = 6.62606957*(10**-34), c = 299792458, k = 1.3806488*(10**-23), T = 2.7):
    return (2*(h*(nu**3))/(c**2))*(1/(np.exp((h*(nu))/(k*T)) - 1))

# Exponential dust function
def dustofL(l):
    return (l/80.0)**-0.42

# Provides the ratio to multiply the dust by when given two frequencies of blackbody
def dustRatio(nu1, nu2):
    return (dustFreqPowLaw(nu1)*blackbodyConvertofNu(nu2))/(dustFreqPowLaw(nu2)*blackbodyConvertofNu(nu1))
