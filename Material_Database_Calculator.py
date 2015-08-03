# -*- coding: utf-8 -*-
# A Thermal Conduction Calculator

from Material_Database_Data import materials
from Material_Database_Functions import *
# import Thermal_Conduction_Calculator_functions as fn  # (already imported in _data)

'''
This program is a Thermal Conduction Calculator
Units are in MKS (meters, kilograms, seconds), aka SI units
lamda should be in units of W/(m*K)
'''

# ------------------------------------------
# Finding Heat Load
'''
Steps:
 1. Reference the correct material from the materials dictionary.
    ex: materials['Aluminum']['1100']
 2. Select the correct model with .thermCond['modelname]
    make sure it has a .function
 2. call heat_load() with inputs: material, modelname, minT, maxT.
'''

# 100 Quad-Twist
heat_load(materials['Saved Objects']['General Materials']['Quad-Twist']['36 AWG'], 'fit', bounds=[4, 300], length=0.55, dL=0.1*0.55)


# ------------------------------------------
# Example section for referencing this code

# Examples:
# test = Mat_Class('test', {'matProp': 'thermCond', 'name': 'testfit', 'eq': lambda x: 3*(x**2), 'eqr': [0, 10e4], 'd_eq': lambda x: 0.1*x}, fullname='test material')
# heat_load(test, 'testfit', [0, 100], A=5e-3, dA=5e-5, L=1e-4, dL=5e-6)

'''
Things to Do:
 Balsa wood does not work
    - make a post_function field
 add plotting lamda ?
 add more data points to Saphire
 The eternal struggle to debug
'''
