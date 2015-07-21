# -*- coding: utf-8 -*-
# A Thermal Conduction Calculator

from Material_Database_Data import *
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

# Black
heat_load(materials['Saved Objects']['Johanna Dewar']['Vespel Legs']['Black'], 'fit', bounds=[0.25, 0.45], area=4.84444268598568e-6, dA=1.99491133502962e-8, length=21.1e-3, dL=0.1e-3)

# brown
heat_load(materials['Saved Objects']['Johanna Dewar']['Vespel Legs']['Brown'], 'fit', bounds=[0.45, 6.5], area=4.82449357263538e-6, dA=3.98982267005898e-9, length=20.2e-3, dL=0.1e-3)



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
