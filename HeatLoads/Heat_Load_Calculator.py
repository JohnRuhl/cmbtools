# -*- coding: utf-8 -*-
# A Thermal Conduction Calculator

from Heat_Load_Material_Data import *
# import Thermal_Conduction_Calculator_functions as fn  # (already imported in _data)

'''
This program is a Thermal Conduction Calculator
Units are in KgMS (kilograms, meters, seconds), aka SI units
lamda should be in units of W/(m*K)
'''

'''
Things to Do:
 make it so class can update an existing dictionary. Global variables maybe?
 Balsa wood does not work
 NIST's Specific Heat Method as a model option
 finish populating from NIST
 do errors and units  <= very difficult
'''

# ------------------------------------------
# Finding Heat Load
'''
Steps:
 1. Reference the correct material from the materials dictionary.
    ex: materials['Aluminum']['1100']
 2. Select the correct model with .lamda['modelname]
    make sure it has a .function
 2. call heat_load() with inputs: material, modelname, minT, maxT.
 3. print answer
'''

QuadTist = Mat_Class('QuadTist', {'mdl': 'constant', 'eq': lambda T: 38, 'eqr': [1, 300]}, fullname='QuadTist Wire')
answer = heat_load(QuadTist, 'constant', 77, 300, area=0.0000019949, length=0.0508)
print 'answer: {0} +/- {1}\n'.format(answer, 'TBA')

QuadTist = Mat_Class('QuadTist', {'mdl': 'constant', 'eq': lambda T: 15, 'eqr': [1, 300]}, fullname='QuadTist Wire')
answer = heat_load(QuadTist, 'constant', 4, 77, area=0.0000019949, length=0.0508)
print 'answer: {0} +/- {1}\n'.format(answer, 'TBA')

# ------------------------------------------
# Reference section for referencing this code

# Examples:
# materials['Aluminum']['1100'].add_area(0.005)
# materials['Aluminum']['1100'].add_length(2)
# answer = heat_load(materials['Aluminum']['1100'], 'log10', 4, 300)
# print 'answer: {0} +/- {1}\n'.format(answer, 'TBA')

# answer = heat_load(materials['Aluminum']['3003F'], 'log10', 4, 300, area=4, length=2)
# print 'answer: {0} +/- {1}\n'.format(answer, 'TBA')

# QuadTist = Mat_Class('QuadTist', {'mdl': 'constant', 'eq': lambda T: 38, 'eqr': [1, 300]}, fullname='QuadTist Wire')
# answer = heat_load(QuadTist, 'constant', 77, 300, area=0.0000019949, length=0.0508)
# print 'answer: {0} +/- {1}\n'.format(answer, 'TBA')

# Just some references
# print aluminum_dict['5083'].lamda['log10'].coeffs
# print matDict['Aluminum']['5083'].lamda['log10'].coeffs
# print matDict['Aluminum']['name']
