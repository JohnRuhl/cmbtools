# -*- coding: utf-8 -*-
# A Thermal Conduction Calculator

import scipy.integrate as integrate
import numpy as np

'''
This program is a Thermal Conduction Calculator
Units are in KgMS (kilograms, meters, seconds), aka SI units
lamda should be in units of W/(m*K)
'''

'''
Notes to self:
 make it so class can update an existing dictionary. Global variables maybe?'''

# ------------------------------------------
# Classes and Functions


# Material Class
class New_Mat(object):
    '''This creates a new material for thermal conduction calculations'''
    def __init__(self, name, *args, **kwargs):
        self.name = name
        # Lamda = Thermal Conductivity given as argument, kwarg, & und
        if len(args) > 0:                            # lamda
            self.lamda = self.Lamda(args[0])
        elif 'lamda' in kwargs:
            self.lamda = self.Lamda(kwargs.get('lamda'))
        else:
            self.lamda = self.Lamda(None)
        # cross sectional area given as argument, kwarg, & und
        if len(args) > 1:                            # cross sectional area
            self.A = args[1]
        elif 'area' in kwargs:
            self.A = kwargs.get('area')
        else:
            self.A = None
        # Length given as argument, kwarg, & und
        if len(args) > 2:                            # length of wire/tube
            self.L = args[2]
        elif 'length' in kwargs:
            self.L = kwargs.get('length')
        else:
            self.L = None
        if 'UNS' in kwargs:                          # scientific designation / call code
            self.UNS = kwargs.get('UNS')
        if 'fullname' in kwargs:                     # full name
            self.fullname = kwargs.get('fullname')
        if 'datarange' in kwargs:                    # data range should be in form [min max]
            self.datarange = kwargs.get('fullname')

    class Lamda(object):
        '''does lamda. forms: [constant/equation],
                              [equation, [coeffs]]
            units are W/(m*K)'''
        def __init__(self, arg):
            if arg == None:
                self.equation = None
                self.coeffs = None
                self.function = None
            elif len(arg) == 1:
                self.equation = arg[0]
                self.coeffs = None
                self.function = lambda T: arg[0](T)
                test = self.function(10)
                print test
            else:
                self.equation = arg[0]
                self.coeffs = arg[1]
                self.function = lambda T: arg[0](T, arg[1])

        # if 'dict' in kwargs:                 # make it so class can update an existing dictionary
        #    key = kwargs.get('dict')
        #    kwargs.get('dict') = make_dict(self, dict=kwargs.get('dict'))


# Lamda Models

# Log10 Log10 Curve Fit
def lamda_log10_model(T, coeffs=[]):
    '''this function makes a lamda based on NIST's best fit. note: log base 10
    y=10^(a+b*log(T)+c*(log(T))^2+d*(log(T))^3...
    http://cryogenics.nist.gov/MPropsMAY/materialproperties.htm
    '''
    if len(coeffs) == 0:
        lamda = None
        return lamda
    else:
        temp = 1
        for i, v in enumerate(coeffs):
            temp += v*np.power(np.log10(T), i)
        lamda = np.power(10, temp)
        return lamda


# Other Model
def lamda_test_model(T):
    lamda = T
    return lamda


# Integral Function
def lamda_integral(function, *args):
    '''This function will do the integral, from T1 to T2 of lamda dT.
    If only one T is given, it will go from 0 to T'''
    if len(args) == 1:
        integral = integrate.quad(function, 0.0, args[0])
    elif len(args) == 2:
        integral = integrate.quad(function, args[0], args[1])
    return integral


# Dictionary Creation
def make_dict(*args, **kwargs):
    '''This function will either make a dictionary or update an existing one'''
    if 'dict' in kwargs:
        new_dict = kwargs.get('dict')
    else:
        new_dict = {}
    if 'name' in kwargs:
        new_dict['name'] = kwargs.get('name')
    for i in range(len(args)):
        if isinstance(args[i], dict):
            if 'name' in args[i]:
                new_dict.update({args[i]['name']: args[i]})
            elif 'subname' in kwargs:
                new_dict.update({kwargs.get('subname'): args[i]})
        else:
            new_dict.update({args[i].name: args[i]})
    return new_dict


# ------------------------------------------
# Data
'''
currently, the lamda data needs to be enclosed in [].
if a nondefined function is required, use "lambda x: expression".
'''

# Aluminum Data
aluminum1100 = New_Mat('1100', [lamda_log10_model, [23.39172, -148.5733, 422.1917, -653.6664, 607.0402, -346.152, 118.4276, -22.2781, 1.770187]], UNS='A91100', datarange=[4, 300], fullname='Aluminum 1100')
aluminum3003F = New_Mat('3003F', [lamda_log10_model, [0.63736, -1.1437, 7.4624, -12.6905, 11.9165, -6.18721, 1.63939, -0.172667]], UNS='A93003', datarange=[4, 300], fullname='Aluminum 3003F')
aluminum5083 = New_Mat('5083', [lamda_log10_model, [-0.90933, 5.751, -11.112, 13.612, -9.3977, 3.6873, -0.77295, 0.067336]], UNS='A95083', datarange=[4, 300], fullname='Aluminum 5083')
aluminum6061T6 = New_Mat('6061-T6', [lamda_log10_model, [0.07918, 1.0957, -0.07277, 0.08084, 0.02803, -0.09464, 0.04179, -0.00571, 2.96344]], UNS='A96061', datarange=[4, 300], fullname='Aluminum 6061-T6')
aluminum6063T5 = New_Mat('6063-T5', [lamda_log10_model, [22.401433, -141.13433, 394.95461, -601.15377, 547.83202, -305.99691, 102.38656, -18.810237, 1.4576882]], UNS='A96063', datarange=[4, 300], fullname='Aluminum 6063-T5')
# Aluminum Dictionary
aluminum_dict = make_dict(aluminum1100, aluminum3003F, aluminum5083, aluminum6061T6, aluminum6063T5, name='Aluminum')

# Balsa Data
balsa6ro = New_Mat('6 lb/ft^3', [lamda_log10_model, [4172.447, -11309.97, 12745.09, -7647.584, 2577.309, 2577.309, -462.538, 34.5351]])
balsa11ro = New_Mat('11 lb/ft^3', [lamda_log10_model, [4815.4, -12969.63, 14520.76, -8654.164, 2895.712, -515.7272, 38.19218]])
# Balsa Dictionary
balsa_dict = make_dict(balsa6ro, balsa11ro, name='Balsa')


# Master Dictionary
materials = make_dict(aluminum_dict, balsa_dict, name='Master Dictionary')


# ------------------------------------------
# Finding Integral
'''
Steps:
 1. Reference the correct matrial from the materials dictionary
    ex: materials['Aluminum']['1100'].lamda.function
 2. call lamda_integral() with inputs: function, range. 
    range should either be 'max' (assumed min of 0) or 'min, max'
    ex: ..., 50) or ...,5 , 50) 
'''

answer = lamda_integral(materials['Aluminum']['1100'].lamda.function, 4, 300)
print 'answer: {0} +/- {1}\n'.format(answer[0], answer[1])


# ------------------------------------------
# Reference section for referencing this code

# print aluminum_dict['5083'].lamda.coeffs
# print matDict['Aluminum']['5083'].lamda.coeffs
# print matDict['Aluminum']['name']
