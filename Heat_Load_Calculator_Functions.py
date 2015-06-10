# -*- coding: utf-8 -*-
# Thermal Conduction Calculator Functions

import scipy.integrate as integrate
import numpy as np
import types

'''
   make lamda a dictionary with name lamda. each entry should be a fit model (ex log10).
   each entry is a class with .equation, .coeffs, .function
   b/c right now, it only does one model
'''


class Mat_Class(object):
    ''' Creates a new material for thermal conduction calculations
        must give a name
        can give lamda, area, length, UNS, fullname, and info
        - lamda:
            - must be a dict
            - must have a modelname with key 'mdl'
            - must have an equation with key 'eq'
            - may have coeffs with key 'cfs'
            - must have equation range with key 'eqr'
            ex: {'mdl': 'modelname', 'eq': equation, 'cfs': [1,3,5], 'eqr': [4,200]}
        - area & length must be numbers
        - UNS, fullname, and info must be strings
    '''
    def __init__(self, name, *args, **kwargs):
        super(Mat_Class, self).__init__()
        self.name = name
        self.lamda = make_dict(name='lamda')

        # Lamda. fill out later
        if len(args) > 0:                            # lamda
            d = args[0]
            if 'cfs' in d:
                self.add_model(d['mdl'], d['eq'], d['cfs'], eqrange=d['eqr'], function=True)
            else:
                self.add_model(d['mdl'], d['eq'], eqrange=d['eqr'], function=True)
        elif 'lamda' in kwargs:
            d = kwargs.get('lamda')
            if 'cfs' in d:
                self.add_model(d['mdl'], d['eq'], d['cfs'], eqrange=d['eqr'], function=True)
            else:
                self.add_model(d['mdl'], d['eq'], eqrange=d['eqr'], function=True)

        # cross sectional area given as argument, kwarg, & und
        if len(args) > 1:                            # cross sectional area
            self.area = args[1]
        elif 'area' in kwargs:
            self.area = kwargs.get('area')
        else:
            self.area = None

        # Length given as argument, kwarg, & und
        if len(args) > 2:                            # length of wire/tube
            self.length = args[2]
        elif 'length' in kwargs:
            self.length = kwargs.get('length')
        else:
            self.length = None

        if 'UNS' in kwargs:                          # scientific designation / call code
            self.UNS = kwargs.get('UNS')
        else:
            self.UNS = None

        if 'fullname' in kwargs:                     # full name
            self.fullname = kwargs.get('fullname')
        else:
            self.fullname = None

        if 'info' in kwargs:                         # info is any note. should be str
            self.info = kwargs.get('datarange')
        else:
            self.info = None

    # def add_lamda_complete(self, modelname, equation, eqrange):
    #     self.lamda[modelname] = self.Make_Lamda_Complete(equation, eqrange)

    def add_area(self, *args, **kwargs):
        if 'area' in kwargs:
            self.area = kwargs.get('area')
        elif len(args) > 0:
            self.area = args[0]
        else:
            print 'need to pass an argument'

    def add_length(self, *args, **kwargs):
        if 'length' in kwargs:
            self.length = kwargs.get('length')
        elif len(args) > 0:
            self.length = args[0]
        else:
            print 'need to pass an argument'

    def add_model(self, *args, **kwargs):
        ''' makes material.lamda[modelname] as a class.
            needs: modelname (as 1st argument or kwarg)
            can be given any of: equation, coeffs, & eqrange
                - as either arg or kwarg, except eqrange
                - eqrange must be kwarg
                - preserve order, coeffs never before eq when both args
            if 'function' in kwargs given, then function will be made with *only* given inputs
                - needs an equation
        '''
        if 'modelname' in kwargs:                             # modelname kwarg
            modelname = kwargs.get('modelname')
        else:                                                 # modelname args[0]
            modelname = args[0]
        self.lamda[modelname] = self.Add_Lamda_Model(modelname)
        if 'equation' or 'coeffs' in kwargs:                  # equation or coeff kwarg
            if 'equation' in kwargs:
                equation = kwargs.get('equation')
                self.add_equation(modelname, equation)
            if 'coeffs' in kwargs:
                coeffs = kwargs.get('coeffs')
                self.add_coeffs(modelname, coeffs)
        if len(args) > 0:                                     # equation or coeff arg
            if isinstance(args[0], types.FunctionType):
                equation = args[0]
                self.add_equation(modelname, equation)
            elif isinstance(args[0], list):
                coeffs = args[0]
                self.add_coeffs(modelname, coeffs)
        if len(args) > 1:                                     # equation or coeff arg
            if isinstance(args[1], types.FunctionType):
                equation = args[1]
                self.add_equation(modelname, equation)
            elif isinstance(args[1], list):
                coeffs = args[1]
                self.add_coeffs(modelname, coeffs)
        if len(args) > 2:                                     # coeff arg
            coeffs =  args[2]
            self.add_coeffs(modelname, coeffs)
        if 'eqrange' in kwargs:                               # eqrange kwarg
            self.add_eqrange(modelname, kwargs.get('eqrange'))
        if 'function' in kwargs:
            if 'coeffs' in locals():
                self.add_function(modelname, equation, coeffs)
            else:
                self.add_function(modelname, equation)

    class Add_Lamda_Model(object):
        '''Add_Lamda_Model adds a model to lamda given a modelname'''
        def __init__(self, modelname):
            self.modelname = modelname
            self.equation = None
            self.coeffs = None
            self.eqrange = None
            self.fn = None

    def add_eqrange(self, *args, **kwargs):
        ''' adds an equation range to a model.
            called as material.add_eqrange(modelname, eqrange)  (with either as a kwarg)
        '''
        if 'eqrange' in kwargs:                                      # if eqrange kwarg
            eqrange = kwargs.get('eqrange')
        else:                                                         # if eqrange arg
            if isinstance(args[0], str):                                 # if modelname arg, eq args[1]
                eqrange = args[1]
            else:                                                        # if modelname not arg, eq args[0]
                eqrange = args[0]
        if 'modelname' in kwargs:                                     # if modelname kwarg
            self.lamda[kwargs.get('modelname')].eqrange = eqrange
        else:                                                         # if modelname arg at args[0]
            self.lamda[args[0]].eqrange = eqrange

    def add_equation(self, *args, **kwargs):
        ''' adds an equation to a model.
            called material.add_equation(modelname, equation)  (with either as a kwarg)
            note:  modelname must be a string
        '''
        if 'equation' in kwargs:                                      # if equation kwarg
            equation = kwargs.get('equation')
        else:                                                         # if equation arg
            if isinstance(args[0], str):                                 # if modelname arg, eq args[1]
                equation = args[1]
            else:                                                        # if modelname not arg, eq args[0]
                equation = args[0]
        if 'modelname' in kwargs:                                     # if modelname kwarg
            self.lamda[kwargs.get('modelname')].eq = equation
        elif isinstance(args[0], str):                                # if modelname arg
            self.lamda[args[0]].eq = equation

    def add_coeffs(self, * args, **kwargs):
        ''' adds an coefficients to a model.
            called as material.add_coeffs(modelname, coeffs)  (with either as a kwarg)
        '''
        if 'coeffs' in kwargs:                                      # if coeffs kwarg
            coeffs = kwargs.get('coeffs')
        else:                                                         # if coeffs arg
            if isinstance(args[0], str):                                 # if modelname arg, eq args[1]
                coeffs = args[1]
            else:                                                        # if modelname not arg, eq args[0]
                coeffs = args[0]
        if 'modelname' in kwargs:                                     # if modelname kwarg
            self.lamda[kwargs.get('modelname')].coeffs = coeffs
        elif isinstance(args[0], str):                                # if modelname arg
            self.lamda[args[0]].coeffs = coeffs

    def add_function(self, *args, **kwargs):
        ''' adds a callable function to the model.
            called as material.add_function(modelname, equation, coeffs)
                - any can be a kwarg.
                - if no eq or coeffs given, taken from material.lamda[modelname]
                - modelname must be a str
        '''
        if 'equation' in kwargs:                                       # if equation kwarg
            equation = kwargs.get('equation')
        # option of no equation passed
        if 'coeffs' in kwargs:
            coeffs = kwargs.get('coeffs')
        # option of no coeffs passed

        if len(args) > 0:
            if isinstance(args[0], types.FunctionType):                  # if equation arg at args[0]
                equation = args[0]
            elif isinstance(args[0], list):                                # if coeffs arg at args[0]
                coeffs = args[0]
        if len(args) > 1:
            if isinstance(args[1], types.FunctionType):                  # if equation arg at args[1]
                equation = args[1]
            elif isinstance(args[1], list):                                # if coeffs arg at args[1]
                coeffs = args[1]
        if len(args) > 2:
            if isinstance(args[2], list):                                # if coeffs arg at args[2]
                coeffs = args[2]

        if 'modelname' in kwargs:                                      # if modelname kwarg
            modelname = kwargs.get('modelname')
        else:                                                          # if modelname arg
            modelname = args[0]
        if 'equation' not in locals():                                 # if no equation given
            equation = self.lamda[modelname].eq
        if 'coeffs' not in locals():                                   # if no coeffs given
            coeffs = self.lamda[modelname].coeffs
        if coeffs is None:                                             # if there are no coeffs
            self.lamda[modelname].fn = lambda T: equation(T)
        else:                                                          # if there are coeffs
            self.lamda[modelname].fn = lambda T: equation(T, coeffs)


# Dictionary Creation
def make_dict(*args, **kwargs):
    '''
        This function will either make a dictionary or update an existing one
    '''
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


# Lamda Models

# Log10 Log10 Curve Fit
def lamda_log10_model(T, coeffs=[]):
    ''' this function makes a lamda based on NIST's best fit. note: log base 10
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
def lamda_integral(model, lower, upper):
    ''' This function will do the integral, from T1 to T2 of lamda dT, and divide by the temp range
    '''
    if upper > model.eqrange[1] or lower < model.eqrange[0]:
        if upper > model.eqrange[1] and lower < model.eqrange[0]:
            print 'Error: max too large and min too small.\
                   \n {0}input range was: {1}{2}:{3}\
                   \n {0}max range is:    {1}{4}:{5}'\
                   .format(' '*6, ' '*10, lower, upper, model.eqrange[0], model.eqrange[1])
        else:
            if upper > model.eqrange[1]:
                print 'Error: upper bound too large.\n       user input was:    {0}\
                       \n       max is:            {1}'.format(upper, model.eqrange[1])
            if lower < model.eqrange[0]:
                print 'Error: lower bound too small.\n       user input was: {0}\
                       \n       min is:         {1}'.format(lower, model.eqrange[1])
        return None, None
    else:
        integral, abserror = integrate.quad(model.fn, lower, upper)
        lamda = integral/(upper-lower)
        lamda_error = abserror/(upper-lower)
        return lamda, lamda_error


def heat_load(material, modelname, lowerT, upperT, **kwargs):
    ''' Calculates heat load with equation:   lamda * area * deltaT / length
        average lamda is found with integral, from T1 to T2 of lamda dT, and divide by the temp range
    '''
    if 'area' in kwargs:
        area = kwargs.get('area')
    else:
        area = material.area
    if 'length' in kwargs:
        length = kwargs.get('length')
    else:
        length = material.length
    lamda, lamda_error = lamda_integral(material.lamda[modelname], lowerT, upperT)
    print 'average lamda over {}:{} : {}'.format(lowerT, upperT, lamda)
    answer = np.divide(np.multiply(np.multiply(lamda, area), (upperT-lowerT)), length)
    return answer

#test = Mat_Class('test', area=2, length=4)
#print test.lamda
#test.add_model('log10')
#print test.lamda['log10'].modelname
#test.add_eqrange('log10', [30, 40])
#print test.lamda['log10'].eqrange
