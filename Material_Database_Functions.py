# -*- coding: utf-8 -*-
# Thermal Conduction Calculator Functions

import numpy as np
import scipy.integrate as integrate
from types import FunctionType, NoneType
from scipy.interpolate import UnivariateSpline


class Mat_Class(object):
    ''' Creates a new material for thermal conduction calculations
        must give a name
        can give a material property (ex thermCond), area, length, UNS, fullname, and info
        - all model properties such as thermCond must be args
        - all non-model properties such as dimensions, and designations must be kwargs
        - Model Properties:
            - must be an arg
            - must be a dict
            - must have a material property types with key 'matProp'
            - must have a modelname with key 'name'
            - may have an equation with key 'eq'
            - may have eqinput with key 'eqinput'
              - if eqinput is coeffs: [coeffs]
              - if eqinput is data: [[temps], [data]]
              - if eqinput is mixed: [[coeffs], [[temps],[data]],[[coeffs_range],[data_range]]]
            - may have equation range with key 'eqr'
            - may have equation error with key 'egerror' or 'd_eq'
            ex: {'name': 'modelname', 'eq': equation, 'eqinput': [1,3,5], 'eqr': [4,200]}

        - area & length must be numbers
        - UNS, fullname, and info must be strings
    '''
    def __init__(self, name, *args, **kwargs):
        super(Mat_Class, self).__init__()
        self.name = name
        self.fullname = name

        # thermCond.
        for i in range(len(args)):
            Dict = args[i]
            self.add_model(Dict)

        # Non-Model Properties
        # cross sectional area given as argument, kwarg, & und
        if 'area' in kwargs:                       # cross sectional area
            self.area = kwargs.get('area')
        else:
            self.area = None
        if 'd_area' in kwargs:
            self.area_error = kwargs.get('d_area')
        else:
            self.area_error = None
        # Length given as argument, kwarg, & und
        if 'length' in kwargs:                       # length of wire/tube
            self.length = kwargs.get('length')
        else:
            self.length = None
        if 'd_length' in kwargs:
            self.length_error = kwargs.get('d_length')
        else:
            self.length_error = None
        if 'UNS' in kwargs:                          # scientific designation / call code
            self.UNS = kwargs.get('UNS')
        else:
            self.UNS = None
        if 'fullname' in kwargs:                     # full name
            self.fullname = kwargs.get('fullname')
        if 'info' in kwargs:                         # info is any note. should be str
            self.info = kwargs.get('datarange')
        else:
            self.info = None

    def add_area(self, area=None, **kwargs):
        'adds area to a given model. if no area passed, adds as None'
        if 'area' in kwargs:
            self.area = kwargs.get('area')
        else:
            self.area = area

    def add_area_error(self, d_area=None, **kwargs):
        'adds error to the area. if no area passed, adds as None'
        if 'd_area' in kwargs:
            self.area_error = kwargs.get('d_area')
        else:
            self.d_area = d_area

    def add_length(self, length=None, **kwargs):
        'adds length to a given model. if no area passed, adds as None'
        if 'length' in kwargs:
            self.length = kwargs.get('length')
        else:
            self.length = length

    def add_length_error(self, d_length=None, **kwargs):
        'adds length to a given model. if no length passed, adds as None'
        if 'd_length' in kwargs:
            self.length_error = kwargs.get('d_length')
        else:
            self.d_length = d_length

    def add_info(self, info=None, **kwargs):
        'adds info to a given model. if no info passed, adds as None'
        if 'info' in kwargs:
            self.info = kwargs.get('info')
        else:
            self.info = info

    def add_UNS(self, UNS=None, **kwargs):
        'adds UNS to a given model. if no UNS passed, adds as None'
        if 'UNS' in kwargs:
            self.UNS = kwargs.get('UNS')
        else:
            self.UNS = UNS

    def add_fullname(self, fullname=None, **kwargs):
        'adds fullname to a given material. if no fullname passed, adds as None'
        if 'fullname' in kwargs:
            self.fullname = kwargs.get('fullname')
        else:
            self.fullname = fullname

    def add_anything(self, **kwargs):
        if 'area' in kwargs:
            self.add_area(kwargs.get('area'))
        if 'd_area' in kwargs:
            self.add_area_error(kwargs.get('d_area'))
        if 'length' in kwargs:
            self.add_length(kwargs.get('length'))
        if 'd_length' in kwargs:
            self.add_length_error(kwargs.get('d_length'))
        if 'info' in kwargs:
            self.add_info(kwargs.get('info'))
        if 'UNS' in kwargs:
            self.add_UNS(kwargs.get('UNS'))
        if 'fullname' in kwargs:
            self.add_fullname(kwargs.get('fullname'))
        if 'model' in kwargs:
            self.add_model(kwargs.get('model'))
        if 'eqrange' or 'eqr' in kwargs:
            if 'eqrange' in kwargs:
                eqrange = kwargs.get('eqrange')
            elif 'eqr' in kwargs:
                eqrange = kwargs.get('eqr')
            matProp = kwargs.get('matProp')
            modelname = kwargs.get('name')
            self.add_eqrange(matProp, modelname, eqrange)
        if 'equation' or 'eq' in kwargs:
            if 'equation' in kwargs:
                equation = kwargs.get('equation')
            elif 'eq' in kwargs:
                equation = kwargs.get('eq')
            matProp = kwargs.get('matProp')
            modelname = kwargs.get('name')
            self.add_equation(matProp, modelname, equation)
        if 'eqinput' or 'eqin' in kwargs:
            if 'eqinput' in kwargs:
                eqinput = kwargs.get('eqinput')
            elif 'eqin' in kwargs:
                eqinput = kwargs.get('eqin')
            matProp = kwargs.get('matProp')
            modelname = kwargs.get('name')
            self.add_eqinput(matProp, modelname, eqinput)
        if 'eqerror' or 'd_eq' in kwargs:
            if 'eqerror' in kwargs:
                eqerror = kwargs.get('eqerror')
            elif 'd_eq' in kwargs:
                eqerror = kwargs.get('d_eq')
            matProp = kwargs.get('matProp')
            modelname = kwargs.get('name')
            self.add_eqerror(matProp, modelname, eqerror)

    def add_model(self, Dict, **kwargs):
        ''' makes material.thermCond[modelname] as a class.
            needs: matProp and modelname (as kwarg or in main dict)
            can be given any of: equation, eqinput, & eqrange
                - in main dict or in kwarg
            if 'function' in dict or kwargs, then function will be made with *only* given inputs
                - needs an equation
            equation must be a function (lambda functions work)
            ex: add_model({'matProp': 'thermCond', 'name'= 'NIST', 'eq': lambda x: x, 'eqinput': [1, 23]})
            add_model supports: equation as eq or equation in dict or kwarg
                                eqinput as eqinput or eqinput in dict or kwarg
                                eqrange as eqr or eqrange in dict or kwarg
            kwarg has highest priority, then arg, then dict entries
        '''

        # detects if model already exists, and just adds to it, or..
            # getattr(self, Dict['model'])
        # add model
            # setattr(self, Dict['model'], make_dict(name=Dict['model']))

        if not isinstance(Dict, (dict, NoneType)):
            if Dict is None:
                Dict = {}
            else:
                print 'Error: dict needs to be type dict or None'
                return

        if 'matProp' in kwargs:
            matProp_name = kwargs.get('matProp')
        elif 'matProp' in Dict:
            matProp_name = Dict['matProp']
        else:          #CHANGE THIS SO can be material[X].thermCond.add_model()  # is this issue resolved?
            print 'Error: need to pass a matProp'
            return

        if 'name' in kwargs:
            modelname = kwargs.get('name')
        elif 'name' in Dict:
            modelname = Dict['name']
        else:
            print 'Error: need to pass a modelname'

        if not isinstance(matProp_name, str):
            print 'Error: matProp must be str'
            return
        elif not isinstance(modelname, str):
            print 'Error: modelname must be str'
            return
        else:
            try:
                matProp = getattr(self, matProp_name)
            except:
                setattr(self, matProp_name, make_dict(name=matProp_name))
                matProp = getattr(self, matProp_name)
            matProp[modelname] = self.Add_thermCond_Model(modelname)

        # Equation
        if 'equation' in kwargs:
            equation = kwargs.get('equation')
        elif 'eq' in kwargs:
            equation = kwargs.get('eq')
        elif 'equation' in Dict:
            equation = Dict['equation']
        elif 'eq' in Dict:
            equation = Dict['eq']
        else:
            equation = None
        if isinstance(equation, (FunctionType, NoneType)):
            self.add_equation(matProp_name, modelname, equation)
        else:
            print 'Error: equation nees to be type function or None'
            return

        # eqinput
        if 'eqinput' in kwargs:
            eqinput = kwargs.get('eqinput')
        elif 'eqi' in kwargs:
            eqinput = kwargs.get('eqi')
        elif 'eqinput' in Dict:
            eqinput = Dict['eqinput']
        elif 'eqi' in Dict:
            eqinput = Dict['eqi']
        else:
            eqinput = None
        if isinstance(eqinput, (list, NoneType)):
            self.add_eqinput(matProp_name, modelname, eqinput)
        else:
            print 'Error: eqinput nees to be type list or None'
            return

        # Equation Range
        if 'eqrange' in kwargs:
            eqrange = kwargs.get('eqrange')
        elif 'eqr' in kwargs:
            eqrange = kwargs.get('eqr')
        elif 'eqrange' in Dict:
            eqrange = Dict['eqrange']
        elif 'eqr' in Dict:
            eqrange = Dict['eqr']
        else:
            eqrange = None
        if isinstance(eqrange, (list, NoneType)):
            self.add_eqrange(matProp_name, modelname, eqrange)
        else:
            print 'Error: eqrange nees to be type list or None'
            return

        # Eqerror
        if 'eqerror' in kwargs:
            eqerror = kwargs.get('eqerror')
        elif 'd_eq' in kwargs:
            eqerror = kwargs.get('d_eq')
        elif 'eqerror' in Dict:
            eqerror = Dict['error']
        elif 'd_eq' in Dict:
            eqerror = Dict['d_eq']
        else:
            eqerror = None
        if isinstance(eqerror, (int, float, list, FunctionType, NoneType)):
            self.add_eqerror(matProp_name, modelname, eqerror)
        else:
            print 'Error: eqerror nees to be type number, list, function or None'
            return

    class Add_thermCond_Model(object):
        '''Add_thermCond_Model adds a model to thermCond given a modelname'''
        def __init__(self, modelname):
            self.modelname = modelname
            self.eq = None
            self.eqinput = None
            self.eqrange = None
            self.eqerror = None

    def add_eqrange(self, *args, **kwargs):
        ''' adds an equation range to a specified model.
            called as material.add_eqrange(property, modelname, eqrange)  (with either as a kwarg)
            Cannot give something as both an arg and a kwarg!
        '''
        arg_count = 0
        if 'matProp' in kwargs:
            matProp = getattr(self, kwargs.get('matProp'))
        else:
            matProp = getattr(self, args[arg_count])
            arg_count += 1
        if not isinstance(matProp, dict):
            print 'matProp needs to be a dict'
            return
        if 'name' in kwargs:
            modelname = kwargs.get('name')
        else:
            modelname = args[arg_count]
            arg_count += 1
        if not isinstance(modelname, str):
            print 'modelname needs to be a str'
            # return
        if 'eqrange' in kwargs:
            eqrange = kwargs.get('eqrange')
        else:
            eqrange = args[arg_count]
        matProp[modelname].eqrange = eqrange

    def add_equation(self, *args, **kwargs):
        ''' adds an equation to a model.
            called material.add_equation(property, modelname, equation)  (with either as a kwarg)
            note:  modelname must be a string
            Cannot give something as both an arg and a kwarg!
        '''
        arg_count = 0
        if 'matProp' in kwargs:
            matProp = getattr(self, kwargs.get('matProp'))
        else:
            matProp = getattr(self, args[arg_count])
            arg_count += 1
        if not isinstance(matProp, dict):
            print 'matProp needs to be a dict'
            return
        if 'name' in kwargs:
            modelname = kwargs.get('name')
        else:
            modelname = args[arg_count]
            arg_count += 1
        if not isinstance(modelname, str):
            print 'modelname needs to be a str'
            return
        if 'equation' in kwargs:
            equation = kwargs.get('equation')
        elif 'eq' in kwargs:
            equation = kwargs.get('eq')
        else:
            equation = args[arg_count]
        matProp[modelname].eq = equation

    def add_eqinput(self, *args, **kwargs):
        ''' adds equation input to a model.
            called as material.add_eqinput(property, modelname, eqinput)  (with either as a kwarg)
            Cannot give something as both an arg and a kwarg!
        '''
        arg_count = 0
        if 'matProp' in kwargs:
            matProp = getattr(self, kwargs.get('matProp'))
        else:
            matProp = getattr(self, args[arg_count])
            arg_count += 1
        if not isinstance(matProp, dict):
            print 'matProp needs to be a dict'
            return
        if 'name' in kwargs:
            modelname = kwargs.get('name')
        else:
            modelname = args[arg_count]
            arg_count += 1
        if not isinstance(modelname, str):
            print 'modelname needs to be a str'
            # return
        if 'eqinput' in kwargs:
            eqinput = kwargs.get('eqinput')
        elif 'eqin' in kwargs:
            eqinput = kwargs.get('eqin')
        else:
            eqinput = args[arg_count]
        try:                                                    # eqinput is raw data for a fit
            eqinput = np.array(eqinput[0][:], eqinput[1][:])
        except:                                                 # eqinput is coeffs for an equation
            pass
        matProp[modelname].eqinput = eqinput

    def add_eqerror(self, *args, **kwargs):
        ''' adds an equation error to a model.
            called as material.add_eqerror(property, modelname, eqerror)  (with either as a kwarg)
            Cannot give something as both an arg and a kwarg!
        '''
        arg_count = 0
        if 'matProp' in kwargs:
            matProp = getattr(self, kwargs.get('matProp'))
        else:
            matProp = getattr(self, args[arg_count])
            arg_count += 1
        if not isinstance(matProp, dict):
            print 'matProp needs to be a dict'
            return
        if 'name' in kwargs:
            modelname = kwargs.get('name')
        else:
            modelname = args[arg_count]
            arg_count += 1
        if not isinstance(modelname, str):
            print 'modelname needs to be a str'
            return
        if 'eqerror' in kwargs:
            eqerror = kwargs.get('eqerror')
        elif 'd_eq' in kwargs:
            eqerror = kwargs.get('d_eq')
        else:
            eqerror = args[arg_count]
        matProp[modelname].eqerror = eqerror


# Dictionary Creation
def make_dict(*args, **kwargs):
    '''
        This function will either make a dictionary or update an existing one
    '''
    if 'dict' in kwargs:                                      # gets an existing dictionary
        new_dict = kwargs.get('dict')
    else:                                                     # makes a dictionary
        new_dict = {}
    if 'name' in kwargs:                                      # names the dictionary
        new_dict['name'] = kwargs.get('name')
    for i in range(len(args)):                                # gets an existing dictionary
        if isinstance(args[i], dict):
            if 'name' in args[i]:                             # adds a dictionary to the dictionary
                new_dict.update({args[i]['name']: args[i]})
            else:
                print 'operation failed'
                return
        else:                                                 # adds a class to the dictionary
            new_dict.update({args[i].name: args[i]})
    return new_dict


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


####################################################################################
# THERMAL CONDUCTIVITY

def heat_load(*args, **kwargs):
    ''' Calculates heat load with equation:   thermCond * area * deltaT / length
        average thermCond is found with integral, from T1 to T2 of thermCond dT, and divide by the temp range
        Must be given (material, modelname, bounds) in that order
            - any can be a kwarg, and so given with all the other kwargs
            - Do not give something as both an arg and a kwarg!
        May be given area and/or length as kwargs.
            - can be called as either area or A, length or L
        May be given error in area and length as kwargs.
            - can be called as either area_error or dA, length_error or dL
            - can be an upper error and lower error in form [lower, upper]
    '''
    arg_count = 0

    if 'material' in kwargs:
        material = kwargs.get('material')
    else:
        material = args[0]
        arg_count += 1
    if 'name' in kwargs:
        modelname = kwargs.get('name')
    else:
        modelname = args[arg_count]
        arg_count += 1
    if 'bounds' in kwargs:
        bounds = kwargs.get('bounds')
    else:
        bounds = args[arg_count]
        arg_count += 1
    # Area
    if 'area' in kwargs:
        area = kwargs.get('area')
    elif 'A' in kwargs:
        area = kwargs.get('A')
    else:
        area = material.area
    # Area Error
    if 'area_error' in kwargs:
        d_area = kwargs.get('area_error')
    elif 'd_area' in kwargs:
        d_area = kwargs.get('d_area')
    elif 'dA' in kwargs:
        d_area = kwargs.get('dA')
    else:
        d_area = material.area_error
    # length
    if 'length' in kwargs:
        length = kwargs.get('length')
    elif 'L' in kwargs:
        length = kwargs.get('L')
    else:
        length = material.length
    # length error
    if 'length_error' in kwargs:
        d_length = kwargs.get('length_error')
    elif 'd_length' in kwargs:
        d_length = kwargs.get('d_length')
    elif 'dL' in kwargs:
        d_length = kwargs.get('dL')
    else:
        d_length = material.length_error

    if material is None or modelname is None or bounds is None:
        print 'Error: material, modelname, &/or bounds can"t be None'
        return

    print '\nHeat Load for {0}:'.format(material.fullname)

    model = material.thermCond[modelname]
    integral, abserror = thermCond_integral(model, bounds)
    thermCond = np.multiply(np.divide(area, length), integral)

    d_thermCond = thermCond_error(model, bounds, integral, abserror, area, d_area, length, d_length)

    print 'average thermCond over {0}:{1} : {2}'.format(bounds[0], bounds[1], integral/(bounds[1]-bounds[0]))
    if d_thermCond[0] == d_thermCond[1]:
        print 'Heat Load: {0} +/- {1}\n'.format(thermCond, d_thermCond[0])
    else:
        print 'Heat Load: {0} +{1}/-{2}\n'.format(thermCond, d_thermCond[0], d_thermCond[1])
    return thermCond


# Integral Function
def thermCond_integral(model, bounds):
    '''
        This function will do the integral, from T1 to T2 of thermCond dT, and divide by the temp range
    '''
    try:
        if model.eqrange is None:
            print 'Warning: no preset equation range'
        if bounds[1] > model.eqrange[1] or bounds[0] < model.eqrange[0]:
            if bounds[1] > model.eqrange[1]:
                print 'Warning: upper bound too large.\n       user input was:    {0}\
                       \n       max is:            {1}'.format(bounds[1], model.eqrange[1])
            if bounds[0] < model.eqrange[0]:
                print 'Warning: lower bound too small.\n       user input was: {0}\
                       \n       min is:         {1}'.format(bounds[0], model.eqrange[0])
            print '... Calculating anyway\n'
    except:
        print 'warning: no eqrange'
    if model.eqinput is None:
        integral, abserror = integrate.quad(model.eq, bounds[0], bounds[1])
    else:
        integral, abserror = integrate.quad(model.eq, bounds[0], bounds[1], args=(model.eqinput))
    return integral, abserror


def thermCond_error(model, bounds, integral, abserror, area, d_area, length, d_length):
    '''
        This fucntion will numerically find the error in the heat_load from T1 to T2.
    '''
    # prepping d_area and d_length
    if d_area is None:
        d_area = [0]
    elif isinstance(d_area, (int, float)):
        d_area = [d_area]
    if d_length is None:
        d_length = [0]
    elif isinstance(d_length, (int, float)):
        d_length = [d_length]
    # prepping eqerror
    try:                             # eqerror a list
        eqerror = model.eqerror[:]
    except:                          # eqerror not a list
        if model.eqerror is None:
            print 'Warning: no eqerror. Assuming 0'
            eqerror = [0]
        else:
            eqerror = [model.eqerror]

    # making lists to iterate through
    temp_errors = [0, 0]
    errors = [d_area, d_length, [abserror]*2]
    error_constants = [np.divide(integral, length),
                       np.divide(np.multiply(area, integral), np.square(length)),
                       np.divide(area, length)]

    # comparing eqerror to abserror
    for i, val in enumerate(eqerror):
        if isinstance(val, FunctionType):       # eqerror[i] is a f(model.fn)
            temp = val(integral)
        elif isinstance(val, (int, float)):     # eqerror[i] is constant
            temp, uncertainty = integrate.quad(lambda T: val,
                                               bounds[0], bounds[1])
            temp += uncertainty
        else:
            return
        if temp > errors[2][i]:  # comparing to abserror to see which is bigger
            errors[2][i] = temp
    if i == 0:
        errors[2][1] = errors[2][0]
    # adding the squares of the error
    for i, val in enumerate(errors):
        for j, wal in enumerate(val):
            temp_errors[j] += np.square(np.multiply(wal, error_constants[i]))
        if j == 0:
            temp_errors[1] += np.square(np.multiply(wal, error_constants[i]))
    # taking the square root to get complete quadrature
    answer = np.sqrt(temp_errors)
    return answer


# ----------------------------------------------------------------------------------
# thermCond Models


# NIST Log10 Curve Fit
def thermCond_NIST_model(T, coeffs=[]):
    ''' this function makes a thermCond based on NIST's best fit. note: log base 10
        y=10^(a+b*log(T)+c*(log(T))^2+d*(log(T))^3...
        http://cryogenics.nist.gov/MPropsMAY/materialproperties.htm
        Note: right now, [] or None is = to T
    '''
    if coeffs is None:
        thermCond = T
        return thermCond
    elif len(coeffs) == 0:
        thermCond = T
        return thermCond
    else:
        logY = 0
        for i, v in enumerate(coeffs):
            logY += v*np.power(np.log10(T), i)
        thermCond = np.power(10, logY)
        return thermCond


# log10(NIST Log10 Curve Fit)
def thermCond_log10NIST_model(T, coeffs=[]):
    '''
        this returns log10(thermCond_log10_model)
    '''
    if len(coeffs) == 0:
        thermCond = None
        return thermCond
    else:
        logY = 1
        for i, v in enumerate(coeffs):
            logY += v*np.power(np.log10(T), i)
        log10_thermCond = logY
        return log10_thermCond


# Fitting Model
def thermCond_fit_UnivariateSpline_model(T, data=[]):
    # This is a fitting model for creating an equation when given a list of x,y values
    x = data[0]
    y = data[1]
    if len(x) <= 5:
        order = len(x)-1
    else:
        order = 5
    fit = UnivariateSpline(x, y, k=order)
    thermCond = fit.__call__(T)
    return thermCond


# Mixed NIST and UnivariateSpline
def thermCond_NIST_and_UnivariateSpline_model(T, eqinput=[]):
    coeffs = eqinput[0]                    # coeffs for NIST model
    data = eqinput[1]                      # data for Spline model (3 not included)
    cfs_range = eqinput[2][0]
    fit_range = eqinput[2][1]

    N = thermCond_NIST_model(T, coeffs)
    if T < cfs_range[0]:
        NWeight = np.exp(T-cfs_range[0])
    elif T > cfs_range[1]:
        NWeight = np.exp(cfs_range[0]-T)
    else:
        NWeight = 1
    F = thermCond_fit_UnivariateSpline_model(T, data)
    if T < fit_range[0]:
        FWeight = np.exp(T-fit_range[0])
    elif T > fit_range[1]:
        FWeight = np.exp(fit_range[0]-T)
    else:
        FWeight = 1
    thermCond = (N*NWeight + F*FWeight)/(NWeight+FWeight)
    return thermCond
