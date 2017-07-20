# -*- coding: utf-8 -*-
# ModelEquations

'''
All 'master' realizationEquations are designated with Signal and must only accept 2 inputs
This file holds all the functions for making the theoretical CMB Models
'''

import numpy as np
from Functions import npFloat


##########################################################################################
                                    # Total Signal

def totalSignal(inputs, a0, a1, a2, a3, a4):
    ''' totalSignal((lList, frequency), parameters)
        Creates the total theoretical signal
        Calls the BModeSignal and dustSignal functions
        Adds them
    '''
    # creating parameters
    params_sq, params_sin = parameterSplit([a0, a1, a2, a3, a4])
    # getting the inputs
    inputs_sq, inputs_sin = inputSplit(inputs)  # *** include freqs ***
    # getting functions
    poly = polySignal(inputs_sq, *params_sq)
    sine = sineSignal(inputs_sin, *params_sin)
    # adding into one signal
    signal = np.array(np.add(poly, sine))
    return signal


def totalError(inputs, a0, a1, a2, a3, a4):
    params_sq, params_sin = parameterSplit([a0, a1, a2, a3, a4])
    inputs_sq, inputs_sin = inputSplit(inputs)
    polyerror = polyError(inputs_sq, *params_sq)
    sineerror = sineError()

    error = np.array(polyerror + sineerror)
    return error


##########################################################################################
                                    # PolySignal

def polySignal(inputs, a0, a1, a2):
    ''' BModeSignal(BMode, parameters)
    '''
    parameters = [a0, a1, a2]
    lList, freq = inputs  # just making sure
    # making BMode model
    polysignal = np.sum([param*np.power(lList, power) for param, power in zip(parameters, range(len(parameters))[::-1])], axis=0)
    return np.array(polysignal)

def polyError(inputs, a0, a1, a2):
    return 0.02*polySignal(inputs, a0, a1, a2)

##########################################################################################
                                    # Sin Signal

def sineSignal(inputs, a0, a1):
    lList, freq = inputs
    # getting parameter
    ampl, period = a0, a1  # , dustT
    # making dust model
    sinesignal = ampl*np.sin(period*lList)*freq**1.59
    return np.array(sinesignal)

def sineError(*args):
    return 2.

##########################################################################################
                                    # Parameters

def parameterSplit(parameters):
    ''' This function splits the list of parameters into the BMode and Dust sections
    '''
    # poly function
    try:  # checking if already packaged
        len(parameters[0])
    except:  # just a list, so package
        params_sq = parameters[:-2]
    else:  # already packaged
        params_sq = parameters[0]
    finally:
        params_sq = np.array(params_sq)

    # Sin function
    try:
        len(parameters[1])
    except:
        params_sin = parameters[-2:]
    else:
        params_sin = parameters[1]
    finally:
        params_sin = np.array(params_sin)

    return params_sq, params_sin

def inputSplit(inputs):
    # poly function
    inputs_sq = inputs[:]   # *** right now is all inputs ***

    # sine function
    inputs_sin = inputs[:]   # *** right now is all inputs ***

    return inputs_sq, inputs_sin


def makeParamInputs(params, output="all"):
    params = npFloat(list(params))
    params_sq, params_sin = parameterSplit(params)

    if output == "all":
        return params, params_sq, params_sin
    else:
        return params


# x = np.arange(10)
# parameters = [1, 1, 1, 1, 1]
# print(totalSignal(x, *parameters))
