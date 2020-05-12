# -*- coding: utf-8 -*-
# ModelEquations

'''
All 'master' ModelEquations are designated with Signal and must only accept 2 inputs
This file holds all the functions for making the theoretical CMB Models
'''

import numpy as np
from Functions import npFloat


##########################################################################################
                                    # Total Signal

def totalSignal(inputs, a0, a1, a2, a3, a4, a5):
    ''' totalSignal(parameters, (lList, frequency))
        Creates the total theoretical signal
        Calls the BModeSignal and dustSignal functions
        Adds them
    '''
    # creating parameters
    params_sq, params_sin = parameterSplit([a0, a1, a2, a3, a4, a5])
    # getting the inputs
    inputs_sq, inputs_sin = inputSplit(inputs)  # *** include freqs ***
    # getting functions
    poly = polySignal(inputs_sq, *params_sq)
    sine = sineSignal(inputs_sin, *params_sin)
    # adding into one signal
    signal = np.array(poly + sine)
    return signal



##########################################################################################
                                    # PolySignal

def polySignal(inputs, a0, a1, a2):
    parameters = [a0, a1, a2]
    lList, freq = inputs[:]
    # making BMode model
    PolySignal = np.sum([param*np.power(lList, power) for param, power in zip(parameters, range(len(parameters))[::-1])], axis=0)
    return np.array(PolySignal)


##########################################################################################
                                    # Sin Signal

def sineSignal(inputs, a0, a1, a2):
    lList, freq = inputs[:]
    # getting parameter
    ampl, period, phase = a0, a1, a2
    # making dust model
    cossignal = ampl * np.cos(period * lList + phase)
    return np.array(cossignal)


##########################################################################################
                                    # Parameters

def parameterSplit(parameters):
    ''' This function splits the list of parameters into the BMode and Dust sections
    '''
    # Poly
    try:  # checking if already packaged
        len(parameters[0])
    except:  # just a list, so package
        params_sq = parameters[:3]
    else:  # already packaged
        params_sq = parameters[0]
    finally:
        params_sq = np.array(params_sq)

    # Sin function
    try:
        len(parameters[1])
    except:
        params_sin = parameters[-3:]
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
# parameters = [1, 1, 1, 1, 1, -np.pi/2]
# print(polySignal(x, *parameters[:3]))
