# -*- coding: utf-8 -*-
# All the Functions

import numpy as np
from fractions import gcd


##########################################################################################
                                    # Miscellaneous Functions
def floor(number):
    ''' Returns the floor of a number
    '''
    return int(number)


def ceil(number):
    ''' Returns the ceiling of a number
    '''
    return (number-(number % (number/abs(number)))+1)


def npFloat(arg):
    ''' makes a numpyFloat64 array from a normal list
        numpy must be imported as np
    '''
    for i, val in enumerate(arg):
        arg[i] = float(val)
    arg = np.array(arg)
    return arg


def floatList(arg):
    ''' makes a numpyFloat64 array of arrays
        numpy must be imported as np
    '''
    outList = [0]*len(arg)
    for i, val in enumerate(arg):
        outList[i] = (npFloat(val))
    return outList


def concatenate(*args):
    ''' concatenates a list of 1-dimensional, horizontal arrays
        the output is a numpyFloat64 array
    '''
    args = floatList(args)
    newList = np.concatenate(args)
    return newList


def GCD(*numbers):
    """ Greatest common divisor of given integers
        from https://gist.github.com/endolith/114336
    """
    return reduce(gcd, numbers)


def LCM(numberTuple):
    """ Lowest common multiple
        from https://gist.github.com/endolith/114336
        accepts one list in the form of a tuple
    """
    def lcm(a, b):
        return (a * b) // gcd(a, b)
    return reduce(lcm, numberTuple, 1)


def closest_square(number):
    ''' returns the closest square and its sqrt
    '''
    sqrt = np.sqrt(number)
    closest = int(np.round(sqrt, 0))
    square = int(np.square(closest))
    return square, closest


def lists2Lengths(listOLists):
    ''' returns a list of the lengths of a list of lists
    '''
    listLengths = []
    for params in listOLists:
        try:
            length = len(params)
        except TypeError:
            length = 1
        listLengths.append(length)
    return listLengths


def array_nest(argument, numpy=False):
    if numpy is True:
        return np.array([argument])
    else:
        return [argument]




##########################################################################################
                                    # Lists and Binning
# makes list from upper and lower bounds, and a step
def makeList(start, stop, step=1):
    '''can make a list with float inputs
       can't make a list with upper number > ~6.5e7
    '''
    List = []
    number = start
    if isinstance(step, (int, float)):
        if step < 0:
            step = -step
        if start > stop:
            step = -step
    # make it so it can accept a list here
    if start < stop:
        while number < stop:
            List.append(number)
            number += step
    elif start > stop:
        while number > stop:
            List.append(number)
            number += step
    return np.array(List)


# bins data
def binData(data, xData, xStep=[20]):
    ''' averages the given data into bins of the given size
        It only works on 1 dimensional lists
        Note: unless xStep 1 number, assumes it works for data and xData.
    '''
    if isinstance(xStep, (int, float)):
        xStep = [xStep]
        for i in range((len(xData)-(len(xData) % xStep[0]))/xStep[0]):        # yields the floor of the division of (len(nData))/xStep
            xStep.append(xStep[0])
    elif len(xStep) == 1:                                                     # makes xStep a list
        for i in range((len(xData)-(len(xData) % xStep[0]))/xStep[0]):        # yields the floor of the division of (len(nData))/xStep
            xStep.append(xStep[0])
    outList = []
    index, binnum = 0, 0                                                      # index steps through every datapoint. binnum is current bin number
    for i in range((len(xData)-(len(xData) % xStep[binnum]))/xStep[binnum]):  # steps thru the bins
        temp = []
        for j in range(xStep[binnum]):                                        # steps thru data in each bin
            # print len(data), index
            temp.append(data[index])                                          # appends data to temp list
            index += 1                                                        # index always increases, never reset
        binnum += 1                                                           # steps through the xStep for that bin
        outList.append(np.mean(temp))                                         # appends average of current temp list to the outList
    outList = np.array(outList)
    return outList


def bin2Data(dataList, xData, xStep=[20]):
    '''averages the given lists of data into lists of bins of the given size by calling binData
       Note: unless xStep 1 number, assumes it works for data and xData.
       dataLists are in the form [[data1, data2, data3], [data1, data2, data3]], where each dataX is a list
    '''
    outList = []  # the binned list
    for i, group in enumerate(dataList):  # iterating through the groups in the list of lists of data
        outList.append([])
        for j, data in enumerate(group):  # iterating through the lists of data
            outList[i].append(binData(data, xData, xStep))  # binning the data
    return outList


# bins dictionary of data
def binDict(dictionary, xData, xStep=[20]):
    '''averages the data in the dictionary into bins of the given size
       Note: unless xStep 1 number, assumes it works for data and xData.
       Dictionary must only be composed of data to be binned
    '''
    outDict = {}
    for data, key in enumerate(dictionary):
        outDict[key] = binData(data, xData, xStep)
    return outDict


# Centers of Bins
def binCenter(lBins):
    ''' Automatically finds the centers of each bin and bins the data accordingly
    '''
    centers = []
    temp = [lBins[0]]
    for i, v in enumerate(lBins[1:]):
        temp.append(v)
        centers.append(0.5*(temp[i]+temp[i+1]))
    return np.array(centers)


def addDictData(*args):
    ''' both dictionaries must have the same keys!
    '''
    answer = args[0]
    for dictionary in args[1:]:
        for key, value in dictionary.iteritems():
            answer[key] += value
    return answer


def parameterSplit(parameters):
    ''' This function splits the list of parameters into the BMode and Dust sections
    '''
    # BMode
    try:
        len(parameters[0])
    except:
        params_BM = [parameters[0]]
    else:
        params_BM = parameters[0]
    finally:
        params_BM = np.array(params_BM)
    # Dust
    try:
        len(parameters[1:])
    except:
        params_D = [parameters[1:]]
    else:
        params_D = parameters[1:]
    finally:
        params_D = np.array(params_D)
    return params_BM, params_D


##########################################################################################
                                    # Generating Data
def dudata(data, pct):
    ''' makes dud data
    '''
    noisydata = [np.random.normal(T, abs(T*pct)) for T in data]
    return np.array(noisydata)


def noisydata(data, error):
    noisydata = []
    for i, datalist in enumerate(data):
        noisydata.append([np.random.normal(m, err) for m, err in zip(datalist, error[i])])
    return noisydata


def error(data, pct=0.02, mtd=1):
    if mtd == 1:
        error = [pct*max(data) for d in data]
    elif mtd == 2:
        error = [pct*np.square(max(data))/data[int(d)] for d in data]
    elif mtd == 3:
        error = [pct*data[int(d)] for d in data]
    # else:         # make it so that it can take a user given function
    error = np.array(error)
    return error


# Data extractor
def extractData(filename, columnNumber=3):
    # Extracts data from a file with a known number of columns of data
    data = np.genfromtxt(filename)
    outputList = []
    for row in data:
        outputList.append(row[columnNumber])  # ouputs only 1 column of data
    return np.array(outputList)


##########################################################################################
                                    # Updating Data
def update_fitdata(*args):
    '''
        arguments must be a tuple with the object and the modelname
        ex ((Dust, 'bin'), (Dust, 'raw'))
    '''
    for argument in args:  # iterates through arguments
        for val in argument[1:]:  # iterates through which things to update
            argument[0].add_fitdata(val, None)
    # *** Change to use an instance initial argument, so that it can update any instance, not just fitdtata ***


def update_d_fitdata(*args):
    '''
        arguments must be a tuple with the object and the modelname
        ex ((Dust, 'bin'), (Dust, 'raw'))
    '''
    for argument in args:
        for val in argument[1:]:
            argument[0].add_d_fitdata(val, None)


# def data(function, func_input, iterator, mtd=1, ind_var=0):
#     ''' Generates Data by evaluation a function over a an iterator
#         function: the function to be evaluated
#             must have only one input. Input must be a tuple
#             independent variable must be the first input
#         func_input: a tuple containing the inputs to the function
#             must contain the independent variable
#             independent variable is assumed to be first in tuple unless otherwise stated in "ind_var"
#         iterator: the list over which the function is evaluated
#         mtd: selects which method to use
#         ind_var: specifies location of independent variable in func_input
#     '''
#     if mtd == 1:
#         data = [function(func_input) for func_input[ind_var] in iterator]
#     else:
#         print 'Error: section in development'
#         return
#     return data
