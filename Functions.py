# -*- coding: utf-8 -*-
# All the Functions

import numpy as np
from fractions import gcd


##########################################################################################
                                    # Miscellaneous Functions
def floor(number):
    ''' Returns the floor of a number
    '''
    return (number-(number % (number/abs(number))))


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
    if isinstance(xStep, (int, float)) or len(xStep) == 1:
        if isinstance(xStep, (int, float)):
            temp, xStep = xStep, [xStep]
        elif len(xStep) == 1:
            temp = xStep[0]
        for i in range(int(len(xData)/temp)-1):  # yields the floor of the division of (len(nData))/xStep
            xStep.append(temp)
        if len(xData) % temp != 0:  # if xStep doesn't fit evenly into len(xData)
            xStep.append(len(xData) % temp)

    outList = []
    index = 0  # index steps through every datapoint. binnum is current bin number
    for i, Bin in enumerate(xStep):
        temp = np.mean(data[index:(index+Bin)])
        outList.append(temp)
        index += Bin
    outList = np.array(outList)
    return outList, xStep


def bin2Data(dataList, xData, xStep=[20]):
    '''averages the given lists of data into lists of bins of the given size by calling binData
       Note: unless xStep 1 number, assumes it works for data and xData.
       dataLists are in the form [[data1, data2, data3], [data1, data2, data3]], where each dataX is a list
    '''
    outList = []  # the binned list
    for i, group in enumerate(dataList):  # iterating through the groups in the list of lists of data
        outList.append([])
        for j, data in enumerate(group):  # iterating through the lists of data
            outList[i].append(binData(data, xData, xStep)[0])  # binning the data
    return outList


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

def noisyDataList(datalist, error):
    ''' (datalist, error)
        This function adds noise to the data, which is a list of lists
        the error can be a list
    '''
    noisydata = []
    for i, data in enumerate(datalist):
        noisydata.append([np.random.normal(m, err) for m, err in zip(data, error[i])])
    return noisydata


def noisyData(datalist, error):
    ''' (datalist, error)
        This function adds noise to the data, which is a list of lists
        the error can be a list
    '''
    noisydata = [np.random.normal(m, err) for m, err in zip(datalist, error)]
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

def update(method, *args):
    '''
        Give the method as a string
        arguments must be a tuple with the object and the modelname
        ex ((Dust, 'bin'), (Dust, 'raw'))
    '''
    for argument in args:
        for val in argument[1:]:
            getattr(getattr(getattr(argument[0], val), method), 'add_{}'.format(method))(None)
