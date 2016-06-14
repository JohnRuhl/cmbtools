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
    '''
    arg = [float(val) for val in arg]
    return np.array(arg)


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




##########################################################################################
                                    # Lists and Binning


# BINNING FUNCTIONS HERE

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]


def splitList(l, chunksize, ndarray=True):
    if ndarray is True:
        return np.array(list(chunks(l, chunksize)))
    else:
        return list(chunks(l, chunksize))

##########################################################################################
                                    # Generating Data
def noisyData(datalist, error):
    ''' (datalist, error)
        This function adds noise to the data, which is a list
        the error can be a list

        TO DO: Add options
    '''
    noisydata = np.random.normal(datalist, error)
    return noisydata


# Data extractor
def extractData(filename, columnNumber=3):
    # Extracts data from a file with a known number of columns of data
    data = np.genfromtxt(filename)
    outputList = []
    for row in data:
        outputList.append(row[columnNumber])  # ouputs only 1 column of data
    return np.array(outputList)
