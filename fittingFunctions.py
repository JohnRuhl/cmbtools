# -*- coding: utf-8 -*-
# All the Fitting Functions

import numpy as np


def linearFit(xList, yList, degree=1):
    # Line of Best Fit (using Least Squares Regression):
    '''
    Description: A function that takes in two lists of equal length, corresponding
    to the x and y coords of a graph, and returns a list containing
    the intercept, slope and r^2 of the line of best fit to that graph
    '''
    coeffs = np.polyfit(xList, yList, degree)
    b = coeffs[0]  # highest degree
    a = coeffs[1]  # lowest degree
    # calculating correlation coefficient. r = b *(stdev yList / stdev xList)
    r = b * (np.std(yList)/np.std(xList))
    # calculating r-square value
    rSquare = np.square(r)
    return [a, b, rSquare]


def matrixFit(functionList, fitTo, errorFitTo):
    # Generates normalized A vector using general least squares regression fit
    '''
    Description: A function that takes in 3 lists, corresponding to a composite list
    of usable functions (in the form of y-coordinates) and the function (with errors)
    to which the formter composite list should be fit to. Returns a list containing
    the fit coefficients for the composite list of functions.
    '''
    # Creates the A matrix for use in determining the constants
    columns = len(functionList)
    rows = len(errorFitTo)
    # initialize the matrix with float values of 0
    matrixA = np.matrix([[0.0 for i in functionList] for i in errorFitTo])
    # initialize the list used as temporary storage for the row values
    tempList = np.empty(columns)

    for i in range(rows):
        for j in range(columns):
            s = errorFitTo[j]
            tempList[j] = functionList[j][i]/s
        matrixA[i] = tempList

    # initialization and assignment of vector b
    vectorB = np.empty(len(fitTo))
    # fills empty vector with measured/error
    for y, s, i in zip(fitTo, errorFitTo, range(len(fitTo))):
        vectorB[i] = y/s
    # ((A transpose) dot (A))
    tempA = np.dot(matrixA.T, matrixA)
    # inverse of a matrix
    tempB = tempA.I
    # Dot Product of matrix and b vector
    tempC = np.dot(tempB, matrixA.T)
    # Turns result into python array
    finalTemp = np.dot(tempC, vectorB)
    fitCoeff = np.array(finalTemp)[0]

    return fitCoeff


def iterMatrixFit(functionList, fitTo, errorFitTo, iterate=10**2):
    # Generates normalized A vector using general least squares regression fit, repeatedly
    '''
    Description: does the MatrixFit repeatedly to see if there is any spread. There shouldn't be
    '''
    aList = []
    bList = []
    for i in range(iterate):
        # applying lines of best fit
        temp = matrixFit(functionList, fitTo, errorFitTo)
        # creating an aList, bList, and rSquareList for later analysis
        aList.append(temp[0])
        bList.append(temp[1])
    aMean = np.mean(aList, dtype=np.float64)
    bMean = np.mean(bList, dtype=np.float64)
    aStd = np.std(aList)
    bStd = np.std(bList)
    fitCoeff = [aMean, bMean]
    errorFitCoeff = [aStd, bStd]
    return fitCoeff, errorFitCoeff
