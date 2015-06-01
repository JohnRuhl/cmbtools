# -*- coding: utf-8 -*-
# Generalized Least Squares Regression Fit
import numpy as np

# Generates normalized A vector using general least squares regression fit
def matrixFit(functionList, yMeasured, errorYList):
    
    # Creates the A matrix for use in determining the constants
    columns = len(functionList)
    rows = len(errorYList)
    # initialize the matrix with float values of 0
    matrixA = np.matrix([[0.0 for i in functionList] for i in errorYList])
    # initialize the list used as temporary storage for the row values
    tempList = np.empty(columns)
    
    for i in range(rows):
        for j in range(columns):
            s = errorYList[j]
            tempList[j] = functionList[j][i]/s
        matrixA[i] = tempList
    
    # initialization and assignment of vector b
    vectorB = np.empty(len(yMeasured))
    
    for y, s, i in zip(yMeasured, errorYList, range(len(yMeasured))):
        vectorB[i] = y/s
    
    # ((A transpose) dot (A))
    tempA = np.dot(matrixA.T, matrixA)
    
    # inverse of a matrix
    tempB = tempA.I
    
    # Dot Product of matrix and b vector
    tempC = np.dot(matrixA.T, vectorB)
    
    # Turns result into python array
    finalTemp = np.dot(tempC, tempB)
    vectorA = np.array(finalTemp)[0]
    
    return vectorA
    
