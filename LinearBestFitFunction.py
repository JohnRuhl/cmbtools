# -*- coding: utf-8 -*-
# Line of Best Fit (using Least Squares Regression):

import numpy as np

'''
    A function that takes in two lists of equal length corresponding
    the x-coords and y-coords of a graph and returns the a list containing
    the intercept, slope and r^2 of the line of best fit to that
    graph respectively
'''

def linearFit(xList, yList):
    # finds number of values in xList and yList
    Nx = len(xList)
    Ny = len(yList)
    
    # error if xlist and ylist are different lengths
    if Nx != Ny:
        print "Error: incompatible lists"
        raise Exception("Incompatible List Lengths")
    
    #the degree =1 of the polynomial equation 
    ''' 
        given x-values, y-values, and relationship degree, 
        polyfit returns list of polynomial coefficients 
        starting with coefficient with the highest power 
    '''
    p = np.polyfit(xList, yList, 1)  
    b = p[0] # highest degree
    a = p[1] #lowest degree

    # calculating correlation coefficient. r = b *(stdev yList / stdev xList)
    r = b / (np.std(yList)/np.std(xList))
    
    # calculating r-square value
    rSquare = np.square(r)
    
    return [a, b, rSquare]
