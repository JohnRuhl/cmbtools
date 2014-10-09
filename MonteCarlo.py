# -*- coding: utf-8 -*-
# Monte Carlo

import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

#try:
import LinearBestFitFunction as fit
# trying to import from website
# except:
#    import site as fit
#    site.addsitedir("https://raw.githubusercontent.com/Case-Research-P-N-S/Band-Optimizing-for-Foreground-Subtraction/master/LinearBestFitFunction.py")

def monteCarloGen(xList, aTrue, bTrue, std):   
  
    # ---------------------------------------------------
    # Theoretical Values for a and b (slope and intercept)

    aList = []
    bList = []
    rSquareList = []
    
    # ------------------------------------------------
    # Creating Best Fit Line, storing, and reiterating.
    
    # assuming x's are constant, may later generate x's from a normal distribution
    # how many times the Monte Carlo is repeated
    iterations = 1000
    for i in range(iterations):
        # generating yList
        yList = [aTrue + bTrue*x + np.random.normal(0, std) for x in xList]
        # applying lines of best fit
        results = fit.linearFit(xList, yList)
        # creating an aList, bList, and rSquareList for later analysis
        aList.append(results[0])
        bList.append(results[1])
        rSquareList.append(results[2])
    
    histoDataA = np.histogram(aList)
    histoDataB = np.histogram(bList)
    
    
    plt.figure(0)
    plotFit = fit.linearFit(aList, bList)
    plotList = [x for x in np.arange(min(aList), max(aList), 0.1)]
    plt.plot(aList, bList, 'bo')
    plt.plot(plotList, [plotFit[0] + plotFit[1]*x for x in plotList])
    
    plt.figure(1)
    plt.bar(histoDataB[1][0:-1], histoDataB[0])
    
    plt.figure(2)
    plt.bar(histoDataA[1][0:-1], histoDataA[0])
