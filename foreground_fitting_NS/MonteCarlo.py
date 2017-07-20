# -*- coding: utf-8 -*-
# Monte Carlo

import numpy as np
import fittingFunctions as fit

'''
should a plot of RSquare from the linearFit be made?
'''


def MonteCarlo(TList, measured, std, iterate=10**4):
    aList = []
    bList = []
    for i in range(iterate):
        # generating yList
        yList = [np.random.normal(m, err) for m, err in zip(measured, std)]
        # applying lines of best fit
        fitCoeff = fit.matrixFit(TList, yList, std)
        # creating an aList, bList, and rSquareList for later analysis
        aList.append(fitCoeff[0])
        bList.append(fitCoeff[1])

    # Making the histogram bins
    if iterate < 100:
        binNumber = int(iterate/10)
    else:
        binNumber = 100
    histoDataA = np.histogram(aList, binNumber)
    histoDataB = np.histogram(bList, binNumber)

    plotFit = fit.linearFit(aList, bList)
    aMean = np.mean(aList, dtype=np.float16)
    bMean = np.mean(bList, dtype=np.float16)

    return histoDataA, histoDataB, aMean, bMean, plotFit, aList, bList
