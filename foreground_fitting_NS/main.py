# -*- coding: utf-8 -*-
# Main execution file

import numpy as np
from dataExtractor import extractData  # this extracts the data from a txtfile containing the Lamda Data
import fittingFunctions as fit         # best fit functions
import Equations as eqn                # file with most of the used equations
import plottingFunctions as pltfn      # handles all plotting
from MonteCarlo import MonteCarlo      # runs the monte carlo
from Classes import dataClass          # all the classes

'''
Document Info
- Function of L
'''

'''
Notes to Self
the best fit function sometimes doesn't work great if the errors are wierd
change nu0 to a real value?
the MonteCarlo function first does a plt.hist then passes that to a plt.bar. Change that.
    make the bins be the same as in lBins(histoDataA, histoDataB, aMean, bMean
move the non-Equation functions out of Equations
dust temperature is never called should it be?
'''


# ------------------------------------------
# Constants

# making a dictionary with all the angle l values as the x-axis
lDict = dict(lMin=40, lMax=400, lStep=20)  # minimum l value used, maxiumm l value used, l-step value for binning,
lDict.update(lList=[l for l in range(lDict['lMin'], lDict['lMax'])]) # x axis for all plots
lDict.update(lBins=range(lDict['lMin'], lDict['lMax']+lDict['lStep'], lDict['lStep'])) # list of bins
lDict.update(lBinCent=eqn.binCenter(lDict['lBins'])) # finding bin centers

# Reference dictionary# Reference dictionary of constants of constants
const = dict(nu1=90.0*(10**9), nu2=150.0*(10**9), nu0=1.0, h=6.62606957*(10**-34),
             c=299792458, k=1.3806488*(10**-23), TVac=2.7, TDust=19.6)
             # frequency 1, frequency 2, reference frequency, Planck's constant,
             # speed of light, Boltzmann constant, vacuum temp, dust temp
const.update(list=[const['h'], const['c'], const['k'], const['TVac'], const['TDust']])


# ------------------------------------------
# Creating functions

# ----- Function Lists -----
# Theory data
Dust = dataClass('Dust', [eqn.dustRatio(const)*eqn.dust(l) for l in lDict['lList']])       # Dust of l
BMode = dataClass('BMode', extractData("LAMDA Data")[lDict['lMin']:lDict['lMax']])         # BB(l) extracted from file
Theory = dataClass('Theory', (Dust.data+BMode.data))                                       # the theoretical curve
# Measured data
Measured = dataClass('Measured', [np.random.normal(T, T/10) for T in Theory.data])         # the added noise is fake data until recieve real data
Measured.error = eqn.error(Measured.data, pct=0.2, mtd=1)                                  # the error in the measured data
# Best fit. Not for important use since not binned. referenced as "r-"
rTheoryList = [Dust.data, BMode.data]                                                      # concatenation of Dust and BB
rFitCoeff = fit.matrixFit(rTheoryList, Measured.data, Measured.error)                      # getting the best fit coefficients of theory to measured data
Dust.fitCoeff, BMode.fitCoeff = rFitCoeff[:]
rBestFit = dataClass('Best Fit', (Dust.data*rFitCoeff[0]+BMode.data*rFitCoeff[1]))         # best fit to measured by [coeffs * theory]

# -----  Binned Function Lists -----
# Theory data
DustBin = dataClass('DustBin', eqn.binData(Dust.data, lDict['lList']))                     # binned Dust
DustBin.error = eqn.error(DustBin.data, pct=0.2, mtd=1)                                    # error in binned Dust
BModeBin = dataClass('BModeBin', eqn.binData(BMode.data, lDict['lList']))                  # binned BB
BModeBin.error = eqn.error(BModeBin.data, pct=0.2, mtd=1)                                  # error in binned BB
TheoryBin = dataClass('TheoryBin', (DustBin.data+BModeBin.data))                           # the theoretical curve
TheoryBin.error = eqn.error(TheoryBin.data, pct=0.1, mtd=1)                                # error in theoretical curve      # Important for MonteCarlo
# Measured data
MeasuredBin = dataClass('MeasuredBin', eqn.binData(Measured.data, lDict['lList']))         # binned measured data
MeasuredBin.error = eqn.error(MeasuredBin.data, pct=0.2, mtd=1)                            # error in binned measured data
# Best Fit
theoryList = [DustBin.data, BModeBin.data]                                                 # concatenation of binned Dust and binned BB
fitCoeff = fit.matrixFit(theoryList, MeasuredBin.data, MeasuredBin.error)                  # best fit of Measured Bin to
DustBin.fitCoeff, BModeBin.fitCoeff = fitCoeff[:]
bestFit = dataClass('Best Fit', (DustBin.data*fitCoeff[0]+BModeBin.data*fitCoeff[1]))      # best fit to MeasuredBin by [coeffs * theory]
bestFit.error = eqn.error(bestFit.data, pct=0.2, mtd=1)                                    # error in best fit


# ------------------------------------------
# Plotting Stuff

# Fit Plots
pltfn.plotScatter(lDict['lList'], Measured, BMode, Dust, rBestFit, Theory)
pltfn.plotErrorbar(lDict['lBinCent'], MeasuredBin, BModeBin, DustBin, bestFit, TheoryBin)

# MonteCarlo
histoDataA, histoDataB, aMean, bMean, plotFit, aList, bList = MonteCarlo(theoryList, TheoryBin.data, TheoryBin.error, iterate=10**4)   # Monte Carlo of _________
pltfn.plotHisto(histoDataA, histoDataB, aMean, bMean, fitCoeff)
pltfn.plotCorrelation(aList, bList, plotFit, aMean, bMean)

# seeing some outputs
print "fitCoeff {}".format(fitCoeff)
print "aCoeff - aMean = {}".format(DustBin.fitCoeff - aMean)
print "bCoeff - bMean = {}".format(BModeBin.fitCoeff - bMean)
