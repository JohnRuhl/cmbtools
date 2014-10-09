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


# ---------------------------------------------
# Finding a standard deviation value for the Monte Carlo 
# it's probably a function. don't know yet
std = 1


# ---------------------------------------------------
# Theoretical Values for a and b (slope and intercept)

aTrue = 1
bTrue = 1

aList = []
bList = []
rSquareList = []

# ------------------------------------------------
# Creating Best Fit Line, storing, and reiterating.


xList = [1,2,3,4,5]


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
plt.plot(aList, bList, 'bo')

plt.figure(1)
plt.bar(histoDataB[1][0:-1], histoDataB[0])

plt.figure(2)
plt.bar(histoDataA[1][0:-1], histoDataA[0])
