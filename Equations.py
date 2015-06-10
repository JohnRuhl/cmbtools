# -*- coding: utf-8 -*-
# Equations

import numpy as np


def dustFreqPowLaw(nu, nu0):
    # The power law used in some dust equations (mostly a filler function)
    return (nu/nu0)**1.59


def blackbodyConvertofNu(nu, list):
    # A Plancks Law frequency function to convert the Black Body equation to the right units
    # nu, const = [h, c, k, TVac, TDust]
    # TDust=19.6, h=6.62606957*(10**-34), c=299792458, k=1.3806488*(10**-23), TVac = 2.7
    h, c, k, TVac, TDust = list[:]
    return 2*(h**2)*(nu**4)*np.exp((h*nu)/(k*TVac))/(k*(TVac**2)*(c**2)*((np.exp((h*nu)/(k*TVac)) - 1)**2))


def blackbody(nu, list):                                                     # Why is blackbody function never called?
    # nu, h=6.62606957*(10**-34), c=299792458, k=1.3806488*(10**-23), T=2.7
    # Another Plankcs Law equation
    h, c, k, TVac, TDust = list[:]
    return (2*(h*(nu**3))/(c**2))*(1/(np.exp((h*(nu))/(k*TVac)) - 1))


def dust(l):
    # Exponential dust function
    return (l/80.0)**(-0.42)


def dustRatio(const):
    # Provides the ratio to multiply the dust by when given two frequencies of blackbody
    # This is one of the two "outermost" functions that is called
    # const stands for a list of constants given in main.py
    return (dustFreqPowLaw(const['nu1'], const['nu0'])*blackbodyConvertofNu(const['nu2'], const['list'])) \
           / (dustFreqPowLaw(const['nu2'], const['nu0'])*blackbodyConvertofNu(const['nu1'], const['list']))


def binData(data, xData, xStep=[20]):
    # averages the given data into bins of the given size

    if len(xStep) == 1:                                                       # makes xStep a list
        for i in range((len(xData)-(len(xData) % xStep[0]))/xStep[0]):        # yields the floor of the division of (len(nData))/xStep
            xStep.append(xStep[0])
    outList = []
    index, binnum = 0, 0                                                      # index steps through every datapoint. binnum is current bin number
    for i in range((len(xData)-(len(xData) % xStep[binnum]))/xStep[binnum]):  # steps throu the bins
        temp = []
        for j in range(xStep[binnum]):          # steps thru data in each bin
            temp.append(data[index])            # appends data to temp list
            index += 1                          # index always increases, never reset
        binnum += 1                             # steps through the xStep for that bin
        outList.append(np.mean(temp))           # appends average of current temp list to the outList
    return outList


def binCenter(lBins):
    # finds the center of each bin
    out = []
    temp = [lBins[0]]
    for i, v in enumerate(lBins[1:]):
        temp.append(v)
        out.append(0.5*(temp[i]+temp[i+1]))
    return out


def error(data, pct=0.2, mtd=1):
    if mtd == 1:
        error = [pct*max(data) for d in data]
    elif mtd == 2:
        error = [pct*np.square(max(data))/data[int(d)] for d in data]
    elif mtd == 3:
        error = [pct*data[int(d)] for d in data]
    # else:         # make it so that it can take a user given function
    return error
