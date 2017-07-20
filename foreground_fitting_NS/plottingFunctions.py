# -*- coding: utf-8 -*-
# Plotting Function

import matplotlib.pyplot as plt
import numpy as np

'''
make these more general purpose functions
'''

def plotScatter(xAxis, Measured, BMode, Dust, bestFit, Theory):
    # Makes a scatter plot
    fig = plt.figure()
    ax1 = plt.subplot(2, 1, 1)
    measured, = ax1.plot(xAxis, Measured.data, label='Measured')
    Fit, = ax1.plot(xAxis, bestFit.data, 'k', label='Best Fit')
    plt.title("Measured and Best Fit Function")
    handles, labels = ax1.get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper right')

    ax2 = plt.subplot(2, 1, 2)
    bmode, = ax2.plot(xAxis, [B*BMode.fitCoeff for B in BMode.data], 'g', label='BMode')
    dust, = ax2.plot(xAxis, [D*Dust.fitCoeff for D in Dust.data], 'r', label='Dust')
    bestfit, = ax2.plot(xAxis, bestFit.data, 'k', label='Best Fit')
    theory, = ax2.plot(xAxis, Theory.data, 'b', label='Theory')
    plt.title(str(Dust.fitCoeff)+"*Dust & "+str(BMode.fitCoeff)+"*BMode")
    handles, labels = ax2.get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper right')

    plt.savefig("lineplots.png")
    plt.close(fig)


def plotErrorbar(lBinCent, MeasuredBin, BModeBin, dustBin, bestFit, TheoryBin):
    # Makes a scatter plot
    fig = plt.figure()
    # Measured data and best fit
    ax1 = plt.subplot(2, 1, 1)
    ax1.errorbar(lBinCent, MeasuredBin.data, yerr=MeasuredBin.error, ls='o', label='Measured')  # plots the mean of each bin at the center of each bin
    ax1.errorbar(lBinCent, bestFit.data, yerr=bestFit.error, ls='--', label='Best Fit')  # yMplot[-1][0].set_linestyle('--')
    plt.title("Measured and bestFit")
    plt.xticks(lBinCent)
    handles, labels = ax1.get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper right')

    # compares best fit to adjusted values of theory
    ax2 = plt.subplot(2, 1, 2)
    ax2.errorbar(lBinCent, bestFit.data, yerr=bestFit.error, ls='--', label='Best Fit')
    ax2.errorbar(lBinCent, [B*BModeBin.fitCoeff for B in BModeBin.data], yerr=BModeBin.error, ls='--', label='B Mode')
    ax2.errorbar(lBinCent, [D*dustBin.fitCoeff for D in dustBin.data], yerr=dustBin.error, ls='--', label='Dust')
    ax2.errorbar(lBinCent, TheoryBin.data, yerr=TheoryBin.error, ls='--', label='Theory')
    handles, labels = ax2.get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper right')

    plt.savefig("binplots.png")
    plt.close(fig)


def plotCorrelation(aList, bList, plotFit, aMean, bMean):
    # Makes a correlation plot between _________________
    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    # Correlation plot
    correlplot, = ax.plot(aList, bList, 'bo', label='b by a List')
    # best fit plot
    plotList = np.arange(min(aList), max(aList), 0.1)
    bestfit, = ax.plot(plotList, [plotFit[0] + plotFit[1]*x for x in plotList], 'r-', label='bestfit {:2.0f}'.format(plotFit[2]))
    # mean plots
    bmean = ax.axhline(y=bMean, c='k', linestyle='dashed', linewidth=1, label='b mean')
    amean = ax.axvline(x=aMean, c='k', linestyle='dashed', linewidth=1, label='a mean')

    plt.title('bList by aList')
    plt.xlabel('aList')
    plt.ylabel('bList')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper right')
    plt.savefig('aList by bList.png')
    plt.close(fig)


def plotHisto(histoDataA, histoDataB, aMean, bMean, fitCoeff):
    # Makes a histogram of the binned lists done in the correlation plot using histogram data made in MonteCarlo
    fig = plt.figure()
    # A List Subplot
    ax1 = plt.subplot(2, 1, 1)
    ahist = ax1.bar(histoDataA[1][0:-1], histoDataA[0], histoDataA[1][1]-histoDataA[1][0], label='A Hist')
    amean = ax1.axvline(x=aMean, c='g', linestyle='dashed', label='A Mean')
    aCoeff = ax1.axvline(x=fitCoeff[0], c='r', linestyle='dashed', label='A Coeff')

    plt.title('aList')
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels, loc='upper right')

    # B List Subplot
    ax2 = plt.subplot(2, 1, 2)
    bhist = ax2.bar(histoDataB[1][0:-1], histoDataB[0], histoDataB[1][1]-histoDataB[1][0], label='B Hist')
    bmean = ax2.axvline(x=bMean, c='g', linestyle='dashed', label='B Mean')
    bCoeff = ax2.axvline(x=fitCoeff[1], c='r', linestyle='dashed', label='B Coeff')
    plt.title('bList')
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles, labels, loc='upper right')

    plt.savefig('aList, bList histograms.png')
    plt.close(fig)
