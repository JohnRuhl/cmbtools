# -*- coding: utf-8 -*-
# Plotting Function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from Functions import LCM, lists2Lengths
from fittingFunctions import polyfit

'''
make these more general purpose functions
'''


# class MakePlot(object):
#     """docstring for makePlot"""
#     def __init__(self, *args):
#         super(makePlot, self).__init__()
#         self.args = args[:]

#     # THIS PART IS IN PROGRESS
#     def makePlot(self, fig=None, index=1):
#         if fig is None:
#             fig = plt.figure()
#             ax = plt.subplot(111)
#         else:
#             ax = fig.add_subplot(2,1,index)
#         name = 'curve{}'.format(1)
#         name = ax.errorbar(self.xval, self.yval, yerr=self.error, ls='-', label=name)
#         handles, labels = ax.get_legend_handles_labels()
#         ax.legend(handles, labels, loc='upper right')
#         return fig

def makePlot(xAxis=None, *args, **kw):
    fig, ax = plt.subplots(1, 1)
    title = ''
    for i, val in enumerate(args):
        val._plot(fig, ax, xAxis)
        title += '{}, '.format(val.name)
    if 'title' in kw:
        title = kw.get('title')
    plt.title(title)
    return fig, ax


def plotList(*args):
    plotlist = []
    for i, val in enumerate(args):
        plotlist.append(val)
    return plotlist

# def ComparePlots(*args):
#     fig = plt.figure()
#     for i, val in enumerate(args):
#         ax = plt.subplot(len(args), 1, i+1)
#         ax = val
#     plt.savefig("lineplots2.png")
#     plt.close(fig)


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
    bmode, = ax2.plot(xAxis, [B*BMode.params for B in BMode.data], 'g', label='BMode')
    dust, = ax2.plot(xAxis, [D*Dust.params for D in Dust.data], 'r', label='Dust')
    bestfit, = ax2.plot(xAxis, bestFit.data, 'k', label='Best Fit')
    theory, = ax2.plot(xAxis, Theory.data, 'b', label='Theory')
    plt.title(str(Dust.params)+"*Dust & "+str(BMode.params)+"*BMode")
    handles, labels = ax2.get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper right')

    plt.savefig("lineplots.png")
    plt.close(fig)


def plotErrorbar(lBinCent, MeasuredBin, BModeBin, dustBin, bestFit, TheoryBin):
    # Makes a scatter plot
    fig = plt.figure()
    # Measured data and best fit
    ax1 = plt.subplot(2, 1, 1)
    ax1.errorbar(lBinCent, MeasuredBin.data, yerr=MeasuredBin.d_data, ls='o', label='Measured')  # plots the mean of each bin at the center of each bin
    ax1.errorbar(lBinCent, bestFit.data, yerr=bestFit.d_data, ls='--', label='Best Fit')  # yMplot[-1][0].set_linestyle('--')
    plt.title("Measured and bestFit")
    handles, labels = ax1.get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper right')

    # compares best fit to adjusted values of theory
    ax2 = plt.subplot(2, 1, 2)
    ax2.errorbar(lBinCent, bestFit.fitdata, yerr=bestFit.d_fitdata, ls='--', label='Best Fit')
    ax2.errorbar(lBinCent, BModeBin.fitdata, yerr=BModeBin.d_fitdata, ls='--', label='B Mode')
    ax2.errorbar(lBinCent, dustBin.fitdata, yerr=dustBin.d_fitdata, ls='--', label='Dust')
    ax2.errorbar(lBinCent, TheoryBin.fitdata, yerr=TheoryBin.d_fitdata, ls='--', label='Theory')
    handles, labels = ax2.get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper right')

    plt.savefig("binplots.png")
    plt.close(fig)



##########################################################################################
                                    # Monte-Carlo Plotters

def plotCorrelation(MCData):
    numParams = len(MCData.params)
    numRows, numCols = numParams-1, 2*(numParams-1)  # number of columns and rows
    plots_per_row = [numParams-1-x for x in range(numParams-1)]  # number of plots per row

    fig = plt.figure()
    gs = GridSpec(numRows, numCols)
    for row, numPlots in enumerate(plots_per_row):  # rows as indes and number of plots for the columns
        plotStart, plotEnd = 2*row, 2*row+2   # where the plot starts
        for col in range(row, numPlots+row):
            # making axes
            ax = plt.subplot(gs[row, plotStart:plotEnd])  # placing axes
            plotStart, plotEnd = plotEnd, plotEnd+2  # moving next axes over a column

            # finding x and y data
            xdata = MCData.results[:, row]  # selecting xdata
            ydata = MCData.results[:, col+1]  # selecting ydata

            # correlation plot
            correlplot, = ax.plot(xdata, ydata, 'bo', label='correlation')

            # Best Fit Plot
            # plotList = np.arange(min(xdata), max(xdata), 0.1)
            # fitCoeffs = polyfit(xdata, ydata, 1)
            # bestfit, = ax.plot(plotList, [fitCoeffs[0] + fitCoeffs[1]*x for x in plotList], 'r-', label='bestfit')

            # mean plots
            param1mean = ax.axvline(x=np.mean(xdata), c='k', linestyle='dashed', linewidth=1, label='param{} mean: {:.2f}'.format(row, np.mean(xdata)))
            param2mean = ax.axhline(y=np.mean(ydata), c='k', linestyle='dashed', linewidth=1, label='param{} mean: {:.2f}'.format(col+1, np.mean(ydata)))

            # plot properties
            plt.title('parameter {},{}'.format(row+1, col+2))
            # plt.xlabel('parameter {}'.format(row+1))
            # plt.ylabel('parameter {}'.format(col+2))
            # handles, labels = ax.get_legend_handles_labels()
            # ax.legend(handles, labels, loc='upper right')
            ax.set_xticklabels([])  # turns off x-tick labels
            ax.set_yticklabels([])  # turns off y-tick labels

    # figure properties
    plt.suptitle("Correlation Plots")
    plt.savefig('MC Parameter Correlations.png')
    plt.close(fig)


def plotHisto(MCData, paramsList, name='Histograms.png'):
    ''' Makes histograms from the MonteCarlo data
        paramsList is a list of the params list for each object (like BMode & Dust)
            needs to be in same order as params
    '''
    # Making the histogram bins
    numIterations = len(MCData.results[:,0])
    if numIterations < 100:
        binNumber = int(numIterations/10)
    else:
        binNumber = 100

    # Making the Plot
    fig = plt.figure()
    numRows, numCols = len(paramsList), LCM(tuple(lists2Lengths(paramsList)))  # finding number of rows and columns
    gs = GridSpec(numRows, numCols)
    paramIndex = 0  # which parameter at.  ** works until params packaged differently **
    for row, params in enumerate(paramsList):  # iteraturubg through the parameter lists
        colIndex, colStep = 0, (numCols / len(params))  # spacing the histogram plots correctly
        for j, val in enumerate(params):  # iterating through the parameters in each parameter list
            # initializing the subplot
            newColIndex = colIndex+colStep  # finding the end index of the plot
            ax = plt.subplot(gs[row, colIndex:newColIndex])  # making the subplot
            colIndex = newColIndex  # updating column index

            # making the plot data
            histoData = np.histogram(MCData.results[:, paramIndex], binNumber)  # ** when MCData packaged by params( not 1 2-D-list) will need different index **
            histogram  = ax.bar(histoData[1][0:-1], histoData[0], histoData[1][1]-histoData[1][0], label='Hist')
            mean = ax.axvline(x=MCData.params[paramIndex], c='g', linestyle='dashed', label='Mean')
            parameter = ax.axvline(x=paramsList[row][j], c='r', linestyle='dashed', label='Parameter')
            paramIndex += 1  # which parameter at.  ** works until params packaged differently **

            # subplot properties
            plt.title('parameter {},{}'.format(row, j))
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, loc='upper right')
    plt.suptitle("MonteCarlo Histograms")
    plt.savefig(name)
    plt.close(fig)


def plotMCMeasured(MCData):
    fig, ax = plt.subplots(1, 1)
    for measured in MCData.measured_log[1:]:
        ax.plot(MCData.measured.xaxis, measured, 'b--', alpha=0.5)
    ax.plot(MCData.measured.xaxis, MCData.measured.data, 'r-')
    plt.title("MC 'Measured'")
    fig.savefig('MC Measured.png')
    plt.close(fig)
