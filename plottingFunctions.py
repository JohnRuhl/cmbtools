# -*- coding: utf-8 -*-
# Plotting Function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from Functions import LCM, lists2Lengths, parameterSplit
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


def plotErrorbar(xaxis, MeasuredBin, BModeBin, dustBin, bestFit, TheoryBin, index):
    # Makes a scatter plot
    fig = plt.figure()
    # Measured data and best fit
    ax1 = plt.subplot(2, 1, 1)
    ax1.errorbar(xaxis, MeasuredBin.fitdata[0][index], yerr=MeasuredBin.d_fitdata[0][index], ls='o', label='Measured{}'.format(index))  # plots the mean of each bin at the center of each bin
    ax1.errorbar(xaxis, bestFit.fitdata[-1][index], yerr=None, ls='--', label='Best Fit{}'.format(index))  # yMplot[-1][0].set_linestyle('--')
    plt.title("Measured and bestFit")
    handles, labels = ax1.get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper right')

    # compares best fit to adjusted values of theory
    ax2 = plt.subplot(2, 1, 2)
    ax2.errorbar(xaxis, bestFit.fitdata[-1][index], yerr=bestFit.d_fitdata[-1][index], ls='--', label='Best Fit{}'.format(index))
    ax2.errorbar(xaxis, BModeBin.fitdata[-1][index], yerr=None, ls='--', label='B Mode{}'.format(index))
    ax2.errorbar(xaxis, dustBin.fitdata[-1][index], yerr=None, ls='--', label='Dust{}'.format(index))
    ax2.errorbar(xaxis, TheoryBin.fitdata[-1][index], yerr=None, ls='--', label='Theory{}'.format(index))
    handles, labels = ax2.get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper right')

    plt.savefig("binplots{}.png".format(index))
    plt.close(fig)



##############################################################################################
                                    # Monte-Carlo Plotters

def plotCorrelation(MCData, paramsList, name='MC Parameter Correlations.png'):
    # Finding plot properties
    numParams = len(MCData.params)  # number of parameters
    numRows, numCols = numParams, numParams  # number of columns and rows
    plots_per_row = [numParams-x for x in range(numParams)]  # number of plots per row
    numIterations = len(MCData.results[:, 0])  # number of data points
    if numIterations < 100:  # packaging data points into bins
        binNumber = int(numIterations/10)
    else:
        binNumber = 100

    fig = plt.figure()
    gs = GridSpec(numRows, numCols)
    for row, numPlots in enumerate(plots_per_row):  # rows as indices and number of plots for the columns
        plotStart, plotEnd = row, row+1   # determines the column start and stop point
        for col in range(row, numPlots+row):  # iterates through the parameters
            # making axes
            ax = plt.subplot(gs[row, plotStart:plotEnd])  # placing axes
            plotStart, plotEnd = plotEnd, plotEnd+1  # moving next axes over a column

            # finding x and y data
            xdata = MCData.results[:, row]  # selecting xdata
            ydata = MCData.results[:, col]  # selecting ydata

            # correlation plot
            if row == col:
                histoData = np.histogram(MCData.results[:, row], binNumber)  # *** when MCData packaged by params( not 1 2-D-list) will need different index ***
                histogram  = ax.bar(histoData[1][0:-1], histoData[0], histoData[1][1]-histoData[1][0], label='Hist')
                mean = ax.axvline(x=MCData.params[row], c='g', linestyle='dashed', label='Mean')
                # get right indices
                splitparams = parameterSplit(MCData.params)
                index1, index2 = 0, 0
                for i in range(len(splitparams)):
                    for j in splitparams[i]:
                        if row == (index1+index2):
                            parameter = ax.axvline(x=paramsList[index1][index2], c='r', linestyle='dashed', label='Parameter')
                            plt.title('Parameter {},{}'.format(index1, index2))  # *** FIX THIS ***
                            break
                        index2 += 1
                    index2 = 0
                    index1 += 1
            else:
                correlplot, = ax.plot(xdata, ydata, 'bo', label='correlation')
                # Best Fit Plot
                plotList = np.arange(min(xdata), max(xdata), 0.1)
                fitCoeffs = polyfit(xdata, ydata, 1)
                bestfit, = ax.plot(plotList, [fitCoeffs[0] + fitCoeffs[1]*x for x in plotList], 'r-', label='bestfit')
                # mean plots
                param1mean = ax.axvline(x=np.mean(xdata), c='k', linestyle='dashed', linewidth=1, label='param{} mean: {:.2f}'.format(row, np.mean(xdata)))
                param2mean = ax.axhline(y=np.mean(ydata), c='k', linestyle='dashed', linewidth=1, label='param{} mean: {:.2f}'.format(col+1, np.mean(ydata)))
                plt.title('parameter {},{}'.format(col, row))

            # plot properties
            # plt.xlabel('parameter {}'.format(row+1))
            # plt.ylabel('parameter {}'.format(col+2))
            # handles, labels = ax.get_legend_handles_labels()
            # ax.legend(handles, labels, loc='upper right')
            ax.set_xticklabels([])  # turns off x-tick labels
            ax.set_yticklabels([])  # turns off y-tick labels

    # figure properties
    plt.suptitle("Correlation Plots")
    plt.savefig(name)
    plt.close(fig)


def plotHisto(MCData, paramsList, name='Histograms.png'):
    ''' Makes histograms from the MonteCarlo data
        paramsList is a list of the params list for each object (like BMode & Dust)
            needs to be in same order as params
    '''
    # Making the histogram bins
    numIterations = len(MCData.results[:, 0])
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
            histoData = np.histogram(MCData.results[:, paramIndex], binNumber)  # *** when MCData packaged by params( not 1 2-D-list) will need different index ***
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


def plotMeasured(Measured, lDict):
    fig = plt.figure()
    for index, data in enumerate(Measured.raw.fitdata[0]):
        plt.plot(lDict['lList'], data, label='Measured{}.raw'.format(index))
        plt.plot(lDict['lBinCent'], Measured.bin.fitdata[0][index], label='Measured{}.bin'.format(index))
    plt.legend()
    plt.title('Measured')
    fig.savefig('Measured.png')

def plotBMode(BMode, lDict):
    fig = plt.figure()
    for index, data in enumerate(BMode.raw.data[0]):
        plt.plot(lDict['lList'], data, label='Raw')
        plt.plot(lDict['lList'], BMode.raw.fitdata[-1][index], label='Fit')
        plt.plot(lDict['lBinCent'], BMode.bin.data[0][index], label='bin')
        plt.plot(lDict['lBinCent'], BMode.bin.fitdata[-1][index], label='binFit')
    plt.legend()
    plt.title('BMode')
    fig.savefig('BMode.png')


def plotMeasuredTheory(Measured, Theory, lDict, number):
    fig = plt.figure()
    for index, data in enumerate(Measured.bin.fitdata[-1]):
        plt.plot(lDict['lBinCent'], data, label='Measured{}.raw'.format(index))
        plt.plot(lDict['lBinCent'], Theory.bin.fitdata[-1][index], label='Theory{}.bin'.format(index))
    plt.legend()
    plt.title('Theory vs Measured')
    fig.savefig('TheoryMeasured{}.png'.format(number))


def plotMeasuredTheoryFit(Measured, Theory, Fit, lDict, number):
    fig = plt.figure()
    for index, data in enumerate(Measured.bin.fitdata[-1]):
        plt.plot(lDict['lBinCent'], data, label='Measured{}.raw'.format(index))
        plt.plot(lDict['lBinCent'], Theory.bin.fitdata[-1][index], label='Theory{}.bin'.format(index))
        plt.plot(lDict['lBinCent'], Fit.bin.fitdata[-1][index], label='Fit{}.bin'.format(index))
        # print Theory.bin.fitdata[-1][index]-Fit.bin.fitdata[-1][index]
    plt.legend()
    plt.title('Theory, Measured, Fit')
    fig.savefig('TheoryMeasuredFit{}.png'.format(number))
