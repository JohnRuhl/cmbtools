# -*- coding: utf-8 -*-
# Main execution file

'''
Document Info
- Function of L
'''

import matplotlib.pyplot as plt
import fittingFunctions as fit         # best fit functions
import Equations2 as eqn               # file with most of the used equations
import Functions as fn                 # file with most of the functions
import plottingFunctions as pltfn      # handles all plotting
# import numpy as np
from Classes import DataClass          # all the classes
from Dictionaries import make_lDict, make_Const, make_frequencies  # lDict and constants



##########################################################################################
                            # Files, Constants & Parameters

# Files
lamda_data_filename = 'LAMDA Data'
measured_bin_filename = 'results.txt'

# making a dictionary with all the angle l values as the x-axis
lDict = make_lDict(lamda_data_filename)

# Reference dictionary of constants
const = make_Const()

# Reference dictionary of frequencies
freqs = make_frequencies(90e9, 150e9)

# Best Fit Parameter Inputs
params_BM = fn.npFloat([1])
params_D = fn.npFloat([1, 1.59, 19.6])
parameters = fn.concatenate(params_BM, params_D)



##########################################################################################
                                      # Data

# ----------------------------------------
            # Theory Data
# BMode
BMode = DataClass('BMode',       {'model': 'raw',
                                  'equation': eqn.BModeSignal,
                                  'eqinput': fn.extractData(lamda_data_filename, 3)[lDict['refMin']:lDict['refMax']], 'params': params_BM,
                                  'error': lambda x: fn.error(x, pct=0.02, mtd=1),
                                  'xaxis': lDict['lList']})
BMode.add_model(                 {'model': 'bin',
                                  'equation': eqn.BModeSignal,  # is this the equation?
                                  'eqinput': fn.binData(BMode.raw.data, lDict['lList'], lDict['lStep']), 'params': params_BM,
                                  'error': lambda x: fn.error(x, pct=0.02, mtd=1),
                                  'xaxis': lDict['lBinCent']})

# Dust
Dust = DataClass('Dust',         {'model': 'raw',
                                  'equation': eqn.dustSignal,
                                  'eqinput': (lDict['lList'], const, freqs['nu1']), 'params': params_D,
                                  'error': lambda x: fn.error(x, pct=0.02, mtd=1),
                                  'xaxis': lDict['lList']},
                                  {'model': 'raw2',
                                  'equation': eqn.dustSignal,
                                  'eqinput': (lDict['lList'], const, freqs['nu2']), 'params': params_D,
                                  'error': lambda x: fn.error(x, pct=0.02, mtd=1),
                                  'xaxis': lDict['lList']})
Dust.add_model(                  {'model': 'bin',
                                  'equation': eqn.dustSignal,  # is this the equation or fn.binData?
                                  'eqinput': (lDict['lBinCent'], const, freqs['nu1']), 'params': params_D,
                                  'data': fn.binData(Dust.raw.data, lDict['lList'], lDict['lStep']),  # binning data & re-evaluating are < 1% different
                                  'error': lambda x: fn.error(x, pct=0.02, mtd=1),
                                  'xaxis': lDict['lBinCent']},
                                 {'model': 'bin2',
                                  'equation': eqn.dustSignal,  # is this the equation or fn.binData?
                                  'eqinput': (lDict['lBinCent'], const, freqs['nu2']), 'params': params_D,
                                  'data': fn.binData(Dust.raw2.data, lDict['lList'], lDict['lStep']),  # binning data & re-evaluating are < 1% different
                                  'error': lambda x: fn.error(x, pct=0.02, mtd=1),
                                  'xaxis': lDict['lBinCent']})

# Theory
Theory = DataClass('Theory',     {'model': 'raw',
                                  'equation': eqn.totalSignal,
                                  'eqinput': (lDict['lList'], const, freqs['nu1'], BMode.raw.data), 'params': parameters,
                                  'd_data': BMode.raw.d_data+Dust.raw.d_data,  # ** switch to addDictData(BMode.raw.d_data, Dust.raw.d_data) **
                                  'xaxis': lDict['lList']},
                                 {'model': 'raw2',
                                  'equation': eqn.totalSignal,
                                  'eqinput': (lDict['lList'], const, freqs['nu2'], BMode.raw.data), 'params': parameters,
                                  'd_data': BMode.raw.d_data+Dust.raw2.d_data,  # ** switch to addDictData(BMode.raw.d_data, Dust.raw.d_data) **
                                  'xaxis': lDict['lList']},
                                 {'model': 'bin',
                                  'equation': eqn.totalSignal,
                                  'eqinput': (lDict['lBinCent'], const, freqs['nu1'], BMode.bin.data), 'params': parameters,
                                  'd_data': BMode.bin.d_data+Dust.bin.d_data,  # ** switch to addDictData(BMode.bin.d_data, Dust.bin.d_data) **
                                  'xaxis': lDict['lBinCent']},
                                 {'model': 'bin2',
                                  'equation': eqn.totalSignal,
                                  'eqinput': (lDict['lBinCent'], const, freqs['nu2'], BMode.bin.data), 'params': parameters,
                                  'd_data': BMode.bin.d_data+Dust.bin2.d_data,  # ** switch to addDictData(BMode.bin.d_data, Dust.bin.d_data) **
                                  'xaxis': lDict['lBinCent']})

# ----------------------------------------
            # Measured Data
# Measured
Measured = DataClass('Measured', {'model': 'raw',
                                  'data': fn.dudata(Theory.raw.data, 0.1),
                                  'error': fn.extractData(measured_bin_filename, 2)[lDict['refMin']:lDict['refMax']],
                                  'xaxis': lDict['lList']},
                                  {'model': 'raw2',
                                  'data': fn.dudata(Theory.raw2.data, 0.1),
                                  'error': fn.extractData(measured_bin_filename, 2)[lDict['refMin']:lDict['refMax']],
                                  'xaxis': lDict['lList']})
Measured.add_model(              {'model': 'bin',
                                  'data': fn.binData(Measured.raw.data, lDict['lList'], lDict['lStep']),
                                  'error': lambda x: fn.error(x, pct=0.02, mtd=1),  # *** replace with real errors
                                  'xaxis': lDict['lBinCent']},
                                  {'model': 'bin2',
                                  'data': fn.binData(Measured.raw2.data, lDict['lList'], lDict['lStep']),
                                  'error': lambda x: fn.error(x, pct=0.02, mtd=1),  # *** replace with real errors
                                  'xaxis': lDict['lBinCent']})

# ----------------------------------------
            # Processed Data
# Fit Coefficients
ChiSqFit = fit.ChiSqOpt(Measured.bin, Theory.bin.params, BMode.bin, Dust.bin)
Theory.bin.parameters, BMode.bin.params, Dust.bin.params = ChiSqFit.fit('combo', method='nelder-mead')  # ** retire outputs and put in a fn.updata_params?**
print '\nBest Fit Parameters:', Theory.bin.parameters, '\n'

# Best Fit
# BestFit = DataClass('Best Fit', {'model': 'raw', 'data': (Dust.raw.data*Dust.raw.params + BMode.raw.data*BMode.raw.params)})
BestFit = DataClass('Best Fit')
BestFit.add_model({'model': 'bin',  # different than Theory since parameters changed
                   'equation': eqn.totalSignal,
                   'eqinput': (lDict['lBinCent'], const, freqs['nu1'], BMode.bin.data), 'params': Theory.bin.parameters,
                   'error': lambda x: fn.error(x, pct=0.02, mtd=1),
                   'xaxis': lDict['lBinCent']})

# Temporary:
BMode.raw.params, Dust.raw.params = BMode.bin.params, Dust.bin.params

# Updating Fit Data & d_Fit Data
fn.update_fitdata((BMode, 'raw', 'bin'), (Dust, 'raw', 'bin'), (Theory, 'raw', 'bin'))
fn.update_d_fitdata((BMode, 'raw', 'bin'), (Dust, 'raw', 'bin'), (Theory, 'raw', 'bin'))



##########################################################################################
                                    # Monte Carlo
# performing MC
MCData = fit.MCChiSqFit(Theory.bin, BMode.bin, Dust.bin, iterate=10, method='nelder-mead')  # Monte Carlo of _________
# plotting histograms and parameter correlations
pltfn.plotHisto(MCData, [BMode.bin.params, Dust.bin.params])
pltfn.plotCorrelation(MCData)
# pltfn.plotMCMeasured(MCData)



##########################################################################################
                                    # Outputs

# ----------------------------------------
            # Plotting Stuff
# Fit Plots:
pltfn.plotErrorbar(lDict['lBinCent'], Measured.bin, BMode.bin, Dust.bin, BestFit.bin, Theory.bin)
# pltfn.plotScatter(lDict['lList'], Measured.raw, BMode.raw, Dust.raw, BestFit.raw, Theory.raw)

# Quick Plots
fig = plt.figure()
plt.plot(lDict['lList'], BMode.raw.data, label='Raw')
plt.plot(lDict['lList'], BMode.raw.fitdata, label='Fit')
plt.legend()
plt.title('BMode')
fig.savefig('BMode.png')

fig = plt.figure()
plt.plot(lDict['lList'], Measured.raw.data, label='Measured.raw')
plt.plot(lDict['lBinCent'], Measured.bin.data, label='Measured.bin')
plt.legend()
plt.title('Measured')
fig.savefig('Measured.png')

# ----------------------------------------
            # seeing some outputs
print "\nparams {}".format(MCData.params)
print "BMode Ampl - Ampl Mean = {}".format(BMode.bin.params[0] - MCData.params[0])
print "Dust Ampl - Ampl Mean = {}".format(Dust.bin.params[0] - MCData.params[2])

# New-Style Plots
# figure, axis = pltfn.makePlot(lDict['lList'], Measured.raw, BestFit.raw)
# figure2, axis2 = pltfn.makePlot(lDict['lList'], Measured.raw, BMode.raw, Dust.raw, BestFit.raw, Theory.raw)
# pltfn.ComparePlots(axis, axis2)


''' NEED TO DO:
- change dust equation to something else
- Make better best fitter for the correlation plots
- Switch to lmfit
- ADD PRIORS (BOUNDS )
- Make multiple frequencies:
  - all the things (data, d_data, fitdata, d_fitdata) which need to use other class objects can't update right with the dict
  - switch binData to binDict

    Maybe Do:
- make all params be in form [[BMode[:]], [Dust[:]]] rather than [BMode[:], Dust[:]]
- Then wouldn't need paramsList in plotHisto. replace with MCData.parameters
'''
