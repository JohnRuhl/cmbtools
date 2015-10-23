# -*- coding: utf-8 -*-
# Main execution file

'''
Document Info
- Function of L
'''

import numpy as np
import fittingFunctions as fit         # best fit functions
import Equations as eqn               # file with most of the used equations
import Functions as fn                 # file with most of the functions
import plottingFunctions as pltfn      # handles all plotting
import Dictionaries as dct             # lDict and constants
import lmfit as lmfit
# import numpy as np
from Classes import ModelClass          # all the classes



###############################################################################################
                            # Files, Constants & Parameters

# Files
lamda_data_filename = 'LAMDA Data'
measured_bin_filename = 'results.txt'

# making a dictionary with all the l values as the x-axis
lDict = dct.make_lDict(lamda_data_filename)

# Reference dictionary of constants
const = dct.make_Const()

# Reference dictionary of frequencies
freqDict, freqs = dct.make_frequencies(90e9, 140e9, 160e9, 220e9)
ModelClass.freqs = freqs  # this adds the frequencies to the ModelClass

# Best Fit Parameter Inputs
params_BM = fn.npFloat([1])
params_D = fn.npFloat([1, 1.59, 19.6])
parameters = fn.concatenate(params_BM, params_D)



###############################################################################################
                                      # Data

#########################################
          # Theoretical Models
# ----------------------------
          # BMode Model
BMode_rawdata = fn.extractData(lamda_data_filename, 3)[lDict['refMin']:lDict['refMax']]
BMode_bindata = fn.binData(BMode_rawdata, lDict['lList'], lDict['lStep'])
BMode = ModelClass('BMode',       {'model': 'raw',
                                  'equation': eqn.BModeSignal,
                                  'eqinput': BMode_rawdata, 'params': params_BM,
                                  'xaxis': lDict['lList']})
BMode.add_model(                 {'model': 'bin',
                                  'equation': eqn.BModeSignal,  # is this the equation?
                                  'eqinput': BMode_bindata, 'params': params_BM,
                                  'xaxis': lDict['lBinCent']})

# ----------------------------
           # Dust Model
Dust = ModelClass('Dust',         {'model': 'raw',
                                  'equation': eqn.dustSignal,
                                  'eqinput': [lDict['lList'], const, None], 'params': params_D,
                                  'xaxis': lDict['lList']})
Dust.add_model(                  {'model': 'bin',
                                  'equation': eqn.dustSignal,  # is this the equation or fn.binData?
                                  'eqinput': [lDict['lBinCent'], const, None], 'params': params_D,
                                  # 'data': fn.bin2Data(Dust.raw.data, lDict['lList'], lDict['lStep']),  # binning data & re-evaluating are < 1% different
                                  'xaxis': lDict['lBinCent']})

# print Dust.bin.data[0], Dust.bin.data[0], Dust.add_fitdata('bin', None), Dust.finalize_fitdata('bin'), Dust.bin.fitdata[0][:,0]

# Theory
Theory = ModelClass('Theory',     {'model': 'raw',
                                  'equation': eqn.totalSignal,
                                  'eqinput': [lDict['lList'], const, None, BMode_rawdata], 'params': parameters,
                                  'xaxis': lDict['lList']},
                                 {'model': 'bin',
                                  'equation': eqn.totalSignal,
                                  'eqinput': [lDict['lBinCent'], const, None, BMode_bindata], 'params': parameters,
                                  'xaxis': lDict['lBinCent']})


#########################################
          # Observed Data
# Measured
measured_error = fn.extractData(measured_bin_filename, 2)  #[lDict['refMin']:lDict['refMax']]
Measured = ModelClass('Measured', {'model': 'raw',
                                  'data': fn.dudata(Theory.raw.data, 0.1),
                                  'error': lambda x: fn.error(x, pct=0.02, mtd=1),  # *** replace with real errors ***
                                  'xaxis': lDict['lList']})
Measured.add_model(              {'model': 'bin',
                                  'data': fn.bin2Data(Measured.raw.data, lDict['lList'], lDict['lStep']),
                                  'error': Measured.raw.error,
                                  'xaxis': lDict['lBinCent']})


#########################################
           # Processed Data
# Fit Coefficients
ChiSqFit = fit.ChiSqOpt(Measured.bin, Theory.bin.params, BMode.bin, Dust.bin)
Theory.bin.params, BMode.bin.params, Dust.bin.params = ChiSqFit.fit('combo', method='nelder-mead', tol=1e-4)  # *** retire outputs and put in a fn.updata_params? ***
print '\nBest Fit Parameters:', Theory.bin.params, '\n'
print np.array(ChiSqFit.chisq[-1])

# Best Fit
BestFit = ModelClass('Best Fit')
BestFit.add_model(               {'model': 'bin',  # different than Theory since parameters changed
                                  'equation': eqn.totalSignal,
                                  'eqinput': [lDict['lBinCent'], const, None, BMode_bindata], 'params': Theory.bin.params[:],
                                  'error': lambda x: fn.error(x, pct=0.02, mtd=1),
                                  'xaxis': lDict['lBinCent']})

# Temporary:
Theory.raw.params, BMode.raw.params, Dust.raw.params = Theory.bin.params, BMode.bin.params, Dust.bin.params

# Updating Fit Data & d_Fit Data
fn.update_fitdata((BMode, 'raw', 'bin'), (Dust, 'raw', 'bin'), (Theory, 'raw', 'bin'))
fn.update_d_fitdata((BMode, 'raw', 'bin'), (Dust, 'raw', 'bin'), (Theory, 'raw', 'bin'))


# ----------------------------------------
            # Plotting Stuff
# Fit Plots:
for index in range(len(freqs)):
  pltfn.plotErrorbar(lDict['lBinCent'], Measured.bin, BMode.bin, Dust.bin, BestFit.bin, Theory.bin, index)
pltfn.plotMeasured(Measured, lDict)
pltfn.plotBMode(BMode, lDict)

###############################################################################################
                                    # Monte Carlo
# performing MC
MCData = fit.MCChiSqFit(Theory.bin, Theory.bin.params, BMode.bin, Dust.bin, iterate=20, method='nelder-mead')
# plotting histograms and parameter correlations
pltfn.plotHisto(MCData, [BMode.bin.params, Dust.bin.params])
pltfn.plotCorrelation(MCData, [BMode.bin.params, Dust.bin.params])
# pltfn.plotMCMeasured(MCData)



###############################################################################################
                                    # Outputs

# ----------------------------------------
            # seeing some outputs
print "\nparams {}".format(MCData.params)
print "Difference between parameters", Theory.bin.params-MCData.params
print "BMode Ampl - Ampl Mean = {}".format(BMode.bin.params[0] - MCData.params[0])
print "Dust Ampl - Ampl Mean = {}".format(Dust.bin.params[0] - MCData.params[1])

# New-Style Plots
# figure, axis = pltfn.makePlot(lDict['lList'], Measured.raw, BestFit.raw)
# figure2, axis2 = pltfn.makePlot(lDict['lList'], Measured.raw, BMode.raw, Dust.raw, BestFit.raw, Theory.raw)
# pltfn.ComparePlots(axis, axis2)


''' NEED TO DO:
- change dust equation to something else
- Make better best fitter for the correlation plots
- Switch to lmfit
- ADD PRIORS (BOUNDS )
- Switch structure so that each add_model has all the add_data, etc... inide it.
  - new inheritance structure
- split up stuff into mini-classes
- Set up a suite of self-plotters
'''
