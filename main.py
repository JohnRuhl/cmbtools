# -*- coding: utf-8 -*-
# Main execution file

'''
Document Info
- Function of L
'''

import numpy as np
import fittingFunctions as fit         # best fit functions
import ModelEquations as meqn               # file with most of the used ModelEquations
import Functions as fn                 # file with most of the functions
import plottingFunctions as pltfn      # handles all plotting
import Dictionaries as dct             # lDict and constants
# import lmfit as lmfit
from Classes import ModelClass          # all the classes



###############################################################################################
                            # Files, Constants & Parameters

# Files
lamda_data_filename = 'LAMDA Data'
BB_errors = 'BB_errors.txt'

# making a dictionary with all the l values as the x-axis
lDict = dct.make_lDict(lamda_data_filename, lMin=50, lMax=500)

# Reference dictionary of constants
const = dct.make_Const()

# Reference dictionary of frequencies
freqDict, freqs = dct.make_frequencies(90.e9, 150.e9, 212.e9, 242e9)
ModelClass.freqs = freqs  # this adds the frequencies to the ModelClass

# Best Fit Parameter Inputs
R = float(.01)
params_BM = fn.npFloat([1.])
params_D = fn.npFloat([1., 1.59])  # , 19.6   # *** Add priors ***
parameters = fn.concatenate(params_BM, params_D)



###############################################################################################
                                      # Data

#########################################
          # Theoretical Models
# ----------------------------
          # BMode Model
BMode_rawdata = fn.extractData(lamda_data_filename, 3)[lDict['refMin']:lDict['refMax']] # *** I don't think this is R=1 ***
BMode_bindata, _ = fn.binData(BMode_rawdata, lDict['lList'], lDict['lStep'])
BMode = ModelClass('BMode',      {'model': 'raw',
                                  'equation': meqn.BModeSignal,
                                  'eqinput': BMode_rawdata, 'params': params_BM,
                                  'xaxis': lDict['lList']})
BMode.add_model(                 {'model': 'bin',
                                  'equation': meqn.BModeSignal,  # is this the equation?
                                  'eqinput': BMode_bindata, 'params': params_BM,
                                  'xaxis': lDict['lBinCent']},
                                 {'model': 'realize',
                                  'equation': meqn.BModeSignal,  # is this the equation?
                                  'eqinput': BMode_bindata, 'params': params_BM })
# print BMode.info # BMode._info.overwrite('10') # print BMode._info.info # print BMode.info

# ----------------------------
           # Dust Model
Dust = ModelClass('Dust',        {'model': 'raw',
                                  'equation': meqn.dustSignal,
                                  'eqinput': [lDict['lList'], const, None], 'params': params_D })
Dust.add_model(                  {'model': 'bin',
                                  'equation': meqn.dustSignal,  # is this the equation or fn.binData?
                                  'eqinput': [lDict['lBinCent'], const, None], 'params': params_D },
                                 {'model': 'realize',
                                  'equation': meqn.dustSignal,  # is this the equation or fn.binData?
                                  'eqinput': [lDict['lBinCent'], const, None], 'params': params_D })

# print Dust.bin.data[0], Dust.bin.data[0], Dust.add_fitdata('bin', None), Dust.finalize_fitdata('bin'), Dust.bin.fitdata[0][:,0]

# Theory
realize_error = [fn.extractData(BB_errors, i) for i in range(len(freqs))]  # *** Build a better extractor ***
Theory = ModelClass('Theory',    {'model': 'raw',
                                  'equation': meqn.totalSignal,
                                  'eqinput': [lDict['lList'], const, None, BMode_rawdata], 'params': parameters },
                                 {'model': 'bin',
                                  'equation': meqn.totalSignal,
                                  'eqinput': [lDict['lBinCent'], const, None, BMode_bindata], 'params': parameters },
                                 {'model': 'realize',
                                  'equation': meqn.totalSignal,
                                  'eqinput': [lDict['lBinCent'], const, None, (BMode_bindata*float(R))], 'params': parameters,  # *** CHANGE THIS. can't have R just here ***
                                  'd_data': realize_error})  # *** CHANGE to 'error', but need to change how d_data calls error***


###############################################################################################
                                    # Monte Carlo
# performing MC
MCData = fit.MCChiSqFitClass(Theory.realize, Theory.realize.params, BMode.bin, Dust.bin, iterate=1000, method='nelder-mead')
MCData.runMC()
# plotting histograms and parameter correlations
pltfn.plotHisto(MCData, [[R], Dust.bin.params])
pltfn.plotCorrelation(MCData, [[R], Dust.bin.params])
# pltfn.plotMCMeasured(MCData)

''' NEED TO UPDATE THE params?
    NEED TO GET more info out of the Monte Carlo
'''
# # Updating Fit Data & d_Fit Data
# fn.update_fitdata((BMode, 'raw', 'bin'), (Dust, 'raw', 'bin'), (Theory, 'raw', 'bin'))
# fn.update_d_fitdata((BMode, 'raw', 'bin'), (Dust, 'raw', 'bin'), (Theory, 'raw', 'bin'))


###############################################################################################
                                    # Outputs

# ----------------------------------------
            # seeing some outputs
print "\nparams {}".format(MCData.params)
# print "Difference between parameters", Theory.bin.params-MCData.params
print 'chisq:', MCData.chisq
# print "BMode Ampl - Ampl Mean = {}".format(BMode.bin.params[0] - MCData.params[0])
# print "Dust Ampl - Ampl Mean = {}".format(Dust.bin.params[0] - MCData.params[1])




############################################################
          # Realization of a Model for Graphing
# Realization
Theory.add_model(                {'model': 'realize1',  # *** FIX ***
                                  'data': fn.noisyDataList(Theory.realize.fitdata[0], Theory.realize.d_fitdata[-1]),
                                  'fitdata': fn.noisyDataList(Theory.realize.fitdata[-1], Theory.realize.d_fitdata[-1]),
                                  'params': Theory.bin.params[:],
                                  'd_data': realize_error })
# This is techically in measured_log with the error from Theory.realize.d_fitdata

# Fit Coefficients
ChiSqFit = fit.ChiSqOpt(Theory.realize, Theory.realize1.params, BMode.realize, Dust.realize)
Theory.realize1.params = ChiSqFit.fit('list', method='nelder-mead', tol=1e-4)  # *** retire outputs and put in a fn.update_params? ***
print '\nBest Fit Parameters:', ChiSqFit.params
print np.array(ChiSqFit.chisq)

# Best Fit
Theory.add_model(                {'model': 'grapher',
                                  'equation': meqn.totalSignal,
                                  'eqinput': [lDict['lBinCent'], const, None, BMode_bindata*float(R)], 'params': ChiSqFit.params[:]})

# ----------------------------------------
            # Plotting Stuff
# Fit Plots:
for index in range(len(freqs)):
    pltfn.plotErrorbar(lDict['lBinCent'], Theory.realize1, BMode.bin, Dust.bin, Theory.grapher, Theory.bin, index, freqs)
pltfn.plotMeasured(Theory.realize, lDict)
pltfn.plotBMode(BMode, lDict)



''' NEED TO DO:
- Make better best fitter for the correlation plots
- Switch to lmfit
  - ADD PRIORS (BOUNDS )
- Set up a suite of self-plotters
'''
