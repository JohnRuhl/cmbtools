# -*- coding: utf-8 -*-
# Main execution file

'''
Document Info
- Function of L
'''

import numpy as np
import matplotlib.pyplot as plt
import FittingFunctions as fit         # best fit functions
import ModelEquations as meqn          # file with most of the used ModelEquations
import RealizationEquations as reqn    # file with most of the used realization equations
import Functions as fn                 # file with most of the functions
# import plottingFunctions as pltfn      # handles all plotting
import Dictionaries as dct             # lDict and constants
# import lmfit as lmfit
from Classes import ModelClass          # all the classesf


###############################################################################################
                            # Files, Constants & Parameters

# Files
# lamda_data_filename = 'LAMDA Data'
# BB_errors = 'BB_errors.txt'

# making a dictionary with all the l values as the x-axis
lDict = dct.make_lDict(list(range(int(1e3))), lMin=50, lMax=100)
# Reference dictionary of constants
const = dct.make_Const()
# Reference dictionary of frequencies
freqDict, freqs = dct.make_frequencies(90., 150., 212., 242)
ModelClass.freqs = freqs  # adds freqs to ModelClass

# Model Equation Parameter Inputs
mparams, mparams_sq, mparams_sin = meqn.makeParamInputs([1., 2, 4, 80, 10, -np.pi/2])

# Realization Equation Parameter Inputs
rparams, rparams_sq, rparams_sin = reqn.makeParamInputs([2., 3., 5., 100., 10.])


###############################################################################################
                                      # Data

#########################################
          # Theoretical Models
# ----------------------------
          # Poly Model
Poly = ModelClass("poly", {"model": "model", "equation": {"func": meqn.polySignal, "inputs": (lDict["lList"], "freqs"), "params": mparams_sq}})
Poly.add_model({"model": "rlz", "equation": {"func": reqn.polySignal, "inputs": (lDict["lList"], "freqs"), "params": rparams_sq},
                "error": {"func": reqn.polyError, "inputs": (lDict["lList"], "freqs"), "params": rparams_sq}})

# ----------------------------
          # Sin Model
Sin = ModelClass("sin", {"model": "model", "equation": {"func": meqn.sineSignal, "inputs": (lDict["lList"], "freqs"), "params": mparams_sin}},
                 {"model": "rlz", "equation": {"func": reqn.sineSignal, "inputs": (lDict["lList"], "freqs"), "params": rparams_sin}})
                 # "error": {"func": lambda x, *args: 2., "inputs": None, "params": (None,)}})

# ----------------------------
          # Theory Model
Theory = ModelClass("sin", {"model": "model", "equation": {"func": meqn.totalSignal, "inputs": (lDict["lList"], "freqs"), "params": mparams}})
Theory.add_model({"model": "rlz", "equation": {"func": reqn.totalSignal, "inputs": (lDict["lList"], "freqs"), "params": rparams},
                    "error": {"func": reqn.totalError, "inputs": (lDict["lList"], "freqs"), "params": rparams}})

# ###############################################################################################
#                                     # Fitting
ChiData = fit.MonteCarloClass(Theory.rlz, Theory.model.eqn.params, Poly.model, Sin.model)

x = [fitwith.eqn.inputs for fitwith in ChiData.fitwith]
print(x)
print(ChiData.freqmodel(x, *ChiData.params))

p, C, full = ChiData.runMC(iterate=1e3, print_out=True)
# # ChiData.finalizeOpt()
# # print(ChiData.rlz.eqn.evaln)


# # ###############################################################################################
# #                                     # Outputs

# fig = plt.figure()
# ax = fig.add_subplot(211, title=r"$title$", xlabel=r"$xlabel$", ylabel=r"$ylabel$")
# ax.errorbar(lDict["lList"], ChiData.measured.evaln, yerr=ChiData.measured.d_evaln, fmt="c", label="rlz")
# [ax.plot(lDict["lList"], ChiData.fitwith[i].evaln, "k-.", label="model{}".format(i)) for i in range(len(ChiData.fitwith))]
# ax.plot(lDict["lList"], np.sum([ChiData.fitwith[i].evaln for i in range(len(ChiData.fitwith))], axis=0), "go", label="cumul")
# ax.plot(lDict["lList"], Theory.model.evaln, "g--", label="model")
# ax.minorticks_on()
# ax.legend(loc="best", numpoints=1)

# ax2 = fig.add_subplot(212, title=r"$title$", xlabel=r"$xlabel$", ylabel=r"$ylabel$")
# [ax2.plot(lDict["lList"], np.abs((ChiData.fitwith[i].eqn.evaln-ChiData.measured.evaln)/ChiData.measured.evaln), "r--", label="model{}".format(i)) for i in range(len(ChiData.fitwith))]
# # ax2.plot(lDict["lList"], np.abs((ChiData.fitwith.eqn.func(ChiData.fitwith.eqn.inputs, *ChiData.params)-ChiData.measured.evaln)/ChiData.measured.evaln), "r--", label="model")

# plt.savefig("outputs.png")
