# -*- coding: utf-8 -*-
# All the Dictionaries

import numpy as np
import Functions as fn                # file with most of the functions


##########################################################################################
                                    # lDict

def make_lDict(data=None, **kw):
    ''' making a dictionary with all the angle l values as the x-axis
        Startpoint: either filename str or data list
        data list must be a list or tuple. will be taken as [int(x) for x in data]
        if no inputs other than startpoint is given, then comes w/ some prepopulated values
    '''

    # data List
    if isinstance(data, str):
        dataList = np.array([int(x) for x in fn.extractData(data, 0)])
    elif isinstance(data, (list, np.ndarray, tuple)):
        dataList = np.array([int(x) for x in data])
    else:  # if no dataList is given in any recognizable format
        print('Error: need a filename or dataList')
        return

    # ---- Base ----
    # base steps
    if 'baseStep' in kw:
        baseStep = kw.get('baseStep')
    else:
        baseStep = dataList[1]-dataList[0]  # fix this
    # data file's minimum l value
    if 'dataMin' in kw:
        dataMin = kw.get('dataMin')
    else:
        dataMin = dataList[0]
    # data file's maximum l value
    if 'dataMax' in kw:
        dataMax = kw.get('dataMax')
    else:
        dataMax = dataList[-1]

    # l-step for binning
    if 'lStep' in kw:
        lStep = kw.get('lStep')
    else:
        lStep = 20*baseStep
    # min l value used
    if 'lMin' in kw:
        lMin = kw.get('lMin')
    else:
        lMin = 50
    # floor of lMin, incase lMin is not an integer
    if 'refMin' in kw:
        refMin = kw.get('refMin')
    else:
        refMin = lMin-dataMin
    # max l value used (including)
    if 'lMax' in kw:
        lMax = kw.get('lMax')
    else:
        lMax = 500
    # floor of lMax, incase lMax is not an integer
    if 'refMax' in kw:
        refMax = kw.get('refMax')
    else:
        refMax = lMax-dataMin

    # x-axis for everything
    if 'lList' in kw:
        lList = np.array(kw.get('lList'))
    else:
        lList = np.array(dataList[refMin:refMax])

    # # reference list of bins
    # if 'refBins' in kw:
    #     refBins = np.array(kw.get('refBins'))
    # else:
    #     refBins = fn.makeList(refMin, refMax+lStep, lStep)
    # # list of bins
    # if 'lBins' in kw:
    #     lBins = np.array(kw.get('lBins'))
    # else:
    #     lBins = fn.makeList(lMin, lMax+lStep, lStep)

    # # reference bin centers
    # if 'refBinCent' in kw:
    #     refBinCent = np.array(kw.get('refBinCent'))
    # else:
    #     refBinCent = fn.binCenter(refBins)
    # # l bin centers
    # if 'lBinCent' in kw:
    #     lBinCent = np.array(kw.get('lBinCent'))
    # else:
    #     lBinCent = fn.binCenter(lBins)

    lDict = {
             # steps
             'baseStep': baseStep,
             'lStep': [lStep],  # *** NEED List for binning
             # minimums
             'dataMin': dataMin,
             'lMin': lMin,
             'refMin': refMin,
             # maximums
             'lMax': lMax,
             'dataMax': dataMax,
             'refMax': refMax,
             # list
             'dataList': dataList,
             'lList': lList,
             # binned list
             # 'refBins': refBins,
             # 'lBins': lBins,
             # # centers of binned list
             # 'refBinCent': refBinCent,
             # 'lBinCent': lBinCent
            }
    return lDict


##########################################################################################
                                    # Constants

def make_Const(**kw):
    ''' Reference dictionary of constants
        If no inputs given, comes prepopulated with standard values
    '''
    # if 'nu1' in kw:
    #     nu1 = np.float64(kw.get('nu1'))
    # else:
    #     nu1 = np.float64(90.0*(10**9))
    # if 'nu2' in kw:
    #     nu2 = np.float64(kw.get('nu2'))
    # else:
    #     nu2 = np.float64(150.0*(10**9))
    if 'nu0' in kw:
        nu0 = np.float64(kw.get('nu0'))
    else:
        nu0 = np.float64(220.)
    if 'h' in kw:
        h = np.float64(kw.get('h'))
    else:
        h = np.float64(6.62606957*(10**-34))
    if 'c' in kw:
        c = np.float64(kw.get('c'))
    else:
        c = np.float64(299792458.)
    if 'k' in kw:
        k = np.float64(kw.get('k'))
    else:
        k = np.float64(1.3806488*(10**-23))
    if 'Tcmb' in kw:
        Tcmb = np.float64(kw.get('Tcmb'))
    else:
        Tcmb = np.float64(2.7)
    if 'TDust' in kw:
        TDust = np.float64(kw.get('TDust'))
    else:
        TDust = np.float64(19.6)

    constants = {
             # 'nu1': nu1,      # frequency 1
             # 'nu2': nu2,      # frequency 2
             'nu0': nu0,      # reference frequency
             'h': h,          # Planck's constant
             'c': c,          # speed of light
             'k': k,          # Boltzmann constant
             'Tcmb': Tcmb,    # vacuum temp
             'TDust': TDust,  # dust temp
             'List': [h, c, k, Tcmb, TDust]
            }
    return constants



##########################################################################################
                                    # Frequencies
def make_frequencies(*args, **kw):
    ''' Reference dictionary of constants
        If no inputs given, comes prepopulated with standard values
    '''
    freqs_Dict = {}
    freqs = []
    for index, value in enumerate(args):
        freqs_Dict['nu{}'.format(index+1)] = value
        freqs.append(value)
    for value, key in enumerate(kw):
        freqs_Dict[key] = value
        freqs.append(value)
    return freqs_Dict, freqs
