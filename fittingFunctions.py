# -*- coding: utf-8 -*-
# All the Fitting Functions

import numpy as np
import scipy.optimize as opt
from Functions import parameterSplit, noisyDataList, update

# ----- Fit Finders -----


# Line of Best Fit (using Least Squares Regression):
def polyfit(xList, yList, degree=1, rsqr=False):
    ''' Line of Best Fit (using Least Squares Regression):
    Description: A function that takes in two lists of equal length, corresponding
    to the x and y coords of a graph, and returns a list containing
    the intercept, slope and r^2 of the line of best fit to that graph
    '''
    coeffs = np.polyfit(xList, yList, degree)
    b = coeffs[0]  # highest degree
    a = coeffs[1]  # lowest degree
    if rsqr is True:
        # calculating correlation coefficient. r = b *(stdev yList / stdev xList)
        r = b * (np.std(yList)/np.std(xList))
        # calculating r-square value
        rSquare = np.square(r)
        return [a, b, rSquare]
    else:
        return [a, b]


# ---------- Chi Square Optimizer ----------
class ChiSqOpt(object):
    """docstring for ChiSqOpt
    """
    def __init__(self, Measured, parameters, *args):
        ''' takes all of the input necessary to do a Chi-Sq  optimization
        '''
        self.measured = Measured
        self.measured_log = [Measured.fitdata[-1]]
        self.freqs = Measured.freqs

        self.params = parameters  # current best fit parameters
        self.params_log = []  # empty list of resulting parameter fits
        self.chisq = []
        self.chisq_log = []  # empty list for the chi-square

        self.args = args  # all the models with which to fit



    # Theory Constructions of the Data for fits
    def residual(self, index, freq, *params):
        '''
            needs parameters as input since they will be changed by the fitter
            make sure the parameters are in the same order as the objects (Bmode parameters with BMode equation)
        '''
        # finding the model
        model = 0
        for i, val in enumerate(self.args):  # iterates through the arguments
            eqinput = val.eqinput[:]  # getting the eqinput for the equation, the [:] prevents linking
            if None in eqinput:  # finding if a freq needs to be given
                none_index = eqinput.index(None)  # finding where freq needs to be
                eqinput[none_index] = freq  # giving eqinput the correct freq
            model += val.equation(eqinput, params[i])  # evaluating the equation

        # finding the residual
        residual = self.measured_log[-1][index] - model  # calling the last one since new ones are made by the Monte Carlo
        return residual

    def chisqfc(self, parameters):
        ''' Calculates the sum of the chisq of all the bands'''
        # preparing inputs
        params_BM, params_d = parameterSplit(parameters)

        chisq_sum = 0
        chisq_log = []
        for index, freq in enumerate(self.freqs):  # iterates through the frequency bands
            # finding the residual
            residual = self.residual(index, freq, params_BM, params_d)
            # finding the error
            d_data = self.measured.d_fitdata[-1][index]
            try:  # checking if d_data exists
                residual/d_data
            except TypeError:  # else assuming 2% error
                print "uh oh"
            # finding the chi-square
            chisq = np.sum((residual/d_data)**2)
            # recording the chi-square
            chisq_log.append(chisq)
            chisq_sum += chisq
        self.chisq_log.append(np.array(chisq_log))
        return chisq_sum

    def fit(self, output='list', **kw):
        ''' needs scipy.optimize imported as opt
        '''
        # selection of optimizaiton method
        if 'method' in kw and 'options' in kw and 'tol' in kw:
            method, options, tol = kw.get('method'), kw.get('options'), kw.get('tol')
        else:
            method, options, tol = self.optimization_options(kw)
        # bounds = ((0, None), (0, None), (0, None), (0, 1e3))

        result = opt.minimize(self.chisqfc, self.params, method=method, options=options, tol=tol)  # bounds=bounds
        self.params = np.array(result.x)
        self.params_log.append(result.x)  # appending params_log to self.params_log
        self.chisq = self.chisq_log[-1]
        self.sigma = self.sigma_error()

        if output == 'list':  # fit will return [params_BM[:], params_d[:]] as a unified list
            return self.params
        elif output == 'split':  # fit will return [[params_BM], [params_d]]
            return parameterSplit(self.params)
        elif output == 'combo':
            split = parameterSplit(self.params)
            return self.params, split[0], split[1]
        elif output == 'none':
            return

    def sigma_error(self):
        sigma = np.std(self.params_log, axis=0, dtype=np.float64)
        return sigma

    @staticmethod
    def optimization_options(kw):
        if 'method' in kw:
            method = kw.get('method')
        else:
            method = 'nelder-mead'  # choose a better method as default?
        if method is None:
            method = 'nelder-mead'
        # giving a full dict of options
        if 'options' in kw:
            options = kw.get('options')
            if options is None:
                options = {'maxiter': 1e5, 'maxfev': 1e5}  # pre-populating options
                if 'maxiter' in kw:  # maximum number of iterations
                    options['maxiter'] = kw.get('maxiter')
                if 'maxfev' in kw: # maximum number of function evaluations
                    options['maxfev'] = kw.get('maxfev')
        else:  # giving options individually
            options = {'maxiter': 1e5, 'maxfev': 1e5, 'xtol': 1e-1}  # pre-populating options
            if 'maxiter' in kw:  # maximum number of iterations
                options['maxiter'] = kw.get('maxiter')
            if 'maxfev' in kw: # maximum number of function evaluations
                options['maxfev'] = kw.get('maxfev')
            if 'xtol' in kw:
                options['xtol'] = kw.get('xtol')

        if 'tol' in kw:
            tol = kw.get('method')
        else:
            tol = None

        return method, options, tol

    def finalize(self):
        self.params_log = np.array(self.params_log)
        self.params = np.mean(self.params_log, axis=0)  # *** CHANGE THIS ***
        self.measured_log = np.array(self.measured_log)
        self.chisq = np.array(self.chisq)
        self.chisq_log = np.array(self.chisq_log)
        self.sigma = np.array(self.sigma)


#######################################################################################
#                                 Monte-Carlo

class MCChiSqFitClass(ChiSqOpt):
    """docstring for MCChiSqFitClass"""
    def __init__(self, Measured, parameters, *args, **kw):
        if 'iterate' in kw:  # checking that the iterator is passed
            self.iterate = kw.get('iterate')
        else:  # no iterator was passed
            self.iterate = 10**4
            # args = (iterate,) + args  # appending the misinterpreted arg back into *args
        super(MCChiSqFitClass, self).__init__(Measured, parameters, *args)  # inheriting the methods from ChiSqOpt

        if 'method' in kw:
            self.method = kw.get('method')
        else:
            self.method = None

    # def fit()
    # def chisqfc()
    # def residual
    # def optimization_options()
    # def finalize

    def runMC(self, **kw):
        print '\nChi-Sq Monte Carlo Iteration #:'
        for i in range(int(self.iterate)):
            print i+1,
            # generating a realization
            data = self.measured.fitdata[-1]  # getting most recent fit
            realization = noisyDataList(data, self.measured.d_fitdata[-1])
            self.measured_log.append(realization)

            params = self.fit(method=self.method)  # Doing the best fit
            self.measured.params = params  # updating the parameters

        # Finalizing Fit
        self.finalize()
        print '\n parameters: {}\nsigma: {}'.format(self.params, self.sigma)
