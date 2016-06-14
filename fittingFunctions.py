# -*- coding: utf-8 -*-
# All the Fitting Functions

import numpy as np
import scipy.optimize as opt
import scipy.special as sf
import copy as copy
import Classes as clss
from Functions import noisyData

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


class Params(object):
        """docstring for Params"""
        def __init__(self, params, p0, fitwith):
            self.params = params
            self.p0 = p0
            self.fitwith = fitwith
            self.lenparams = self.begs = self.ends = None
            self.setparamprops()

        def setparamprops(self, lenparams=None, begs=None, ends=None):
            """sets lenparams -- list of len of params in fitwith, begs -- starting indices in *params for fitwith,
            and ends -- starting indices in *params for fitwith.
            can be over-written by a kwarg
            """
            if lenparams is None:
                self.lenparams = [len(fitwith.eqn.params) for fitwith in self.fitwith]  # list of len of params
            else:
                self.lenparams = lenparams
            if begs is None:
                self.begs = np.cumsum(np.insert(self.lenparams[:-1], 0, 0))  # starting indices
            else:
                self.begs = begs
            if ends is None:
                self.ends = np.cumsum(self.lenparams)  # ending indices
            else:
                self.ends = ends

        def splitParams(self, *params, update_params=False):
            if params is None:  # just need to split current params
                return self.splitParams_fast((self.begs, self.ends), *self.params)
            else:  # splitting different params
                lenparams = [len(x.eqn.params) for x in self.fitwith]
                begs = np.cumsum(np.insert(lenparams[:-1], 0, 0))  # starting indices
                ends = np.cumsum(lenparams)  # ending indices
                if update_params is True:
                    setparamprops(lenparams=lenparams, begs=begs, ends=ends)
            return self.splitParams_fast((begs, ends), *params)

        def splitParams_fast(self, indices, *params):
            """indices = (begs, ends), begs = [0, 3, 5]"""
            pList = [params[int(i):int(j)] for i, j in zip(indices[0], indices[1])]
            return pList



# ---------- Chi Square Optimizer ----------
class OptimizeClass(object):
    """docstring for OptimizeClass
    """
    def __init__(self, measured, p0, *fitwith):  # *** measured, parameters, *fitwith***
        ''' takes all of the input necessary to do a Chi-Sq optimization
        '''
        self.freqs = clss.ModelClass.freqs  # *** or get from measured.freqs?***
        self.measured = measured

        self.fitwith = fitwith  # all the models with which to fit

        self.chisq = None
        self.cov = None

        self.Params = Params(p0, p0, self.fitwith)

    # -----------------
    # Params shortcuts
    @property
    def params(self):
        return self.Params.params   # return self.fitwith.eqn.params
    @params.setter
    def params(self, params):
        self.Params.params = params  # self.fitwith.eqn.params = params

    def update_fitwith_params(self, *params, indices=None, reval=False):
            if indices is None:
                p = self.Params.splitParams(*params)
            else:  # Should be faster
                p = self.Params.splitParams_fast(indices, *params)
            for i, fitwith in enumerate(self.fitwith):
                fitwith.eqn.update_params(p[i], reval=reval)

    # Composite model of all fitwiths
    def ymodel(self, inputs, *params):
        pList = self.Params.splitParams_fast((self.Params.begs, self.Params.ends), *params)
        models = [np.concatenate(fitwith.eqn.freqfunc(x, *p)) for fitwith, x, p in zip(self.fitwith, inputs, pList)]
        return np.sum(models, axis=0)  # summing all the fitwiths

    def freqmodel(self, inputs, *params):
        pList = self.Params.splitParams_fast((self.Params.begs, self.Params.ends), *params)
        models = np.array([fitwith.eqn.func(x, *p) for fitwith, x, p in zip(self.fitwith, inputs, pList)])
        return np.sum(models, axis=0)  # summing all the fitwiths

    # CURVE FIT
    def curve_fit(self, full_output=False, print_out=False, reval=False):
        """only handles one function model. Doesn't handle errors which are dependant on the parameters"""
        # Set-Up
        x = [fitwith.eqn.inputs for fitwith in self.fitwith]
        y = np.concatenate(self.measured.evaln)
        sigma_y = np.concatenate(self.measured.d_evaln)

        # print(len(self.ymodel(x, *self.Params.p0)))

        # Evaluating
        (p, C) = opt.curve_fit(self.ymodel, x, y, sigma=sigma_y, p0=self.Params.p0)

        # Updating
        self.params = p
        self.update_fitwith_params(*p, reval=reval)  # *** should be doing this? ***
        self.chisq = np.sum((y - self.ymodel(x, *p))**2 / sigma_y**2)
        self.cov = C

        if full_output + print_out > 0:  # either fulloutput or printout is True
            full = {}  # Full output dictionary
            full["params"] = p
            full["cov"] = C
            full["chisq"] = self.chisq
            full["dof"] = len(y) - len(p)
            full["Q"] = sf.gammaincc(0.5 * full["dof"], 0.5 * full["chisq"])

            if print_out is True:
                print("Best fit:")
                for i, val in enumerate(p):
                    print("a{} = {} +/- {}".format(i, val, np.sqrt(C[i, i])))
                print("chisq = {} \nndof = {} \ngoodness of fit = {}".format(full["chisq"], full["dof"], full["Q"]))
            if full_output is True:
                return p, C, full
        return p, C

    def finalizeOpt(self):
        for i, fitwith in enumerate(self.fitwith):
            fitwith.eqn.reval()


class MonteCarloClass(OptimizeClass):
    """docstring for MonteCarloClass"""
    def __init__(self, measured, p0, *fitwith, **kw):
        super(MonteCarloClass, self).__init__(measured, p0, *fitwith)

        self.measured_log = [measured.eqn.evaln.copy()]   # *** update to get newest, better shallow copy ***
        self.chisq_log = []
        # keeping track of Params
        self.Params.params_log = []
        self.cov_log = []
        self.sigma = None

        # ModelClass Realization for consistent packaging
        self.rlz= None

        if 'iterate' in kw:  # checking that the iterator is passed
            self.iterate = kw.get('iterate')
        else:  # no iterator was passed
            self.iterate = 10**4

        # if 'method' in kw:
        #     self.method = kw.get('method')
        # else:
        #     self.method = None

    def runMC(self, **kw):
        if "iterate" in kw:
            iterate = kw.get("iterate")
        else:
            iterate = self.iterate

        print('\nChi-Sq Monte Carlo Iteration #:')
        for i in range(int(iterate)):
            print(i+1, end=", ")
            # generating a realization
            data = self.measured.evaln  # *** update to get newest ***
            realization = noisyData(data, self.measured.d_evaln)
            self.measured_log.append(realization)  # tracking the realization
            params, cov, full = self.curve_fit(full_output=True)

            self.Params.params_log.append(params)  # *** or self.params ***
            self.chisq_log.append(self.chisq)  # *** or self.chisq ***
            self.cov_log.append(cov)

        # Finalizing Fit
        self.finalizeMC(**kw)

        return self.params, self.cov, full  # this full is only the last

    def finalizeMC(self, **kw):
        self.measured_log = np.array(self.measured_log)
        self.chisq_log = np.array(self.chisq_log)

        self.cov_log = np.array(self.cov_log)
        self.cov = np.mean(self.cov_log, axis=0)

        self.Params.params_log = np.array(self.Params.params_log)
        self.sigma = np.std(self.Params.params_log, axis=0, dtype=np.float64)
        self.params = np.mean(self.Params.params_log, axis=0)
        self.update_fitwith_params(*self.params)
        super(MonteCarloClass, self).finalizeOpt()  # finalize from OptimizeClass

        # print("\n", np.array([fitwith.eqn.inputs for fitwith in self.fitwith]))

        # self.rlz = clss.ModelClass("MC Fit", {"model": "rlz", "equation": {"func": self.freqmodel,
        #                                                                    "inputs": np.array([fitwith.eqn.inputs for fitwith in self.fitwith]),
        #                                                                    "params": np.array(self.params)}})
                  # "error": {"func": reqn.polyError, "inputs": lDict["lList"], "params": rparams_sq}})  # *** figure out error ***
        try:
            print_out = kw["print_out"]
        except Exception:
            pass
        else:
            if print_out is True:
                print("Monte Carlo fit:")
                for i, val in enumerate(self.params):
                    print("a{} = {} +/- {}".format(i, val, np.sqrt(self.cov[i, i])))
                print("sigma = {} ".format(self.sigma))
