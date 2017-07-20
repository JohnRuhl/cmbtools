# -*- coding: utf-8 -*-
# All the Classes

import numpy as np
import scipy.optimize as opt
import copy as copy
# import matplotlib.pyplot as plt
# from types import FunctionType, NoneType

'''
Possible Name Changes:
- eqinput to eqin
- erinput to erin

'''


class BaseMethods(object):
    """These are methods all the classes can use"""
    def __init__(self, name):
        self._name = name
        # self.freqs = ModelClass.freqs

    def check(self, arg, name, checktype, raise_error=False):
        """ This checks if the argument is the right type"""
        if isinstance(arg, checktype):
            return True
        else:
            if raise_error is False:
                return False
            else:
                raise TypeError("""{0} in {1} is not the right type.
                {1}.{0} is type {2}\n""".format(name, self.name, type(arg)))

    # def setandcheck(self, arg, name, checking):
    #     """ This checks if the argument is the right type and will set it if it is"""
    #     if self.check(name, arg, checking):
    #         self.arg = arg
    #     else:
    #         print self.name, type(arg)
    #         raise TypeError('{} is not the right type'.format(name))


class ModelClass(BaseMethods):
    """All the objects: dust, CMB, synchotron are instantiated as a model.
    """
    freqs = None

    def __init__(self, name, *args, **kw):
        super(ModelClass, self).__init__(name)
        self.name = name

        # ----------------------
        # Non-Model Properties

        # full name
        if 'fullname' in kw:
            self.fullname = kw.get('fullname')
        else:
            self.fullname = name
        # info
        if 'info' in kw:  # info is any note. should be str
            self.info = kw.get('info')
        else:
            self.info = None

        # ----------------------
        # Model
        for arg in args:
            self.add_model(arg)

    def add_model(self, *Dicts):
        """ adds data model(s)
        ex: {"model":name, "equation":{"func": func, "params": params, "inputs": inputs},
             "error": {"func": func, "params": params, "inputs": inputs}}
        Inputs:
            - name
            - frequencies evaluated at
            - dictionary of the information
        * Note: all arrays must be np.array()

        The model must have an equation function
          func(parameters, args=extra input).
          Input as {"equation": {"func": func, "params": params, "eqinput": eqinput}}
            - equation function
            - equation parameters (can be None)
            - optional extra input.
        """

        for Dict in Dicts:  # gets each model in the list of models
            if not isinstance(Dict, dict):  # check dict, else error
                print(self.name, type(Dict))
                raise TypeError('model is not a dict')
            else:  # it's a dictionary with the model
                if 'model' in Dict:  # checking there is a modelname
                    model_name = Dict['model']
                else:  # no modelname is a problem
                    raise TypeError('Error: need to pass a model name')

                if isinstance(model_name, str):  # checking the model name is a string
                    setattr(self, model_name, ModelMaker(model_name, self.freqs, Dict))
                else:  # modelname must be a string
                    print(self.name, type(model_name))
                    raise TypeError('model is not the right type')


class ModelMaker(object):
    """This makes a specific model.
    Inputs:
        - name
        - frequencies evaluated at
        - dictionary of the information
    * Note: all arrays must be np.arra()

    The model must have an equation function
      func(parameters, args=extra input).
      Input as {"equation": {"func": func, "params": params, "eqinput": eqinput}}
        - equation function
        - equation parameters (can be None)
        - optional extra input.
    """

    def __init__(self, name, freqs, Dict):
        super(ModelMaker, self).__init__()
        self.name = name
        self.freqs = freqs

        # ----------------------
        # The Model
        """the equation is a function for evaluating the model.
        It is of the form func(parameters, args=extra input).
        Input as {"equation": {"func": func, "params": params, "inputs": inputs}}
        parameter array must be np.array() """
        try:
            equation = Dict['equation']
        except:
            KeyError('Error: need to pass an equation')
        else:
            self.eqn = Equation(equation, freqs)

        # ----------------------
        # The Error
        """the error is a function for evaluating the error in the model.
        It is of the form func(parameters, args=extra input).
        Input as {"error": {"func": func, "params": params, "inputs": inputs}}
        parameter array must be np.array() """
        try:
            error = Dict['error']
        except:
            error = None
        self.error = Error(error, freqs)


    # ----------------------
    # Faster Access to Evaln
    # Evaln
    @property
    def evaln(self):
        return self.eqn.evaln
    @evaln.setter
    def evaln(self, evaln):
        self.eqn.evaln = evaln

    # D_Evaln
    @property
    def d_evaln(self):
        return self.error.d_evaln
    @d_evaln.setter
    def d_evaln(self, d_evaln):
        self.error.d_evaln = d_evaln



###############################################################################
class AbstractModelMethods(BaseMethods):
    """docstring for AbstractModelMethods"""
    def __init__(self, Dict, freqs):
        self.func = self.inputs = self.params = self.evaln = None
        self.freqs = freqs

        if Dict is not None:
            # ----------------------
            # Equation Function
            try:
                func = Dict['func']
            except:
                KeyError('Error: need to pass an equation function')
            else:
                self.func = func

            # ----------------------
            # Equation Input
            try:
                inputs = Dict['inputs']
            except:
                inputs = None
            finally:
                self.inputs = np.array(inputs)

            # ----------------------
            # Equation Parameters
            try:
                params = Dict['params']
            except:
                KeyError('Error: need to pass equation parameters')
            else:
                self.params = params
                self.initparams = params
        else:
            self.func = lambda inputs, *params: None
            self.inputs = np.array(["freqs", ])  # *** Should be ["Freqs"]
            self.params = self.initparams = (None,)

    def update_params(self, params, reval=False):
        self.params = params
        if reval is True:
            return self.reval()

    # def freqiter(self, inputs):
    #     inputsList = np.zeros_like(inputs)
    #     ind = (inputs == "freqs")
    #     for i, freq in enumerate(self.freqs):
    #         inputs[ind] = freq
    #         inputsList[i] = inputs
    #     inputs[ind] = "freqs"
    #     return inputsList

    def freqfunc(self, inputs, *params):
        evaln = [0]*len(self.freqs)  # each freq is an elt
        ind = np.where(inputs == "freqs")[0][0]  # index of freqs  # *** FIX ***
        for i, freq in enumerate(self.freqs):  # frequencies
            inputs[ind] = freq  # substitute freq
            evaln[i] = self.func(self.inputs, *params)
        inputs[ind] = "freqs"  # resetting inputs
        return evaln



class Equation(AbstractModelMethods, object):
    """the equation is a function for evaluating the model.
        Form func(parameters, args=extra inputs).
        Input as {"equation": {"func": func, "params": params, "inputs": inputs}}
        parameter array must be np.array() """
    def __init__(self, Dict, freqs):
        super(Equation, self).__init__(Dict, freqs)
        # ----------------------
        # Evaluated Equation
        self.evaln = self.freqfunc(self.inputs, *self.params)

    def reval(self, params=None, update_params=True):
        if params is not None:  # passed new params
            if update_params is False:  # updating params
                return self.freqfunc(self.inputs, *params)
            self.params = params
        self.evaln = self.freqfunc(self.inputs, *self.params)
        return self.evaln


class Error(AbstractModelMethods, object):
    """the error is a function for evaluating the error in the model.
        Form func(parameters, args=extra input).
        Input as {"error": {"func": func, "params": params, "erinput": inputs}}
        parameter array must be np.array() """
    def __init__(self, Dict, freqs):
        super(Error, self).__init__(Dict, freqs)
        # ----------------------
        # Evaluated Error
        self.d_evaln = self.freqfunc(self.inputs, *self.params)

    def reval(self, params=None, update_params=True):
        if params is not None:  # passed new params
            if update_params is False:  # updating params
                return self.freqfunc(self.inputs, *params)
            self.params = params
        self.d_evaln = self.freqfunc(self.inputs, *self.params)
        return self.d_evaln
