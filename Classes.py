# -*- coding: utf-8 -*-
# All the Classes

import numpy as np
import matplotlib.pyplot as plt
from types import FunctionType, NoneType

'''
PROBLEM with _model. Can't get it to access all the data in DataClass
'''


class DataClass(object):
    ''' must give a name
        can give a model (ex raw), fullname, and info
        - all model properties must be args
        - all non-model properties must be kw
        - Model Properties:
            - must be an arg
            - must be a dict
            - must have a material property types with key 'model'
            - may have data with key 'data
            - may have params with key 'params'
            - may have error with key 'error'
            ex: {'model': 'modelname', 'data': [data], 'error': 2}
        - fullname, and info must be strings
    '''

    def __init__(self, name, *args, **kw):
        super(DataClass, self).__init__()
        self.name = name
        self.fullname = name
        # Models
        for i, Dict in enumerate(args):
            self.add_model(Dict)

        # Non-Model Properties
        if 'fullname' in kw:                     # full name
            self.fullname = kw.get('fullname')
        if 'info' in kw:                         # info is any note. should be str
            self.info = kw.get('info')
        else:
            self.info = None

    def add_info(self, info=None):
        'adds info to a given model. if no info passed, adds as None'
        self.info = info

    def add_fullname(self, fullname=None):
        'adds fullname to a given model. if no fullname passed, adds as None'
        self.fullname = fullname

    def add_model(self, *Dicts):
        ''' adds a data model
            needs: modelname in main dict
            can be given any of (in dict): data, params, & error
            data must be a list
            ex: add_model({'matProp': 'thermCond', 'name':'NIST', 'eq': lambda x: x, 'eqinput': [1, 23]})
        '''

        for Dict in Dicts:
            if not isinstance(Dict, (dict, NoneType)):
                if Dict is None:
                    Dict = {}
                else:
                    print 'Error: dict needs to be type dict or None'
                    return

            # Model
            if 'model' in Dict:
                model_name = Dict['model']
            else:
                print 'Error: need to pass a model name'
                return
            if isinstance(model_name, str):
                setattr(self, model_name, _model(self.name, model_name))  # *** temporary solution? ***
            else:
                print self.name, type(model_name)
                raise TypeError('model is not the right type')

            # Equation
            if 'equation' in Dict:
                equation = Dict['equation']
            else:
                equation = None
            self.add_equation(model_name, equation)

            # eqinput
            if 'eqinput' in Dict:
                eqinput = Dict['eqinput']
            else:
                eqinput = None
            self.add_eqinput(model_name, eqinput)  # *** prob need some ifinstance() check ***

            # params
            if 'params' in Dict:
                parameters = Dict['params']
            else:
                parameters = None
            self.add_params(model_name, parameters)

            # Data
            if 'data' in Dict:
                data = Dict['data']
            else:
                try:
                    data = equation(eqinput, parameters)
                except:
                    data = None
            self.add_data(model_name, data)

            # Error
            # needs to be after eq, eqinput, params, and data since it can use those to calculate the error
            if 'error' in Dict:
                error = Dict['error']
            else:
                error = None
            self.add_error(model_name, error)

            # Error in Data
            if 'd_data' in Dict:
                d_data = Dict['d_data']
            else:
                d_data = None
            self.add_d_data(model_name, d_data)

            # Fit Data
            if 'fitdata' in Dict:
                fitdata = Dict['fitdata']
            else:
                fitdata = None
            self.add_fitdata(model_name, fitdata)

            # Error in Fit Data
            if 'd_fitdata' in Dict:
                d_fitdata = Dict['d_fitdata']
            else:
                d_fitdata = None
            self.add_d_fitdata(model_name, d_fitdata)

            # xaxis
            if 'xaxis' in Dict:
                xaxis = Dict['xaxis']
            else:
                xaxis = None
            self.add_xaxis(model_name, xaxis)

    def add_equation(self, model=None, equation=None):
        ''' add an equation to a model.
            called material.add_equation(property, modelname, equation)
            note: modelname must be a string
        '''
        if isinstance(model, str):                              # check if model string
            model = getattr(self, model)
        else:                                                   # raise error if model not string
            print self.name, type(model)
            raise TypeError('model not a string')

        if isinstance(equation, (FunctionType, NoneType, list, np.ndarray)):
            model.equation = equation
        else:                                                   # raise error if not right type
            print self.name, type(equation)
            raise TypeError('equation is not the right type')

    def add_eqinput(self, model=None, eqinput=None):
        ''' adds equation input to a model.
            called as material.add_params(property, modelname, params)  (with either as a kwarg)
            Cannot give something as both an arg and a kwarg!
        '''
        if isinstance(model, str):                              # check if model string
            model = getattr(self, model)
        else:                                                   # raise error if model not string
            print self.name, type(model)
            raise TypeError('model not a string')

        if isinstance(eqinput, (NoneType, list, np.ndarray, tuple)):
            model.eqinput = eqinput
        else:                                                   # raise error if not right type
            print self.name, type(eqinput)
            raise TypeError('eqinput is not the right type')

    def add_params(self, model=None, params=None):
        ''' adds equation input parameters to a model.
            called as material.add_params(property, modelname, params)  (with either as a kwarg)
            Cannot give something as both an arg and a kwarg!
        '''
        if isinstance(model, str):                              # check if model string
            model = getattr(self, model)
        else:
            print self.name, type(model)
            raise TypeError('model not a string')

        if isinstance(params, (NoneType, list, np.ndarray, dict)):
            model.params = params
        else:                                                   # raise error if not right type
            print self.name, type(params)
            raise TypeError('params is not the right type')

    def add_error(self, model=None, error=None):
        ''' adds error equation to a model.
            Cannot give something as both an arg and a kwarg!
        '''
        if isinstance(model, str):                              # check if model string
            model = getattr(self, model)
        else:
            print self.name, type(model)
            raise TypeError('model not a string')

        if isinstance(error, (int, float, list, np.ndarray, FunctionType, NoneType)):
            model.error = error
        else:                                                   # raise error if not right type
            print self.name, type(error)
            raise TypeError('error is not the right type')

    def add_data(self, model=None, data=None):
        ''' adds the data to a model.
            called material.add_data(property, modelname, equation)  (with either as a kwarg)
            note:  modelname must be a string
            Cannot give something as both an arg and a kwarg!
        '''
        if isinstance(model, str):                              # check if model string
            model = getattr(self, model)
        else:                                                   # raise error if model not string
            print self.name, type(model)
            raise TypeError('model not a string')

        if isinstance(data, (list, np.ndarray, NoneType)):
            if data is None:
                try:                                            # assessing if data from equation
                    data = model.equation(model.eqinput, model.params)
                except:                                         # guess the equation had a problem
                    pass                                        # data stays as None
            data = np.array(data)
            model.data = data
        else:                                                   # raise error if not right type
            print self.name, type(data)
            raise TypeError('data is not the right type')

    def add_d_data(self, model=None, d_data=None):
        ''' adds error to the data.
            adding d_data as None makes it first try to add as model.error(model.data)
            not giving a d_data makes it first try to add as model.error(model.data) then as model.error
        '''
        if isinstance(model, str):                              # check if model string
            model = getattr(self, model)
        else:                                                   # raise error if model not string
            print self.name, type(model)
            raise TypeError('model not a string')

        if isinstance(d_data, (int, float, list, np.ndarray, NoneType)):  # if same for all freqs
            if d_data is None:
                try:                                            # evaulate for the error in the data
                    d_data = np.array(model.error(model.data))
                except TypeError:                               # already list, so add error in data
                    d_data = np.array(model.error)
            else:                                               # d_data was given
                d_data = np.array(d_data)
            model.d_data = d_data
        else:                                                   # raise error if not right type
            print self.name, type(d_data)
            raise TypeError('d_data is not the right type')

    def add_fitdata(self, model=None, fitdata=None):
        ''' adds the data, evaulated with the best fit parameters, to a model.
            called material.add_fitdata(property, modelname, equation)
            note:  modelname must be a string
        '''
        if isinstance(model, str):                              # check if model string
            model = getattr(self, model)
        else:                                                   # raise error if model not string
            print self.name, type(model)
            raise TypeError('model not a string')
        if isinstance(fitdata, (list, np.ndarray, NoneType, dict)):
            if fitdata is None:
                try:                                            # evaulate equation with parameters
                    fitdata = np.array(model.equation(model.eqinput, model.params))
                except:                                         # already list, so add fitdata
                    print model.name, model.modelname, 'fitdata added as None'
            else:                                               # fitdata was given
                fitdata = np.array(fitdata)
            model.fitdata = fitdata
        else:                                                   # raise error if not right type
            print self.name, type(fitdata)
            raise TypeError('fitdata is not the right type')

    def add_d_fitdata(self, model=None, d_fitdata=None):
        ''' adds error to fitdata.
        '''
        if isinstance(model, str):                              # check if model string
            model = getattr(self, model)
        else:                                                   # raise error if model not string
            print self.name, type(model)
            raise TypeError('model not a string')

        if isinstance(d_fitdata, (list, np.ndarray, NoneType)):
            if d_fitdata is None:
                try:                                            # evaulate equation with parameters
                    d_fitdata = np.array(model.error(model.fitdata))
                except:                                         # already list, so add error in data
                    d_fitdata = model.d_data
                    print model.name, model.modelname, 'd_fitdata added as d_data'
            else:                                               # d_fitdata already given
                d_fitdata = np.array(d_fitdata)
            model.d_fitdata = d_fitdata
        else:                                                   # raise error if not right type
            print self.name, type(d_fitdata)
            raise TypeError('d_fitdata is not the right type')

    def add_xaxis(self, model=None, xaxis=None):
        ''' adds an axis to a model.
            called material.add_axis(modelname, equation)  (with either as a kwarg)
            note:  modelname must be a string
            Cannot give something as both an arg and a kwarg!
        '''
        if isinstance(model, str):                              # check if model string
            model = getattr(self, model)
        else:                                                   # raise error if model not string
            print self.name, type(model)
            raise TypeError('model not a string')
        if isinstance(xaxis, (list, np.ndarray, NoneType)):
            model.xaxis = xaxis
        else:                                                   # raise error if not right type
            print self.name, type(xaxis)
            raise TypeError('xaxis is not the right type')


class _model(DataClass):
    ''''''
    def __init__(self, name, modelname):                       # this solution is temporary
        self.name = name                                       # this solution is temporary
        self.modelname = modelname
        # data as equation and inputs and error
        self.equation = None
        self.eqinput = None
        self.params = None
        self.error = None
        # original data as a list for easy access
        self.data = None
        self.d_data = None
        # data evaulated with best fit parameters, as a list for easy access
        self.fitdata = None
        self.d_fitdata = None
        # plotting axis
        self.xaxis = None

    def _plot(self, fig=None, ax=111, xaxis=None, **kw):  # come up with better name
        ''' creates a plot of self.
            can create a figure not passing any args
            can add to an existing figure by passing a figure
                can specify which axis with a second arg as the number
            can add to an existing subplots py passing the fig and the axis
        '''
        if 'name' in kw:
            name = kw.get('name')
        else:
            name = self.name  # this solution is temporary
        if fig is None:                     # no figure given
            fig = plt.figure()
            ax = plt.subplot(ax)
        elif isinstance(ax, (int, float)):  # figure given
            ax = fig.add_subplot(ax)
        else:                               # figure and axis give
            pass
        if xaxis is None:
            xaxis = self.xaxis
        name = ax.errorbar(xaxis, self.data, yerr=self.error, ls='-', label=name)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='upper right')
        return fig



# equation_dict = {}
# for key in self.frequencies:                        # iterate thru freqs and assign same value
#     equation_dict[key] = equation
