# -*- coding: utf-8 -*-
# All the Classes

import numpy as np
import matplotlib.pyplot as plt
from types import FunctionType, NoneType

'''
'''


########################################################################################################################

class BaseMethods(object):
    # Equality
    def __eq__(self, other):
        """Checks equality and inequality. Credit to Stack Exchange"""
        if type(other) is type(self.arg):
            return True
        return False

    def __ne__(self, other):
        """Checks equality and inequality. Credit to Stack Exchange"""
        return not self.__eq__(other)

    def __str__(self):
        return str(self.arg)

    # Slices
    def __getitem__(self, index):
        return self.arg[index]

    def __setitem__(self, index, value):
        self.arg[index] = value

    def __len__(self):
        return len(self.arg)

    # Math
    # def __add__(self, other):
    #     return self.arg + other
    # def __sub__(self, other):
    #     return self.arg - other
    # def __mul__(self, other):
    #     returns self.arg * other
    # def __floordiv__(self, other):
    #     return self.arg / other
    # def __mod__(self, other):
    #     return self.arg % other

    def get(self):  # *** Should it be __get__ ? ***
        ''' returns the isinstance.'''
        return self.arg

    def set(self, arg):
        self.arg = arg

    def index(self, value):
        return self.arg.index(value)

    def show(self):
        print self.arg

    def setandcheck(self, arg, name, checking):
        if isinstance(arg, checking):
            self.arg = arg
        else:  # raise error if not right type
            print self.name, type(arg)
            raise TypeError('{} is not the right type'.format(name))


class FunctionMethods(BaseMethods, object):
    """docstring for FunctionMethods"""
    def __call__(self, *args):
        return self.arg.__call__(*args)

    def callable(self):
        return self.arg.__call__


class PropertyMethods(BaseMethods, object):
    """contains some basic homemade methods"""
    def __call__(self):
        ''' returns the isinstance. meant to be overwritten by a local __call__'''
        return self.arg


# Holds some methods for strings
class StringClass(BaseMethods, object):
    def __init__(self, arg):
        self.arg = arg

    def __add__(self, other):
        if isinstance(self.arg, str):
            return self.arg+str(other)
        else:
            return str(other)

    def show(self):
        print self.arg.__str__()

    def append(self, arg):
        self.arg = self.arg+str(arg)
        return self.__str__()


class PropertyClass(PropertyMethods, object):
    """docstring for PropertyClass"""
    def __init__(self, arg):
        super(PropertyClass, self).__init__()
        self.arg = arg


########################################################################################################################


class ModelClass(object):
    ''' must pass a name as a str
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
    freqs = PropertyClass(None)  # the frequency bands. Shared across all ModelClass objects

    def __init__(self, name, *args, **kw):
        super(ModelClass, self).__init__()
        self.name = StringClass(name)

        # Non-Model Properties
        if 'fullname' in kw:  # full name
            self.fullname = StringClass(kw.get('fullname'))
        else:
            self.fullname = StringClass(name)
        if 'info' in kw:  # info is any note. should be str
            self.info = StringClass(kw.get('info'))
        else:
            self.info = StringClass(None)

        # Models
        # self.models = dict()
        for arg in args:
            self.add_model(arg)

    def get(self, arg):
        ''' General get argument
        '''
        return self.__dict__[arg]

    def show(self, arg):
        ''' General print argument
        '''
        print arg, ':  ', self.__dict__[arg]

    def add_model(self, *Dicts):
        ''' adds a data model
            needs: modelname in main dict
            can be given any of (in dict): data, params, & error
            data must be a list
            ex: add_model({'matProp': 'thermCond', 'name':'NIST', 'eq': lambda x: x, 'eqinput': [1, 23]})
        '''

        for Dict in Dicts:  # gets each model in the list of models
            if not isinstance(Dict, dict):  # checking it's a dict, else raise an error.
                print self.name, type()
                raise TypeError('model is not a dict')
            else:
                if 'model' in Dict:  # checking there is a modelname
                    model_name = Dict['model']
                else:
                    raise TypeError('Error: need to pass a model name')

                if isinstance(model_name, str):  # checking the model name is a string
                    setattr(self, model_name, Model(model_name, self.name, self.freqs, Dict))  # *** temporary solution? ***
                else:
                    print self.name, type(model_name)
                    raise TypeError('model is not the right type')


class Model(object):
    ''''''
    def __init__(self, name, upname, freqs, Dict):
        super(Model, self).__init__()
        self.freqs = freqs  # *** Doesn't work ***
        self.name = StringClass(name)  # *** Temporary, need to link to main name ***
        self.upname = upname
        # data as equation and inputs and error
        # Equation
        try:
            equation = Dict['equation']
        except:
            equation = None
        finally:
            self.equation = Equation(equation)

        # eqinput
        try:
            eqinput = Dict['eqinput']
        except:
            eqinput = None
        finally:
            self.eqinput = Eqinput(eqinput)  # *** prob need some ifinstance() check ***

        # params
        try:
            parameters = Dict['params']
        except:
            parameters = None
        finally:
            self.params = Parameters(parameters)

        # Error
        try:
            error = Dict['error']
        except:
            error = None
        self.error = Error(error)

        # Data
        try:
            data = Dict['data']
        except:
            data = None
        self.data = Data(data, self.freqs, self.equation, self.params, self.eqinput)

        # Error in Data
        try:
            d_data = Dict['d_data']
        except:
            d_data = None
        self.d_data = D_Data(d_data, self.freqs, self.data, self.error)

        # Fit Data
        try:
            fitdata = Dict['fitdata']
        except:
            fitdata = None
        self.fitdata = FitData(fitdata, self.freqs, self.equation, self.eqinput, self.params, self.data)

        # Error in Fit Data
        try:
            d_fitdata = Dict['d_fitdata']
        except:
            d_fitdata = None
        self.d_fitdata = D_Fitdata(d_fitdata, self.freqs, self.fitdata, self.error, self.d_data)



class Equation(FunctionMethods, object):
    """docstring for Equation"""
    def __init__(self, arg):
        super(Equation, self).__init__()
        self.set(arg)

    def set(self, arg=None):
        ''' add an equation to a model.'''
        self.setandcheck(arg, 'equation', (FunctionType, NoneType, list, np.ndarray))


class Eqinput(PropertyMethods, object):
    """adds equation input to a model equation."""
    def __init__(self, arg):
        super(Eqinput, self).__init__()
        self.set(arg)

    def set(self, arg=None):
        self.setandcheck(arg, 'eqinput', (NoneType, list, np.ndarray))


class Parameters(PropertyMethods, object):
    """adds equation input parameters to a model equation."""
    def __init__(self, arg):
        super(Parameters, self).__init__()
        self.set(arg)

    def set(self, arg=None):
        self.setandcheck(arg, 'params', (NoneType, list, np.ndarray, dict))


class Error(FunctionMethods, object):
    """docstring for Error"""
    def __init__(self, arg):
        super(Error, self).__init__()
        self.set(arg)

    def set(self, arg):
        self.setandcheck(arg, 'error', (int, float, list, np.ndarray, FunctionType, NoneType))


class Data(PropertyMethods, object):
    """docstring for Data"""
    def __init__(self, arg, freqs, equation, params, eqinput):
        super(Data, self).__init__()
        self.freqs = freqs
        self.equation = equation
        self.params = params
        self.eqinput = eqinput
        self.add_data(arg)

    def add_data(self, data=None):
        ''' adds the data to a model.
            data must be given as None or [[data]], where data is a list.
        '''
        # print self.eqinput[1]

        # Adding the Data to the Model
        if isinstance(data, (list, np.ndarray, NoneType)):
            # Getting all past fitdata
            if data is None:
                eqinput = self.eqinput  # the input to the model's equation
                if eqinput is not None:  # checking if the model actually has an equation, otherwise it is measured data
                    if None in eqinput:  # There is a list of frequencies to iterate through
                        data = [[data]*len(self.freqs)]  # making iterable blank data
                        none_index = eqinput.index(None)  # finding where to iterate the eqinput
                        for i, freq in enumerate(self.freqs):
                            eqinput[none_index] = freq
                            data[0][i] = self.equation(eqinput, self.params)
                            eqinput[none_index] = None  # IMPORTANT resetting back to None
                    else:  # The data is frequency invariant
                        try:  # eqinput contains the data in the final form. Never true for auto-ModelEquations
                            data = self.equation(eqinput, self.params)
                            data[0][0, 0]
                        except:  # the data needs to be copied into each frequency band
                            data = [[data]*len(self.freqs)]  # data is going to be replaced
                            for i in range(len(self.freqs)):  # repeating data for each frequency
                                data[0][i] = self.equation(eqinput, self.params)
            else:  # If data is given
                try:
                    if len(data[0]) == len(self.freqs):  # data in final form
                        pass
                    else:
                        data = [[data] * len(self.freqs)]
                except TypeError:
                    data = [[data] * len(self.freqs)]
            self.arg = np.array(data)
        else:  # raise error if not right type
            print self.name, type(data)
            raise TypeError('data is not the right type')


class D_Data(PropertyMethods, object):
    """docstring for D_Data"""
    def __init__(self, arg, freqs, data, error):
        super(D_Data, self).__init__()
        self.freqs = freqs
        self.data = data
        self.error = error
        self.add_d_data(arg)

    def add_d_data(self, d_data=None):
        ''' adds error to the data.
            adding d_data as None makes it first try to add as model.error(model.data)
            not giving a d_data makes it first try to add as model.error(model.data) then as model.error
        '''
        # Adding the Error to the Data
        if isinstance(d_data, (int, float, list, np.ndarray, NoneType)):  # if same for all freqs
            if d_data == None:
                d_data = [[d_data]*len(self.freqs)]  # making iterable blank d_data
                for i, freq in enumerate(self.freqs):  # iterating over the frequencies
                    try:  # evaulate for the error in the data
                        d_data[0][i] = np.array(self.error(self.data[0][i]))
                    except:  # already list, so add error in data
                        d_data[0][i] = np.array(self.error)
            else:
                if len(d_data[0]) == len(self.freqs):  # if d_data already right length
                    pass
                else:  # only 1 d_data is given, so it is repeated
                    d_data = [np.array(d_data) for i in range(len(self.freqs))]
            self.arg = np.array(d_data)
        else:  # raise error if not right type
            print self.name, type(d_data)
            raise TypeError('d_data is not the right type')


class FitData(PropertyMethods, object):
    """docstring for D_FitData"""
    def __init__(self, arg, freqs, equation, eqinput, params, data):
        super(FitData, self).__init__()
        self.arg = arg
        self.freqs = freqs
        self.equation = equation
        self.eqinput = eqinput
        self.params = params
        self.data = data
        self.add_fitdata(arg)

    def add_fitdata(self, newfitdata=None):
        ''' adds the data, evaulated with the best fit parameters, to a model.
            called material.add_fitdata(property, modelname, equation)
            note:  modelname must be a string
        '''
        # Adding the Fitdata
        if isinstance(newfitdata, (list, np.ndarray, NoneType)):
            # Getting all past fitdata
            if self.arg == None:  # add_arg has never been called
                self.arg = []  # starting with a blank list to add to
            else:  # newfitdata already has entries
                pass  # so don't need to do anything

            if newfitdata is None:
                eqinput = self.eqinput  # the input to the self's equation
                if eqinput is not None:  # checking if the self actually has an equation, otherwise it is measured data
                    if None in eqinput:  # There is a list of frequencies to iterate through
                        newfitdata = [newfitdata]*len(self.freqs)  # making iterable blank data
                        none_index = eqinput.index(None)  # finding where to iterate the eqinput
                        for i, freq in enumerate(self.freqs):
                            eqinput[none_index] = freq
                            newfitdata[i] = self.equation(eqinput, self.params)
                            eqinput[none_index] = None  # Important resetting back to None
                    else:  # The data is frequency invariant
                        print
                        try:  # eqinput contains the data in the final form. This is never auto-called
                            newfitdata = self.equation(eqinput, self.params)  # data generated. Is it in final form?
                            newfitdata[0][0, 0]  # only works if the data is all np arrayed up  # *** IS THIS TRUE FOR newfitdata? ***
                        except:  # the data needs to be copied into each frequency band
                            newfitdata = [newfitdata]*len(self.freqs)  # newfitdata is going to be replaced
                            for i in range(len(self.freqs)):  # repeating data for each frequency
                                newfitdata[i] = self.equation(eqinput, self.params)
                        else:  # if data was in right final form
                            pass  # *** MAKE SURE IT IS IN RIGHT FORM ***
                else:  # only called if data was given and fitdata called w/ nothing given. So copy data to fitdata
                    newfitdata = self.data[0]
            else:  # newfitdata was given  # *** THIS SECTION NEEDS WORK ***
                if len(newfitdata) == len(self.freqs):  # fitdata in final form
                    for i, val in enumerate(newfitdata):  # checking if subsections are right datatype
                        if isinstance(val, np.ndarray):
                            pass
                        else:
                            newfitdata[i] = np.array(val)
                else:  # fitdata needs to be repeated
                    newfitdata = [np.array(newfitdata) for i in self.freqs]
            self.arg.append(newfitdata)
        else:  # raise error if not right type
            raise TypeError('newfitdata is not the right type')

    def finalize_fitdata(self, model):  # *** CHANGE TO FITDATA CALL FUNCTION doing np.array(fitdata) ***
        # Getting the Model
        if isinstance(model, str):  # check if model string
            model = getattr(self, model)  # getting the model
        else:  # raise error if model not string
            print self.name, type(model)
            raise TypeError('model not a string')
        model.fitdata = np.array(model.fitdata)


class D_Fitdata(PropertyMethods, object):
    """docstring for D_Fitdata"""
    def __init__(self, arg, freqs, fitdata, error, d_data):
        super(D_Fitdata, self).__init__()
        self.arg = arg
        self.freqs = freqs
        self.fitdata = fitdata
        self.error = error
        self.d_data = d_data
        self.add_d_fitdata(arg)

    def add_d_fitdata(self, newd_fitdata=None):
        ''' adds error to fitdata.
        '''
        # Adding the Error to the Data
        if isinstance(newd_fitdata, (int, float, list, np.ndarray, NoneType)):
            # Getting all past fitdata
            arg = getattr(self, 'arg')
            if arg is None:  # add_fitdata has never been called
                self.arg = []  # starting with a blank list to add to
            else:  # newfitdata already has entries
                pass  # so don't need to do anything

            if newd_fitdata is None:
                newd_fitdata = [newd_fitdata]*len(self.freqs)  # making iterable blank newd_fitdata
                for i, freq in enumerate(self.freqs):  # iterating over the frequencies
                    try:  # evaulate for the error in the data
                        newd_fitdata[i] = self.error(self.fitdata[-1][i])
                    except TypeError:  # already list, so add error in data
                        if not isinstance(self.error, NoneType):
                            newd_fitdata[i] = np.array(self.error)
                        else:
                            newd_fitdata[i] = np.array(self.d_data[0][i])
                    except AttributeError:
                        newd_fitdata[i] = np.array(self.d_data[0][i])
            else:
                if len(newd_fitdata) == len(self.freqs):  # if newd_fitdata already right length
                    pass
                else:  # only 1 newd_fitdata is given, so it is repeated
                    newd_fitdata = [np.array(newd_fitdata) for i in range(len(self.freqs))]
            self.arg.append(newd_fitdata)
        else:  # raise error if not right type
            print self.name, type(newd_fitdata)
            raise TypeError('newd_fitdata is not the right type')

    def finalize_d_fitdata(self):  # *** CHANGE TO FITDATA CALL FUNCTION doing np.array(fitdata) ***
        self.arg = np.array(self.arg)


# from ModelEquations import BModeSignal
# ModelClass.freqs = np.array([5.0])
# test = ModelClass('name', {'model': 'raw', 'equation': BModeSignal, 'params': np.array([1.]), 'eqinput': range(20)})
# print test.raw.data
# test.raw.data.show()
# test.raw.data.other.set([2])
# print test.raw.data.other
# print test.name, test.raw.name
# test.raw.name.set('name2')
# print test.name, test.raw.name
# print test.raw.name.append('2')
# print '1',  test.raw.equation(4)
# test.raw.equation.set(lambda x: x**2+3)
# print test.raw.equation(4)
# print test.raw.equation.__call__(3)
# test2 = test.raw.equation.callable()
# print test2(2)
# test4 = test.raw.equation
# print test.raw.equation
# print test.raw.params
# test.raw.params.set([2])
# print test.raw.params
# test.raw.params.show()

# test.raw.upname.set('upname')
# test.raw.upname.show()
# test.name.show()

# if test.raw.equation == (lambda x: x**2+2):
#     print 'y'

# test.freqs = 'y'
# print test.raw.freqs
