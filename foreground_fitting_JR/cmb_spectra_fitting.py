# -*- coding: utf-8 -*-
"""
Created 14 Nov 2014
J. Ruhl

Functions to read in and fit cmb spectra to cmb and foreground models

"""
def read_binned_Dell(filename)
   # Data file should be in format
   #   photon_freq(GHz)   ell   D_ell    sigma_Dell

   return data

def galactic_dust(nu_GHz,ell)
   # given a scalar value of photon frequency = nu (in GHz)
   # and a vector list of values of ell,
   # return a vector list of amplitudes of D_ell for galactic dust from Planck model,
   # "pretty good" area on sky.

   return D_dust

def cmb_Dell_template(ell)
    
    return 

def model_cmbdust(index,A_cmb,A_dust)
    return A_cmb*cmb_Dell(ell) + A_dust*dust_Dell(nu,ell)

def fit_data_to_models
import numpy as np
from scipy.optimize import curve_fit

>>> def func(x, a, b, c):
    ...     return a * np.exp(-b * x) + c
    >>>
    >>> xdata = np.linspace(0, 4, 50)
    >>> y = func(xdata, 2.5, 1.3, 0.5)
    >>> ydata = y + 0.2 * np.random.normal(size=len(xdata))
    >>>
    >>> popt, pcov = curve_fit(func, xdata, ydata)


# ---------------------------------------
# Following Garcia's numerical methods book
# generate A matrix, start with all zeros

A = numpy.zeros((len(ell),n_params));

# Generate first column first
j = 1
# "naiive" way
for i in 1:N_points
    A(i,j) =  1/sigma(i);
    end
    % better way
    % A(:,j) = 1./sigma;

    j=2;
    % "naiive" way
    for i = 1:N_points
        A(i,j) = T(i)./sigma(i);
        end
        % better way
        % A(:,j) = T./sigma;


        % ---------------------------------------
        % generate b vector...
        b = R_meas./sigma;
        b = b.';  % make it a column vector

        At = A.';

        C = inv(At*A)

        a = C*(At*b)


        %  Error on a resistance at T = 270K
        var_R_naiive = 1*C(1,1) + 270.^2*C(2,2);
        var_R = var_R_naiive + 2*1*270*C(1,2);

        sigma_R_naiive = sqrt(var_R_naiive)
        sigma_R = sqrt(var_R)

