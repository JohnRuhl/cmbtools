import numpy as np
import pylab

pylab.ion()

c = 3e8
mu0 = 4*np.pi*1e-7
eps0 = 8.85e-12

rho = 2.21e-8  # gold
#rho = 1.59e-8  # silver
#rho = 78.e-8  # stainless steel

sigma = 1/rho
nu = 100e9
#nu = 5e15
w = 2*np.pi*nu

A1 = w*np.sqrt(eps0*mu0/2)
B1 = np.sqrt(1. + (sigma/(eps0*w))**2)
k = A1*np.sqrt(B1 + 1)
kappa = A1*np.sqrt(B1 - 1)

kc = k + 1j*kappa

betac = (c/w)*kc

Rfac = (1-betac)/(1+betac)
Tfac = 2./(1+betac)

R = np.abs(Rfac*np.conj(Rfac))
T = np.abs(Tfac*np.conj(Tfac))

print('A1 = {0:5.4e}'.format(A1))
print('B1 = {0:5.4e}'.format(B1))
print('k = {0:5.4e}'.format(k))
print('beta = {0:5.4e}'.format(betac))
print('R = {0:5.4f}'.format(R))
print('T = {0:5.4f}'.format(T))



