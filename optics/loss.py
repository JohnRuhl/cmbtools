import numpy as np
import matplotlib.pyplot as plt

import bbody

plt.ion()

# eps = eps_r + j*eps_i
# tan_delta = (w*eps_i + sigma)/eps_r
# P = P0*exp(-delta*k*z)
# k = w*sqrt(mu*eps_r)

# global physical constants
c = 3.0e8    # m/s
e0 = 8.85*10**-12
u0 = 4*np.pi*10**-7
n0 = 1.000
epsmu = np.sqrt(e0/u0)
Y0 = epsmu*n0

def delta_from_eps(eps):
    delta = np.arctan2(eps_i,eps_r)
    return delta

def eps_from_tandelta(tandelta,n):
    eps_r = n*n
    eps_i = eps_r*tandelta
    eps_tot = eps_r + 1j*eps_i
    return eps_tot

def T_loss(nu,eps,length):
    w = 2*np.pi*nu
    eps_r = eps.real
    eps_i = eps.imag
    n = np.sqrt(eps_r)
    kprime = eps_i/(2*n)
    alpha = 2*kprime*w/c
    trans = np.exp(-alpha*length)
    return trans

nu_GHz = np.arange(100., 400.)
nu = 1e9*nu_GHz

# Stycast 2850FT (300K, 300GHz, lamb)
#n = 2.28
#thickness = 0.001*25.4/1000. # 1mil, in m
#tandelta = 275.e-4   # lamb indicates this is freq dep, but this is at 300GHz.
#eps_tot = eps_from_tandelta(tandelta,n)
#loss = T_loss(nu,eps_tot,thickness)

## HDPE
n = 1.53
thickness = (1/16.)*25.4/1000. # 
#tandelta = 3.e-4   # lamb indicates this is freq dep, but this is at 300GHz.
tandelta = 5.6e-4   # lamb indicates this is freq dep, but this is at 300GHz.
eps_tot = eps_from_tandelta(tandelta,n)
loss = T_loss(nu,eps_tot,thickness)

# Stycast 2850FT 300K (Halpern Gush etal)
#n = 2.2
##a = 7.0e-3
##b = 2.2
#a = 2.5e-2
#b = 2.2
#alpha = a*(nu/30e9)**b
#b = 2.2
#alpha = a*(nu/30e9)**b
#alpha = alpha*100.
#loss = np.exp(-alpha*thickness)


# Cirlex
# Lau etal
# https://arxiv.org/pdf/astro-ph/0701091.pdf
#eps_r = 3.37*e0
## eps_i = eps_prefactor* (nu/150GHz)**eps_i_beta
## These are 300K values
#eps_i_prefactor = 0.3    #Table 1, one sample 0.27, another 0.37.
#eps_i_beta = 0.52        # Table 1
#thickness = 2.*0.121e-3  # 2 layers of that thickness cirlex
#eps_tot = eps_r + eps_i



plt.figure(1)
plt.clf()
plt.plot(nu_GHz, loss)
plt.grid()

plt.show()



