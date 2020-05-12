# tau_func.py net detectors days fsky ellmin ellmax
#  This version takes four arguments to populate the experiment description.
#  expt = [freq_GHz ndet NET_uKcmb days_integration fwhm_arcmin sky_area_sqdeg]


import numpy as np
from knox_func import knox_func
import pylab

pylab.ion()

# camb output files are in D_ell = T_0^2*ell*(ell+1) * C_ell / (2*pi)  (uK^2)
#file1 = '/Users/ruhl/camb_runs/Standard_PlanckTau0.055_r=0.0/camb_25839877_lensedtotcls.dat'
file1 = '/Users/ruhl/camb_runs/Standard_PlanckTau0.055_r=0.0/camb_66227559_lensedtotcls.dat'
tau1 = 0.055
file2 = '/Users/ruhl/camb_runs/Standard_Tau0.060_r=0.0/camb_76267983_lensedtotcls.dat'
tau2 = 0.060
#file2 = '/Users/ruhl/camb_runs/Standard_Tau0.056_r=0.0/camb_10406913_lensedtotcls.dat'
Dls1 = np.loadtxt(file1)
Dls2 = np.loadtxt(file2)

# 5K, 10K, 20K, 40K, 100K, 1000K
# P depth: 67, 47, 33, 24, 15, 4.7
# Our baseline experiment
net = 80.
ndet = 1000
ndays = 1000.
fsky = 0.4
ellmin = 4
ellmax = 40

# set up which variable we're going to vary
# We're going to call it xvar, and set the range appropriately as we call the for loop, 
# and put it in the right places in the definition of expt .
#
sigma_tau_vec = np.array([])
xvariable_vec = np.array([])
xlabelstr = 'ell_min'
titlestr = '(net,fsky,lmax) = (80,0.4,40)'
ii = 0
for xvar in np.arange(4,30):
    expt={'title': 'expt varying ellmin',
        'freq' : 150,
        'ndet' : ndet,     # this is the number of polarization sensitive bolos;  
        'net': net,  # per detector at 150GHz
        'days': ndays,  # keep map depth constant
        'fwhm': 30.,    # arcmin
        'fsky': fsky,
        'sqdeg': fsky*42000,
        'deltaell': 2,
        'topbin': ellmax,
        'lowestell': xvar,
        }

    # Do the knox calculation on the second model.  
    # The returned power spectrum values are binned, and they 
    #   are C_ell not D_ell.
    data2 = knox_func(expt,Dls2)
    Cls2 = data2['Cls']
    l = Cls2[:,0]
    f = l*(l+1)/(2*np.pi)  # conversion back to D_ell...
    EE2_binned=f*Cls2[:,2]
    #cov_Cls2 = data2['cov_Cls']

    # Do the knox calculation on first model.  The returned power spectrum values are binned, and they 
    #   are C_ell not D_ell.
    data1 = knox_func(expt,Dls1)
    Cls1 = data1['Cls']
    EE1_binned=f*Cls1[:,2]

    cov_Cls1 = data1['cov_Cls']
    #cov_TT = f*cov_Cls1[:,1]
    cov_EE = f*f*cov_Cls1[:,2]
    cov_BB = f*f*cov_Cls1[:,3]
    #cov_TE = f*cov_Cls1[:,4]

    # binned curve
    dEE_binned = EE2_binned - EE1_binned
    dEE_binned_error = np.sqrt(cov_EE)

    s_to_n = dEE_binned/dEE_binned_error
    S2N = np.sqrt(np.sum(s_to_n*s_to_n))
    sigma_tau_vec = np.append(sigma_tau_vec,(tau2-tau1)/S2N)
    xvariable_vec = np.append(xvariable_vec,xvar)
    ii = ii+1

pylab.figure(1)
pylab.plot(xvariable_vec,sigma_tau_vec)
pylab.ylabel('sigma_tau')
pylab.xlabel(xlabelstr)
pylab.title(titlestr)
pylab.grid()
pylab.ylim([0,0.015])

pylab.show()
