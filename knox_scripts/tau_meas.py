#  Originally in matlab by J.Ruhl, converted to Python by S.Wen Nov. 2014
#  Evaluate knox-formula estimates for tau-measuring experiments.
#
#  Define best "money" band for each experiment;  
#  expt = [freq_GHz ndet NET_uKcmb days_integration fwhm_arcmin sky_area_sqdeg]

import numpy
from knox_func import knox_func
import matplotlib.pyplot as plt
import pylab
pylab.ion()

# Planck 150GHz
planck150={'title': 'Planck, 150GHz, 1 year',
        'freq' : 150,
        'ndet' : 8,           # this is the number of polarization sensitive bolos;  4 more unpol.
        'net': 62,  # per detector at 150GHz
        'days': .24*365*(14/12),  # 14 months, to match bluebook,
        'fwhm': 7.1,          # from bluebook
        'fsky': 0.24,
        'sqdeg': 0.24*42000,
        'deltaell': 2,
        'topbin': 20,
        'lowestell': 4,
        }
        
# Planck 90GHz is similar:  NET = 50, Ndet=8, 9.5arcmin beams.
planck90={'title': 'Planck, 90GHz, 1 year',
        'freq' : 90,
        'ndet' : 8,           # this is the number of polarization sensitive bolos;  4 more unpol.
        'net': 50,  # per detector at 150GHz
        'days': 365*(14/12),  # 14 months, to match bluebook,
        'fwhm': 9.5,          # from bluebook
        'fsky': 0.8,
        'sqdeg': 0.8*42000,
        'deltaell': 2,
        'topbin': 20,
        'lowestell': 4,
        }


# taumachine, 1000 detectors at 150GHz, Wanaka, 30 day flight
taumachine={
        'title': 'TauMachine, 1 flight, 150GHz only',
        'freq' : 150,   #GHz
        'net': 80,  # uKrtsec_cmb
        'days': 0.5*10,  #ULDB, nights
        'fwhm': 30,  #arcmin   
        'ndet': 1000, # 
        'fsky': 0.4,
        'sqdeg': 0.4*41253.,  #sky area coverage
        'deltaell': 2,
        'topbin': 30,
        'lowestell': 3,
        }

# for calculating s/n
ell_max = 20
ell_min = 10 

        
####################################################################
# set which experiment, basefile.
expt = taumachine
#expt = planck150
#expt = planck90
print('Experiment: ' + expt['title'])
print('Sky fraction: ' + str(expt['fsky']))

# camb output files are in D_ell = T_0^2*ell*(ell+1) * C_ell / (2*pi)  (uK^2)
#file1 = '/Users/ruhl/camb_runs/Standard_PlanckTau0.055_r=0.0/camb_25839877_lensedtotcls.dat'
file1 = '/Users/ruhl/camb_runs/Standard_PlanckTau0.055_r=0.0/camb_66227559_lensedtotcls.dat'
file2 = '/Users/ruhl/camb_runs/Standard_Tau0.060_r=0.0/camb_76267983_lensedtotcls.dat'
#file2 = '/Users/ruhl/camb_runs/Standard_Tau0.056_r=0.0/camb_10406913_lensedtotcls.dat'
Dls1 = numpy.loadtxt(file1)
Dls2 = numpy.loadtxt(file2)

# These are ell by ell (ie no binning) for plotting curves
ell1 = Dls1[:,0]
EE1 = Dls1[:,2]
EE2 = Dls2[:,2]
delta_EE = EE2 - EE1

# Do the knox calculation on the second model.  The returned power spectrum values are binned, and they 
#   are C_ell not D_ell.
data2 = knox_func(expt,Dls2)
Cls2 = data2['Cls']
l = Cls2[:,0]
f = l*(l+1)/(2*numpy.pi)  # conversion back to D_ell...
EE2_binned=f*Cls2[:,2]
#cov_Cls2 = data2['cov_Cls']

# Do the knox calculation on first model.  The returned power spectrum values are binned, and they 
#   are C_ell not D_ell.
data1 = knox_func(expt,Dls1)
Cls1 = data1['Cls']
EE1_binned=f*Cls1[:,2]

labels = ['TT','EE','BB','TE']
cov_Cls1 = data1['cov_Cls']
#cov_TT = f*cov_Cls1[:,1]
cov_EE = f*f*cov_Cls1[:,2]
cov_BB = f*f*cov_Cls1[:,3]
#cov_TE = f*cov_Cls1[:,4]


# binned curve
dEE_binned = EE2_binned - EE1_binned
dEE_binned_error = numpy.sqrt(cov_EE)

# plots
pylab.figure(1)
ax = plt.subplot(211)
plt.errorbar(l, EE1_binned, xerr=0, yerr=numpy.sqrt(cov_EE), linewidth=1.5)
plt.errorbar(l, dEE_binned, xerr=0, yerr=numpy.sqrt(cov_EE), linewidth=1.5)

pylab.plot(l,numpy.sqrt(cov_BB),'.',label='BB error')
pylab.plot(l,numpy.sqrt(cov_EE),'.',label='EE error')
#pylab.plot(ell1, EE1, ell1, EE2 )   # difference in EE for delta_tau = 0.005
#pylab.plot(ell1, numpy.abs(EE2- EE1),label='delta EE' )   # difference in EE for delta_tau = 0.005
ax.set_yscale("log", nonposy='clip')
ax.set_ylim(0.0001,0.3)
pylab.legend()
ax.set_title(expt['title'], fontsize=16)
ax.set_xlabel('$l$', fontsize=16)
ax.set_ylabel(r'$D_l$', fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.xaxis.grid(True,'major')
ax.yaxis.grid(True,'major')
ax.yaxis.grid(True,'minor')
pylab.xlim([0,50])

# smooth curves
ell_low = ell1[ell_min:ell_max-1]
EE1_low = EE1[ell_min:ell_max-1]
EE2_low = EE2[ell_min:ell_max-1]
dEE = numpy.abs(EE2_low - EE1_low)

# calculate the S/N on the *binned* power spectra, 

# calculate binstarts the same way as in knox_func.py

s_to_n = dEE_binned/dEE_binned_error
S2N = numpy.sqrt(numpy.sum(s_to_n*s_to_n))
print('Total signal to noise = ' + str(S2N))
sigma_tau = 0.005/S2N
print('Naiive sigma_tau = ' + str(sigma_tau))
pylab.subplot(212)
pylab.plot(l,s_to_n)
pylab.xlim([0,50])
pylab.grid()
pylab.ylabel('Signal to noise, per ell')
string1 = 'S2N total = ' + str(S2N)
pylab.text(30,0.4,string1)

pylab.show()


pylab.show()
