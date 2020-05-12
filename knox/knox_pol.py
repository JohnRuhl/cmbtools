#  Originally in matlab by J.Ruhl, converted to Python by S.Wen Nov. 2014
#  Evaluate knox-formula estimates for various experiments
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
        'days': 365*(14/12),  # 14 months, to match bluebook,
        'fwhm': 7.1,          # from bluebook
        'sqdeg': 0.8*42000,
        'deltaell': 10,
        'lowestell': 2,
        }
        
# Planck 90GHz is similar:  NET = 50, Ndet=8, 9.5arcmin beams.
planck90={'title': 'Planck, 90GHz, 1 year',
        'freq' : 90,
        'ndet' : 8,           # this is the number of polarization sensitive bolos;  4 more unpol.
        'net': 50,  # per detector at 150GHz
        'days': 365*(14/12),  # 14 months, to match bluebook,
        'fwhm': 9.5,          # from bluebook
        'sqdeg': 0.8*42000,
        'deltaell': 10,
        'lowestell': 2,
        }


# spider, 3 telescopes at 150GHz, McMurdo
spider={'title': 'Spider, 3 telescopes, 1 McMurdo flight',
        'freq' : 150,   #GHz
        'net': 130,  # uKrtsec_cmb
        'days': 0.7*15,  #days in an LDB flight
        'fwhm': 30,  #arcmin   
        'ndet': int(0.5*3*(4*2*64)), # 3 telescopes, 50% yield
        'sqdeg': 4000.,  #sky area coverage
        'deltaell': 50,
        'topbin': 425,
        'lowestell': 25,
        }

# taumachine, 1000 detectors at 150GHz, Wanaka, 30 day flight
taumachine={'title': 'TauMachine, 1 flight, 150GHz only',
        'freq' : 150,   #GHz
        'net': 80,  # uKrtsec_cmb
        'days': 0.5*10,  #ULDB, nights
        'fwhm': 30,  #arcmin   
        'ndet': 1000, # 
        'sqdeg': 20000.,  #sky area coverage
        'deltaell': 1,
        'topbin': 100,
        'lowestell': 4,
        }

# keck at 3 telescopes, 1 year
keck1={'title': 'Keck, 1 year, 3 telescopes',
       'freq' : 150,
        'net': 450,  # per detector at 150GHz
        'days': 0.3*3*365,
        'fwhm': 30,
        'ndet': 3*4*2*64, 
        'sqdeg': 1000,
        'deltaell': 25,
        'topbin': 500,
        'lowestell': 50,
        }
# keck at 3 years
keck3 = keck1;
keck3={'title': 'Keck, 3 years, 3 telescopes'}

# sptsz
sptsz={'title': 'SPT-SZ, 3 years, 2000sqdeg',
       'freq' : 150,
        'net': 350,  # per detector at 150GHz
        'days': 0.3*3*365,
        'fwhm': 1,
        'ndet': 650, 
        'sqdeg': 2000,
        'deltaell': 200,
        'topbin': 4000,
        'dl' : 50,
        'lowestell': 600,
        }
        
# sptpol
sptpol = sptsz;
sptpol={'title': 'SPTpol 150GHz, 30 sq deg, 3 years',
        'net': 448/2,  # per detector at 150GHz
     #   'net': 394,  # per detector at 90GHz
        'days': 0.25*3*365,
        'fwhm': 1,
        'ndet': (7*84)*2, # at 150GHz
        'sqdeg': 30.0,
        'deltaell': 100,
        'topbin': 4000,
        'lowestell': 50,
        }

# sptpol2011
sptpol2011 = sptpol;
sptpol2011={'title': 'SPTpol2011',
        'net': 394,  
        'days': 0.25*3*365,
        'fwhm': 1.6,
        'ndet': 2*192, # 90's only
        'sqdeg': 100,
        'deltaell': 100,
        'lowestell': 50,
        }
        
# sptpol2
sptpol2 = sptpol;
sptpol2={'title': 'SPTpol2',
        'net': 300,  # per detector at 150GHz
        'ndet': 2*sptpol['ndet'],
        'sqdeg': 2000,
        'deltaell': 50,
        'lowestell': 50,
        }
        
# Polar1
polar1={'title': 'Polar 1, 1000sqdeg',
        'net': 350,  
        'days': 0.3*3*365,
        'fwhm': 5,
        'ndet': 2000, 
        'sqdeg': 1000,
        'deltaell': 50,
        'lowestell': 50,
        }
       
# Polar10
polar10 = polar1;
polar10={'title': 'Polar10, 2000sqdeg',
        'days': 10*polar1['days'],
        'sqdeg': 2000,
        }
        
####################################################################
# set which experiment, basefile.
expt = taumachine
#basefile = 'camb_24109479_lensedtotcls.dat'  # lensed with r= 0.01
basefile = '/Users/ruhl/camb_runs/Standard_PlanckTau0.055_r=0.0/camb_25839877_lensedtotcls.dat'  # lensed, r=0.0, planck tau

baseCls = numpy.loadtxt(basefile)  # T_0^2*ell*(ell+1) * C_ell / (2*pi)  (uK^2)

data = knox_func(expt,baseCls)

Cls = data['Cls']
cov_Cls = data['cov_Cls']

# save data to a file;  overwrites.
# These are the actual Cls (not l(l+1)Cl ) and cov on those.  
# (ie not D_l's)
#numpy.save('covCls.dat', cov_Cls)
#numpy.save('Cls.dat', Cls)

l = cov_Cls[:,0]
labels = ['TT','EE','BB','TE']
cov_TT = cov_Cls[:,1]
cov_EE = cov_Cls[:,2]
cov_BB = cov_Cls[:,3]
cov_TE = cov_Cls[:,4]

f = l*(l+1)/(2*numpy.pi)  # conversion back to D_ell...

# plots
pylab.figure(1)
pylab.clf()
ax = plt.subplot(111)
for XX in range(4):
    YY= XX+1
    plt.errorbar(l, numpy.abs(f*Cls[:,YY]), xerr=0, yerr=f*numpy.sqrt(cov_Cls[:,YY]), linewidth=1.5,label=labels[XX])

pylab.plot(l,f*numpy.sqrt(cov_BB))

ax.set_yscale("log", nonposy='clip')
ax.set_ylim(0.0001,100)
pylab.legend()
ax.set_title(expt['title'], fontsize=16)
ax.set_xlabel('$l$', fontsize=16)
ax.set_ylabel(r'$C_l^{XX}$', fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.xaxis.grid(True,'major')
ax.yaxis.grid(True,'major')
ax.yaxis.grid(True,'minor')

pylab.show()
