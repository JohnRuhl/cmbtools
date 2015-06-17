# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:44:59 2015

@author: sam
"""
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
        
def plot_knox():
    import matplotlib.pyplot as plt
    import numpy
    expt = sptpol
    basefile = 'camb_24109479_lensedtotcls.dat';  # lensed with r= 0.01
    baseCls = numpy.loadtxt(basefile);  # T_0^2*ell*(ell+1) * C_ell / (2*pi)  (uK^2)
    ell=(baseCls[:,0]).astype(int)
    factor=ell*((ell + 1)) / (2 *numpy.pi)
    C_TT=baseCls[:,1] / factor
    C_EE=baseCls[:,2] / factor
    C_BB=baseCls[:,3] / factor
    C_TE=baseCls[:,4] / factor
    ax1 = plt.subplot(221)
    ax1.plot(ell, C_TT)
    ax1.set_title("C_TT")
    ax2 = plt.subplot(222)
    ax2.plot(ell, C_EE)    
    ax2.set_title('C_EE')
    ax3 = plt.subplot(223)
    ax3.plot(ell, C_BB)
    ax3.set_title('C_BB')
    ax4 = plt.subplot(224)
    ax4.plot(ell, C_TE)
    ax4.set_title('C_TE')
    plt.show()