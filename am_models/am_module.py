# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 09:15:30 2015

This module is to be used with the atmospheric modeling program am, available from 
Harvard-Smithsonian Center for Astrophysics. Upon importing the module, you will be
requested to give the name of the am output file. The module will then load the 
data from this file, and the functions inside the module will draw from the data.
Most functions also require that you give the name of a file that contains data 
specific to the instrumentation used for the purpose of noise calculation. For 
example, the file 'SPT90.txt' contains parameters of the 90 GHz detector of the 
South Pole Telescope. Using it with a function like NET_loop or noise_calc returns
calculations specific to that detector system.

@author: Samuel Musilli, szm5@case.edu
"""
import os
os.system
global c, h, k, T_cmb, numpy, T_RJ, I
import numpy, operator, pylab, scipy, math, bolo_module
    #physical constants in SI units
c = float(299792456)
h = 6.626068e-34
k = 1.3806503e-23
T_cmb = 2.725 #Kelvin
file = input('please give an am output file name: ')
f = open(file, 'r')
lines = f.readlines()
f.close()
nu = []
T_RJ = []
I = []
for line in lines: 
    p = line.split()
    nu.append(float(p[0]))
    T_RJ.append(float(p[1]))
    I.append(float(p[2]))
nu = numpy.array(nu)
T_RJ = numpy.array(T_RJ)
I = numpy.array(I)
bump = nu[1] - nu[0]
data_start = nu[0]

def t_mean():
    return numpy.mean(T_RJ)
    
def NET_loop(name, band_shift, command):
    """Plots the NET as a function of band edge shift. 'name' is to be a filename of
    some experiment, arranged how the examples like SPT90.txt and SPT150.txt are arranged.
    'band_shift' is to be the amount (in GHz) you would like to shift the edge of the 
    observing frequency range. The code is written in such a way that a positive value
    will always increase band size - this is regardless of whether or not you're 
    changing the upper or lower edge of the band. 'command' is to be a string, either
    'upper' or 'lower' depending on which band edge you want to shift."""
    band_shift = float(band_shift)
    bump = band_shift/20.
    NET = numpy.empty([20, 1])
    shift = numpy.empty([20, 1])
    T_rj_mean = numpy.empty([20, 1])
    if command == 'lower':    
        for t in range(20):
            NET[t], T_rj_mean[t] = noise_calc(name, t*bump, command)
            shift[t] = t*bump
            if t == 0:
                label = nom
            print "iteration "+str(t+1)+" of 20 complete"
    elif command == 'upper':
        for t in range(20):
            NET[t], T_rj_mean[t] = noise_calc(name, t*bump, command)
            shift[t] = t*bump
            if t == 0:
                label = nom
            print "iteration "+str(t+1)+" of 20 complete"
    pylab.subplot(2, 1, 1)
    if command == 'lower':
        if band_shift > 0:
            pylab.plot((nom-0.5*bandwidth)-shift, numpy.fliplr(NET), 'b')
        else:
            pylab.plot((nom-0.5*bandwidth)-shift, NET, 'b')
    elif command == 'upper':
        pylab.plot(shift+(nom+0.5*bandwidth), NET, 'b')
    pylab.ylabel('NET')
    pylab.xlabel('observing frequency in GHz')
    pylab.title(command+' edge shift for '+str(label)+'GHz band')
    pylab.subplot(2, 1, 2)
    if command == 'lower':
        if band_shift > 0:
            pylab.plot((nom-0.5*bandwidth)-shift, T_rj_mean, 'b')
        else:
            pylab.plot((nom-0.5*bandwidth)-shift, T_rj_mean, 'b')
    elif command == 'upper':
        pylab.plot(shift+(nom+0.5*bandwidth), T_rj_mean, 'b')
    pylab.ylabel('T_rj_mean')
    pylab.xlabel('observing frequency in GHz')
    
def plot_rj(lower, upper):
    """Plots the Rayleigh-Jeans temperature as a function of frequency."""
    pylab.plot(nu[int((lower-data_start)/bump):int((upper-data_start)/bump)]\
    ,T_RJ[int((lower-data_start)/bump):int((upper-data_start)/bump)])
    pylab.xlabel('frequency in GHz')
    pylab.ylabel('T_RJ')
    pylab.title('Rayleigh-Jeans temp as a function of frequency')
    
def waterplot(file1, file2, atm_temp): 
    """receives two am output files, as well as the atmospheric temperature. It then 
    plots the difference between the atmospheric power loading of the two output files
    and plots that difference as a function of frequency. Also plots the difference
    in Rayleigh-Jeans temperature as a function of frequency."""
    f = open(file1, 'r')
    lines = f.readlines()
    f.close()
    nu1 = []
    T_RJ1 = []
    I1 = []
    for line in lines: 
        p = line.split()
        nu1.append(float(p[0]))
        T_RJ1.append(float(p[1]))
        I1.append(float(p[2]))
    nu1 = numpy.array(nu1)
    T_RJ1 = numpy.array(T_RJ1)
    I1 = numpy.array(I1)
    fi = open(file2, 'r')
    lines = fi.readlines()
    fi.close()
    #nu2 = []
    T_RJ2 = []
    I2 = []
    bump1 = nu1[1] - nu1[0]
    for line in lines: 
        p = line.split()
        #nu2.append(float(p[0]))
        T_RJ2.append(float(p[1]))
        I2.append(float(p[2]))
    graph_nu = []
    T_RJ2 = numpy.array(T_RJ2)
    I2 = numpy.array(I2)
    j = len(nu1)
    P_opt1 = []
    P_opt2 = []
    for i in range(int(j*bump)): #P_opt = eta*k*bandwidth*T_rj_mean. This creates unit bandwidth
        P_opt1.append(numpy.longdouble(k*T_RJ1[int(i/bump)]**2/atm_temp))
        P_opt2.append(numpy.longdouble(k*T_RJ2[int(i/bump)]**2/atm_temp))
        graph_nu.append(nu1[int(i/bump)])
    P_opt1 = numpy.array(P_opt1)
    P_opt2 = numpy.array(P_opt2)
    pylab.subplot(2,1,1)    
    pylab.plot(graph_nu, P_opt1-P_opt2)
    pylab.ylabel('power difference in watts')    
    pylab.xlabel('frequency, in GHz')
    pylab.title(file1+' - '+ file2 + ' atmospheric loading vs. frequency')
    pylab.subplot(2,1,2)
    pylab.plot(nu1, T_RJ1-T_RJ2)
    pylab.ylabel('T_RJ difference')    
    pylab.xlabel('frequency, in GHz')
    pylab.title('T_RJ difference by frequency')
    
def noise_calc(name, bandshift, command):
"""this function takes the name of a data file which contains the parameters of an experiment, 
the amount you want to edge of a band to be shifted (which band that is must be described in the
textfile), and whether that shift is for the upper or the lower band (the argument 'command' 
should be, as a string, either 'upper' or 'lower'). a positive bandshift number will always 
give a larger band. That is, if you give a bandshift of 3 GHz and the lower edge is at 80 GHz, 
it will then go to 77 GHz, NOT 83 GHz."""
    global nu, nom, bandwidth
    bandshift = float(bandshift)    
    c = 2.99792456e8
    h = 6.626e-34
    k = 1.38e-23
    T_cmb = 2.725
    import imp, numpy
    f = open(name)
    exp = imp.load_source('experiment', 'r', f)
    f.close()
    data = exp.data
    datasrc = exp.datasrc
    nom = data['nom_band']
    bandwidth = data['band_width']
    low = nom - 0.5*bandwidth
    lower = round(int(low/bump)*bump,1)
    up = nom + 0.5*bandwidth
    upper = round(int(up/bump)*bump,1)
    lower_n = nom - (0.5*bandwidth + bandshift)
    lower_new = round(int(lower_n/bump)*bump,1)
    upper_n = nom + (0.5*bandwidth + bandshift)
    upper_new = round(int(upper_n/bump)*bump,1)
    if command == 'upper':
        upper_i = numpy.where(nu == upper_new)[0][0]
        lower_i = numpy.where(nu == lower)[0][0]
    if command == 'lower':
        upper_i = numpy.where(nu == upper)[0][0]
        lower_i = numpy.where(nu == lower_new)[0][0]
    nuGHz = nu[lower_i:upper_i]*1.e9
    R0 = data['R_bolo']
    RL = data['R_load']
    tau_0 = data['tau_0']
    tau_el = data['tau_electronics']
    L = data['L']
    Tbolo = data['T_bolo']
    Tbase = data['T_base']
    beta = data['beta']
    alpha = data['alpha']
    eta_to_sky = data['eta']
    Qtot = 0
    NEP2_photon_total = 0
    AOmega = data['Nmodes']*c**2/nuGHz**2
    x = h*nuGHz/(k*T_cmb)
    prefactor = data['tau']*eta_to_sky*data['Npol']*data['Nmodes']
    intr = AOmega*(prefactor*h**2*(numpy.power(nuGHz,4)))*numpy.exp(x)\
/(k*(c**2)*(T_cmb**2)*(numpy.exp(x)-1.0)**2)
    dPdTcmb = numpy.trapz(intr,nuGHz)
    for i in range(len(nuGHz)-1):    
        loopnu = numpy.array([nuGHz[i], nuGHz[i+1]])        
        eta = T_RJ[lower_i+i]/datasrc[14]['T']
        P_opt1 = T_RJ[lower_i+i]*eta*k*(bandwidth+bandshift)
        P_opt2 = T_RJ[lower_i+i+1]*(T_RJ[lower_i+i+1]/datasrc[14]['T'])*k*(bandwidth+bandshift)
        y = h*loopnu/(k*datasrc[14]['T'])        
        n = eta*eta_to_sky*numpy.array([1, 1])*data['tau']/(numpy.exp(y)-1)
        NEP2_photon = numpy.trapz(2*h**2*loopnu**2*2*data['Nmodes']*data['Npol']*n, loopnu)
        Q = numpy.trapz(numpy.array([P_opt1, P_opt2]), loopnu)
        NEP2_photon_total += NEP2_photon        
        Qtot += Q
    eta_tot = 1.
    band = numpy.ones(len(nuGHz),float)
    for j in range(len(datasrc)):
        if datasrc[j]['name'] != 'atm':
            x = h*nuGHz/(k*datasrc[j]['T'])
            n = datasrc[j]['eps']*eta_tot*band*data['tau']\
            /(numpy.exp(x)-1)
            integrand = 2*h**2*nuGHz**2*2*data['Nmodes']*data['Npol']*n
            NEP2_photon = numpy.trapz(integrand, nuGHz)
            NEP2_photon_total += NEP2_photon
            eta_tot = eta_tot*(1-datasrc[j]['eps'])
    n = data['n']
    w = 2*math.pi*numpy.arange(1.0, 10.0, 0.01)
    W = 2*Qtot    
    Lg = W*alpha/(2*numpy.mean(T_RJ[lower_i:upper_i])**2/datasrc[15]['T']*k*(bandwidth+bandshift)*Tbolo)
    P_el = W
    I_0 = numpy.sqrt(P_el/R0)
    tau_i = numpy.longdouble(tau_0/(1-Lg))
    B = 0.5*numpy.sqrt((1/tau_el - 1/tau_i)**2 - 4*R0/L*(Lg*(2+beta)/tau_0))
    A = 1/(2*tau_el)+1/(2*tau_i)
    tau_p = numpy.longdouble(1./(A+B))
    tau_m = numpy.longdouble(1./(A-B))
    T1 = -1./(I_0*R0)
    T2 = 1./(2.+beta) 
    T3 = (1.-tau_p/tau_i)/(1.+1j*w*tau_p)
    T4 = (1.-tau_m/tau_i)/(1.+1j*w*tau_m)
    s_w = T1*T2*T3*T4
    Si_tes = 4*k*Tbolo*R0*I_0**2*(1.+w**2*tau_0**2)*abs(s_w)**2/Lg**2
    T_L = 4.2    
    Si_L = 4*k*T_L*I_0**2*RL*(Lg-1)**2*(1.+w**2*tau_i**2)*abs(s_w)**2/Lg**2
    b = n-1
    T = numpy.arange(Tbase, Tbolo, (Tbolo-Tbase)/1000.)
    F_Tbolo_Tbase = sum((T*(T**b)/(Tbolo*(Tbolo**b)))**2)/sum((T**b)/(Tbolo**b))
    NEP_phonon = numpy.sqrt(4*k*Tbolo**2*numpy.mean(T_RJ[lower_i:upper_i])**2/\
datasrc[15]['T']*k*(bandwidth+bandshift)*F_Tbolo_Tbase)
    JohnsonNEI = numpy.longdouble(numpy.sqrt(Si_tes+Si_L))
    NEP_Johnson = JohnsonNEI/abs(s_w)
    NEP_total = numpy.sqrt(NEP2_photon_total+NEP_phonon**2+numpy.mean(NEP_Johnson)**2)
    NET = NEP_total/dPdTcmb/numpy.sqrt(2.)
    return NET, numpy.mean(T_RJ[lower_i:upper_i])
    
    
    
    
