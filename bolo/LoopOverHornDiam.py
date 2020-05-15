import matplotlib.pyplot as plt
import numpy as np

import optical_calcs
import DetectorDefs
import tes_calcs

# CDT
# NET_cmb = [270, 238, 309, 331] uKrtsec
# at [85, 95, 145, 155] GHz, using
# using roughly 4 tel at 85/145, 4 tel at 95/155 ?  Maybe 5 at each.
# Ndet_tot = [17000, 18000, 21000, 21000]
# So say Ndet = 19500 in 4 or 5 150GHz telescopes
# So Ndet = 3900 or 4874 in 5 or 4 telescopes.
CDT_npix_4 = 4874/2.
CDT_npix_5 = 3900/2.
# So CDT mapping speed is about 4874/320uKrtsec^2 = 0.0475, or 3900/320uKrtsec^2 = 0.038 s/uK^2.

# Bicep numbers: B2=bicep2, BA == bicep array, B3 = Bicep3
# f/#: B2= 2.0, BA=1.6, B3 = 1.6
# pixel size:BA=5.2mm, B3 = 7.8mm@90GHz
# Ndet at 150GHz: B3=2560@90GHz, BA = 7776
#  I think bicep array is about 300 pixels per tile, so 600 detectors, so 12*600 = 7200 det/fpu for 150GHz.
#  That's a 55cm aperture.  
#

# Simons numbers
#  5.6mm pixel spacing
#  From Aikito, 18.98deg, f/1.4 or so.
#  -2.71dB at 90GHz, -6.29dB at 150GHz.  P(theta)= 0.54, 0.23 at stop edge)
#    similar but a little different for sinuous lenslets.
# 

DO150 = True
#DO150 = False

c = 3e8 #m/s

plt.ion()

half_angle_deg = 18.98
half_angle_rad = half_angle_deg*np.pi/180.
fnumber = 1/(2*np.tan(half_angle_rad))
primary_diameter = 0.46   # meters
FParea = 1200*1e-4  # based on SO, m^2

#-------------
lambda_horn = 0.002
inter_horn_spacing = 0.0005  # 0.5mm, very aggressive!
d_horn_vec = np.arange(0.2*lambda_horn, 5.5*lambda_horn, 0.1*lambda_horn)

# Run first freq band
if DO150:
    detector, optics = DetectorDefs.S4_150()
    legendstring = '150 GHz'
else:
    detector, optics = DetectorDefs.S4_90()
    legendstring = '90 GHz'

nu_center = (detector['nu'][0] + detector['nu'][-1])/2.
lambda_center = c/nu_center

# Figure out how many f*lambda these horns are
flambda_horn = d_horn_vec/(fnumber*lambda_center)

# calculate relative NET_cmb's
#eps_vec = np.arange(0.01,0.7,0.01)
eps_vec = np.array([])
dBedge_vec = np.array([])
theta_fwhm_vec = np.array([])
netcmb_tot_vec = np.array([])
netcmb_j_vec = np.array([])
netcmb_phonon_vec = np.array([])
netcmb_photon_vec = np.array([])
Density = np.array([])
Npixels = np.array([])
MappingSpeed = np.array([])

for d_horn in d_horn_vec:
    # calculate spillover
    w0 = d_horn/3.125     # Phil Mauskopf thinks this is ugly because one is a diam, the other a radius.
    theta0 = np.arctan(lambda_center/(np.pi*w0))
    theta_fwhm_feed = 1.18*theta0
    sigma = theta_fwhm_feed/2.355
    #  or
    #sigma = np.arctan(lambda_center/(np.pi*w0)) / 2.
    Pedge = np.exp(-half_angle_rad**2/(2*sigma**2))
    dBedge = np.abs(10*np.log10(Pedge))
    theta_fwhm = (1.02 + 0.0135*dBedge)*lambda_center/primary_diameter
    eps = Pedge # true for a 2D gaussian
    for item in optics:
        if item['name'] == 'Cold stop':
            item['eps'] = eps
    theta_fwhm_vec = np.append(theta_fwhm_vec,theta_fwhm)
    dBedge_vec = np.append(dBedge_vec,dBedge)
    eps_vec = np.append(eps_vec,eps)
    detector, optics = optical_calcs.optical_calcs(detector,optics)
    detector = tes_calcs.tes_calcs(detector)

    density_1 = 1/d_horn**2
    net_cmb_tot_1 = detector['NoiseTotal']['NET_CMB']/np.sqrt(2)  # convert from uKrtHz to uK/rtsec
    netcmb_tot_vec = np.append(netcmb_tot_vec,net_cmb_tot_1)
    # These are in uK/rtHz
    #netcmb_j_vec = np.append(netcmb_j_vec,detector['Johnson']['NET_CMB'])
    #netcmb_phonon_vec = np.append(netcmb_phonon_vec,detector['Phonon']['NET_CMB'])
    #netcmb_photon_vec = np.append(netcmb_photon_vec,detector['NET_photon_total_cmb'])

    # hex close packed, but allowing for lam/2 extra around each horn.
    density_1 = 0.9069*((d_horn/(d_horn+inter_horn_spacing))**2)*1/d_horn**2  # hex close packed
    Density = np.append(Density, density_1)
    npixels_1 = density_1*FParea
    Npixels = np.append(Npixels,npixels_1)
    MappingSpeed = np.append(MappingSpeed,2*npixels_1/net_cmb_tot_1**2)


# Calculate relative detector density for each element of eps_vec 
# How:  given the angle to the stop, find the diameter D in n*f*lambda of a corrugated feed.
# Density is proportional to 1/D^2.
# Detector count is proportional to density.
# Mapping speed is proportional to detector_count/NET^2

titlestring = 'T_cs = 0.1K, T_stop = 1K, D_primary = 46cm, theta_half_stop = 18.98deg'
plt.subplot(321)
plt.plot(1000*d_horn_vec,1e6*netcmb_tot_vec,label=legendstring)
if DO150:
    CDTNET = 0*d_horn_vec + 320.0  #uKrtsec
    plt.plot(1000*d_horn_vec,CDTNET, label='CDT150')
else:
    CDTNET = 0*d_horn_vec + 255  #uKrtsec
    plt.plot(1000*d_horn_vec,CDTNET,'--', label='CDT90')
#plt.plot(eps_vec,1e6*netcmb_j_vec,label='johnson')
#plt.plot(eps_vec,1e6*netcmb_phonon_vec,label='phonon')
#plt.plot(eps_vec,1e6*netcmb_photon_vec,label='photon')
plt.legend(loc='lower right')
plt.ylim(ymin=0)
plt.grid(b=True)
plt.title(titlestring)
plt.ylabel('NET_cmb (uKrtsec)')
#
plt.subplot(322)
plt.plot(1000*d_horn_vec, theta_fwhm_vec*180./np.pi*60.,label=legendstring)
plt.ylabel('FWHM (arcmin)')
plt.legend(loc='lower right')
plt.grid(b=True)
#
plt.subplot(323)
plt.plot(1000*d_horn_vec,Npixels,label=legendstring)
if DO150:
    plt.plot(1000*d_horn_vec, 0*d_horn_vec+CDT_npix_4, label='CDT_4tel')
    plt.plot(1000*d_horn_vec, 0*d_horn_vec+CDT_npix_5, label='CDT_5tel')
plt.ylim(ymin=0)
plt.ylabel('Pixels per telescope')
plt.legend(loc='lower right')
plt.grid(b=True)
#
plt.subplot(324)
plt.plot(1000*d_horn_vec, eps_vec,label=legendstring)
plt.ylabel('Spillover fraction at cold stop')
plt.legend(loc='lower right')
plt.grid(b=True)
#
plt.subplot(325)
if DO150:
    plt.plot(1000*d_horn_vec, MappingSpeed/1e12,label=legendstring)
    CDTMS_5 = 0*d_horn_vec + 0.038   # Ndet/uKrtsec^2 = [sec/uK^2]
    CDTMS_4 = 0*d_horn_vec + 0.0475   # Ndet/uKrtsec^2 = [sec/uK^2]
    plt.plot(1000*d_horn_vec, CDTMS_4, label = 'CDT-150_4tel')
    plt.plot(1000*d_horn_vec, CDTMS_5, label = 'CDT-150_5tel')
    plt.xlabel('Horn Diameter (mm)')
else:
    plt.plot(1000*d_horn_vec, MappingSpeed/1e12,'--',label=legendstring)
    CDTMS_4 = 0*d_horn_vec + 0.075   # Ndet/uKrtsec^2 = [sec/uK^2]
    CDTMS_5 = 0*d_horn_vec + 0.060   # Ndet/uKrtsec^2 = [sec/uK^2]
    plt.plot(1000*d_horn_vec, CDTMS_4,'--',label = 'CDT-90_4tel')
    plt.plot(1000*d_horn_vec, CDTMS_5,'--', label = 'CDT-90_5tel')
plt.ylim(ymin=0)
plt.ylabel('Mapping Speed (s/uK^2)')
plt.legend(loc='lower right')
plt.grid(b=True)
#
plt.subplot(326)
plt.plot(1000*d_horn_vec, dBedge_vec,label=legendstring)
plt.xlabel('Horn Diameter (mm)')
plt.ylabel('Edge taper (dB)')
plt.ylim(ymin=0)
plt.legend()
plt.grid(b=True)

plt.figure(2)
plt.plot(d_horn_vec,flambda_horn)

plt.show()

    





