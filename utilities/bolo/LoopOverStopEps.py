import matplotlib.pyplot as plt
import numpy as np

import optical_calcs
import DetectorDefs
import tes_calcs

plt.ion()


detector, optics = DetectorDefs.generic150()

half_angle_deg = 17.5 
half_angle_rad = half_angle_deg*np.pi/180.

# calculate relative NET_cmb's
eps_vec = np.arange(0.01,0.7,0.01)
netcmb_tot_vec = np.array([])
netcmb_j_vec = np.array([])
netcmb_phonon_vec = np.array([])
netcmb_photon_vec = np.array([])
Density = np.array([])
MappingSpeed = np.array([])
for eps in eps_vec:
    for item in optics:
        if item['name'] == 'Cold stop':
            item['eps'] = eps
    detector, optics = optical_calcs.optical_calcs(detector,optics)
    detector = tes_calcs.tes_calcs(detector)
    netcmb_tot_vec = np.append(netcmb_tot_vec,detector['NoiseTotal']['NET_CMB'])
    netcmb_j_vec = np.append(netcmb_j_vec,detector['Johnson']['NET_CMB'])
    netcmb_phonon_vec = np.append(netcmb_phonon_vec,detector['Phonon']['NET_CMB'])
    netcmb_photon_vec = np.append(netcmb_photon_vec,detector['NET_photon_total_cmb'])

    beam_power_at_edge = eps_vec
    beam_sigma = half_angle_rad*np.sqrt(-1./(2*np.log(beam_power_at_edge)))
    theta_waist_rad = beam_sigma*2
    w0_over_lambda = 1/(np.pi*np.tan(theta_waist_rad))
    horn_diam_over_lambda = 3.125*w0_over_lambda
    Density = np.append(Density, 1/horn_diam_over_lambda**2)
    MappingSpeed = np.append(MappingSpeed,Density/netcmb_tot_vec**2)

    w0_overlambda_lowerband


# Calculate relative detector density for each element of eps_vec 
# How:  given the angle to the stop, find the diameter D in n*f*lambda of a corrugated feed.
# Density is proportional to 1/D^2.
# Detector count is proportional to density.
# Mapping speed is proportional to detector_count/NET^2


titlestring = 'T_cs = 0.1K, T_stop = 1K'
plt.figure(1)
plt.clf()
plt.subplot(311)
plt.plot(eps_vec,1e6*netcmb_tot_vec,label='total')
plt.plot(eps_vec,1e6*netcmb_j_vec,label='johnson')
plt.plot(eps_vec,1e6*netcmb_phonon_vec,label='phonon')
plt.plot(eps_vec,1e6*netcmb_photon_vec,label='photon')
plt.legend()
plt.grid()
plt.title(titlestring)
plt.ylabel('NET_cmb (uK/rtHz)')
#
plt.subplot(312)
plt.plot(eps_vec,Density/Density[0],label='Det. Density')
plt.ylim(ymin=0)
plt.ylabel('Relative number per unit area')
plt.legend()
plt.grid()
#
plt.subplot(313)
plt.plot(eps_vec, MappingSpeed/MappingSpeed[0],label='Mapping Speed')
plt.ylim(ymin=0)
plt.ylabel('Relative Mapping Speed')
plt.xlabel('Spillover fraction at cold stop')
plt.legend()
plt.grid()

plt.show()



