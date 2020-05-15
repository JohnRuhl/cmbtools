import matplotlib.pyplot as plt
import numpy as np

import optical_calcs
import DetectorDefs
import tes_calcs

plt.ion()

detector, optics = DetectorDefs.S4_SAT_155()
detector, optics = optical_calcs.optical_calcs(detector, optics)

T_cs_vec = np.arange(0.020,0.4,0.001)

plt.figure(1)
plt.clf()

# 280 GHz LAT
detector['PsatFactor'] = 3.0
detector['Qtot'] = 10.1e-12   # optical power
detector['W'] = detector['PsatFactor']*detector['Qtot']

NEP_johnson = np.array([])
NEP_phonon  = np.array([])

for T_coldstage in T_cs_vec:
    detector['T_base'] = T_coldstage
    detector['T_bolo'] = 1.6*T_coldstage

    detector = tes_calcs.tes_calcs(detector)
    NEP_johnson = np.append(NEP_johnson, detector['Johnson']['NEP'])
    NEP_phonon = np.append(NEP_phonon, detector['Phonon']['NEP'])

NEP_photon = 0*NEP_phonon + 91e-18
NEP_tot = np.sqrt(NEP_photon**2 + NEP_phonon**2 + NEP_johnson**2)

ax=plt.subplot(311)
plt.plot(T_cs_vec,NEP_photon*1e18, label='Photon')
plt.plot(T_cs_vec,NEP_phonon*1e18, label='Phonon')
plt.plot(T_cs_vec,NEP_johnson*1e18,label='Johnson')
plt.plot(T_cs_vec,NEP_tot*1e18,label='Total')
plt.text(0.05, 0.75, '280GHz LAT, P_sat = 30.3pW', transform=ax.transAxes)
plt.xlabel('T_base')
plt.ylabel('NEP [aW/rtHz]')
plt.grid()
plt.legend(loc='center right')
#plt.title('280GHz SAT, P_sat = 30.3pW, T_bolo = 1.6*T_base')


# 93 GHz LAT
detector['PsatFactor'] = 3.0
detector['Qtot'] = 1.18e-12   # optical power
detector['W'] = detector['PsatFactor']*detector['Qtot']

NEP_johnson = np.array([])
NEP_phonon  = np.array([])

for T_coldstage in T_cs_vec:
    detector['T_base'] = T_coldstage
    detector['T_bolo'] = 1.6*T_coldstage

    detector = tes_calcs.tes_calcs(detector)
    NEP_johnson = np.append(NEP_johnson, detector['Johnson']['NEP'])
    NEP_phonon = np.append(NEP_phonon, detector['Phonon']['NEP'])

NEP_photon = 0*NEP_phonon + 15.1e-18
NEP_tot = np.sqrt(NEP_photon**2 + NEP_phonon**2 + NEP_johnson**2)

ax = plt.subplot(312)
plt.plot(T_cs_vec,NEP_photon*1e18, label='Photon')
plt.plot(T_cs_vec,NEP_phonon*1e18, label='Phonon')
plt.plot(T_cs_vec,NEP_johnson*1e18,label='Johnson')
plt.plot(T_cs_vec,NEP_tot*1e18,label='Total')
plt.text(0.05, 0.85, '93GHz LAT, P_sat = 3.5pW', transform=ax.transAxes)
plt.xlabel('T_base')
plt.ylabel('NEP [aW/rtHz]')
plt.grid()
plt.legend(loc='center right')


# 30 GHz LAT
detector['PsatFactor'] = 3.0
detector['Qtot'] = 0.20e-12   # optical power
detector['W'] = detector['PsatFactor']*detector['Qtot']

NEP_johnson = np.array([])
NEP_phonon  = np.array([])

for T_coldstage in T_cs_vec:
    detector['T_base'] = T_coldstage
    detector['T_bolo'] = 1.6*T_coldstage

    detector = tes_calcs.tes_calcs(detector)
    NEP_johnson = np.append(NEP_johnson, detector['Johnson']['NEP'])
    NEP_phonon = np.append(NEP_phonon, detector['Phonon']['NEP'])

NEP_photon = 0*NEP_phonon + 4.0e-18
NEP_tot = np.sqrt(NEP_photon**2 + NEP_phonon**2 + NEP_johnson**2)

ax = plt.subplot(313)
plt.plot(T_cs_vec,NEP_photon*1e18, label='Photon')
plt.plot(T_cs_vec,NEP_phonon*1e18, label='Phonon')
plt.plot(T_cs_vec,NEP_johnson*1e18,label='Johnson')
plt.plot(T_cs_vec,NEP_tot*1e18,label='Total')
plt.text(0.05, 0.85, '30GHz LAT, P_sat = 0.6pW', transform=ax.transAxes)
plt.xlabel('T_base')
plt.ylabel('NEP [aW/rtHz]')
plt.grid()
plt.legend(loc='center right')
#plt.title('20GHz LAT, P_sat = 0.7pW, T_bolo = 1.6*T_base')



plt.show()
