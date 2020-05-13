import matplotlib.pyplot as plt
import numpy as np

import optical_calcs
import DetectorDefs
import tes_calcs

plt.ion()


NET_tots = {}
bfacs = [1.0, 0.0]
for bfac in bfacs:
    detector, optics = DetectorDefs.S4_SAT_85()
    detector['BoseFactor'] = bfac
    detector['NEI_squid'] = 5.e-14 # a factor of 100 lower to take it out of contention.
    detector, optics = optical_calcs.optical_calcs(detector, optics)

    ap_eff_vec = np.arange(0.75,0.96,0.01)
    ap_eps_vec = 1.0 - ap_eff_vec

    blabel = str(bfac)
    NET_tots[blabel] = np.array([])
    for eps in ap_eps_vec:
        for thing in optics:
            if thing['name'] == '5K lenses':
                thing['eps'] = eps
        detector, optics = optical_calcs.optical_calcs(detector, optics)
        detector = tes_calcs.tes_calcs(detector)
        NET_tots[blabel] = np.append(NET_tots[blabel],detector['NoiseTotal']['NET_CMB']) 

plt.figure(1)
plt.clf()
plt.subplot(211)
plt.title('S4_SAT_85')
for bfac in bfacs:
    blabel = str(bfac)
    plt.plot(ap_eff_vec,NET_tots[blabel]/np.sqrt(2), label = blabel)
plt.legend()
plt.xlabel('Aperture efficiency')
plt.ylabel('NET (uKrtsec)')
plt.grid()

plt.subplot(212)
for bfac in bfacs:
    blabel = str(bfac)
    plt.plot(ap_eff_vec,NET_tots[blabel]/NET_tots[blabel][10], label = blabel)
plt.legend()
plt.xlabel('Aperture efficiency')
plt.ylabel('Re-normalized NET')
plt.grid()





plt.show()
