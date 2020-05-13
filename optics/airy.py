import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp

plt.ion()

lam = 500e-9  # wavelength in meters
k = 2*np.pi/lam

###### Play with these to change pattern
a = 4*lam

theta_rayleigh = 1.22*lam/(2*a) * 180./np.pi
print('Rayleigh: ',theta_rayleigh)

theta_max = 89.
offset_deg = 8.7
#offset_deg = theta_rayleigh
offset= offset_deg*(np.pi/180.)

theta_deg = np.linspace(-theta_max,theta_max,1000)
theta = (np.pi/180.)*theta_deg

alpha = (k*a)*np.sin(theta)
alpha2 = (k*a)*np.sin(theta+offset)

I1 = (sp.jv(1,alpha)/alpha)**2
I2 = (sp.jv(1,alpha2)/alpha2)**2

plt.figure(1)
plt.plot(theta_deg,I1)
#plt.plot(theta_deg,I2)
#plt.plot(theta_deg,I1+I2)
plt.grid()

plt.show()
