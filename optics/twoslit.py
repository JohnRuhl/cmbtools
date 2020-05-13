import matplotlib.pyplot as plt
import numpy as np

plt.ion()

lam = 500e-9  # wavelength in meters
k = 2*np.pi/lam

###### Play with these to change pattern
b = 4*lam    # slit width
a = 50.*lam   # slit separation


theta_max = 89.

theta_deg = np.linspace(-theta_max,theta_max,10000)
theta = (np.pi/180.)*theta_deg

beta = (k*b/2)*np.sin(theta)
alpha = (k*a/2)*np.sin(theta)

I = (np.sin(beta)**2 / beta**2) * np.cos(alpha)**2

plt.figure(1)
plt.plot(theta_deg,I)

plt.show()
