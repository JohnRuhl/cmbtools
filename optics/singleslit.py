import numpy as np
import matplotlib.pyplot as plt

plt.ion()
plt.clf()

# single slit diffraction pattern, a sinc function

a = 10.e-6  # slit width in meters
D = 1 # meter to screen
lam = 500e-9 # wavelength in meters

# index of refraction of medium
n= 1.

# angles over which to calculate
theta_max= np.pi/16.
theta = np.linspace(-theta_max, theta_max, 1000)
theta_deg = (180./np.pi)*theta
x = D * np.tan(theta)


beta = np.pi*(a/lam)*np.sin(theta)
I = (np.sin(beta)/beta)**2
I = I/I.max()

#plt.plot(theta_deg,I, theta_deg,I_2) #plot sinc^2 as well
plt.plot(1000.*x,I) 
#plt.ylim(bottom=0)
plt.xlabel('x (mm)')
plt.ylabel('I')

plt.show()
