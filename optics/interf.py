import numpy as np
import matplotlib.pyplot as plt

plt.ion()
plt.clf()

# Number of elements
N = 700.

# length of line
L_over_lambda =  .03     # play with >1 vs <1

# spacing between elements
d_over_lambda = L_over_lambda/N

# index of refraction of medium
n= 1.

# angles over which to calculate
theta_max= np.pi/2.
theta = np.linspace(-theta_max, theta_max, 1000)
theta_deg = (180./np.pi)*theta

# phase 
delta = 2*np.pi * n*d_over_lambda *np.sin(theta)

# irradiance formula
I = np.sin(N*delta/2.)**2 / np.sin(delta/2.)**2
I = I/I.max()  # normalize it for convenience

beta = (np.pi*L_over_lambda)*np.sin(theta)
I_2 = (np.sin(beta)/beta)**2
I_2 = I_2/I_2.max()

#plt.plot(theta_deg,I, theta_deg,I_2) #plot sinc^2 as well
plt.plot(theta_deg,I_2) 
#plt.ylim(bottom=0)
plt.xlabel('Theta (deg)')
plt.ylabel('I')

plt.show()
