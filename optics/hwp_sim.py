import pylab 
import numpy as np
from sympy import Matrix
from sympy import Symbol

pylab.ion()

def rotmat(theta):
    c = np.cos(theta)
    s = np.sin(theta)
    Mtheta = np.array([ [1,0,0,0],[0,c,s,0],[0,-s,c,0],[0,0,0,1] ])
    return Mtheta;

def getMtot(Mdet,Mhwp,xi,theta,psi):
    Mtheta = rotmat(2*theta_hwp)
    Mntheta = rotmat(-2*theta_hwp)
    Mxi = rotmat(2*xi_det)
    Mpsi = rotmat(2*psi)
    Mtot = np.matmul(Mdet,Mxi)
    Mtot = np.matmul(Mtot,Mntheta)
    Mtot = np.matmul(Mtot,Mhwp)
    Mtot = np.matmul(Mtot,Mtheta)
    Mtot = np.matmul(Mtot,Mpsi)
    return Mtot;

deg2rad = np.pi/180.
rad2deg = 1/deg2rad

# HWP matrix in its own basis
T = 1.00
c = -1.00
rho = 0.00
s = 0.00
Mhwp = np.array([ [T, rho, 0, 0], [rho,T,0,0],[0,0,c,-s],[0,0,s,c] ] )

# Detector, with partial polarizer
eta = 1.00
delta = 0.00
x = (eta*eta + delta*delta)/2.
y =  (eta*eta - delta*delta)/2.
Mdet = np.array([ [x,y,0,0],[y,x,0,0],[0,0,eta*delta,0],[0,0,0,eta*delta] ] )

theta_hwp = 0.0*deg2rad     # HWP rotation angle
psi_inst = 0.0*deg2rad      # Instrument rotation angle
xi_det = 0.0*deg2rad    # Detector rotation angle

psi_inst = deg2rad*np.arange(0,360)
M_II = np.array([])
M_IQ = np.array([])
M_IU = np.array([])
M_IV = np.array([])
for psi in psi_inst:
    Mtot = getMtot(Mdet,Mhwp,xi_det,theta_hwp,psi)
    M_II = np.append(M_II,Mtot[0][0])
    M_IQ = np.append(M_IQ,Mtot[0][1])
    M_IU = np.append(M_IU,Mtot[0][2])
    M_IV = np.append(M_IV,Mtot[0][3])


pylab.figure(1)
pylab.clf()
pylab.subplot(4,1,1)
pylab.plot(psi_inst,M_II)
pylab.subplot(4,1,2)
pylab.plot(psi_inst,M_IQ)
pylab.subplot(4,1,3)
pylab.plot(psi_inst,M_IU)
pylab.subplot(4,1,4)
pylab.plot(psi_inst,M_IV)






