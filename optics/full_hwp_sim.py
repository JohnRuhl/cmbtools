import pylab 
import numpy as np
from sympy import Matrix
from sympy import Symbol

def rotmat(theta):
    c = np.cos(theta)
    s = np.sin(theta)
    Mtheta = np.array([ [1.0,0.0,0.0,0.0],[0.0,c,s,0.0],[0.0,-s,c,0.0],[0.0,0.0,0.0,1.0] ])
    return Mtheta;

def getMtot(Mdet,Mhwp,xi,theta,psi):
    Mtheta = rotmat(2.0*theta)
    Mntheta = rotmat(-2.0*theta)
    Mxi = rotmat(2.0*xi_det)
    Mpsi = rotmat(2.0*psi)
    Mtot = np.matmul(Mdet,Mxi)
    Mtot = np.matmul(Mtot,Mntheta)
    Mtot = np.matmul(Mtot,Mhwp)
    Mtot = np.matmul(Mtot,Mtheta)
    Mtot = np.matmul(Mtot,Mpsi)
    return Mtot;

deg2rad = np.pi/180.
rad2deg = 1/deg2rad

# HWP matrix in its own basis, using numbers for X1 for now
T = 0.975#1.00
c = -0.944#-1.00
rho = 0.01#0.00 #made something up because I don't have a number for this
s = -0.088#0.00
Mhwp = np.array([ [T, rho, 0.0, 0.0], [rho,T,0.0,0.0],[0.0,0.0,c,-s],[0.0,0.0,s,c] ] )

# Detector, with partial polarizer
eta = np.sqrt(0.95)#1.00 .995
delta = np.sqrt(0.05)#0.00 .005
x = (eta*eta + delta*delta)/2.
y =  (eta*eta - delta*delta)/2.
Mdet = np.array([ [x,y,0.0,0.0],[y,x,0.0,0.0],[0.0,0.0,eta*delta,0.0],[0.0,0.0,0.0,eta*delta] ] )

theta_hwp = deg2rad*np.arange(0,360)     # HWP rotation angle
xi_det = 0.0*deg2rad    # Detector rotation angle, leave this fixed at zero because it's degenerate with the other two angles
psi_inst = deg2rad*np.arange(0,360)


all_gamma=np.zeros((len(psi_inst)))
all_epsilon=np.zeros((len(psi_inst)))

all_four_theta=np.zeros((len(psi_inst)))
all_two_theta=np.zeros((len(psi_inst)))
all_dc=np.zeros((len(psi_inst)))

M_II_all = np.array([])
M_IQ_all = np.array([])
M_IU_all = np.array([])
M_IV_all = np.array([])
det_all = np.zeros((len(psi_inst), len(theta_hwp)))
fft_det_all = np.zeros((len(psi_inst), len(theta_hwp)))

for ind, psi in enumerate(psi_inst):
    #print 'Ind: '+str(ind) +' Psi:'+str(np.rad2deg(psi))
    #S_in=[1, np.cos(2.0*psi), np.sin(2.0*psi), 0]#This is the Stokes vector [I, Q, U, V].  It gets multiplied by M_tot to get the timestreams, and it depends on the instrument angle
    S_in=[1, 1, 0 , 0]#This is the Stokes vector [I, Q, U, V].  It gets multiplied by M_tot to get the timestreams, and it depends on the instrument angle
    M_II = np.array([])
    M_IQ = np.array([])
    M_IU = np.array([])
    M_IV = np.array([])
    det=np.zeros((len(theta_hwp)))
    for theta_ind, theta in enumerate(theta_hwp):
        Mtot = getMtot(Mdet,Mhwp,xi_det,theta,psi)
        M_II = np.append(M_II,Mtot[0][0])
        M_IQ = np.append(M_IQ,Mtot[0][1])
        M_IU = np.append(M_IU,Mtot[0][2])
        M_IV = np.append(M_IV,Mtot[0][3])
        det[theta_ind]=Mtot[0][0]*S_in[0]+Mtot[0][1]*S_in[1]+Mtot[0][2]*S_in[2]+Mtot[0][3]*S_in[3]#det is a function of HWP angle
    #Store stuff in case it's useful
    M_II_all=np.append(M_II_all, M_II)
    M_IQ_all=np.append(M_IQ_all, M_IQ)
    M_IU_all=np.append(M_IU_all, M_IU)
    M_IV_all=np.append(M_IV_all, M_IV)
    det_all[ind, :]=det
    #Take the FFT to get the signal components (DC, 2 theta, and 4 theta)
    fft_det=(4.0/len(det))*np.fft.fft(det)
    fft_det_all[ind, :]=fft_det
    four_theta=np.abs(fft_det[4])
    two_theta=np.abs(fft_det[2])
    dc=np.abs(fft_det[0]/2.0)
    all_four_theta[ind]=four_theta
    all_two_theta[ind]=two_theta
    all_dc[ind]=dc
    #Compute epsilon
    epsilon=(dc-four_theta)/(dc+four_theta)
    all_epsilon[ind]=epsilon
    #Compute gamma
    all_gamma[ind]=(1.0-epsilon)/(1.0+epsilon)


pylab.close('all')
pylab.figure(1)
pylab.plot(rad2deg*psi_inst, all_gamma)
pylab.xlabel('Instrument Angle ($\Psi_{inst}$)')
pylab.ylabel('Polarization Efficiency ($\gamma$)')
pylab.title('Polarization Efficiency vs. Instrument Angle for a Real HWP (X1)')

pylab.figure(2)
pylab.plot(rad2deg*theta_hwp, np.transpose(det_all))
pylab.xlabel('HWP Angle ($\theta_{HWP}$)')
pylab.ylabel('Timestream Signal')
pylab.title('Timestream Signals vs. HWP angle for all Instrument Angles')

pylab.figure(3)
pylab.plot(rad2deg*psi_inst, all_four_theta)
pylab.xlabel('Instrument Angle ($\Psi_{inst}$)')
pylab.ylabel('Four Theta Detector Signal')
pylab.title('Four Theta Detector Signal for a Real HWP (X1)')

pylab.figure(4)
pylab.plot(rad2deg*psi_inst, all_two_theta)
pylab.xlabel('Instrument Angle ($\Psi_{inst}$)')
pylab.ylabel('Two Theta Detector Signal')
pylab.title('Two Theta Detector Signal for a Real HWP (X1)')


#pylab.figure(1)
#pylab.clf()
#pylab.subplot(4,1,1)
#pylab.plot(psi_inst,M_II)
#pylab.subplot(4,1,2)
#pylab.plot(psi_inst,M_IQ)
#pylab.subplot(4,1,3)
#pylab.plot(psi_inst,M_IU)
#pylab.subplot(4,1,4)
#pylab.plot(psi_inst,M_IV)

pylab.show()






