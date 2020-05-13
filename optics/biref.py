import pylab 
import numpy as np

pylab.ion()

def rotmat(theta):
    ''' return a 2x2 rotation matrix, for angle theta '''
    c = np.cos(theta)
    s = np.sin(theta)
    Mtheta = np.array([ [c,-s],[s,c] ])
    return Mtheta;

def getMpock(Vpock,theta_pock):
    ''' 
    return a retarder matrix for the pockels cell, given rotation angle of axes of pockels cell,
    which normally should be 45deg
    '''
    Vpock_max = 5000. # volts
    phi0 = 0  # phase at zero volts
    phimax = np.pi  #phase at maximum voltage
    phipock = phi0 + phimax*(Vpock/Vpock_max)
    M_pock = np.array([ [1, 0], [0, np.exp(1j*phipock)] ])
    M_pock = np.matmul(rotmat(-theta_pock),np.matmul(M_pock,rotmat(theta_pock)))
    return M_pock

def getMsamp(phisamp):
    '''
    Return a retarder matrix for the given phase delay, for the sample.
    '''
    M_sample = np.array([ [1, 0], [0, np.exp(1j*phisamp)] ])
    return M_sample

def getMtot(M_pol, M_sample, M_pock, M_analyzer):
    '''
    Return the total Jones matrix for all optical elements in the chain.
    '''
    Mtot = M_pol
    Mtot = np.matmul(Mtot,M_sample)
    Mtot = np.matmul(Mtot,M_pock)
    Mtot = np.matmul(Mtot,M_analyzer)
    return Mtot;

# define angle of first polarizer after laser here
M_horizpol = np.array([ [1,0],[0,0] ])
theta_pol = -45*np.pi/180.
M_pol = np.matmul(rotmat(-theta_pol),np.matmul(M_horizpol,rotmat(theta_pol)))

# define angle of analyzer polarizer here
theta_analyzer = 45*np.pi/180.
M_analyzer = np.matmul(rotmat(-theta_analyzer),np.matmul(M_horizpol,rotmat(theta_analyzer)))

# laser has to put out some light parallel to first polarizer... don't worry about normalization
E_in = np.array([[1],[0]])

# loop over sample birefringence phase delay between two axes
phisamps = np.arange(0,np.pi,0.3)
pylab.figure(1)
pylab.clf()
for phisamp in phisamps:
    M_sample = getMsamp(phisamp)
    theta_pock= 0
    Vpock = np.arange(0,5000,100)
    I = np.zeros(Vpock.size)
    ii=0
    # loop over pockel cell voltage 
    for V in Vpock:
        M_pock = getMpock(V,0*np.pi/180)
        M_tot = getMtot(M_pol,M_sample,M_pock,M_analyzer)
        # Output E field is M_tot*E_in
        E_out = np.matmul(M_tot,E_in)
        # calculate the intensity at the diode, which is the dot product of E fields.
        intensity = np.asscalar(np.real(np.dot(E_out.T,np.conj(E_out))))
        I[ii] = intensity
        ii = ii + 1
    # plot Intensity vs Pockel cell voltage
    sss = str(phisamp) + ' rad'
    pylab.plot(Vpock,I,label=sss)

pylab.legend()
pylab.xlabel('Pockels cell voltage')
pylab.ylabel('Intensity of light at diode')
pylab.title('Intensity vs. V_pockels, for different sample delta_phi')


