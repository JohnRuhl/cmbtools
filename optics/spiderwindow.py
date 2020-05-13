import pylab
import numpy as np
# SlabPlus1LayerAR.py- Program to calculate reflectance of e and o axes 
# of sapphire with one AR coating layer.
#

pylab.ion()

#  Set physical constants
c = 3.0e8    # m/s
e0 = 8.85*10**-12
u0 = 4*np.pi*10**-7
n0 = 1.000
epsmu = np.sqrt(e0/u0)
Y0 = epsmu*n0

# Set instrument parameters
band_center_ghz = 280. # GHz
lam_center = c/(band_center_ghz*1e9)  # meters, design center of AR coat

# Set plotting range
nu_vector_ghz = np.arange(220,320,.1)
nu_vector = nu_vector_ghz * 1e9

# properties of material and AR coat
n_slab = 1.52           # index of the main material slab
d_slab = 3.175e-3       # thickness  of the main slab
n_ar = 1.22            # index of the ar_coat
d_ar = 0.125e-3  # meters
#n_ar_ideal = np.sqrt(n_slab)            # index of the ar_coat
#d_ar_ideal = (lam_center/4.)/n_ar_ideal  #thickness of ideal AR coat
#opt_thickness_ideal = n_ar_ideal*d_ar_ideal

#n_ar = 1.0
#d_ar = opt_thickness_ideal/n_ar

# define admittances
Y_slab = epsmu*n_slab
Y_ar = epsmu*n_ar

Refl = np.empty(0)
Trans = np.empty(0)
# loop over frequencies
for nu in nu_vector:
    #  Set k*h for each layer and axis
    factor1 = 2*np.pi*(nu/c)
    kh_ar = n_ar*d_ar*factor1
    kh_slab = n_slab*d_slab*factor1

    # Define M matrices for each layer
    A = np.cos(kh_ar)
    B = np.sin(kh_ar)
    M_ar = np.array([[A, 1j*B/(Y_ar)], [1j*(Y_ar)*B, A ]])

    A = np.cos(kh_slab)
    B = np.sin(kh_slab)
    M_slab = np.array([[A, 1j*B/Y_slab],[1j*Y_slab*B, A]])

    #  Multiply themj to get the total M matrix
    M_tot = np.matmul(M_ar,M_slab) 
    M_tot = np.matmul(M_tot,M_ar)

    # Find the reflection coefficients
    denom = (Y0*M_tot[0,0] + 
            Y0*Y0*M_tot[0,1] + 
            M_tot[1,0] + 
            Y0*M_tot[1,1])
    
    r = (Y0*M_tot[0,0] + 
            Y0*Y0*M_tot[0,1] - 
            M_tot[1,0] - 
            Y0*M_tot[1,1])/denom

    t = 2*Y0/denom


    Refl = np.append(Refl, r.real**2 + r.imag**2 )
    Trans = np.append(Trans, t.real**2 + t.imag**2) 
    
data = np.loadtxt('UHMWPEWindowLoss280.txt',comments='#',delimiter=',')
nu_meas_ghz = data[:,0]
trans_meas = data[:,1]/100.
refl_meas = data[:,2]/100.


pylab.figure(1)
pylab.subplot(211)
pylab.plot(nu_vector_ghz, Trans)
pylab.plot(nu_meas_ghz,trans_meas, '.')
pylab.ylabel('Transmission')
pylab.xlabel('Frequency  (GHz)')

pylab.subplot(212)
pylab.plot(nu_vector_ghz, Refl)
pylab.plot(nu_meas_ghz,refl_meas, '.')
pylab.ylabel('Reflection')
pylab.xlabel('Frequency  (GHz)')


