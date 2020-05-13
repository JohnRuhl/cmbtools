import pylab
import numpy as np
# HechtModels.py- Program to calculate reflectance and transmission
# of simple (non-birefringent) slabs
#

#  Set physical constants
c = 3.0e8    # m/s
e0 = 8.85*10**-12
u0 = 4*np.pi*10**-7
n0 = 1.000
epsmu = np.sqrt(e0/u0)
# Admittances
Y0 = epsmu*n0


def TRchisq(nu_vector, T_meas, R_meas, n_vector, d_vector):
    T, R = Find_TR(nu_vector, n_vector, d_vector)
    chisq_T = sum(T_meas - T)**2
    chisq_R = sum(R_meas - R)**2
    chisq = chisq_T + chisqR
    return chisq

def Find_TR(nu_vector, n_vector, d_vector):
    # nu : frequency at which to calculate T, R... in Hz.
    # n_vector:  contains n's for each slab.
    # d_vector:  contains thicknesses (d) of each slab, in meters.
    # returns Refl, Trans at frequency nu

    Y_vector = epsmu*n_vector

    Refl = np.empty(0)
    Trans = np.empty(0)
    for nu in nu_vector:
        #  Set k*h for each layer and axis
        factor1 = 2*np.pi*(nu/c)

        # Define M matrices for each layer, and use them to construct full thing.
        M_tot = np.identity(2)
        for layer in np.arange(n_vector.size):
            kh = n_vector[layer]*d_vector[layer]*factor1
            A = np.cos(kh)
            B = np.sin(kh)
            M_thislayer = np.array([[A, 1j*B/(Y_vector[layer])], 
                                    [1j*(Y_vector[layer])*B, A ]])
            #  Multiply themj to get the total M matrix
            M_tot = np.matmul(M_tot,M_thislayer)
            # END loop over layers to find M_tot at this frequency

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

        Refl = np.append(Refl, r.real**2 + r.imag**2)
        Trans = np.append(Trans, t.real**2 + t.imag**2 )
    
    return Trans, Refl
