import numpy as nm

def calculate_hwp_mueller_matrix(use_preset_thickness_and_index='150GHz', \
                                 cooled=True, \
                                 n_s=False,n_f=False, \
                                 n_ar=False,n_gap=False, \
                                 d_hwp=False,d_ar_near=False,d_ar_far=False, \
                                 d_gap_near=False, d_gap_far=False, \
                                 freq=nm.linspace(135,165,31)):
    # SAB
    #
    # I've been reading "Microwave Engineering" by Pozar to sort of get my bearings for KIDs stuff
    # Both because Spider's python code could use it, and also as a good exercise,
    # this code calculates the HWP mueller matrix vs frequency
    # but using equations just from that textbook and my thesis.
    
    # first determine if we are using presets for the indices of refraction and the thicknesses
    # input "False" here, and none of this preset stuff will trigger, and the values entered
    # as function arguments will be used instead
    if use_preset_thickness_and_index:
        # then see if we are cold or not
        if cooled:
            # set sapphire indices of refraction, 5 K measured values
            n_s=3.336
            n_f=3.019

            # set the preset thicknesses and materials for each band
            if use_preset_thickness_and_index=='150GHz':
                n_ar=1.935 # index of refraction of the cirlex AR, 5 K estimate
                d_ar_near=0.254 # mm thickness
                d_ar_far=0.254
                n_gap=1.51 # index of refraction of the bond layer, Lamb measured value
                d_gap_near=0.006 # mm thickness
                d_gap_far=0.006
                d_hwp=3.160 # mm thickness

            if use_preset_thickness_and_index=='95GHz':
                n_ar=1.951 # index of refraction of the fused quartz AR, Lamb measured value
                d_ar_near=0.427 # mm thickness
                d_ar_far=0.427
                n_gap=1.0 # index of refraction of the air gap
                d_gap_near=0.000 # mm thickness - consider a nominal 10 micron air gap
                d_gap_far=0.000
                d_hwp=4.930 # mm thickness

        if not cooled:
            # set sapphire indices of refraction, 295 K measured values
            n_s=3.3736
            n_f=3.0385

            # set the preset thicknesses and materials for each band
            if use_preset_thickness_and_index=='150GHz':
                n_ar=1.935 # index of refraction of the cirlex AR, 5 K estimate
                d_ar_near=0.254 # mm thickness
                d_ar_far=0.254
                n_gap=1.51 # index of refraction of the bond layer, Lamb measured value
                d_gap_near=0.006 # mm thickness
                d_gap_far=0.006
                d_hwp=3.160 # mm thickness

            if use_preset_thickness_and_index=='95GHz':
                n_ar=1.951 # index of refraction of the fused quartz AR, Lamb measured value
                d_ar_near=0.427 # mm thickness
                d_ar_far=0.427
                n_gap=1.0 # index of refraction of the air gap
                d_gap_near=0.000 # mm thickness - consider a nominal 10 micron air gap
                d_gap_far=0.000
                d_hwp=4.930 # mm thickness

    # convert frequency from GHz
    f_Hz = freq*1.0e9
    # convert the thicknesses to meters
    d_hwp = d_hwp*1.0e-3
    d_ar_near = d_ar_near*1.0e-3
    d_ar_far = d_ar_far*1.0e-3
    d_gap_near = d_gap_near*1.0e-3
    d_gap_far = d_gap_far*1.0e-3
    # convert the indices to impedances
    c_light = 299792458.0 #m/s, exact speed of light
    Z0 = 119.9169832*nm.pi #ohms, exact impedance of free space, from wikipedia
    Z_s = Z0/n_s
    Z_f = Z0/n_f
    Z_ar = Z0/n_ar
    Z_gap = Z0/n_gap

    # pre-allocate the arrays to store the four transmission mueller matrix elements vs freqnency
    # equation 8.8 in my thesis
    T = nm.zeros(len(freq))
    rho = nm.zeros(len(freq))
    c = nm.zeros(len(freq))
    s = nm.zeros(len(freq))
    # and also the reflection mueller matrix elements too
    # equation 9.6 in my thesis
    R = nm.zeros(len(freq))
    tau = nm.zeros(len(freq))
    h = nm.zeros(len(freq))
    q = nm.zeros(len(freq))
    # loop over frequencies
    # this can't be vectorized because a seperate set of matrix multiplications
    # is required to simulate each individual frequency
    for i in xrange(len(freq)):
        # set up ABCD matrices for the AR coat and gap from table 4.1 in Pozar
        bl_ar=((2.0*nm.pi*f_Hz[i])/c_light)*(Z0/Z_ar)*d_ar_near
        ABCD_ar_near=nm.matrix([[nm.cos(bl_ar),1j*Z_ar*nm.sin(bl_ar)],[(1j/Z_ar)*nm.sin(bl_ar),nm.cos(bl_ar)]])
        bl_ar=((2.0*nm.pi*f_Hz[i])/c_light)*(Z0/Z_ar)*d_ar_far
        ABCD_ar_far=nm.matrix([[nm.cos(bl_ar),1j*Z_ar*nm.sin(bl_ar)],[(1j/Z_ar)*nm.sin(bl_ar),nm.cos(bl_ar)]])
        bl_gap=((2.0*nm.pi*f_Hz[i])/c_light)*(Z0/Z_gap)*d_gap_near
        ABCD_gap_near=nm.matrix([[nm.cos(bl_gap),1j*Z_gap*nm.sin(bl_gap)],[(1j/Z_gap)*nm.sin(bl_gap),nm.cos(bl_gap)]])
        bl_gap=((2.0*nm.pi*f_Hz[i])/c_light)*(Z0/Z_gap)*d_gap_far
        ABCD_gap_far=nm.matrix([[nm.cos(bl_gap),1j*Z_gap*nm.sin(bl_gap)],[(1j/Z_gap)*nm.sin(bl_gap),nm.cos(bl_gap)]])

        
        # set up two ABCD matrices, one for each index of the sapphire
        bl_s=((2.0*nm.pi*f_Hz[i])/c_light)*(Z0/Z_s)*d_hwp
        ABCD_s=nm.matrix([[nm.cos(bl_s),1j*Z_s*nm.sin(bl_s)],[(1j/Z_s)*nm.sin(bl_s),nm.cos(bl_s)]])
        bl_f=((2.0*nm.pi*f_Hz[i])/c_light)*(Z0/Z_f)*d_hwp
        ABCD_f=nm.matrix([[nm.cos(bl_f),1j*Z_f*nm.sin(bl_f)],[(1j/Z_f)*nm.sin(bl_f),nm.cos(bl_f)]])

        # multiply out the ABCD matrices for the slow index
        # note that since I am using matrices, NOT numpy arrays
        # you can just do normal readable matrix multiplication and it will work
        ABCD_tot_slow = ABCD_ar_near*ABCD_gap_near*ABCD_s*ABCD_gap_far*ABCD_ar_far
        # again for the fast index
        ABCD_tot_fast = ABCD_ar_near*ABCD_gap_near*ABCD_f*ABCD_gap_far*ABCD_ar_far

        # use table 4.2 in Pozar to convert these ABCD matrices to scattering matrices
        # and extract the relevant terms
        # for slow-axis co-polar thru transmission, we want s21 from the slow axis scattering matrix
        A = ABCD_tot_slow[0,0]; B = ABCD_tot_slow[0,1]; C = ABCD_tot_slow[1,0]; D = ABCD_tot_slow[1,1]
        t_slow = 2.0 / (A + B/Z0 + C*Z0 + D)
        # for co-polar reflection, we want s11
        r_slow = (A + B/Z0 - C*Z0 - D)/(A + B/Z0 + C*Z0 + D)
        # for fast-axis co-polar thru transmission, we want s21 from the fast axis scattering matrix
        A = ABCD_tot_fast[0,0]; B = ABCD_tot_fast[0,1]; C = ABCD_tot_fast[1,0]; D = ABCD_tot_fast[1,1]
        t_fast = 2.0 / (A + B/Z0 + C*Z0 + D)
        # for co-polar reflection, we want s11
        r_fast = (A + B/Z0 - C*Z0 - D)/(A + B/Z0 + C*Z0 + D)

        # now build up jones matrices from these results, following equations 8.9 and 8.10 in my thesis
        # to fill in the entries in equation 8.3 in my thesis
        J_ret = nm.matrix([[t_slow,0.0],[0.0,t_fast]])
        # same for reflections
        J_ref = nm.matrix([[r_slow,0.0],[0.0,r_fast]])

        # now convert from the Jones matrix to a Mueller matrix using equations 8.4 and 8.5 in my thesis
        # really, we just want certain selected matrix elements
        sigma1 = nm.matrix([[1.0,0.0],[0.0, 1.0]])
        sigma2 = nm.matrix([[1.0,0.0],[0.0,-1.0]])
        sigma3 = nm.matrix([[0.0,1.0],[1.0, 0.0]])
        sigma4 = nm.matrix([[0.0,-1j],[ 1j, 0.0]])
        # the matrix elements for transmission are
        T[i]   = 0.5*nm.trace( nm.real( sigma1*J_ret*sigma1*J_ret.H ) )
        rho[i] = 0.5*nm.trace( nm.real( sigma1*J_ret*sigma2*J_ret.H ) )
        c[i]   = 0.5*nm.trace( nm.real( sigma3*J_ret*sigma3*J_ret.H ) )
        s[i]   = 0.5*nm.trace( nm.real( sigma4*J_ret*sigma3*J_ret.H ) )
        # and for reflection
        R[i]   = 0.5*nm.trace( nm.real( sigma1*J_ref*sigma1*J_ref.H ) )
        tau[i] = 0.5*nm.trace( nm.real( sigma2*J_ref*sigma2*J_ref.H ) )
        h[i]   = 0.5*nm.trace( nm.real( sigma3*J_ref*sigma3*J_ref.H ) )
        q[i]   = 0.5*nm.trace( nm.real( sigma4*J_ref*sigma3*J_ref.H ) )
        # the quantites already are real numbers
        # using the real-part function is just so the trace function
        # doesn't display a warning that it is taking the real part

    return freq,T,rho,c,s,R,tau,h,q

def calculate_detector_couplings(theta_hwp,psi = 0.0,isA=True,T=0.97683,rho=0.00962,c=-0.95992,s=-0.01622,eta=1.0,delta=0.0):
    # convert theta_hwp from degrees
    theta_hwp = theta_hwp * (nm.pi/180.0)

    # set whether this is an A or B detector
    if isA:
        pass
    else :
        # must be a B detector
        # interchange the eta and delta values
        eta,delta = (delta,eta)

    # calculate the common factors in equations 8.18 and 8.19 in my thesis
    F = -(1.0/4.0)*         (T-c)*nm.sin(4.0*theta_hwp) *(eta**2.0 - delta**2.0) - (1.0/2.0)*rho*nm.sin(2.0*theta_hwp)*(delta**2.0 + eta**2.0)
    G =  (1.0/4.0)*(T + c + (T-c)*nm.cos(4.0*theta_hwp))*(eta**2.0 - delta**2.0) + (1.0/2.0)*rho*nm.cos(2.0*theta_hwp)*(delta**2.0 + eta**2.0)

    # calculate matrix elements in equation 8.17 in my thesis
    M_II = (1.0/2.0) * (T*(eta**2.0 + delta**2.0) + rho*nm.cos(2.0 * theta_hwp)*(eta**2.0 - delta**2.0))
    M_IQ =  F * nm.sin(2.0*psi) + G * nm.cos(2.0*psi)
    M_IU = -F * nm.cos(2.0*psi) + G * nm.sin(2.0*psi)
    M_IV = (1.0/2.0) * s * nm.sin(2.0*theta_hwp)*(eta**2.0 - delta**2.0)

    # return the matrix elements
    return M_II,M_IQ,M_IU,M_IV
