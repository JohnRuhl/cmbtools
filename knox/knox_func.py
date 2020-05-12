#   Originally in matlab by J.Ruhl, converted to Python by S.Wen Nov. 2014
#   Evaluate knox-formula estimates for a given experiment.
#   returns the covariances for C_l bins.

import numpy
def knox_func(expt=None,baseCls=None):
    SecPerDay=24 * 60 * 60.
    ell=(baseCls[:,0]).astype(int)

    # convert from D_ell (camb output) to C_ell (what is needed for the calculations here)
    factor=ell*((ell + 1)) / (2 *numpy.pi)
    C_TT=baseCls[:,1] / factor
    C_EE=baseCls[:,2] / factor
    C_BB=baseCls[:,3] / factor
    C_TE=baseCls[:,4] / factor

#  create binned spectra
#  binX is the index
    lastbinstart=min(expt['topbin'],ell[-1] - 1 - expt['deltaell'])
    binstarts=numpy.arange(expt['lowestell'] - 1,lastbinstart+1,expt['deltaell'])  #this is an indexer
    binends = binstarts+expt['deltaell']-1

    l_b=numpy.zeros((len(binstarts),1))
    C_TT_b=numpy.zeros((len(binstarts),1))
    C_EE_b=numpy.zeros((len(binstarts),1))
    C_BB_b=numpy.zeros((len(binstarts),1))
    C_TE_b=numpy.zeros((len(binstarts),1))

    for ii in range(len(binstarts)):
        l_b[ii,0]=numpy.mean(ell[binstarts[ii]-1:binends[ii]])
        C_TT_b[ii,0]=numpy.mean(C_TT[binstarts[ii]-1:binends[ii]])
        C_EE_b[ii,0]=numpy.mean(C_EE[binstarts[ii]-1:binends[ii]])
        C_BB_b[ii,0]=numpy.mean(C_BB[binstarts[ii]-1:binends[ii]])
        C_TE_b[ii,0]=numpy.mean(C_TE[binstarts[ii]-1:binends[ii]])
    
    Cls=numpy.concatenate((l_b,C_TT_b,C_EE_b,C_BB_b,C_TE_b),1)
    #print(Cls)
    fullsky=41253 #square degrees
    
    fsky=expt['sqdeg'] / fullsky
    sigma_beam=(expt['fwhm']/ (numpy.sqrt(8 * numpy.log(2)))) * (1 / 60) * (numpy.pi / 180)
    sigb2=sigma_beam ** 2
    sig2_T=(expt['net'] / numpy.sqrt(expt['ndet'] * expt['days'] * SecPerDay)) ** 2
    sig2_P=2 * sig2_T  #in uK^2

    #Calculate the weights per unit solid angle
    w_T=(1 / sig2_T) / (fsky * 4 * numpy.pi)  #1/(uK^2 rad^2);
    w_P=(1 / sig2_P) / (fsky * 4 * numpy.pi)
    inv_w_T=1 / w_T
    inv_w_P=1 / w_P

    # Depth in uK-arcmin, for comparison purposes
    M = 1.485e8 # square arcminutes on the sphere
    D_T = numpy.sqrt(inv_w_T*M/(4*numpy.pi))
    D_P = numpy.sqrt(inv_w_P*M/(4*numpy.pi))
    print('T map depth = {0:.2f} uK_cmb-arcmin'.format(D_T))
    print('P map depth = {0:.2f} uK_cmb-arcmin'.format(D_P))

    A=1.0 / (fsky * (2.0 * l_b + 1) * expt['deltaell'])
    wB_T=inv_w_T * numpy.exp((l_b*l_b) * sigb2)    #weight times beam factor
    wB_P=inv_w_P * numpy.exp((l_b*l_b) * sigb2)

    cov_TT=2 * A*((C_TT_b + wB_T) ** 2)
    cov_EE=2 * A*((C_EE_b + wB_P) ** 2)
    cov_BB=2 * A*((C_BB_b + wB_P) ** 2)
    cov_cosvar_TT = 2*A*C_TT_b**2;
    cov_cosvar_EE = 2*A*C_EE_b**2;
    cov_cosvar_BB = 2*A*C_EE_b**2;
    cov_noise_TT =  2*A*wB_T**2;
    cov_noise_EE =  2*A*wB_P**2;
    cov_noise_BB =  2*A*wB_P**2;
    
    cov_TE=A*((C_TE_b ** 2 + (C_TT_b + wB_T)*((C_EE_b + wB_P))))
    cov_cosvar_TE=A*(C_TE_b ** 2 + C_TT_b*C_EE_b)
    cov_noise_TE =A*(wB_T*wB_P)

    cov_Cls=numpy.concatenate((l_b,cov_TT,cov_EE,cov_BB,cov_TE),1)
    cov_cosvar_Cls=numpy.concatenate((l_b,cov_cosvar_TT,cov_cosvar_EE,cov_cosvar_BB,cov_cosvar_TE),1)
    cov_noise_Cls =numpy.concatenate((l_b,cov_noise_TT,cov_noise_EE,cov_noise_BB,cov_noise_TE),1)

    # Note that the returned values are Cl's not Dell's.
    data={'Cls':Cls,'cov_Cls':cov_Cls,'cov_noise_cls':cov_noise_Cls, 'cov_cosvar_Cls':cov_cosvar_Cls}

    return data
