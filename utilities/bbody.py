import numpy 

#-------------------------------------------------
def Bnu(nu,T):  # Calculate the planck blackbody brightness in W/m^2/sr/Hz
    '''
    Planck blackbody brightness function (2 polarizations, in W/m^2/sr/Hz).
    Frequency (nu) must be in Hz.
    Bnu(nu,T) = 2.*(h*nu)*(nu/c)**2 * (1./(numpy.exp(x) - 1.))
    '''
    "Blackbody brightness function (2 polarization, in W/m^2/sr/Hz)"

    # Frequency vector must be in Hz.
    # Physical constants in SI units
    c= 2.99792458e8  
    h = 6.626068e-34
    k = 1.3806503e-23
    x = (h*nu)/(k*T)

    # Factor of 2 on RHS means this is for 2 polarizations
    brightness = 2.*(h*nu)*(nu/c)**2 * (1./(numpy.exp(x) - 1.))

    return brightness


#-------------------------------------------------
def Brj(nu,T):
    '''
    Rayleigh-Jeans brightness function (2 polarization, in W/m^2/sr/Hz).
    Frequency (nu) must be in Hz.
    Brj(nu,T) = 2.*k*T*(nu/c)**2
    '''

    # Frequency vector must be in Hz.
    # Physical constants in SI units
    c= 2.99792458e8  
    h = 6.626068e-34
    k = 1.3806503e-23

    # Factor of 2 on RHS means this is for 2 polarizations
    brightness = 2.*k*T*(nu/c)**2

    return brightness

#-------------------------------------------------
def Trj(B,nu):
    '''
    Rayleigh-Jeans temperature as a function 2-polarizaiton brightness in W/m^2/sr/Hz).
    Frequency (nu) must be in Hz.
    Trj = B(nu,T)/( 2.*k*(nu/c)**2 )
    '''

    # Frequency vector must be in Hz.
    # Physical constants in SI units
    c= 2.99792458e8  
    h = 6.626068e-34
    k = 1.3806503e-23

    # Factor of 2 on RHS means this is for 2 polarizations
    Trj = B/( 2.*k*(nu/c)**2 )

    return Trj

#-------------------------------------------------
def dBnudT(nu,T):
    '''
    Derivative of Bnu(2 polarization, in W/m^2/sr/Hz) with respect to T.
    Frequency (nu) must be in Hz.
    '''

    # Frequency vector must be in Hz.
    # Physical constants in SI units
    c= 2.99792458e8  
    h = 6.626068e-34
    k = 1.3806503e-23
    x = (h*nu)/(k*T)

    # Factor of 2 on RHS means this is for 2 polarizations
    prefac = 2*h**2/(k*c**2)
    expx = numpy.exp(x)
    dBdT = prefac*(nu**4/T**2)*(expx/(expx - 1)**2)

    return dBdT

#-------------------------------------------------
def dBnudTrj(nu,T):
    '''
    Derivative of the Rayleigh Jeans brightness, 
    Brj(2 polarization, in W/m^2/sr/Hz) with respect to T.
    Frequency (nu) must be in Hz.
    '''

    # Frequency vector must be in Hz.
    # Physical constants in SI units
    c= 2.99792458e8  
    k = 1.3806503e-23

    # Factor of 2 on RHS means this is for 2 polarizations
    dBdT = 2.*k*(nu/c)**2

    return dBdT

#-------------------------------------------------
def dTrjdT(nu,T):
    '''
    Derivative of dTrj/dT
    Frequency (nu) must be in Hz.
    '''

    # Frequency vector must be in Hz.
    # Physical constants in SI units
    c= 2.99792458e8
    h = 6.626068e-34
    k = 1.3806503e-23
    x = (h*nu)/(k*T)

    prefac = 2*h**2/(k*c**2)
    expx = numpy.exp(x)
    dBdT = prefac*(nu**4/T**2)*(expx/(expx - 1)**2)

    dBdTrj = 2.*k*(nu/c)**2

    dTrjdT = dBdT/dBdTrj

    return dTrjdT

#-------------------------------------------------

def BBintegrate_1mode(nulow,nuhi,T):
    '''
    Integrate Bnu from nuhi to nulow, return total, for a single mode single polarization detector.
    '''
    c= 2.99792458e8
    dnu = (nuhi-nulow)/1000.
    nuvec = numpy.arange(nulow,nuhi,dnu)

    lambdavec = c/nuvec
    AOmegavec = lambdavec**2   # single mode
    integrand = AOmegavec*Bnu(nuvec,T)

    power = numpy.trapz(integrand,nuvec)

    return power 
#-------------------------------------------------

def BBintegrate_surface(nulow,nuhi,T,AOmega):
        '''
        Integrate Bnu from nuhi to nulow, return total, for a single mode single polarization detector.
        '''
        c= 2.99792458e8
        dnu = (nuhi-nulow)/1000.
        nuvec = numpy.arange(nulow,nuhi,dnu)

        lambdavec = c/nuvec
        integrand = AOmega*Bnu(nuvec,T)

        power = numpy.trapz(integrand,nuvec)

        return power
#-------------------------------------------------


