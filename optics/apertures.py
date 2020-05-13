import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp

# This file defines a square field (aperture field), and tools you 
# can use to set up an E-field apodization.  
# The physical extent of the aperture field (==X, in meters) is defined
# at the top.  
# 
# Taking the FT of the aperture field distribution gives us the map
# of the k_x and k_y (for the aperture in the XY plane).  
# Note that k_x[1]*x must run from 0 to 2*pi across the aperture.
# that is, k_x[n] = n*(2*pi/D).
#
# To scale from our FT output in (k_x,k_y) to angular space (theta_x, theta_y)
# Note that k_tot = 2*pi/lambda, and 
# k_tot^2 = k_x^2 + k_y^2 + k_z^2 = (2*pi/lambda)^2, 
# tan(theta_x) = k_x/k_z, etc.
#

#
c = 3e8 # m/s

def inboxA(screen,x,y,center,fullwidth):
    # screen must be 2D square numpy array.
    # center must be a tuple, (xcenter,ycenter)
    n_elements = aperture.shape[0]
    middle = center
    xmin = middle - fullwidth/2
    xmax = middle + fullwidth/2
    ymin = middle - fullwidth/2
    ymax = middle + fullwidth/2
    for x in range(n_elements):
        for y in range(n_elements):
            if ( (x>xmin) and (x<xmax) and (y>ymin) and (y<ymax)): 
                val = 1
            else:
                val = 0
            screen[x,y] = val
    return screen

def incircleA(screen,x,y,center,radius):
    # screen must be 2D square numpy array.
    # center must be a tuple, (xcenter,ycenter)
    x0 = center[0]
    y0 = center[1]
    rpix = np.sqrt( (x-x0)**2 + (y-y0)**2)
    for x in range(n_elements):
        for y in range(n_elements):
            if ( rpix < radius):
                val = 1
            else:
                val = 0
            screen[x,y] = val
    return screen

def gaussian(screen, center, sigma):
    # screen must be 2D square numpy array.
    # center must be a tuple, (xcenter,ycenter)
    x0 = center[0]
    y0 = center[1]
    n_elements = screen.shape[0]
    for x in range(n_elements):
        for y in range(n_elements):
            rpix = np.sqrt( (x-x0)**2 + (y-y0)**2)
            val = np.exp(-(rpix**2)/(2*sigma**2))
            screen[x,y] = val
    return screen


def xslit(screen,x,y,xval,width):
    for x in range(n_elements):
        for y in range(n_elements):
            if np.abs(x-xval)<width:
                val = 1
            else:
                val =0
            screen[x,y] = val
    return screen

#def doubleslit(x,y):
#    x0 = 64 - 10
#    x1 = 64 + 10
#    width = 3
#    if (np.abs(x0-x)<width) or (np.abs(x1-x)<width):
#        val = 1
#    else:
#        val =0
#    return val
#
#def slitgrid():
#    gridsize = 512
#    mygrid=np.zeros((gridsize,gridsize))
#    x0 = np.arange(0,gridsize,int(gridsize/5.))
#    y0 = x0
#    mygrid[x0,:] = 1
#    mygrid[x0+1,:] = 1
#    mygrid[:,y0] = 1
#    mygrid[:,y0+1] = 1
#    return mygrid

def BeamFromAperture(screen,screensize,wavelength):
    # screen is 2D numpy array containing screen with E-field distribution
    # screensize is the physical length in meters across the (square) screen
    # wavelength is the wavelength of light we're considering
    # 
    # FFT to get kx, ky plane map
    # convert kx,ky map to angular scale map
    # return angular scale beam map (E-field;  square to get power)
    k_tot = 2*np.pi/wavelength
    n_elements = screen.shape[0]
    dx = screensize/n_elements  # spacing between elements in screen grid (meters)
    dk = 2*np.pi/screensize  # spacing between elements in kx,ky grid (radians/meters)
    kvec = dk*np.fft.fftfreq(n_elements)*n_elements # renormalize fftfreq grid 
    #
    # Generate 2D array of kx and ky values for each point in the k_plane
    kx_array = np.fft.fftshift(np.tile(kvec,(kvec.size,1)))
    ky_array = kx_array.transpose()
    thetax_array = np.arcsin(kx_array/k_tot)  # units of radians
    thetay_array = np.arcsin(ky_array/k_tot)  # units of radians

    # Fourier transform gives kx,ky vector amplitudes, corresponding to thetas above.
    E_beam = np.fft.fftshift(np.fft.fft2(screen))
    Beam_dict = {'E_beam':E_beam,'thetax_array':thetax_array,'thetay_array':thetay_array,'kx_array':kx_array,'ky_array':ky_array,'k_tot':k_tot}
    return Beam_dict

def BeamToScreen(screensize, Beam_dict, distance):
    # screen is xy physical screen (eg on aperture plane)
    # Beam_dict specifies angular pattern of E-field, and thetax and thetay arrays
    # distance is meters from beam generation point (ie base of angles)
    #   to screen.
    x = np.linspace(-screensize/2, screensize/2,screen.shape[0])
    y = x
    xx,yy = np.meshgrid(x,y)

    xxbeam = np.ndarray.flatten(np.tan(Beam_dict['thetax_array'])*distance)
    yybeam = np.ndarray.flatten(np.tan(Beam_dict['thetay_array'])*distance)
    points = np.vstack((xxbeam,yybeam)).T
    beam = np.ndarray.flatten(np.abs(Beam_dict['E_beam']))

    # ScreenBeam is an object we will use to do the right thing. 
    #ScreenBeam = interp.griddata((xxbeam,yybeam),beam, kind='linear', copy=True, 
     #          bounds_error=False, fill_value=0)
     #InterpScreen = np.abs(ScreenBeam(x,y))
    # griddata(points, values, (grid_x, grid_y), method='linear')
    # points is a numpy array that is npts by 2 (for x and y)
    # values is a numpy array that is npts long (func value at x,y)
    # grid_x and grid_y are meshgrid grids of x,y data points.
    InterpScreenValues = interp.griddata(points,beam, (xx,yy), method='linear', fill_value=0)

    return InterpScreen


def ApEff(beam,aperture_mask):
    # beam is E-field, ie square to get power pattern
    #  both beam and aperture_mask are on same physical (meters) grid
    power_pattern = beam*np.conj(beam)  
    denominator = np.sum(power_pattern)
    numerator = np.sum(power_pattern*aperture_mask)
    result = numerator/denominator
    return result


# Generate a screen, apodize it, and do stuff
D1 = 0.1  # meters across screen
N = 200
x0 = N/2.
y0 = x0
sigma_meters = 0.01 # meters
sigma_elements = (sigma_meters/D1)*N
screen = np.zeros((N,N))

distance = 3  # meter from first aperture to final screen.

# set up original screen.  
screen = gaussian(screen, (x0,y0),sigma_elements)

# Use that screen to make an angular beam;  D1 is physical extent of screen
Beam_dict = BeamFromAperture(screen, D1, 0.002)

# Project that beam on a second screen, a distance away
D2 = 1.0  # meters
InterpScreen = np.zeros((N,N))
InterpScreen = BeamToScreen(D2, Beam_dict, distance)


# plotting
#fftxcoord = thetax_array*180./np.pi
#fftycoord = thetay_array*180./np.pi

plt.figure(1)
plt.clf()
plt.imshow(screen,origin='lower',interpolation='none',extent=(-D1/2., D1/2., -D1/2., D1/2.))
plt.title('Screen1 (horn)')

plt.figure(2)
plt.clf()
thetamax = (180./np.pi)*Beam_dict['thetax_array'][0][-1]
thetamin = (180./np.pi)*Beam_dict['thetax_array'][0][0]
#
plt.imshow(np.abs(Beam_dict['E_beam']),origin='lower',interpolation='none',extent=[thetamin,thetamax,thetamin,thetamax])
plt.colorbar()
plt.title('Angular beam from horn')

plt.figure(3)
plt.clf()
plt.imshow(InterpScreen,origin='lower',interpolation='none',extent=(-D2/2., D2/2., -D2/2., D2/2.))
plt.title('E-field amplitude on aperture')


#plt.imshow(np.abs(G),origin='lower',interpolation='none')
#plt.colorbar()

plt.show()




