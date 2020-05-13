import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage

plt.ion()

dist = 2000 #mm
center=[512,512]

def radial_profile(data, center):
    sx,sy = data.shape # first determine radii of all pixels
    X, Y = np.ogrid[0:sx, 0:sy]
    r = np.hypot(X - sx/2, Y - sy/2)*xysize
    theta_deg = np.arctan(r/dist)*180/np.pi
    print(theta_deg.max())
    nbins = 300
    rbin = (nbins*theta_deg/theta_deg.max()).astype(np.int)
    thetabins = np.arange(0,theta_deg.max(),theta_deg.max()/nbins)
    radial_mean = ndimage.mean(data, labels=rbin, index=np.arange(0, rbin.max()))
    #radial_mean = ndimage.sum(data, labels=rbin, index=np.arange(0, rbin.max() ))
    #
    return thetabins, radial_mean

#mapp = np.loadtxt('150_input_beam_2m.txt',comments='#')
#mapp = np.loadtxt('band_edge_beams/170_input_beam_2m_highres2.txt',comments='#')

#plt.figure(1)
#plt.clf()
#nmapp = mapp/mapp.max()
#lmapp = np.log10(nmapp)
#plt.imshow(mapp)
#plt.colorbar()
#plt.title('Log10(beammap)')
#

#xysize = 5.3896 #mm 170, highres2

# 7.48mm is patch edge length;  this is 7.8mm pixel pitch, according to lorenzo

xysize = 14.099  #mm 130, reg
mapp1 = np.loadtxt('band_edge_beams/130_input_beam_2m.txt',comments='#')
thetabins1, Ptheta1 = radial_profile(mapp1,center)

xysize = 12.332  #mm 150, reg
mapp2 = np.loadtxt('band_edge_beams/150_input_beam_2m.txt',comments='#')
thetabins2, Ptheta2 = radial_profile(mapp2,center)

xysize = 10.779  #mm 170, reg
mapp3 = np.loadtxt('band_edge_beams/170_input_beam_2m.txt',comments='#')
thetabins3, Ptheta3 = radial_profile(mapp3,center)

theta = np.arange(0,70.,0.1)
Ptheta_150 = np.interp(theta, thetabins2, Ptheta2)
Ptheta = Ptheta_150 + np.interp(theta, thetabins1, Ptheta1)
Ptheta = Ptheta + np.interp(theta, thetabins3, Ptheta3)
Ptheta = Ptheta/Ptheta.max()
Ptheta2 = Ptheta2/Ptheta2.max()


# Scale by hand to average over a frequency band
nu0 = 150 #GHz
Ptheta_scaled = np.zeros(theta.size)
for nu in [130,135,140,145,150,155,160,165,170]:
    theta_1 = thetabins2*nu0/nu
    Ptheta_scaled = Ptheta_scaled + np.interp(theta,theta_1, Ptheta2)
Ptheta_scaled = Ptheta_scaled/Ptheta_scaled.max()


# P(theta) for 150, band avg, and gaussian
plt.figure(2)
plt.clf()
plt.semilogy(thetabins2, Ptheta2, label = '150 patch')
plt.semilogy(theta,Ptheta_scaled, label='scaled band avg')
plt.semilogy(theta,Ptheta, label='(130,150,170) avg')
plt.grid()

data_pro = np.loadtxt('profiled_feeds/155_band_profile_9p4mm.txt',comments='#',unpack=True)
profiled_9p4_150 = data_pro[:][1]
theta_pro = data_pro[:][0]
plt.semilogy(theta_pro,profiled_9p4_150)


lambda_center = 2.0e-3
gauss_horn_diam = 10.5e-3  #m
w0 = gauss_horn_diam/3.125     
theta0 = np.arctan(lambda_center/(np.pi*w0))
theta_fwhm_feed = 1.18*theta0
sigma = theta_fwhm_feed/2.355
theta_rad = theta*np.pi/180.
#P_horn = 0.005*np.exp(-theta**2/(2*25**2)) + 
P_horn = np.exp(-theta_rad**2/(2*sigma**2))
P_horn = P_horn/P_horn.max()
plt.semilogy(theta,P_horn, label='Gauss')
plt.ylim([1e-4,2.0])
plt.legend()
plt.xlim([0,70])
plt.xlabel('Degrees')
plt.ylabel('Radial avg P(theta)')
plt.title('For patch side length = guassian feed diameter')




#Calculate the encircled energy by doing an integral of r^2*P(theta)
#150
Penc_1 = np.cumsum(Ptheta2*np.sin(thetabins2*np.pi/180.))
Penc_1 = Penc_1/Penc_1.max()
# bandavg
Penc_2 = np.cumsum(Ptheta_scaled*np.sin(theta*np.pi/180.))
Penc_2 = Penc_2/Penc_2.max()
# Gaussian feed
Penc_3 = np.cumsum(P_horn*np.sin(theta*np.pi/180.))
Penc_3 = Penc_3/Penc_3.max()


plt.figure(3)
plt.clf()
rp,Ppatch = np.loadtxt('150_input_beam_encircled_2m.txt',comments='#', unpack=True)
plt.plot(np.arctan(rp/2000.)*180/np.pi,Ppatch, label='from zemax')
plt.plot(theta,Penc_3, label='7.48mm Gauss feed + wide Gauss sidelobe')
plt.plot(thetabins2,Penc_1,'--',label='150GHz, 7.48 patch')
plt.plot(theta,Penc_2,'-.',label='band')
plt.legend()
plt.xlim([0,70])
plt.grid()

plt.figure(4)
plt.clf()
plt.plot(thetabins2, Ptheta2, label = '150 Spider patch')
#plt.plot(theta,Ptheta_scaled, label='scaled band avg')
plt.plot(theta,Ptheta, label='(130,150,170) avg')
plt.plot(theta,P_horn, label='Gaussian, 10.5mm')
plt.plot(theta_pro,profiled_9p4_150, label='Profiled, 9.4mm')
plt.legend()
plt.grid()


plt.show()


