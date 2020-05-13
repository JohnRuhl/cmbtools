import pylab
import numpy as np
from scipy.optimize import curve_fit
import HechtModels as hm

import lmfit

data = np.loadtxt('UHMWPEWindowLoss280_cut.txt',comments='#',delimiter=',')
nu_meas_ghz = data[:,0]
trans_meas = data[:,1]/100.
refl_meas = data[:,2]/100.
meas_data = np.column_stack((trans_meas, refl_meas))


def ThreeLayer(nu,n1,d1,n2,d2):
    n_vector = np.array([n1, n2, n1])
    d_vector = np.array([d1, d2, d1])
    T,R = hm.Find_TR(nu, n_vector, d_vector)
    data_out = np.column_stack((T,R))
    return data_out

ThreeLayer_model = lmfit.Model(ThreeLayer)
params = ThreeLayer_model.make_params(n1=1.2, d1=0.125e-3, n2=1.5, d2 = 3.1e-3)

nu = nu_meas_ghz*1e9
T_model = ThreeLayer_model.eval(params, nu=nu)

result = ThreeLayer_model.fit(meas_data, params, nu=nu)
T_bestfit = result.best_fit[:,0]
R_bestfit = result.best_fit[:,1]

pylab.ion()

pylab.figure(1)
pylab.clf()
pylab.subplot(211)
pylab.plot(nu_meas_ghz, T_bestfit)
pylab.plot(nu_meas_ghz,trans_meas, '.')
pylab.ylabel('Transmission')
pylab.xlabel('Frequency  (GHz)')

pylab.subplot(212)
pylab.plot(nu_meas_ghz, R_bestfit)
pylab.plot(nu_meas_ghz,refl_meas, '.')
pylab.ylabel('Reflection')
pylab.xlabel('Frequency  (GHz)')

pylab.show()

