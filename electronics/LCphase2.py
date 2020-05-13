import pylab
import numpy as np

# plot phase of signal vs bias frequency, for variations around w0

pylab.ion()


f0 = 5500.e3
w0 = 2*np.pi*f0
f = np.arange(f0-100., f0+100., 1.)
w = 2*np.pi*f
Rbolo = 2.0 #ohms

# shunt
Rsh = 0.03
Zsh_complex = 0 # Rsh/10.
Lsh = Zsh_complex/(1j*w0)
Zshunt = Rsh + 1j*Zsh_complex

# comb
L = 60e-6
C = 1/(w0*w0*L)
Zcomb = 1j*w*L + 1/(1j*w*C) + Rbolo

#parallel of comb and shunt
Zparallel = 1/(1/Zcomb + 1/Zshunt)

# find currents
V1 = 1.
Rbias = 1e4
Itot = V1/(Rbias + Zparallel)

Vcomb = V1 - Itot*Rbias
I_shunt = Vcomb/Zshunt
I_bolo = Vcomb/Zcomb

phase_Ibolo = np.angle(I_bolo,deg=True)
phase_Ibolo = phase_Ibolo - phase_Ibolo[-1]

pylab.figure(1)
pylab.clf()

pylab.subplot(211)
pylab.plot(f,phase_Ibolo)


pylab.subplot(212)
pylab.plot(f,np.abs(Zcomb))


pylab.show()



