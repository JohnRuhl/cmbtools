import pylab
import numpy as np

pylab.ion()


linelabel = 'df = -200, Xsh = Rsh/10.'

f0 = 2500.e3
w0 = 2*np.pi*f0
f = f0 - 100.
w = 2*np.pi*f

# shunt
Rsh = 0.03
Zsh_complex = Rsh/10.
Lsh = Zsh_complex/(1j*w0)
Zshunt = Rsh + 1j*Zsh_complex

# comb
L = 60e-6
C = 1/(w0*w0*L)
Rbolo = np.arange(0.01,1.0,0.01)
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

pylab.subplot(211)
pylab.plot(Rbolo,phase_Ibolo,label=linelabel)
pylab.legend()
pylab.title('f0 = 2.5MHz, Rsh=0.03Ohms, L = 60uH')

pylab.subplot(212)
pylab.semilogx(Rbolo,Vcomb/np.abs(I_bolo))


pylab.show()



