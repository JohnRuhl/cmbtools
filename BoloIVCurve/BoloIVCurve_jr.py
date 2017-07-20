# -*- coding: utf-8 -*-
"""
Created on Tue Sep 02 12:42:11 2014

@author: sma120
"""
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
#import pickle
import scipy.optimize as sciopt

#import csv

def Findk(n, Tc, Tcs, Psat):
  # finds "k" for bolometer leg thermal conductivity, Psat = k*(Tc**n - Tcs**n)
  k = (Psat)/(Tc**n - Tcs**n)
  return k

def FitRT(Tb, Rb, Tc):
    # Fits data (Tb,Rb)to the function RTFitFunc, 
    # currently implemented as RTFitFunc(1,100,Tc,0)
    # where RTFitFunc is the logistic equation ?? (Scott please correct this.)
    def RTFitFunc(t, a, b, c, d):
       return a/(1+np.exp(-b*(t-c))) + d
       
    coefs, var = sciopt.curve_fit(RTFitFunc, Temps, Resists, p0=(1, 100, Tc, 0))
    T = np.arange(0, np.amax(Temps), 0.001*np.amax(Temps))
    RT = np.array([])
    RT = np.append(RT, RTFitFunc(T, *coefs))
    
    # Eventually take this out, as plotting should occur in the calling script
    plt.figure('R(T) vs T')
    plt.plot(Temps, Resists, 'k.', label = 'R(T) from Data')
    plt.plot(T, RT, 'c.', label = 'R(T) from Fit')
    
    TransitionParameter = coefs[1]  # b in above formula, like a width
    TcFromFit = coefs[2]  # c in above formula, like the Tc
    
    return RT, T, TransitionParameter, TcFromFit

def TestPlot(Voltage_at_min, Temps_data, Temps_data_n, Resists_data, n, k, Tcs, Tc, Current_at_min, Temp_at_min, Resists_from_fit, Temps_from_fit):
        #Vmin, Tb, Tbn, Rb, n, k, Tcs, Tc, Imin, Tmin, RT, T):
   #Plot R(T) and terms in function, to ensure zeros can be found in FindT
   #Also provides a visual check on the R(T) fit

   yr = Voltage_at_min**2/Resists_data
   yt = k*(Temps_data_n - Tcs**n)
   t = np.arange(0,len(Temps_data), 1)
 
   plt.figure('Zero Check and R(T)')
   plt.subplot(211)
   plt.plot(t, yr, 'r')
   plt.plot(t, yt, 'c')
   plt.xlabel('i')
   plt.ylabel('Rterm (Red) and Tterm (Cyan)')
   plt.title('Zero Check')
   
   plt.subplot(212)
   plt.plot(Temps_data, Resists_data, 'k.')
   plt.plot(Temps_from_fit, Resists_from_fit, 'c.')
   plt.plot(Temp_at_min, Voltage_at_min/Current_at_min, 'r.', ms = 10, label = 'Minimum of IV Curve')
   plt.xlabel('T (K)')
   plt.ylabel(r'R ($\Omega$)')
   plt.title('R(T)')
   
   plt.subplots_adjust(hspace=0.5)
   
   """
   plt.figure('R(T) from data vs fit')
   plt.plot(Tb, Rb, 'k.', label='From data')
   plt.plot(Tb, RT, 'c.', label='From fit')
   plt.legend()
   """
   
   plt.show()

def newton(func, x0, fprime=None, args=(), tol=1.48e-8, maxiter=50, fprime2=None):
    
    if tol <= 0:
        raise ValueError("tol too small (%g <= 0)" % tol)
        
    if fprime is not None:
        # Newton-Rapheson method
        # Multiply by 1.0 to convert to floating point.  We don't use float(x0)
        # so it still works if x0 is complex.
        p0 = 1.0 * x0
        fder2 = 0
        
        for iter in range(maxiter):
            myargs = (p0,) + args
            fder = fprime(*myargs)
            
            
            if fder == 0:
                msg = "derivative was zero."
                print msg
                return p0
            fval = func(*myargs)
            
            
            if fprime2 is not None:
                fder2 = fprime2(*myargs)
                
                
            if fder2 == 0:
                # Newton step
                p = p0 - fval / fder
                
            else:
                # Parabolic Halley's method
                discr = fder ** 2 - 2 * fval * fder2
                
                if discr < 0:
                    p = p0 - fder / fder2
                    
                else:
                    p = p0 - 2*fval / (fder + np.sign(fder) * np.sqrt(discr))

            if abs(p - p0) < tol:
                return p
            
            p0 = p
            
            if np.isnan(p0):
                p0 = 0.
            
    else:
        # Secant method
        p0 = x0
        
        if x0 >= 0:
            p1 = x0*(1 + 1e-4) + 1e-4
            
        else:
            p1 = x0*(1 + 1e-4) - 1e-4
            
        q0 = func(*((p0,) + args))
        q1 = func(*((p1,) + args))
        
        for iter in range(maxiter):
            if q1 == q0:
                if p1 != p0:                    
                    msg = "Tolerance of %s reached" % (p1 - p0)
                    print msg                    
                return (p1 + p0)/2.0
                
            else:
                p = p1 - q1*(p1 - p0)/(q1 - q0)
                
            if abs(p - p1) < tol:
                return p
                
            p0 = p1
            q0 = q1
            p1 = p
            q1 = func(*((p1,) + args))
            
    msg = "Failed to converge after %d iterations, value is %s" % (maxiter, p)
    raise RuntimeError(msg)

def FindR(Resists,Temps,Ti):
    # Given data (Resists,Temps), return an R that is appropriate for Ti via linear interpolation or extrapolation
   i = 0
 
   while (i == 0):
      j=0
      
      if Ti < np.amin(Temps): #Looking at a temperature below whats in our Tb array, linear extrapolation
         minarg = np.argmin(Temps)
         slope=(Resists[minarg + 1]-Resists[minarg])/(Temps[minarg + 1]-Temps[minarg])
         R = Resists[minarg] - slope*(Temps[minarg] - Ti)
         i = 1
             
      elif Ti > np.amax(Temps): #Looking at a temp above whats in our Tb array, linear extrapolation
         maxarg = np.argmax(Temps)
         slope = (Resists[maxarg] - Resists[maxarg + 1])/(Temps[maxarg] - Temps[maxarg + 1])
         R = Resists[maxarg] + slope * (Ti - Temps[maxarg])
         i = 1
      
      else:
      
         for j in xrange(len(Temps)): #If the temp matches a Tb value exactly, just grab the R
            if Ti == Temps[j]:
               R = Resists[j]
               i = 1

            elif Ti < Temps[j] and Ti > Temps[j-1]: #If it's between two, interpolate
               slope = (Resists[j]-Resists[j-1])/(Temps[j]-Temps[j-1])
               R = Resists[j-1] + slope * (Ti-Temps[j-1])
               i = 1
            
            elif Ti < Temps[j-1] and Ti > Temps[j]: #Same, but for arrays arranged in the opposite direction
               slope = (Resists[j-1]-Resists[j])/(Temps[j-1]-Temps[j])
               R = Resists[j] + slope * (Ti-Temps[j])
               i = 1
            
   return R
    
def FindT(Voltage, Resists_from_fit, Temps_from_fit, Temps_from_fit_n, n, k, Tcs, itermax, tolerance, dT):

#def FindT(V, Rb, Tb, Tbn, n, k, Tcs, itermax, tolerance, dT):
   """
   Finds the temperature of the bolometer for a specific voltage across the bolometer.
   Uses Newton's method to find where V^2/R = k(Tbolo^n - Tstage^n)
   V = bias voltage across bolometer,
   (Rb,Tb) = bolometer resistance vs temperature relationship (data points)
   Tbn = Tbolo**n 
   n = index of thermal leg conduction, Psat = k*(Tb**n - Tcs**n)
   k = prefactor in that equation
   Tcs = coldstage temp
   """
   
   funcValues = np.array([])

   """   
   Find initial guess to use in Newton's method
   """   
   Rterm = Voltage**2/Resists_from_fit   #  power dissipated in bolo
   Tterm = k*(Temps_from_fit - Tcs**n)  # power required to heat bolo to Tbn
   funcValues = abs(Rterm - Tterm)  # Want this to be zero
   
   Tbmin = np.argmin(funcValues)
   To = Temps_from_fit[Tbmin]
   
   """
   Defines function so that V^2/R = k(T^n - Tcs^n)    
   Uses Newton's Method to find the Temperatures that 
   satisfy this relation.                             
   """
   def func(Ti, Voltage, Resists_from_fit, Temps_from_fit, k, n, Tcs, dT):
      # this calculates the function that we're trying to find the zero of,
      # which means the two sides of the power equation are balanced.
      R1 = FindR(Resists_from_fit, Temps_from_fit, Ti)
      return (Voltage**2/R1)-k*(Ti**n-Tcs**n) # given that T and R, what is the value of the function?
         
   def funcder(Ti, Voltage, Resists_from_fit, Temps_from_fit, k, n, Tcs, dT):
      R1 = FindR(Resists_from_fit, Temps_from_fit, Ti)
      R2 = FindR(Resists_from_fit, Temps_from_fit, Ti+dT)
      dRdT = (R2-R1)/dT
      return -(Voltage**2/R1**2)*dRdT-k*n*Ti**(n-1)

   def funcder2(Ti, Voltage, Resists_from_fit, Temps_from_fit, k, n, Tcs, dT):
      R1 = FindR(Resists_from_fit, Temps_from_fit, Ti)
      R2 = FindR(Resists_from_fit, Temps_from_fit, Ti+dT)
      R3 = FindR(Resists_from_fit, Temps_from_fit, Ti-dT)
      dRdT = (R2-R1)/dT
      dRdT2 = ((R2-R1)/dT - (R1-R3)/dT)/dT
      return (2*Voltage**2/R1**3)*((dRdT)**2)-(Voltage**2/R1**2)*dRdT2-k*n*(n-1)*Ti**(n-2)

   T = newton(func, To, fprime=funcder, args=(Voltage, Resists_from_fit, Temps_from_fit, k, n, Tcs, dT), tol=tolerance, maxiter=itermax, fprime2 = funcder2)

   return T

def CalcPower(V, R):
    Power = V**2/R
    return Power

def CalcCurrent(V, R):
   current = V/R
   return current

def CreateIVPlot(Vmax, minV, delV, R, T, Tbn, n, k, Tcs, itermax, tolerance, dT):
  
   Volt = np.array([])
   Cur = np.array([])
   Temp = np.array([])
   Resist = np.array([])
   
   """
   Step through voltages, find T and R, then find the current.
   """
   Volt = np.arange(minV, Vmax, delV)
      
   for i in xrange(len(Volt)): 
      V = Volt[i]
      Ti = FindT(V, R, Tb, Tbn, n, k, Tcs, itermax, tolerance, dT)
      Temp = np.append(Temp, Ti)
      Ri = FindR(R, T, Ti)
      Resist = np.append(Resist, Ri)
      
   Cur = CalcCurrent(Volt, Resist)
   
   return Volt, Cur, Temp, Resist

def FindIVMin(Volts, Curs, Tc):
    # find where the current is the smallest in the IV curve.
    #  (careful... it's the lowest value, not necessarily where the derivative is zero.)
    minarg = np.argmin(Curs)
    Vmin = Volts[minarg]
    Cmin = Curs[minarg]
    return Vmin, Cmin, minarg
    
def IVCFindR(Voltages, Currents, Volt):
    i = 0
 
    while (i == 0):
       j=0
      
       if Volt < np.amin(Voltages):
           minarg = np.argmin(Voltages)
           slope=(Currents[minarg + 1]-Currents[minarg])/(Voltages[minarg + 1]-Voltages[minarg])
           Cur = Currents[minarg] - slope*(Voltages[minarg] - Volt)
           i = 1
             
       elif Volt > np.amax(Voltages):
           maxarg = np.argmax(Voltages)
           slope = (Currents[maxarg] - Currents[maxarg + 1])/(Voltages[maxarg] - Voltages[maxarg + 1])
           Cur = Currents[maxarg] + slope * (Volt - Voltages[maxarg])
           i = 1
      
       else:
      
           for j in xrange(len(Voltages)):
               if Volt == Voltages[j]:
                   Cur = Currents[j]
                   i = 1

               elif Volt < Voltages[j] and Volt > Voltages[j-1]:
                   slope = (Currents[j]-Currents[j-1])/(Voltages[j]-Voltages[j-1])
                   Cur = Currents[j-1] + slope * (Volt-Voltages[j-1])
                   i = 1
            
               elif Volt < Voltages[j-1] and Volt > Voltages[j]:
                   slope = (Currents[j-1]-Currents[j])/(Voltages[j-1]-Voltages[j])
                   Cur = Currents[j] + slope * (Volt-Voltages[j])
                   i = 1
                   
    Resist = Volt/Cur
    
    return Resist, Cur
    
def FinddRdT(Resistances, Temperatures, dT):
    # ? Scott:  This needs cleaning and comments to document.
    dRdT = np.array([])
    
    T = Temperatures[0]
    R = FindR(RT, Tb, T+dT)  # Scott:  I don't get it; RT, Tb, and T aren't defined yet...
    drdt = (R-Resistances[0])/(dT)
    dRdT = np.append(dRdT, drdt)
    T = Temperatures[1]
    R1 = FindR(RT, Tb, T+dT)
    R2 = FindR(RT, Tb, T-dT)
    drdt = ((Resistances[1]-R2)/dT+(R1-Resistances[1])/dT)/2
    dRdT = np.append(dRdT, drdt)
    
    for i in range(2, len(Resistances)-2):
        T = Temperatures[i]
        R1 = FindR(RT, Tb, T+dT)
        R2 = FindR(RT, Tb, T+2*dT)
        R3 = FindR(RT, Tb, T-dT)
        R4 = FindR(RT, Tb, T- 2*dT)
        drdt = ((R1-Resistances[i])/dT + (R2-Resistances[i])/(2*dT) + (Resistances[i]-R3)/dT + (Resistances[i]-R4)/(2*dT))/4
        dRdT = np.append(dRdT, drdt)
    
    T = Temperatures[-2]
    R1 = FindR(RT, Tb, T-dT)
    R2 = FindR(RT, Tb, T+dT)
    drdt = ((Resistances[-2]-R1)/dT + (R2-Resistances[-2])/dT)/2
    dRdT = np.append(dRdT, drdt)
    T = Temperatures[-1]
    R = FindR(RT, Tb, T-dT)
    drdt = (Resistances[-1]-R)/dT
    dRdT = np.append(dRdT, drdt)
    
    return dRdT

def GetRP(V, I):
    
    R = V/I
    P = V*I

    dR = R[1:]-R[:-1]
    dP = P[1:]-P[:-1]
    dR = np.append(dR, R[-2]-R[-1])
    dP = np.append(dP, P[-2]-P[-1])
    
    dRdP = dR/dP
    
    return R, P, dRdP
    
def FindLG1(Voltages, Resistances, Temperatures, n, k, dT):
    # Given a list of (V_bias, Rbolo, Tbolo) values, find loop gain = P_elec*alpha/(Gdynamic*Tc)
    # The list (V_bias,Rbolo,Tbolo) presumably comes from analyzing an IV curve along with R(T) curve.
    #
    dPdT = n*k*(Temperatures**(n-1)) # == Gdynamic
    P = Voltages**2/Resistances # == P_elec
    dRdT = np.array([])

    dRdT = FinddRdT(Resistances, Temperatures, dT)
    A = (Temperatures/Resistances)*dRdT  # alpha

    LoopGains = (A*P)/(dPdT*Temperatures)
    
    return LoopGains
    
def FindLG2(Resistances, Temperatures, n, k, Tcs, dT):
    #  Given a list of (Rbolo,Tbolo), find the loop gain.
    #  This version calculates the electrical power from the thermal link equation
    #  so this all comes from an R(T) curve, no IV curve information.
    #
    dRdT = np.array([])
    dRdT = FinddRdT(Resistances, Temperatures, dT)
    A = (Temperatures/Resistances)*dRdT  # alpha
    
    G = n*k*Temperatures**(n-1)  # == Gdynamic
    
    P = k*(Temperatures**n - Tcs**n)  # power from the thermal link eqn, not the IV curve
    
    LG = (A*P)/(G*Temperatures)
    
    return LG
    
def FindLG3(Voltages, Resistances, Temperatures, Currents, n, k, dT):
    #  Similar to FindLG1, but uses P_elec = Ibolo*Vbolo rather than Vbolo**2/R_bolo.
    #   Scott: I'm not sure how this is different than FindLG1, 
    #    given that we get R_bolo from I_bolo and Vbolo.
    #
    dRdT = FinddRdT(Resistances, Temperatures, dT)   
    A = (Temperatures/Resistances)*dRdT     # alpha

    G = n*k*Temperatures**(n-1) # == Gdynamic

    P = Voltages*Currents  # power from Ibolo*Vbolo

    LoopGains = (A*P)/(G*Temperatures)
        
    return LoopGains
    
def FindPResistFromIV(I, V, lower_range, upper_range):
    # finds the dynamic R for a given range of the IV curve by fitting for the slope of V vs I.
    #  Most likely use of this is to find the residual resistance (normal) below the transition.
    #
    if len(I) or len(V) < 2:
        PResist = 0
    else:
        coefs = poly.polyfit(V[lower_range:upper_range], I[lower_range:upper_range], 1)
        PResist = 1/coefs[1]  # coefs[1] is the slope = dI/dV, so 1/coefs[1] is dV/dI = R (dynamic)
        
    return PResist
    
def FindLGFromIV(V, I, lower_range, upper_range, Imin):
    # Find the loop gain only from the IV curve, does not use R(T) curve.
    #
    PResist = FindPResistFromIV(I, V, lower_range, upper_range)
    V = V - I*PResist  # correct for residual resistance.
    DataR = V/I
    DataP = V*I

    #Susceptible to noise in I and V !!! Advise only using values from a fited curve!
    dR = DataR[1:] - DataR[:-1]
    dP = DataP[1:] - DataP[:-1]
    dRdP = dR/dP
    dRdP = np.append(dRdP, (DataR[-1]-DataR[-2])/(DataP[-1]-DataP[-2]))

    LG = DataP/DataR * dRdP   # mathematically equivalent to alpha*Pelec/(Gdyn*Tc)

    return LG
    
#Constants
n = 2.7
Tcs = 0.357
Tc = 0.539
itermax = 100000
SampleRate = 100
Vmax = 9e-6
minV = 1e-9
delV = 5e-9
tolerance = 1e-4
dT = 1e-4
Psat = 40e-12
Range = 50

#Multiplicative factors for R and T
RbUnitFactor = 1.35
TbUnitFactor = 1.

"""
#Pickle data files
filename = './RT.pkl'
rt_dic = pickle.load(open(filename,'r'))
filename2 = './IV.pkl'
iv_dic = pickle.load(open(filename2, 'r'))
"""

#Test Function
Tb = np.arange(0, 1.0, .0001)
Rb = 1/(1+np.exp(-100*(Tb-Tc)))
noise = np.random.normal(0,.025,len(Rb))
Rb = Rb + noise


"""
#Export data files as .csv files. Used to help fit data that doesn't seem to want to be fit
w = csv.writer(open("r_frac.csv", "w"))

for key in rt_dic['r_frac']:
    w.writerow([key])
"""
"""
Rb = np.array(rt_dic['r_frac'])*RbUnitFactor
Tb = np.array(rt_dic['t'])*TbUnitFactor
"""

while (len(Rb) > len(Tb)):
    Rb = np.delete(Rb,len(Rb)-1)
    
while (len(Rb) < len(Tb)):
    Tb = np.delete(Tb,len(Tb)-1)

plt.close('all')
#PResist = FindPResistFromIV(iv_dic, 0, 20)
PResist = 0
#Rb = Rb - PResist

k=Findk(n, Tc, Tcs, Psat)
RT, T, TransitionParameter, TcFromFit = FitRT(Tb, Rb, Tc)
print TcFromFit, 'Tc From Fit'
Tn = np.power(T, n)
Tbn = np.power(Tb, n)

Voltages, Currents, Temperatures, Resistances = CreateIVPlot(Vmax, minV, delV, RT, T, Tn, n, k, Tcs, itermax, tolerance, dT)
Vmin, Cmin, argmin = FindIVMin(Voltages, Currents, TcFromFit)
print Vmin, 'Vmin'
print Cmin, 'Cmin'
Ps = Vmin * Cmin
print Ps, 'Psat'
Tmin = Temperatures[argmin]
print Tmin, 'Tmin'

TestPlot(Vmin, Tb, Tbn, Rb, n, k, Tcs, TcFromFit, Cmin, Vmin, RT, T)

LGIV = FindLGFromIV(Voltages, Currents, 0, 0, Cmin)
LG1 = FindLG1(Voltages, Resistances, Temperatures, n, k, dT)
LG2 = FindLG2(Resistances, Temperatures, n, k, Tcs, dT)
LG3 = FindLG3(Voltages, Resistances, Temperatures, Currents, n, k, dT)

print LGIV[argmin], 'Loop Gain from IV Curve at minimum'
print LG1[argmin], 'Loop Gain from Method 1'
print LG2[argmin], 'Loop Gain from Method 2'
print LG3[argmin], 'Loop Gain from Method 3'

plt.figure('Loop Gain vs Voltage')
plt.plot(Voltages, LGIV, 'k.', label = 'LG from IV')
plt.plot(Voltages, LG1, 'c.', label = 'LG from Method 1')
plt.plot(Voltages, LG2, 'g.', label = 'LG from Method 2')
plt.plot(Voltages, LG3, 'y.', label = 'LG from Method 3')
plt.legend()
"""
plt.figure('Loop Gains from different methods')
plt.plot(Voltages, LGIV, 'k.', label = 'LG From IV')
plt.plot(Voltages, LG1, 'c.', label = 'LG from method 1')
plt.plot(Voltages, LG2, 'y.', label = 'LG from method 2')
plt.plot(Voltages, LG3, 'g.', label = 'LG from method 3')
plt.legend()
"""
