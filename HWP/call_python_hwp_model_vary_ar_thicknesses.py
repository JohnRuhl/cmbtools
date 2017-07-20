import python_hwp_model
import numpy as nm

# this code calculates and displays the band-averaged HWP mueller matrix elements
# while varying the thickness of one of the AR coat layers
# leaving the other thickness unchanged

# calculate HWP mueller matrix elements for nominal thicknesses of everything
freq,T,rho,c,s,R,tau,h,q = python_hwp_model.calculate_hwp_mueller_matrix(use_preset_thickness_and_index=False, \
	                                                    freq=nm.linspace(85,105,2100), \
	                                                    n_ar=1.951, \
                										d_ar_near=0.427, \
                										d_ar_far=0.427, \
                										n_gap=1.0, \
                										d_gap_near=0.000, \
                										d_gap_far=0.000, \
                										d_hwp=4.930, \
                										n_s=3.336, \
                										n_f=3.019)

# print the band-averaged mueller matrix elements (from 85 GHz to 105 GHz)
print 'Nominal AR Coat Thicknesses (both are 427 microns):'
print '(Agrees with top row of Table 8.2 in my thesis astro-ph/1402.2591)'
print 'T   = ' + str(nm.mean(T))
print 'rho = ' + str(nm.mean(rho))
print 'c   = ' + str(nm.mean(c))
print 's   = ' + str(nm.mean(s))
print ' '
print 'Reflection = ' + str(nm.mean(R))
print ' '
print ' '
print ' '

# calculate HWP mueller matrix elements for nominal thicknesses of everything
# EXCEPT make one of the AR coats a different thickness
varied_thickness = 0.385
freq,T,rho,c,s,R,tau,h,q = python_hwp_model.calculate_hwp_mueller_matrix(use_preset_thickness_and_index=False, \
	                                                    freq=nm.linspace(85,105,2100), \
	                                                    n_ar=1.951, \
                										d_ar_near=0.427, \
                										d_ar_far=varied_thickness, \
                										n_gap=1.0, \
                										d_gap_near=0.000, \
                										d_gap_far=0.000, \
                										d_hwp=4.930, \
                										n_s=3.336, \
                										n_f=3.019)

# print the band-averaged mueller matrix elements (from 85 GHz to 105 GHz)
print 'Varied One AR Coat Thicknesses: ' + str(varied_thickness * 1000.0) + ' microns'
print 'T   = ' + str(nm.mean(T))
print 'rho = ' + str(nm.mean(rho))
print 'c   = ' + str(nm.mean(c))
print 's   = ' + str(nm.mean(s))
print ' '
print 'Reflection = ' + str(nm.mean(R))
print ' '
print ' '
print ' '


# calculate HWP mueller matrix elements for nominal thicknesses of everything
# EXCEPT make one of the AR coats a different thickness
varied_thickness = 0.458
freq,T,rho,c,s,R,tau,h,q = python_hwp_model.calculate_hwp_mueller_matrix(use_preset_thickness_and_index=False, \
	                                                    freq=nm.linspace(85,105,2100), \
	                                                    n_ar=1.951, \
                										d_ar_near=0.427, \
                										d_ar_far=varied_thickness, \
                										n_gap=1.0, \
                										d_gap_near=0.000, \
                										d_gap_far=0.000, \
                										d_hwp=4.930, \
                										n_s=3.336, \
                										n_f=3.019)

# print the band-averaged mueller matrix elements (from 85 GHz to 105 GHz)
print 'Varied One AR Coat Thicknesses: ' + str(varied_thickness * 1000.0) + ' microns'
print 'T   = ' + str(nm.mean(T))
print 'rho = ' + str(nm.mean(rho))
print 'c   = ' + str(nm.mean(c))
print 's   = ' + str(nm.mean(s))
print ' '
print 'Reflection = ' + str(nm.mean(R))
print ' '
print ' '
print ' '


# calculate HWP mueller matrix elements for nominal thicknesses of everything
# EXCEPT make one of the AR coats a different thickness
varied_thickness = 0.395
freq,T,rho,c,s,R,tau,h,q = python_hwp_model.calculate_hwp_mueller_matrix(use_preset_thickness_and_index=False, \
	                                                    freq=nm.linspace(85,105,2100), \
	                                                    n_ar=1.951, \
                										d_ar_near=0.427, \
                										d_ar_far=varied_thickness, \
                										n_gap=1.0, \
                										d_gap_near=0.000, \
                										d_gap_far=0.000, \
                										d_hwp=4.930, \
                										n_s=3.336, \
                										n_f=3.019)

# print the band-averaged mueller matrix elements (from 85 GHz to 105 GHz)
print 'Varied One AR Coat Thicknesses: ' + str(varied_thickness * 1000.0) + ' microns'
print 'T   = ' + str(nm.mean(T))
print 'rho = ' + str(nm.mean(rho))
print 'c   = ' + str(nm.mean(c))
print 's   = ' + str(nm.mean(s))
print ' '
print 'Reflection = ' + str(nm.mean(R))
print ' '
print ' '
print ' '
