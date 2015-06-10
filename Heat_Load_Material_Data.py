# -*- coding: utf-8 -*-

from Heat_Load_Calculator_Functions import *

# ------------------------------------------
# Data
'''
currently, the lamda data needs to be enclosed in [].
if a nondefined function is required, use "lambda x: expression".

need to make it so that alloys can be accessed from any of their component's dictionaries, not just one
'''

# Aluminum Data
aluminum1100 = Mat_Class('1100', {'mdl': 'log10', 'eq': lamda_log10_model, 'cfs': [23.39172, -148.5733, 422.1917, -653.6664, 607.0402, -346.152, 118.4276, -22.2781, 1.770187], 'eqr': [4, 300]}, UNS='A91100', fullname='Aluminum 1100')
aluminum3003F = Mat_Class('3003F', {'mdl': 'log10', 'eq': lamda_log10_model, 'cfs': [0.63736, -1.1437, 7.4624, -12.6905, 11.9165, -6.18721, 1.63939, -0.172667], 'eqr': [1, 300]}, UNS='A93003', fullname='Aluminum 3003F')
aluminum5083 = Mat_Class('5083', {'mdl': 'log10', 'eq': lamda_log10_model, 'cfs': [-0.90933, 5.751, -11.112, 13.612, -9.3977, 3.6873, -0.77295, 0.067336], 'eqr': [1, 300]}, UNS='A95083', fullname='Aluminum 5083')
aluminum6061T6 = Mat_Class('6061-T6', {'mdl': 'log10', 'eq': lamda_log10_model, 'cfs': [0.07918, 1.0957, -0.07277, 0.08084, 0.02803, -0.09464, 0.04179, -0.00571, 2.96344], 'eqr': [1, 300]}, UNS='A96061', fullname='Aluminum 6061-T6')
aluminum6063T5 = Mat_Class('6063-T5', {'mdl': 'log10', 'eq': lamda_log10_model, 'cfs': [22.401433, -141.13433, 394.95461, -601.15377, 547.83202, -305.99691, 102.38656, -18.810237, 1.4576882], 'eqr': [4, 295]}, UNS='A96063', fullname='Aluminum 6063-T5')
# Aluminum Dictionary
aluminum_dict = make_dict(aluminum1100, aluminum3003F, aluminum5083, aluminum6061T6, aluminum6063T5, name='Aluminum')

# Balsa Data
balsa6ro = Mat_Class('6 lb/ft^3', {'mdl': 'log10', 'eq': lamda_log10_model, 'cfs': [4172.447, -11309.97, 12745.09, -7647.584, 2577.309, 2577.309, -462.538, 34.5351], 'eqr': [70, 300]})
balsa11ro = Mat_Class('11 lb/ft^3', {'mdl': 'log10', 'eq': lamda_log10_model, 'cfs': [4815.4, -12969.63, 14520.76, -8654.164, 2895.712, -515.7272, 38.19218], 'eqr': [75, 300]})
# Balsa Dictionary
balsa_dict = make_dict(balsa6ro, balsa11ro, name='Balsa')

# Beechwood-Phenolic Data
beechwoodgd = Mat_Class('grain direct', {'mdl': 'log10', 'eq': lamda_log10_model, 'cfs': [-1375.11, 3740.69, -4238.465, 2559.333, -868.6067, 157.1018, -11.82957], 'eqr': [80, 300]}, fullname='Beechwood-Phenolic Grain Direction', info='cross-laminate [0/90], grain direction')
beechwoodf = Mat_Class('flatwise', {'mdl': 'log10', 'eq': lamda_log10_model, 'cfs': [1035.33, -2191.85, 1470.505, 39.845, -541.9035, 289.844, -65.2253, 5.59956], 'eqr': [80, 300]}, fullname='Beechwood-Phenolic Flatwise', info='cross-laminate [0/90], flatwise')
# Beechwood-Phenolic Dictionary
beechwood_dict = make_dict(beechwoodgd, beechwoodf, name='Beechwood')

# Beryllium-Copper Data
beryllium_copper = Mat_Class('Beryllium Copper', {'mdl': 'log10', 'eq': lamda_log10_model, 'cfs': [-0.50015, 1.93190, -1.69540, 0.71218, 1.27880, -1.61450, 0.68722, -0.10501], 'eqr': [1, 120]}, fullname='Beryllium Copper')

# Brass Data
brass = Mat_Class('Brass', {'mdl': 'log10', 'eq': lamda_log10_model, 'cfs': [0.021035, -1.01835, 4.54083, -5.03374, 3.20536, -1.12933, 0.174057, -0.0038151], 'eqr': [5, 110]}, UNS='C26000', fullname='Brass')

# Master Dictionary
# Adding in single-item definitions
materials = make_dict(beryllium_copper, brass, name='Master Dictionary')
# Adding in dictionaries
materials = make_dict(aluminum_dict, balsa_dict, beechwood_dict, dict=materials)


# Reference:
# One Liner:
# aluminum3003F = Mat_Class('3003F', {'mdl':'log10','eq':lamda_log10_model,'cfs':[0.63736, -1.1437, 7.4624, -12.6905, 11.9165, -6.18721, 1.63939, -0.172667],'eqr':[1, 300]}, UNS='A93003', fullname='Aluminum 3003F')

# 2 Liner:
# aluminum1100 = Mat_Class('1100', UNS='A91100', fullname='Aluminum 1100')
# aluminum1100.add_model('log10', lamda_log10_model, [23.39172, -148.5733, 422.1917, -653.6664, 607.0402, -346.152, 118.4276, -22.2781, 1.770187], eqrange=[4, 300], function=True)

# Multiple Lines:
# aluminum1100 = Mat_Class('1100', UNS='A91100', fullname='Aluminum 1100')
# aluminum1100.add_model('log10')
# aluminum1100.add_eqrange('log10', [4, 300])
# aluminum1100.add_equation('log10', lamda_log10_model)
# aluminum1100.add_coeffs('log10', [23.39172, -148.5733, 422.1917, -653.6664, 607.0402, -346.152, 118.4276, -22.2781, 1.770187])
# aluminum1100.add_function('log10')
