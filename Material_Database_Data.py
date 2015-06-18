# -*- coding: utf-8 -*-

from Material_Database_Functions import *

# ------------------------------------------
# Data
'''
if a nondefined function is required, use "lambda x: expression".
need to make it so that alloys can be accessed from any of their component's dictionaries, not just one
'''

# Aluminum Data
aluminum1100 = Mat_Class('1100', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [23.39172, -148.5733, 422.1917, -653.6664, 607.0402, -346.152, 118.4276, -22.2781, 1.770187], 'eqr': [4, 300], 'd_eq': lambda T: 0.02*T}, UNS='A91100', fullname='Aluminum 1100')
aluminum3003F = Mat_Class('3003F', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [0.63736, -1.1437, 7.4624, -12.6905, 11.9165, -6.18721, 1.63939, -0.172667], 'eqr': [1, 300], 'd_eq': lambda T: 0.02*T}, UNS='A93003', fullname='Aluminum 3003F')
aluminum5083 = Mat_Class('5083', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [-0.90933, 5.751, -11.112, 13.612, -9.3977, 3.6873, -0.77295, 0.067336], 'eqr': [1, 300], 'd_eq': lambda T: 0.01*T}, UNS='A95083', fullname='Aluminum 5083')
aluminum6061T6 = Mat_Class('6061-T6', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [0.07918, 1.0957, -0.07277, 0.08084, 0.02803, -0.09464, 0.04179, -0.00571, 2.96344], 'eqr': [1, 300], 'd_eq': lambda T: 0.005*T}, UNS='A96061', fullname='Aluminum 6061-T6')
aluminum6063T5 = Mat_Class('6063-T5', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [22.401433, -141.13433, 394.95461, -601.15377, 547.83202, -305.99691, 102.38656, -18.810237, 1.4576882], 'eqr': [4, 295], 'd_eq': lambda T: 0.02*T}, UNS='A96063', fullname='Aluminum 6063-T5')
# Aluminum Dictionary
aluminum_dict = make_dict(aluminum1100, aluminum3003F, aluminum5083, aluminum6061T6, aluminum6063T5, name='Aluminum')

# Balsa Data
balsa6ro = Mat_Class('6 lb/ft^3', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [4172.447, -11309.97, 12745.09, -7647.584, 2577.309, 2577.309, -462.538, 34.5351], 'eqr': [70, 300], 'd_eq': lambda T: 0.015*T}, fullname='Balsa Wood 6 lb/ft^3', info='lamda(T) from NIST website')
balsa11ro = Mat_Class('11 lb/ft^3', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [4815.4, -12969.63, 14520.76, -8654.164, 2895.712, -515.7272, 38.19218], 'eqr': [75, 300], 'd_eq': lambda T: 0.03*T}, fullname='Balsa Wood 11 lb/ft^3', info='lamda(T) from NIST website')
# Balsa Dictionary
balsa_dict = make_dict(balsa6ro, balsa11ro, name='Balsa')

# Beechwood-Phenolic Data
beechwoodgd = Mat_Class('grain direct', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [-1375.11, 3740.69, -4238.465, 2559.333, -868.6067, 157.1018, -11.82957], 'eqr': [80, 300], 'd_eq': lambda T: 0.01*T}, fullname='Beechwood-Phenolic Grain Direction', info='cross-laminate [0/90], grain direction. lamda(T) from NIST website')
beechwoodf = Mat_Class('flatwise', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [1035.33, -2191.85, 1470.505, 39.845, -541.9035, 289.844, -65.2253, 5.59956], 'eqr': [80, 300], 'd_eq': lambda T: 0.015*T}, fullname='Beechwood-Phenolic Flatwise', info='cross-laminate [0/90], flatwise. lamda(T) from NIST website')
# Beechwood-Phenolic Dictionary
beechwood_dict = make_dict(beechwoodgd, beechwoodf, name='Beechwood')

# Beryllium-Copper Data
beryllium_copper = Mat_Class('Beryllium Copper', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [-0.50015, 1.93190, -1.69540, 0.71218, 1.27880, -1.61450, 0.68722, -0.10501], 'eqr': [1, 120], 'd_eq': lambda T: 0.02*T}, fullname='Beryllium Copper', info='lamda(T) from NIST website')

# Brass Data
brass = Mat_Class('Brass', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [0.021035, -1.01835, 4.54083, -5.03374, 3.20536, -1.12933, 0.174057, -0.0038151], 'eqr': [5, 110], 'd_eq': lambda T: 0.015*T}, UNS='C26000', fullname='Brass', info='lamda(T) from NIST website')

# Carbon Fiber Data
carbonfiber = Mat_Class('Carbon Fiber', {'matProp': 'thermCond', 'name': 'fit', 'eq': thermCond_fit_UnivariateSpline_model, 'eqinput': [[0.3, 1.4, 4.2, 10.0, 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 275.0], [0.0018, 0.0118, 0.025, 0.055, 0.15, 0.175, 0.25, 0.3, 0.7, 0.9, 2.0, 3.0, 4.0, 5.0]], 'eqr': [0.3, 275], 'd_eq': None}, fullname='Carbon Fiber', info='Values < 4.2K from Runyan-Jones, 4K-50K estimated from plot in Chen et al., 50K-275K estimated from plot in Radcliffe et al.')

# Copper Data
# -------------------
# OFHC Copper Data
OFHCcopper50 = Mat_Class('RRR50', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [1.8743, -0.41538, -0.6018, 0.13294, 0.26426, -0.0219, -0.051276, 0.0014871, 0.003723], 'eqr': [4, 300], 'd_eq': None}, UNS='C10100/C10200', fullname='OFHC Copper RRR50', info='OFHC: Oxygen Free High Conductivity.  RRR: Residual-Resistance Ratio. lamda(T) from NIST website')
OFHCcopper100 = Mat_Class('RRR100', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [2.2154, -0.47461, -0.88068, 0.13871, 0.29505, -0.02043, -0.04831, 0.001281, 0.003207], 'eqr': [4, 300], 'd_eq': None}, UNS='C10100/C10200', fullname='OFHC Copper RRR100', info='OFHC: Oxygen Free High Conductivity.  RRR: Residual-Resistance Ratio. lamda(T) from NIST website')
OFHCcopper150 = Mat_Class('RRR150', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [2.3797, -0.4918, -0.98615, 0.13942, 0.30475, -0.019713, -0.046897, 0.0011969, 0.0029988], 'eqr': [4, 300], 'd_eq': None}, UNS='C10100/C10200', fullname='OFHC Copper RRR150', info='OFHC: Oxygen Free High Conductivity.  RRR: Residual-Resistance Ratio. lamda(T) from NIST website')
OFHCcopper300 = Mat_Class('RRR300', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [1.357, 0.3981, 2.669, -0.1346, -0.6683, 0.01342, 0.05773, 0.0002147], 'eqr': [4, 300], 'd_eq': None}, UNS='C10100/C10200', fullname='OFHC Copper RRR300', info='OFHC: Oxygen Free High Conductivity.  RRR: Residual-Resistance Ratio. lamda(T) from NIST website')
OFHCcopper500 = Mat_Class('RRR500', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [2.8075, -0.54074, -1.2777, 0.15362, 0.36444, -0.02105, -0.051727, 0.0012226, 0.0030964], 'eqr': [4, 300], 'd_eq': None}, UNS='C10100/C10200', fullname='OFHC Copper RRR500', info='OFHC: Oxygen Free High Conductivity.  RRR: Residual-Resistance Ratio. lamda(T) from NIST website')
# OFCH Copper Dictionary
OFHCcopper_dict = make_dict(OFHCcopper50, OFHCcopper100, OFHCcopper150, OFHCcopper300, OFHCcopper500, name='OFHC')
# Copper Dictionary
copper_dict = make_dict(OFHCcopper_dict, name='Copper')

# G-10 Data
G10 = Mat_Class('G-10', {'matProp': 'thermCond', 'name': 'Normal Direction', 'eq': thermCond_NIST_and_UnivariateSpline_model, 'eqinput': [[-4.1236, 13.788, -26.068, 26.272, -14.663, 4.4954, -0.6905, 0.0397], [[0.3, 1.4, 4.2], [1.64e-3, 2.06e-2, 6.56e-2]], [[10, 300], [0.3, 4.2]]], 'eqr': [0.3, 300], 'd_eq': lambda T: 0.05*T}, fullname='Fiberglass Epoxy')
G10.add_model({'matProp': 'thermCond', 'name': 'Warp Direction', 'eq': thermCond_NIST_model, 'eqinput': [-2.64827, 8.80228, -24.8998, 41.1625, -39.8754, 23.1778, -7.95635, 1.48806, -0.11701], 'eqr': [12, 300], 'd_eq': lambda T: 0.05*T})

# Kapton Data
kapton = Mat_Class('Kapton', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [5.73101, -39.5199, 79.9313, -83.8572, 50.9157, -17.9835, 3.42413, -0.27133], 'eqr': [4, 300], 'd_eq': lambda T: 0.02*T}, fullname='Kapton', info='lamda(T) from NIST website')

# Manganin Data
manganin = Mat_Class('Manganin', {'matProp': 'thermCond', 'name': 'fit', 'eq': thermCond_fit_UnivariateSpline_model, 'eqinput': [[4, 10, 20, 40, 77], [0.44, 1.4, 3.2, 6.8, 11]], 'eqr': [4, 77], 'd_eq': None}, fullname='Manganin', info="lambda(T) from Ekin's book")

# Phosphor-Bronze Data
phosphor_bronze = Mat_Class('Phosphor-Bronze', {'matProp': 'thermCond', 'name': 'fit', 'eq': thermCond_fit_UnivariateSpline_model, 'eqinput': [[4.2, 20, 100, 300], [1.7, 10, 28, 48]], 'eqr': [4.2, 300], 'd_eq': None}, fullname='Phosphor-Bronze', info='lamda(T) from Lakeshore values')

# Saphire Data
saphire = Mat_Class('Saphire', {'matProp': 'thermCond', 'name': 'fit', 'eq': thermCond_fit_UnivariateSpline_model, 'eqinput': [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 36.07, 40, 50, 60, 70, 80, 90, 100, 200, 300], [2.56, 17.89, 48.72, 100.0, 174.3, 267.6, 391.3, 543.2, 709.5, 884.7, 3425, 5463, 5722, 5708, 4773, 2887, 1401, 891.6, 591.3, 423.2, 74.73, 42.90]], 'eqr': [1, 300], 'd_eq': None}, fullname='Saphire', info='lambda(T) from Lakeshore')

# Steel Data
# -------------------
# Stainless Steel Data
stainless304 = Mat_Class('304', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [-1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256, -0.4658, 0.1650, -0.0199], 'eqr': [1, 300], 'd_eq': lambda T: 0.02*T}, UNS='S30400', fullname='stainless Steel 304', info='lamda(T) from NIST website')
stainless304L = Mat_Class('304L', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [-1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256, -0.4658, 0.1650, -0.0199], 'eqr': [1, 300], 'd_eq': lambda T: 0.02*T}, UNS='S30403', fullname='Stainless Steel 304L', info='lamda(T) from NIST website')
stainless310 = Mat_Class('310', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [-0.81907, -2.1967, 9.1059, -13.078, 10.853, -5.1269, 1.2583, -0.12395], 'eqr': [1, 300], 'd_eq': lambda T: 0.02*T}, UNS='S31000', fullname='Stainsess Steel 310', info='lamda(T) from NIST website')
stainless316 = Mat_Class('316', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [-1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256, -0.4658, 0.1650, -0.0199], 'eqr': [1, 300], 'd_eq': lambda T: 0.02*T}, UNS='S31600', fullname='stainless Steel 316', info='lamda(T) from NIST website')
# Stainless Steel Dictionary
stainless_dict = make_dict(stainless304, stainless304L, stainless310, stainless316, name='Stainless')
# Steel Dictionary
steel_dict = make_dict(stainless_dict, name='Steel')

# Teflon Data
teflon = Mat_Class('Teflon', {'matProp': 'thermCond', 'name': 'NIST', 'eq': thermCond_NIST_model, 'eqinput': [2.7380, -30.677, 89.430, -136.99, 124.69, -69.556, 23.320, -4.3135, 0.33829], 'eqr': [4, 300], 'd_eq': None}, fullname='Teflon', info='lamda(T) from NIST website')

# Vespel Data
vespelSP1 = Mat_Class('SP-1', {'matProp': 'thermCond', 'name': 'fit', 'eq': thermCond_fit_UnivariateSpline_model, 'eqinput': [[0.3, 1.4, 4.2], [5.53e-4, 3.21e-3,  9.74e-3]], 'eqr': [0.3, 4.2], 'd_eq': None}, fullname='Vespel SP-1', info='lamda(T) from Runyan-Jones values')
vespelSP22 = Mat_Class('SP-22', {'matProp': 'thermCond', 'name': 'fit', 'eq': thermCond_fit_UnivariateSpline_model, 'eqinput': [[0.3, 1.4, 4.2], [2.17e-4, 2.46e-3, 1.43e-2]], 'eqr': [0.3, 4.2], 'd_eq': None}, fullname='Vespel SP-22', info='lamda(T) from Runyan-Jones values')
# Vespel Dictionary
vespel_dict = make_dict(vespelSP1, vespelSP22, name='Vespel')


##################################################################

# Master Dictionary
materials = make_dict(name='Materials Dictionary')
# Adding in single-item definitions
materials = make_dict(beryllium_copper, brass, carbonfiber, G10, kapton, manganin, phosphor_bronze, saphire, teflon, dict=materials)
# Adding in dictionaries
materials = make_dict(aluminum_dict, balsa_dict, beechwood_dict, copper_dict, steel_dict, vespel_dict, dict=materials)
