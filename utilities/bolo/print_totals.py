
def print_totals(detector):
    # print out the totals
    print('')
    print('nu_min = {0:6.2f} GHz'.format(detector['nuGHZ'].min()))
    print('nu_mean = {0:6.2f} GHz'.format(detector['nuGHZ'].mean()))
    print('nu_max = {0:6.2f} GHz'.format(detector['nuGHZ'].max()))
    print('')
    print('Qtot = {0:6.2e} W'.format(detector['Qtot']))
    print('W = {0:6.2e} W'.format(detector['W']))
    print('T_tes = {0:6.3f} K'.format(detector['T_bolo']))
    print('T_base = {0:6.3f} K'.format(detector['T_base']))
    print('')
    print('Total: ----------------------------------------------')
    print(' P_optical = {0:10.3e}'.format(detector['Qtot']))
    print(' P_total = {0:10.3e}'.format(detector['W']))
    print(' T_RJ_tot = {0:10.1e}'.format(detector['T_RJ_tot']))
    print(' dPdT_RJ = {0:10.3e}'.format(detector['dPdT_RJ']))
    print(' dPdT_cmb = {0:10.3e}'.format(detector['dPdT_cmb']))
    print('')
    #
    print(' I_bolo = {0:3.2e} Amps'.format(detector['I_0']))
    print(' R_bolo = {0:3.2e} Ohms'.format(detector['R_bolo']))
    print(' V_bolo = {0:3.2e} Volts'.format(detector['V_0']))
    print(' S_DC = {0:3.2e} Amps/Watt'.format(detector['Sdc']))
    print('')
    #
    print(' NEP_photon_total = {0:10.3e} '.format(detector['Photon']['NEP']))
    print(' NEP_phonon_total = {0:10.3e} '.format(detector['Phonon']['NEP']))
    print(' NEP_Johnson_total = {0:10.3e} '.format(detector['Johnson']['NEP']))
    print(' NEP_Readout_total = {0:10.3e} '.format(detector['Readout']['NEP']))
    print(' NEP_total = {0:10.3e} '.format(detector['NoiseTotal']['NEP']))
    #
    print(' NEI_photon_total = {0:10.3e} '.format(detector['Photon']['NEI']))
    print(' NEI_phonon_total = {0:10.3e} '.format(detector['Phonon']['NEI']))
    print(' NEI_Johnson_total = {0:10.3e} '.format(detector['Johnson']['NEI']))
    print(' NEI_Readout_total = {0:10.3e} '.format(detector['Readout']['NEI']))
    #
    print(' NET_photon_total_RJ = {0:10.3e} '.format(detector['Photon']['NET_RJ']))
    print(' NET_phonon_total_RJ = {0:10.3e} '.format(detector['Phonon']['NET_RJ']))
    print(' NET_Johnson_total_RJ = {0:10.3e} '.format(detector['Johnson']['NET_RJ']))
    print(' NET_Readout_total_RJ = {0:10.3e} '.format(detector['Readout']['NET_RJ']))
    #
    print(' NET_total_RJ = {0:10.3e} K/sqrt(Hz)'.format(detector['NoiseTotal']['NET_RJ']))
    print(' NET_total_cmb = {0:10.3e} K/sqrt(Hz)'.format(detector['NoiseTotal']['NET_CMB']))
    print('')
    print('-----------------------------------------------------')

