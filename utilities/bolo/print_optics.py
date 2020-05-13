"""
print_optics.py

prints out a table of the properties of the optical elements, after running optical_calcs().

"""

def print_optics(datasrc):
    # Print out a table of values for each source
    print(' \t   Source      ',end="")
    print(' eps    ',end="")
    print(' T_src  ',end="")
    print('eta_to_bolo',end="")
    print(' Q\t',end="")
    #print('    T_RJ',end="")
    print(' NEP_photon')
    #print('NET_photon_RJ NET_photon_cmb')
    for i in range(len(datasrc)):
        print("%20s" %datasrc[i]['name'],end="")
        print("%7.3f" %datasrc[i]['eps'],end="")
        print( "%8.1f" %datasrc[i]['T'],end="")
        print( "%9.3f" %datasrc[i]['eta_to_bolo'],end="")
        print( "%12.3e" %datasrc[i]['Q'],end="")
        #print( "%7.1f" %datasrc[i]['T_RJ'],end="")
        print( "%10.3e" %datasrc[i]['NEP_photon'])
        #print( "%12.3e" %datasrc[i]['NET_photon_RJ'],end="")
        #print( "%12.3e" %datasrc[i]['NET_photon_cmb'])


