#!/usr/bin/python
# ConvertEOSTable.py

""" 
Convert EOS tables from Lorene to SGRID format
"""

__author__ = 'sbernuzzi'
__date__ = '01/2014 02/2020'

import sys
import numpy as np
from scipy import interpolate

def unit_conv_fact( ):
    """
    Return conversion factors
    """
    return dict( Length_cgs = 1.476625038500000e+05,
                 Time_cgs = 4.925490949141889e-06,
                 Mass_cgs = 1.988920000000000e+33,
                 Energy_cgs = 1.787552150093231e+54,
                 Surface_cgs = 2.180421504325126e+10,
                 Volume_cgs = 3.219664987770316e+15,
                 Press_cgs = 5.551981826938919e+38,
                 Mdens_cgs = 6.177412890952259e+17,
                 Edens_cgs = 5.551981826938919e+38,
                 Amom_cgs = 8.804571936403334e+48,
                 Length_fm = 1.476999442301651e+18,
                 Length_km = 1.476625038500000e+00,
                 Density_kgm3 = 6.177412890952258e+20,
                 Energy_MeV = 1.115702344637173e+60,
                 Volume_fm3 = 3.219664987770317e+54,
                 fm_to_cm = 1e-13,
                 barmass_g = 1.66e-24);

def read_lorene_eosfile_header( fname ):
    """ 
    Read EOS file header in Lorene format

    header:

    |#
    |#
    |#
    |#
    |#
    |{} <-- Number of lines
    |#
    |#       n_B [fm^{-3}]  rho [g/cm^3]   p [dyn/cm^2]
    |# 
    |...

    """
    # read header
    with open(fname, 'r') as fh:
        header = []
        for line in fh:
            lstart = line.lstrip()
            # do nothing if it's not the header
            if lstart[0] != '#':
                continue
            # skip empty comments
            if lstart == '#\n':
                continue
            # skip cgs stuff
            if lstart.find('dyn') >=0 or lstart.find('g/cm')>=0:
                continue
            header.append(lstart)
    return header

def read_lorene_eosfile( fname, skipn=9 ):
    """ 
    Read EOS file in Lorene format

    Skip the header:

    |#
    |#
    |#
    |#
    |#
    |{} <-- Number of lines
    |#
    |#       n_B [fm^{-3}]  rho [g/cm^3]   p [dyn/cm^2]
    |# 
    |...

    """
    return np.loadtxt(fname, skiprows=skipn)

def write_lorene_eosfile( fname, data ):
    """
    Write EOS file in Lorene format 
    """
    n = len(data);
    # header
    with open(fname, 'w') as fh:
        fh.write('#\n#\n#\n');
        fh.write('{} <-- Number of lines\n'.format(n));
        fh.write('#\n');
        fh.write('#        n_B [fm^{-3}]  rho [g/cm^3]   p [dyn/cm^2]\n');
        fh.write('#\n');
    # data
    with open(fname,'a') as fh:
        np.savetxt(fh,data, fmt='%d %.16e %.16e %.16e')
    return

def convert_lorene2sgrid( datain ):
    """
    Convert Lorene data to SGRID
    """
    units = unit_conv_fact()
    dataout = np.zeros_like(datain);
    # rho0 = mB * nB
    dataout[:,0] = ( units["barmass_g"] * datain[:,1] / (units["fm_to_cm"]**3) ) / units["Mdens_cgs"]; 
    # epsl  = ene/rho0 -1 
    dataout[:,1] = ( datain[:,2]/units["Mdens_cgs"] ) / dataout[:,0] - 1.  # (*)
    # pres
    dataout[:,2] = datain[:,3] / units["Press_cgs"]
    # reshape and return positive (*)
    dataout = dataout[:,0:3]
    return dataout[dataout[...,1]>0]

if __name__ == '__main__':

    # read input files in Lorene format and convert to SGRID

    nf = len(sys.argv)-1
    print ("Files to process: %d " % nf)
    print ("Filenames: %s " % str(sys.argv[1:]))
    if nf<1:
        print ("Usage:\n "+sys.argv[0]+" <input files>")
        sys.exit();

    for f in sys.argv[1:]:
        header  = read_lorene_eosfile_header(f)
        print(header)
        data = read_lorene_eosfile(f)
        convdata = convert_lorene2sgrid(data)
        name = f.rsplit( ".", 1 )[ 0 ]+".txt"
        with open(name, 'wb') as fh:
            fh.write('# ')
            fh.write(f)
            fh.write(' converted to units with G=c=Msun=1\n')
            for line in header:
                fh.write(line)
            fh.write('#\n')
            fh.write('# rho0=nB*mB           epsl=(rho-rho0)/rho0    P\n')
            fh.write('#\n')
            np.savetxt( fh, convdata, fmt='%.16e %.16e %.16e')
        #print(data[0,0],data[-1,0])
        print('n={} col={}'.format(np.shape(data)[0],np.shape(data)[1]))
        print('{} => {}'.format( f, name ))
