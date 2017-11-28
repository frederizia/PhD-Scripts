#! /usr/bin/env python
# Quick diffusion coefficient calculation given VACF

import argparse
import matplotlib.pyplot as plt
import matplotlib
import sys
import numpy as np
from scipy.integrate import simps

def GetArgs():
    ## Parse command line arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", type=str, 
                        help = "Name of LAMMPS dump file to be read", 
                        required = True)
    parser.add_argument('-u', '--units', required=False, type=str, action='store', default=[0,100],
                        help='metal or real')
    args = parser.parse_args()
    return args

def diffusion(vacf, DT, t_conv, s_conv):

    int_vacf = simps(vacf, dx=DT)

    print 'The integral of the VACF is', int_vacf

    diff = int_vacf/3
    diff_2d = int_vacf/2

    conv = s_conv**2/t_conv

    diff_si = diff*conv
    diff_2d_si = diff_2d*conv

    print 'The self-diffusion coefficient in 3D is', diff_si
    print 'The self-diffusion coefficient in 2D is', diff_2d_si

    return diff

def main():


    args = GetArgs()

    fname = args.filename
    units = args.units

    if units == 'metal':
        time_conv = 1e-12
        space_conv = 1e-10
        dt = 0.0005
    elif units == 'metal2':
        time_conv = 1e-12
        space_conv = 1e-10
        dt = 0.001
    elif units == 'real':
        time_conv = 1e-15
        space_conv = 1e-10
        dt = 2
    elif units == 'lj':
        time_conv = 1
        space_conv = 1
        dt = 0.001
    else:
        print 'Invalid unit.'
        sys.exit(1)

    print 'Calculating the diffusion coefficient for:', fname

    try:
        C_vv_array = np.loadtxt("{}".format(fname))
        times =  C_vv_array[:,0]
        C_vv_ave = C_vv_array[:,1]
        diffusion(C_vv_ave, dt, time_conv, space_conv)
    except IOError:
        print 'File does not exist.'


if __name__ == "__main__":
    sys.path.insert(0,'/home/fred/SCRIPTS/Python')
    main()