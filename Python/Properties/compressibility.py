#!/usr/bin/python

from __future__ import division
import numpy as np
from class_bulk_props import *
import matplotlib.pyplot as plt
import argparse
import sys

def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--method', type=int, required=False, default=0, action='store',
                       help='Fluctuation or slope')
    parser.add_argument('-p', '--pressure', required=False, nargs='+', default=[0,500,1000,2000], action='store',
                       help='Pressures')
    args = parser.parse_args()
    return args

def poly1(xlist,m,b):
    result = []
    for x in xlist:
        result.append(m*x+b)
    return result

def main():
    args    = GetArgs()
    method  = args.method
    press   = args.pressure

    if method == 1:
        print 'Using volume fluctuations at given pressure.'
        
        kB = 1.380648813e-23
        for p in press:
            print '+++++++++++++++++ p = {} ++++++++++++++++++'.format(p)
            T, V, varV, V2 = press_fluct(p)
            comp = varV/(kB*T*V)
            #comp = (V2-V**2)/(kB*T*V)

            print 'The compressibility is:', comp, 'bar^-1'
            print 'The bulk modulus is:', 1e-5*(1/comp)/1e9, 'GPa'

    elif method == 2:
        print 'Using fit at different pressures'
        press_vals = []
        vol_vals = []
        for p in press:
            print 'Reading in p={}....'.format(p)
            vol_val, press_val = mean_vals(p) 
            vol_vals.append(vol_val)
            press_vals.append(press_val)

        m, v0 = np.polyfit(np.array(press_vals), np.array(vol_vals),1)

        comp = (-1/v0)*m
        print 'The compressibility is:', comp/1e5, 'bar^-1'
        print 'The bulk modulus is:', (1/comp)/1e9, 'GPa'

        fitted = poly1(press_vals, m, v0)
        plt.figure()
        plt.plot(press_vals, vol_vals, linestyle='None', marker='D')
        plt.plot(press_vals, fitted)
        plt.show()


    else:
        print 'Method unknown. Try again.'
        sys.exit(1)

	return


if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()
