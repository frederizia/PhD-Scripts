#! /usr/bin/env python
'''Hexagonal and square order parameter script by Maziar Fayas-Torshizi, amended by Frederike Jaeger'''

from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import argparse
from scipy.signal import argrelextrema

def GetArgs():
    ## Parse command line arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', '--width', nargs='+', required=False, default='liquid', action='store',
                       help='Width of channel')
    parser.add_argument('-rho', '--rho', required=False, default='RHO1', action='store',
                       help='Average channel density')
    parser.add_argument('-s', '--start', nargs='+',required=False, default='None', action='store',
                       help='Starting configurations')
    parser.add_argument('-r', '--density', nargs='+',required=False, default=None, action='store',
                       help='Density for constant channel')

    args = parser.parse_args()
    return args


def avg_density(DIR, z_name):

    # open density file
    file = open('Water_Graphene/{}/dens.{}'.format(DIR, z_name), 'r')
    _, __, data = (file.readline(), file.readline(),
                   file.readline())

    RHOeff = float(data.split()[1])
    print '\n\n#---------------------Average density = {} --------------------#'.format(RHOeff)
    
    return 


def peaks(DIR, z_name):

    # read data
    data = np.loadtxt('DATA/{}/rdf_2d_{}.dat'.format(DIR, z_name), delimiter=',')
    x, rdf = data[:,0], data[:,1]
    x, rdf = x[np.where(x>2)], rdf[np.where(x>2)]

    # find peaks
    idx = argrelextrema(rdf, np.greater)
    peaks = x[idx[0]]

    # idx = (np.gradient(np.sign(np.gradient(rdf))) < 0).nonzero()[0] 
    # peaks = x[idx]


    return peaks

def print_peaks(plist):
    print '1st peak:', plist[0]
    print '2nd peak:', plist[1], 'or', plist[1]/plist[0], 'P1'
    print '3rd peak:', plist[2], 'or', plist[2]/plist[0], 'P1'
    print '4th peak:', plist[3], 'or', plist[3]/plist[0], 'P1'
    print '5th peak:', plist[4], 'or', plist[4]/plist[0], 'P1'

    return


def main():


    args = GetArgs()

    DZ = args.width
    DIR = args.rho
    configs = args.start
    dens = args.density



    print '\n\nsqrt(2) = ', np.sqrt(2)
    print 'sqrt(3) = ', np.sqrt(3)
    print 'sqrt(4) = ', np.sqrt(4)
    print 'sqrt(5) = ', np.sqrt(5)
    print 'sqrt(6) = ', np.sqrt(6)
    print 'sqrt(7) = ', np.sqrt(7)
    print 'sqrt(8) = ', np.sqrt(8)


    for dz in DZ:
        for c in configs:
            if dens == None:
                z_name = 'spce_T298_z{}_eps1.0_{}'.format(dz, c)
                print '\n\n#---------------------Reading in dz = {}, DEN = {} --------------------#'.format(dz,c)
                peaks_list = peaks(DIR, z_name)
                print_peaks(peaks_list)
            else:
                for r in dens:
                    z_name = 'spce_T298_z{}_r{}_eps1.0_{}'.format(dz, r, c)
                    avg_density(DIR, z_name)
                    print '\n\n#---------------------Reading in dz = {}, DEN = {}, r = {} --------------------#'.format(dz,c,r)
                    peaks_list = peaks(DIR, z_name)
                    print_peaks(peaks_list)


    return

if __name__ == "__main__":
    sys.path.insert(0,'/home/fred/SCRIPTS/Python')
    from plotting_params import *
    sys.path.insert(0,'/home/fred/SCRIPTS/Other/Internet')
    from paircorrelation_ShockSolution import *
    main()