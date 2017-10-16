#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
import sys


def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--wavenumber', nargs=1, required=False, type=float, default='1', action='store',
                       help='wave number')
    parser.add_argument('-T', '--T', nargs=1, required=False, type=float, default=[300], action='store',
                       help='temperature')
    args = parser.parse_args()
    return args

def main():
    args    = GetArgs()
    k       = args.wavenumber[0] # 1/lambda
    T       = args.T[0]


    # Definitions
    fig_size_long = (14,5)
    fig_size_sq   = (9,7)

    EXT = 'PDF'
    ext = 'pdf'


    markers = ['o', 'D', 's', 'v', '^', 'd', '*']
    colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']

    kB  = 1.38e-23
    h   = 6.626e-34
    c   = 3e8
    R   = 8.314
    Cv  = 28.82 # J/molK
    Cp  = 37.136 # J/molK

    Ptau = 7e-6*101325 # sPa

    kSI = k*1e2

    print 'Temperature:', T, 'K.'

    print 'The wave number is', k, 'cm-1.'

    # vibrational temperature
    Evib = h*c*kSI
    Tvib = Evib/kB
    Erat = Evib/(kB*T)

    print 'The vibrational temperature is', Tvib, 'K.'
    #print Erat

    # probability
    f = np.exp(-Erat)

    print 'The excitation probability is', f*100, '%.'

    # heat capacity
    C_v = R*(Erat)**2*(np.exp(Erat)/(np.exp(Erat)-1)**2)

    print 'The heat capacity is', C_v, 'J/molK.'

    # bulk viscosity contribution
    fA   = (Cp/Cv)-1
    Kint = (C_v/Cv)*Ptau*1e5

    print 'f(A) is 1 or', fA
    print 'The bulk viscosity contribution is', Kint, 'or', fA*Kint, '10^(-5)Pas.'






    '''# Example figure

    fig1 = plt.figure(figsize=fig_size_sq)
    ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])
    ax1.errorbar(XDAT, YDAT, yerr=YERR, marker='o', linestyle='None')
    ax1.set_xlabel('no. of sheets')
    ax1.set_ylabel('$k$ (10$^{-17}$ m$^2$/sPa)')
    fig1.savefig('PLOTS/{}/PHYS_{}.{}'.format(EXT, name, ext))'''

    return

if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()