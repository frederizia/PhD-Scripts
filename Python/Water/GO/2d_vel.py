#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
import sys


def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--sheets', nargs='+', required=False, default='1', action='store',
                       help='Number of sheets')
    parser.add_argument('-dx', '--delx', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-o', '--off', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-f', '--force', required=False, nargs='+', default='Water', action='store')
    parser.add_argument('-png', '--png', nargs=1, type=str, required=False, default='n', action='store')
    args = parser.parse_args()
    return args

def main():
    args    = GetArgs()
    sheets  = args.sheets
    delx    = args.delx
    offset  = args.off
    force   = args.force
    png     = args.png[0]

    # Definitions
    fig_size_long = (14,5)
    fig_size_sq   = (9,7)
    if png == 'n':
        EXT = 'PDF'
        ext = 'pdf'
    else:
        EXT = 'PNG'
        ext = 'png'

    markers = ['o', 'D', 's', 'v', '^', 'd', '*']
    colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']



    # Example figure

    fig1 = plt.figure(figsize=fig_size_sq)
    ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])
    ax1.errorbar(XDAT, YDAT, yerr=YERR, marker='o', linestyle='None')
    ax1.set_xlabel('no. of sheets')
    ax1.set_ylabel('$k$ (10$^{-17}$ m$^2$/sPa)')
    fig1.savefig('PLOTS/{}/PHYS_{}.{}'.format(EXT, name, ext))

    return

if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()