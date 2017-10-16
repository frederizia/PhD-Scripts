#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
from go_tools import *
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

    fig1 = plt.figure(figsize=fig_size_sq)
    ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])
    fig2 = plt.figure(figsize=fig_size_sq)
    ax2  = fig2.add_axes([0.15,0.15,0.75,0.75])


    for n in sheets:
        for o in offset:
            for dx in delx:
                for f in force:
                    x,y,VX,VXi,VY,VYi,V,Vi = vel_prof(n,o,dx,f)
                    #X, Y = np.meshgrid(x,y)
                    Y,X = np.meshgrid(y,x)
                    #print X, Y, len(x), len(y), np.transpose(VX).shape
                    #M = np.hypot(U, V)
                    Qi = ax1.quiver(Y, X, VXi, VYi, Vi, units='x', pivot='tail', scale_units='width',scale=0.01)
                    qki = ax1.quiverkey(Qi, 0.8, 0.95, 0.001, '$10^{-3}$ \AA/ps', labelpos='E',
                                       coordinates='figure')
                    ax1.scatter(Y, X, color='k', s=4)
                    #ax1.plot([np.min(y),np.max(y)],[5,5], c='r', linestyle='--')
                    #ax1.plot([np.min(y),np.max(y)],[10,10], c='r', linestyle='--')
                    #ax1.plot([np.min(y),np.max(y)],[15.6,15.6], c='r', linestyle='--')
                    #ax1.plot([np.min(y),np.max(y)],[22,22], c='r', linestyle='--')
                    ax1.set_ylabel('$x$ (\AA)')
                    ax1.set_xlabel('$y$ (\AA)')


                    Q = ax2.quiver(Y, X, VX, VY, V, units='x', pivot='tail', scale_units='width',scale=0.01)#, width=0.022,scale=1 / 0.15)
                    qk = ax2.quiverkey(Q, 0.8, 0.95, 0.001, '$10^{-3}$ \AA/ps', labelpos='E',
                                       coordinates='figure')
                    ax2.scatter(Y, X, color='k', s=4)
                    #ax2.plot([np.min(y),np.max(y)],[5,5], c='r', linestyle='--')
                    #ax2.plot([np.min(y),np.max(y)],[10,10], c='r', linestyle='--')
                    #ax2.plot([np.min(y),np.max(y)],[15.6,15.6], c='r', linestyle='--')
                    #ax2.plot([np.min(y),np.max(y)],[22,22], c='r', linestyle='--')
                    ax2.set_ylabel('$x$ (\AA)')
                    ax2.set_xlabel('$y$ (\AA)')

                    plt.show()



                fig1.savefig('PLOTS/{}/vxy_init_{}_n{}_o{}_delx{}_f{}.{}'.format(EXT,'CH',n,o,dx,f,ext))
                fig2.savefig('PLOTS/{}/vxy_ave_{}_n{}_o{}_delx{}_f{}.{}'.format(EXT,'CH',n,o,dx,f,ext))


    return

if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()