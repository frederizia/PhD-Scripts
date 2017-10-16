#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
from couette_tools import *
import sys


def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vel', nargs='+', required=False, default='1', action='store',
                       help='Wall velocity')
    parser.add_argument('-png', '--png', nargs=1, type=str, required=False, default='n', action='store')
    args = parser.parse_args()
    return args

def main():
    args    = GetArgs()
    VEL  = args.vel
    png     = args.png[0]

    folder = 'eps1.0/s_21.44/WT'

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
    ls = ['-', '--', ':', ':']

    # Figure initialisation

    # velocity
    fig1 = plt.figure(figsize=fig_size_sq)

    # velcoity gradient
    fig2 = plt.figure(figsize=fig_size_sq)

    # density
    fig3 = plt.figure(figsize=fig_size_sq)

    # isotropic stress
    fig4 = plt.figure(figsize=fig_size_sq)

    # anisotropic stress
    fig5 = plt.figure(figsize=fig_size_sq)

    # shear viscosity
    fig6 = plt.figure(figsize=fig_size_sq)

    for v in VEL:
        f = 'nemd_H21.44_v{}'.format(v)
        try:
            print 'Reading in {}.{}'.format('vel',f)
            Z, RHO, VX = read_velprof(folder, 'vel',f)
            MID, LEFT, RIGHT = mid_point(Z,VX,0.3)
            Z = Z-Z[MID]
            Z_left = Z[LEFT]
            Z_right = Z[RIGHT]
            #Z, RHO, VX = Z[LEFT:RIGHT], RHO[LEFT:RIGHT], VX[LEFT:RIGHT]
            dVX = derivative(Z, VX)
            dZ = Z[:-1]+(Z[1]-Z[0])/2
            label = '$v_x$ = '.format(v)
            print dVX.shape, VX.shape

            fig1.clear()
            ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])
            ax1.plot(Z, VX, label=label, marker='o', fillstyle='none')
            ax1.set_ylabel('$v_x$ (\AA$^2$/ps)')
            ax1.set_xlabel('$z-z_{\mathrm{mid}}$ (\AA)')
            #ax2a.set_ylim(Z_left-0.5, Z_right+0.5)#limits[dz][0], limits[dz][1])
            #ax3.set_ylim(-8, 8)
            #ax3.set_xlim(0, 7)
            fig1.savefig('PLOTS/{}/velprof_{}.{}'.format(EXT,f,ext))

            fig2.clear()
            ax2  = fig2.add_axes([0.15,0.15,0.75,0.75])
            ax2.plot(dZ, dVX, label=label, marker='o', fillstyle='none')
            ax2.set_ylabel('$dv_x/dz$ (\AA/ps)')
            ax2.set_xlabel('$z-z_{\mathrm{mid}}$ (\AA)')
            #ax2a.set_ylim(Z_left-0.5, Z_right+0.5)#limits[dz][0], limits[dz][1])
            #ax3.set_ylim(-8, 8)
            #ax3.set_xlim(0, 7)
            fig2.savefig('PLOTS/{}/velder_{}.{}'.format(EXT,f,ext))
            

            fig3.clear()
            ax3  = fig3.add_axes([0.15,0.15,0.75,0.75])
            ax3.plot(Z, RHO, label=label)
            ax3.set_ylabel('$\\rho$ (g/cm$^3$)')
            ax3.set_xlabel('$z-z_{\mathrm{mid}}$ (\AA)')
            #ax2a.set_ylim(Z_left-0.5, Z_right+0.5)#limits[dz][0], limits[dz][1])
            #ax3.set_ylim(-8, 8)
            #ax3.set_xlim(0, 7)
            fig3.savefig('PLOTS/{}/densprof_{}.{}'.format(EXT,f,ext))


        except IOError:
            print 'File {}.{} does not exist.'.format('vel',f)

        # pressure profile
        try:
            print 'Reading in stress.{}'.format(f)
            Z_P, P, Ptot, Pxy, Pxz, Pyz, Pxx, Pyy, Pzz, delz = stress_prof(folder, 'stress',f)
            MID_P, LEFT_P, RIGHT_P = mid_point(Z_P,P, 0.2)
            Z_P = Z_P-Z_P[MID_P]
            #Z_P, Pxx, Pyy, Pzz, Pxy, Pxz, Pyz = \
            #        Z_P[LEFT_P:RIGHT_P], Pxx[LEFT_P:RIGHT_P], \
            #        Pyy[LEFT_P:RIGHT_P], Pzz[LEFT_P:RIGHT_P], \
            #        Pxy[LEFT_P:RIGHT_P], Pxz[LEFT_P:RIGHT_P], Pyz[LEFT_P:RIGHT_P]

            
            fig4.clear()
            ax4  = fig4.add_axes([0.15,0.15,0.75,0.75])
            ax4.plot(Z_P, Pxx, label='Pxx', linestyle=ls[0])
            ax4.plot(Z_P, Pyy, label='Pyy', linestyle=ls[1])
            ax4.plot(Z_P, Pzz, label='Pzz', linestyle=ls[2])
            ax4.set_ylabel('$P_{\\alpha\\alpha}$ (MPa)')
            ax4.set_xlabel('$z-z_{\mathrm{mid}}$ (\AA)')
            ax4.legend()
            fig4.savefig('PLOTS/{}/Piso_{}.{}'.format(EXT,f,ext))

            fig5.clear()
            ax5  = fig5.add_axes([0.15,0.15,0.75,0.75])
            ax5.plot(Z_P, Pxz, label='Pxz', linestyle=ls[0])
            ax5.plot(Z_P, Pxy, label='Pxy', linestyle=ls[1])
            ax5.plot(Z_P, Pyz, label='Pyz', linestyle=ls[2])
            ax5.set_ylabel('$P_{\\alpha\\beta}$ (MPa)')
            ax5.set_xlabel('$z-z_{\mathrm{mid}}$ (\AA)')
            ax5.legend()
            fig5.savefig('PLOTS/{}/Paniso_{}.{}'.format(EXT,f,ext))

            print Z_P.shape, Pxz.shape, dVX.shape
            fig6.clear()
            ax6  = fig6.add_axes([0.15,0.15,0.75,0.75])
            ax6.plot(Z_P[:-1], -Pxz[:-1]/dVX, label='Pxz', marker='o', fillstyle='none', linestyle='--')
            ax6.set_ylabel('$\eta$')
            ax6.set_xlabel('$z-z_{\mathrm{mid}}$ (\AA)')
            fig6.savefig('PLOTS/{}/shear_{}.{}'.format(EXT,f,ext))



        except IOError:
            print 'File stress.{} does not exist.'.format(f)


    return

if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()