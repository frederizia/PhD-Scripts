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
    parser.add_argument('-n', '--sheets', required=False, default='1', action='store',
                       help='Number of sheets')
    parser.add_argument('-dx', '--delx', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-o', '--off', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-f', '--force', required=False, nargs='+', default='Water', action='store')
    parser.add_argument('-png', '--png', nargs=1, type=str, required=False, default='n', action='store')
    args = parser.parse_args()
    return args

def main():
    args    = GetArgs()
    sheets  = args.sheets[0]
    delx    = args.delx
    offset  = args.off[0]
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



    fig4 = plt.figure(figsize=fig_size_sq)
    ax4  = fig4.add_axes([0.15,0.15,0.75,0.75])
    fig6 = plt.figure(figsize=fig_size_sq)
    ax6  = fig6.add_axes([0.15,0.15,0.75,0.75])

    markers = ['o', 'D', 's', 'v', '^', 'd', '*']


    name_plot = 'n{}_o{}_delx{}_F'.format(sheets,offset,delx)
    name_comb = 'n{}_o{}_delx'.format(sheets,offset)
    count = 0
    for dx in delx:
        vel_max = 0

        print '\n#----------------------- delx = {} ------------------------#\n'.format(dx)
        name_plot = 'n{}_o{}_delx{}_F'.format(sheets,offset,dx)
        name_comb += dx
        delP = []
        delRho = []
        Jz = []

        fig1 = plt.figure(figsize=fig_size_long)
        ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])
        fig2 = plt.figure(figsize=fig_size_long)
        ax2  = fig2.add_axes([0.15,0.15,0.75,0.75])
        fig3 = plt.figure(figsize=fig_size_sq)
        ax3  = fig3.add_axes([0.15,0.15,0.75,0.75])
        fig5 = plt.figure(figsize=fig_size_sq)
        ax5  = fig5.add_axes([0.15,0.15,0.75,0.75])

        for f in force:
            name_plot += f

            flow = FLOW(sheets,offset,dx,f)
            C, P, dP = flow.stress_prof()
            CZ, RHOZ, CX, RHOX, dRho, VELX, VMAX = flow.density_prof('None')
            JZ_MOL = flow.flux()

            delP.append(dP)
            delRho.append(dRho)
            Jz.append(JZ_MOL)

            if VMAX > vel_max:
                vel_max = VMAX

            ax1.plot(C,P, label='F={}'.format(f))
            ax2.plot(CZ,RHOZ, label='F={}'.format(f))
            ax3.plot(CX,RHOX, label='F={}'.format(f))
            #ax5.plot(CX,VELX, label='F={}'.format(f))


        ax4.plot(delRho, delP, marker=markers[count], label='$d_{\mathrm{slit}}$ = %s'%(dx))
        ax5.plot(delP, Jz, marker=markers[count], linestyle='None', label='$d_{\mathrm{slit}}$ = %s'%(dx))
        ax6.plot(delP, Jz, marker=markers[count], linestyle='None', label='$d_{\mathrm{slit}}$ = %s'%(dx))


        ax1.set_xlabel('z (\AA)')
        ax1.set_ylabel('P (MPa)')
        ax1.set_ylim(-100,400)
        ax1.set_xlim(0,100)
        ax1.legend()

        ax2.set_xlabel('z (\AA)')
        ax2.set_ylabel('$\\rho$ (g/cm$^3$)')
        #ax2.set_ylim(-50,300)
        ax2.set_xlim(0,100)
        ax2.legend()

        ax3.set_xlabel('x (\AA)')
        ax3.set_ylabel('$\\rho$ (g/cm$^3$)')
        #ax3.set_ylim(-50,300)
        #ax3.set_xlim(0,100)
        ax3.legend()

        ax4.set_xlabel('$\Delta\\rho$ (g/cm$^3$)')
        ax4.set_ylabel('$\Delta$P (MPa)')
        #ax3.set_ylim(-50,300)
        #ax3.set_xlim(0,100)
        ax4.legend()

        ax5.set_xlabel('$\Delta$P (MPa)')
        ax5.set_ylabel('$J_z$ (10$^3$ mol/m$^2$s)')
        #ax5.set_ylim(0,vel_max)
        ax5.set_xlim(0, 300)
        ax5.legend()

        ax6.set_xlabel('$\Delta$P (MPa)')
        ax6.set_ylabel('$J_z$ (10$^3$ mol/m$^2$s)')
        #ax5.set_ylim(0,vel_max)
        ax6.set_xlim(0, 300)
        ax6.legend()


        fig1.savefig('PLOTS/TESTS/{}/press_{}.{}'.format(EXT, name_plot, ext))
        fig2.savefig('PLOTS/TESTS/{}/rhoz_{}.{}'.format(EXT, name_plot, ext))
        fig3.savefig('PLOTS/TESTS/{}/rhox_{}.{}'.format(EXT, name_plot, ext))
        fig5.savefig('PLOTS/TESTS/{}/Jz_{}.{}'.format(EXT, name_plot, ext))

        count += 1
    fig4.savefig('PLOTS/TESTS/{}/drop_{}.{}'.format(EXT, name_comb, ext))
    fig6.savefig('PLOTS/TESTS/{}/Jz_{}.{}'.format(EXT, name_comb, ext))

    #plt.show()


    

    return

if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()