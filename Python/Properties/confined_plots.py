#! /usr/bin/env python
# Quick diffusion coefficient calculation given VACF

import argparse
import matplotlib.pyplot as plt
import matplotlib
import sys
import numpy as np
from scipy.integrate import simps
from confined_properties import *

def GetArgs():
    ## Parse command line arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', '--width', nargs='+', required=False, default='liquid', action='store',
                       help='Width of channel')
    parser.add_argument('-r', '--region', nargs='+', required=False, default='liquid', action='store',
                       help='Which region? wall, intermediate, bulk')
    parser.add_argument('-u', '--units', required=False, type=str, action='store', default=[0,100],
                        help='metal or real')
    parser.add_argument('-p', '--pre', nargs='+', required=False, type=str, action='store', default='C_vv',
                        help='C_vv or C_vv_y')
    parser.add_argument('-n', '--name', required=False, default='None', action='store',
                       help='Plot name')
    args = parser.parse_args()
    return args

def main():


    args = GetArgs()

    units = args.units
    dZ  = args.width
    region  = args.region
    plot_name = args.name
    pre = args.pre

    if units == 'metal':
        time_conv = 1e-12
        space_conv = 1e-10
        dt = 0.0005
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


    # Definitions
    fig_size_long = (14,5)
    fig_size_sq   = (9,7)
    fig_size_sq2 = (7,9)

    markers = ['o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P','X']
    colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']
    LABEL = {'C_vv': '$xyz$', 'C_vv_y': '$xy$'}


    ls = ['-', '--', ':', ':']

    # ------------------- Initialise figures -----------------------

    # diffusion
    fig1 = plt.figure(figsize=fig_size_sq)
    ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])

    # diff and rho
    fig3 = plt.figure(figsize=fig_size_sq)
    ax3  = fig3.add_axes([0.15,0.15,0.75,0.75])

    # density profile
    fig2 = plt.figure(figsize=fig_size_long)
    ax2  = fig2.add_axes([0.15,0.15,0.75,0.75])

    # density profile sideways
    fig2a = plt.figure(figsize=fig_size_sq2)

    # stress profile
    fig4 = plt.figure(figsize=fig_size_sq)
    ax4  = fig4.add_axes([0.15,0.15,0.75,0.75])

    # stress and density profile
    fig5 = plt.figure(figsize=fig_size_sq)

    # eta
    fig6 = plt.figure(figsize=fig_size_sq)
    


    # Properties independent of dimensionality
    print '\n#-----------------------Analysing density, pressure and viscosities------------------\n#'
    RHOdict = {}
    for dz in dZ:
        count = 1
        count_list, diff_list = [], []
        f = 'spce_T298_z{}_eps1.0'.format(dz)
        print '\n#-----------------------dz = {}------------------\n#'.format(dz)

        # density profile
        try:
            print 'Reading in densprof.{}'.format(f)
            Z, RHO = read_densprof(f)
            MID, LEFT, RIGHT = mid_point(Z,RHO)
            # store data for use in diffusion plots
            RHOdict[dz] = (RHO, LEFT)
            Z = Z-Z[MID]
            Z_left = Z[LEFT]
            Z_right = Z[RIGHT]
            label = '$\Delta z = {}$ \AA'.format(dz)
            ax2.plot(Z, RHO, label=label)

            fig2a.clear()
            ax2a  = fig2a.add_axes([0.15,0.15,0.75,0.75])
            ax2a.plot(RHO,Z)
            ax2a.set_xlabel('$\\rho$ (g/cm$^3$)')
            ax2a.set_ylabel('$z-z_{\mathrm{mid}}$')
            #ax2a.set_ylim(Z_left-0.5, Z_right+0.5)#limits[dz][0], limits[dz][1])
            ax2a.set_ylim(-8, 8)
            ax2a.set_xlim(0, 7)
            fig2a.savefig('PLOTS_C/densprof_sideways_{}.pdf'.format(dz))


        except IOError:
            print 'File densprof.{} does not exist.'.format(f)

        # pressure profile
        try:
            print 'Reading in stress.{}'.format(f)
            Z_P, P, Ptot, Pxy, Pxz, Pyz, delz = stress_prof(f, 'None')
            MID_P, LEFT_P, RIGHT_P = mid_point(Z_P,P)
            Z_P = Z_P-Z_P[MID_P]
            label = '$\Delta z = {}$ \AA'.format(dz)
            ax4.plot(Z_P, P, label=label)

        except IOError:
            print 'File stress.{}_1 does not exist.'.format(f)

        try:
            fig5.clear()
            ax5  = fig5.add_axes([0.15,0.15,0.75,0.75])
            ax5.plot(Z_P, P, c=colours[0])
            ax5.set_ylabel('P (MPa)')
            ax5.set_ylim(0,np.max(P)+0.05*np.max(P))
            ax5.yaxis.label.set_color(colours[0])
            for tl in ax5.get_yticklabels():   
                tl.set_color(colours[0])

            ax5b = ax5.twinx()
            ax5b.plot(Z, RHO, c=colours[6], linestyle='dashed')
            ax5b.set_ylim(0,np.max(RHO)+0.05*np.max(RHO))
            fig5.text(0.975, 0.5, '$\\rho_{\mathrm{wall}}$ (g/cm$^3$)', color=colours[6], ha='center', va='center', fontsize=24,rotation=270)
            for tl in ax5b.get_yticklabels():   
                tl.set_color(colours[6])
            ax5.set_xlabel('$z-z_{\mathrm{mid}}$')
            #ax5.set_xlim(-17.4, 17.4)
            fig5.savefig('PLOTS_C/stressdens_{}.pdf'.format(dz))
        except:
            pass
        # viscosity profile
        '''try:
            print 'Reading in stresseta.{}'.format(f)
            Z_P, P, Ptot, Pxy, Pxz, Pyz, delz = stress_prof(f, 'eta')
            MID_P, LEFT_P, RIGHT_P = mid_point(Z_P,P)
            Z_P = Z_P-Z_P[MID_P]

            # Calculate viscosity
            Pxy_list, Pxz_list, Pyz_list = [], [], []
            print Pxy.shape
            for i in range(Pxy.shape[1]):
                Pxy_list.append(eta(f,Pxy[:-1,i],delz))
                Pxz_list.append(eta(f,Pxz[:-1,i],delz))
                Pyz_list.append(eta(f,Pyz[:-1,i],delz))
            Paniso = np.mean(np.array([Pxy_list, Pxz_list, Pyz_list]), axis=0)
            #ax6.plot(Z_P, Paniso)
            fig6.clear()
            ax6  = fig6.add_axes([0.15,0.15,0.75,0.75])
            
            ax6.plot(Z_P, Pxz_list, linestyle='--', marker=markers[1], label='$xz$')
            ax6.plot(Z_P, Pyz_list, linestyle=':', marker=markers[2], label='$yz$')
            ax6.plot(Z_P, Pxy_list, linestyle='-', marker=markers[0], label='$xy$')
            ax6.set_ylabel('$\eta$')
            ax6.set_xlabel('$z-z_{\mathrm{mid}}$')
            #ax2.set_xlim(Z_left, Z_right)
            #ax6.set_xlim(-17.4, 17.4)
            ax6.legend()
            #ax4.set_ylim(0, 4)
            fig6.savefig('PLOTS_C/SV_z{}.pdf'.format(dz))

        except IOError:
            print 'File stresseta.{}_1 does not exist.'.format(f)'''


    # Diffusion related properties
    print '\n#-----------------------Analysing diffusion coefficient------------------\n#'
    count_p = 0
    diff_bulk = 2.86624087647e-09
    pre_name = ''
    for p in pre:
        pre_name += LABEL[p].split('$')[1]+'_'
        z_list = []
        d_list = []
        rho_list = []
    

        z_name = 'z_'
        for dz in dZ:
            count = 1
            count_list, diff_list = [], []
            f = 'spce_T298_z{}_eps1.0'.format(dz)

            '\n#-----------------------dz = {}------------------\n#'.format(dz)

            z_name += dz

            for r in region:

                # diffusion
                try:
                    print 'Reading in {}_{}.dat'.format(f,r)
                    C_vv_array = np.loadtxt("Water_Graphene/{}_{}_{}.dat".format(p,f,r))
                    times =  C_vv_array[:,0]
                    C_vv_ave = C_vv_array[:,1]
                    diff_tmp = diffusion(C_vv_ave, dt, time_conv, space_conv, p)
                    count_list.append(count)
                    diff_list.append(diff_tmp/diff_bulk)
                    count += 1
                    if r == 'wall' and p == 'C_vv_y':
                        print dz
                        RHOdat = RHOdict[dz]
                        RHO_max = rho_wall_max(RHOdat[0], RHOdat[1])
                        z_list.append(dz)
                        d_list.append(diff_tmp/diff_bulk)
                        rho_list.append(RHO_max)
                        print diff_tmp/diff_bulk, RHO_max



                except IOError:
                    print 'File {}_{}.dat does not exist.'.format(f,r)

            ax1.plot(count_list, diff_list, marker=markers[count_p], label=LABEL[p], markersize=12, linestyle='dashed')

            count_p += 1

    print '\n#-----------------------Final plotting------------------\n#'

    #ax1.set_xlabel('region')
    ax1.set_ylabel('$D/D_{\mathrm{bulk}}$')
    ax1.set_xlim(0.7,3.3)
    ax1.set_xticks(count_list)
    ax1.set_xticklabels(['wall', 'intermediate', 'bulk'])
    ax1.legend()

    ax2.set_ylabel('$\\rho$ (g/cm$^3$)')
    ax2.set_xlabel('$z-z_{\mathrm{mid}}$')
    #ax2.set_xlim(Z_left, Z_right)
    ax2.set_xlim(-18, 18)
    ax2.set_ylim(0, 4)
    #ax2.legend()

    

    ax3.plot(z_list, d_list, c=colours[0], linestyle='dashed', marker='D')
    ax3.set_ylabel('$D/D_{\mathrm{bulk}}$')
    ax3.yaxis.label.set_color(colours[0])
    for tl in ax3.get_yticklabels():   
        tl.set_color(colours[0])
    ax3b = ax3.twinx()
    ax3b.plot(z_list, rho_list, c=colours[6], linestyle='dashed', marker='o')
    #ax3b.set_ylabel('$\\rho_{\mathrm{wall}}$ (g/cm$^3$)', rotation=270)
    #ax3b.yaxis.label.set_color(colours[6])
    fig3.text(0.975, 0.5, '$\\rho_{\mathrm{wall}}$ (g/cm$^3$)', color=colours[6], ha='center', va='center', fontsize=24,rotation=270)
    for tl in ax3b.get_yticklabels():   
        tl.set_color(colours[6])
    ax3.set_xlabel('$\Delta z$ (\AA)')
    #ax3.set_xlabel('$z-z_{\mathrm{mid}}$')
    #ax3.legend()

    ax4.set_ylabel('P (MPa)')
    ax4.set_xlabel('$z-z_{\mathrm{mid}}$')
    #ax2.set_xlim(Z_left, Z_right)
    ax4.set_xlim(-18, 18)
    #ax4.set_ylim(0, 4)

    if len(dZ) == 1:
        try:
            ax5.plot(Z_P, P, c=colours[0])
            ax5.set_ylabel('P (MPa)')
            ax5.yaxis.label.set_color(colours[0])
            for tl in ax5.get_yticklabels():   
                tl.set_color(colours[0])
            ax5b = ax5.twinx()
            ax5b.plot(Z, RHO, c=colours[6], linestyle='dashed')
            fig5.text(0.975, 0.5, '$\\rho_{\mathrm{wall}}$ (g/cm$^3$)', color=colours[6], ha='center', va='center', fontsize=24,rotation=270)
            for tl in ax5b.get_yticklabels():   
                tl.set_color(colours[6])
            ax5.set_xlabel('$z-z_{\mathrm{mid}}$')
            ax5.set_xlim(-17.4, 17.4)
            fig5.savefig('PLOTS_C/stressdens_{}.pdf'.format(z_name))
        except:
            pass

    

    fig2.savefig('PLOTS_C/densprof_{}.pdf'.format(z_name))
    fig4.savefig('PLOTS_C/stressprof_{}.pdf'.format(z_name))
    fig3.savefig('PLOTS_C/diffdens_wall_{}.pdf'.format(z_name))
    fig1.savefig('PLOTS_C/diffusion_{}{}.pdf'.format(pre_name, z_name))
    
    print '\n#-----------------------Done------------------\n#'
    return

if __name__ == "__main__":
    sys.path.insert(0,'/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()