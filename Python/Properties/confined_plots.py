#! /usr/bin/env python
# Quick diffusion coefficient calculation given VACF

import argparse
import matplotlib.pyplot as plt
import matplotlib
import sys
import csv
import numpy as np
from scipy.integrate import simps
from itertools import izip,imap
from confined_tools import *
from matplotlib import cm

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
    parser.add_argument('-s', '--start', nargs='+',required=False, default='None', action='store',
                       help='Starting configurations')
    args = parser.parse_args()
    return args

def main():


    args = GetArgs()

    units = args.units
    dZ  = args.width
    region  = args.region
    plot_name = args.name
    pre = args.pre
    configs = args.start

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

    markers = ['o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P','X', '8', 'p',\
    'o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P','X', '8', 'p']
    colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026',\
    '#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']
    LABEL = {'C_vv': '$xyz$', 'C_vv_y': '$xy$'}
    layers =  {'6':2, '6.5':2, '7':2, '8':2, '8.5':2 ,'9':2, '9.5':3, '10':3, \
    '11':3, '12':3, '13':4, '14':4, '15':4, '16':5, '17':5, '20':6, '22':6, '25':7,\
    '27':8, '30':9, '33':10, '35':10, '40':12}


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

    # pmf
    fig7 = plt.figure(figsize=fig_size_sq)

    # eta wall
    fig9 = plt.figure(figsize=fig_size_sq)
    ax9  = fig9.add_axes([0.15,0.15,0.75,0.75])

    fig9a = plt.figure(figsize=fig_size_sq)
    ax9a  = fig9a.add_axes([0.15,0.15,0.75,0.75])

    # eta total
    fig10 = plt.figure(figsize=fig_size_sq)
    ax10  = fig10.add_axes([0.15,0.15,0.75,0.75])

    fig10a = plt.figure(figsize=fig_size_sq)
    ax10a  = fig10a.add_axes([0.15,0.15,0.75,0.75])

    # eta wall z
    fig11 = plt.figure(figsize=fig_size_sq)
    ax11  = fig11.add_axes([0.15,0.15,0.75,0.75])

    # eta total gk
    fig12 = plt.figure(figsize=fig_size_sq)
    ax12  = fig12.add_axes([0.15,0.15,0.75,0.75])

    fig12a = plt.figure(figsize=fig_size_sq)
    ax12a = fig12a.add_axes([0.15,0.15,0.75,0.75])

    # eta GK vs Ds
    fig13 = plt.figure(figsize=fig_size_sq)
    ax13  = fig13.add_axes([0.15,0.15,0.75,0.75])

    # kappa
    fig14 = plt.figure(figsize=fig_size_sq)
    ax14  = fig14.add_axes([0.15,0.15,0.75,0.75])

    # kappa/eta
    fig15 = plt.figure(figsize=fig_size_sq)
    ax15  = fig15.add_axes([0.15,0.15,0.75,0.75])

    # diffusion
    fig16 = plt.figure(figsize=fig_size_sq)
    ax16  = fig16.add_axes([0.15,0.15,0.75,0.75])

    # height layer
    fig17 = plt.figure(figsize=fig_size_sq)
    ax17  = fig17.add_axes([0.15,0.15,0.75,0.75])

    fig17a = plt.figure(figsize=fig_size_sq)
    ax17a  = fig17a.add_axes([0.15,0.15,0.75,0.75])

    # 2d dens
    #fig8 = plt.figure(figsize=fig_size_sq)
    


    # Properties independent of dimensionality
    
    print '\n#-----------------------Analysing density, pressure and viscosities------------------\n#'
    RHOdict = {}
    Hlist = []
    Z_init = []
    eta_diff_wall, eta_diff_tot, eta_diff_tot_y, eta_tot_gk, eta_tot_gk_acf, eta_xy_gk = [], [], [], [], [], []
    kappa_tot_gk, kappa_eta_tot = [], []
    diff_tot_y, diff_tot, diff_tot_y_msd, diff_tot_y_msd_err = [], [], [], []
    for dz in dZ:
        for c in configs:
            Z_init.append(dz)
            count = 1
            count_list, diff_list = [], []
            f = 'spce_T298_z{}_eps1.0_{}'.format(dz,c)
            print '\n#-----------------------dz = {}, DEN = {}------------------\n#'.format(dz,c)

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
                Hlist.append(Z_right-Z_left)
                height = Z_right-Z_left
                print 'Channel height:', height
                label = '$\Delta z = {}$ \AA'.format(dz)
                #RHOfit, roots = interpolate_derivative(Z, RHO)
                #roots = roots[(roots>Z_left) & (roots < Z_right)]
                #print 'Roots for {} are:'.format(dz), roots


                ax2.plot(Z, RHO, label=label)

                fig2a.clear()
                matplotlib.rcParams.update({'font.size': 30})
                ax2a  = fig2a.add_axes([0.15,0.15,0.75,0.75])
                ax2a.plot(RHO,Z)
                #ax2a.plot(RHOfit,Z, linestyle='dashed')
                ax2a.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax2a.set_ylabel('$z-z_{\mathrm{mid}}$')
                ax2a.set_ylim(Z_left-0.5, Z_right+0.5)#limits[dz][0], limits[dz][1])
                #ax2a.set_ylim(-8, 8)
                ax2a.set_xlim(0, 7)
                fig2a.savefig('PLOTS_C/densprof_sideways_{}_{}.pdf'.format(dz,c),bbox_inches='tight')


            except IOError:
                print 'File densprof.{} does not exist.'.format(f)

            # oxygen and hydrogen densities
            try:
                print 'Reading in denso.{}'.format(f)
                RHO_O, CY, CZ = read_densprof_2d(f,'o')

                den_max = np.max(RHO_O)
                den_min = np.min(RHO_O)
                fig8 = plt.figure(figsize=(fig_size_sq)) 
                ax8  = fig8.add_axes([0.1,0.15,0.8,0.75])
                #fig1.text(0.44, 0.025, '$y$ (\AA)', ha='center', va='center', fontsize=26)
                #ax8.set_xlim(ylo,yhi)
                #ax8.set_ylim(zlo,zhi)
                ax8.set_xlabel('$y$ (\AA)')
                ax8.set_ylabel('$z$ (\AA)')
                #plt.tick_params(pad=7)
                ctest=ax8.contourf(CY, CZ, np.transpose(RHO_O), cmap=cm.RdBu, levels=np.linspace(den_min,den_max,500))
                fig8.colorbar(ctest)
                fig8.savefig('PLOTS_C/denso_2d_z{}.png'.format(dz),bbox_inches='tight')
                fig8.clear()

                print 'Reading in densh.{}'.format(f)
                RHO_H, CY, CZ = read_densprof_2d(f,'h')

                den_max = np.max(RHO_H)
                den_min = np.min(RHO_H)
                fig8 = plt.figure(figsize=(fig_size_sq)) 
                ax8  = fig8.add_axes([0.1,0.15,0.8,0.75])
                #fig1.text(0.44, 0.025, '$y$ (\AA)', ha='center', va='center', fontsize=26)
                #ax8.set_xlim(ylo,yhi)
                #ax8.set_ylim(zlo,zhi)
                ax8.set_xlabel('$y$ (\AA)')
                ax8.set_ylabel('$z$ (\AA)')
                #plt.tick_params(pad=7)
                ctest=ax8.contourf(CY, CZ, np.transpose(RHO_H), cmap=cm.RdBu, levels=np.linspace(den_min,den_max,500))
                fig8.colorbar(ctest)
                fig8.savefig('PLOTS_C/densh_2d_z{}_{}.png'.format(dz,c),bbox_inches='tight')
                fig8.clear()


            except IOError:
                print 'File denso.{} does not exist.'.format(f)

            # pressure profile
            try:
                print 'Reading in stress.{}'.format(f)
                Z_P, P, Ptot, Pxy, Pxz, Pyz, Pxx, Pyy, Pzz, delz = stress_prof(f, 'None')
                MID_P, LEFT_P, RIGHT_P = mid_point(Z_P,P)
                Z_P = Z_P-Z_P[MID_P]
                label = '$\Delta z = {}$ \AA'.format(dz)
                #ax4.plot(Z_P, P, label=label)
                matplotlib.rcParams.update({'font.size': 22})
                ax4.plot(Z_P, Pxx, label='Pxx', linestyle=ls[0])
                ax4.plot(Z_P, Pyy, label='Pyy', linestyle=ls[1])
                ax4.plot(Z_P, Pzz, label='Pzz', linestyle=ls[2])

            except IOError:
                print 'File stress.{} does not exist.'.format(f)

            try:
                fig5.clear()
                ax5  = fig5.add_axes([0.15,0.15,0.75,0.75])
                ax5.plot(Z_P, P, c=colours[0])
                ax5.set_ylabel('P (MPa)')
                ax5.set_ylim(-(np.max(P)+0.05*np.max(P))/2,np.max(P)+0.05*np.max(P))
                ax5.yaxis.label.set_color(colours[0])
                for tl in ax5.get_yticklabels():   
                    tl.set_color(colours[0])

                ax5b = ax5.twinx()
                ax5b.plot(Z, RHO, c=colours[6], linestyle='dashed')
                ax5b.set_ylim(-(np.max(RHO)+0.05*np.max(RHO))/2,np.max(RHO)+0.05*np.max(RHO))
                fig5.text(0.975, 0.5, '$\\rho_{\mathrm{wall}}$ (g/cm$^3$)', color=colours[6], ha='center', va='center', fontsize=24,rotation=270)
                for tl in ax5b.get_yticklabels():   
                    tl.set_color(colours[6])
                ax5.set_xlabel('$z-z_{\mathrm{mid}}$')
                #ax5.set_xlim(-17.4, 17.4)
                fig5.savefig('PLOTS_C/stressdens_{}_{}.pdf'.format(dz,c))
            except:
                pass

            # pmf
            try:
                print 'Reading in pmf.{}'.format(f)
                Z_PMF, PMF = pmf_prof(f)
                MID_PMF = int(len(PMF)/2)
                PMF_avg = np.mean(PMF[MID_PMF-20:MID_PMF+20])
                PMF = PMF-PMF_avg
                Z_PMF = Z_PMF-Z_PMF[MID_PMF]
                PMF_min = np.min(PMF[40:-60])
                PMF_max = np.max(PMF[60:-60])
                dPMF = PMF_max-PMF_min
                print 'The energy barrier for the first layer is:', dPMF
                label = '$\Delta z = {}$ \AA'.format(dz)

                fig7.clear()
                ax7  = fig7.add_axes([0.15,0.15,0.75,0.75])
                ax7.plot(Z_PMF, PMF, c=colours[0], lw=3,linestyle='-')
                ax7.set_ylabel('$\Delta$ PMF (kcal/mol)')
                ax7.set_ylim(-0.75,2)
                ax7.yaxis.label.set_color(colours[0])
                for tl in ax7.get_yticklabels():   
                    tl.set_color(colours[0])
                ax7b = ax7.twinx()
                ax7b.set_ylim(0,3.5)
                ax7b.plot(Z, RHO, c=colours[6], lw=3,linestyle='dashed')
                #ax3b.set_ylabel('$\\rho_{\mathrm{wall}}$ (g/cm$^3$)', rotation=270)
                #ax3b.yaxis.label.set_color(colours[6])
                fig7.text(0.995, 0.5, '$\\rho$ (g/cm$^3$)', color=colours[6], ha='center', va='center', fontsize=24,rotation=270)
                for tl in ax7b.get_yticklabels():   
                    tl.set_color(colours[6])
                ax7.set_xlabel('$\Delta z$ (\AA)')
                ax7.set_xlim(-15,15)
                fig7.savefig('PLOTS_C/pmf_{}.pdf'.format(dz),bbox_inches='tight')

            except IOError:
                print 'File pmf.{}_1 does not exist.'.format(f)

            # viscosity profile from Diffusion
            
            try:

                # wall viscosity
                print 'Reading in C_vv_y_{}_wall.dat'.format(f)
                C_vv_array = np.loadtxt("Water_Graphene/C_vv_y_{}_wall.dat".format(f))
                times =  C_vv_array[:,0]
                C_vv_ave = C_vv_array[:,1]
                diff_wall_tmp = diffusion(C_vv_ave, dt, time_conv, space_conv, 'C_vv_y')
                eta_diff_wall_tmp = eta_diff(diff_wall_tmp)
                eta_diff_wall.append(eta_diff_wall_tmp*1e3)

                # bulk shear viscosity
                print 'Reading in C_vv_{}_1_10000_z0_{}.dat'.format(f,int(dz)-1)
                C_vv_array = np.loadtxt("Water_Graphene/C_vv_{}_1_10000_z0_{}.dat".format(f, int(dz)-1))
                times =  C_vv_array[:,0]
                C_vv_ave = C_vv_array[:,1]
                diff_tot_tmp = diffusion(C_vv_ave, dt, time_conv, space_conv, 'C_vv')
                eta_diff_tot_tmp = eta_diff(diff_tot_tmp)
                eta_diff_tot.append(eta_diff_tot_tmp*1e3)
                diff_tot.append(diff_tot_tmp*1e9)

                # bulk shear viscosity surface
                # diff from vacf
                print 'Reading in C_vv_y_{}_1_10000_z0_{}.dat'.format(f,int(dz)-1)
                C_vv_array = np.loadtxt("Water_Graphene/C_vv_y_{}_1_10000_z0_{}.dat".format(f, int(dz)-1))
                times =  C_vv_array[:,0]
                C_vv_ave = C_vv_array[:,1]
                diff_tot_y_tmp = diffusion(C_vv_ave, dt, time_conv, space_conv, 'C_vv_y')
                # diff from msd
                diff_tot_y_msd_tmp, diff_tot_y_msd_err_tmp = diff_msd(f)
                # eta from diff vacf
                eta_diff_tot_y_tmp = eta_diff(diff_tot_y_tmp)
                eta_diff_tot_y.append(eta_diff_tot_y_tmp*1e3)
                diff_tot_y.append(diff_tot_y_tmp*1e9)
                diff_tot_y_msd.append(diff_tot_y_msd_tmp*1e9)
                diff_tot_y_msd_err.append(diff_tot_y_msd_err_tmp*1e9)


                print 'The viscosities from diffusion are:', eta_diff_wall_tmp*1e3, eta_diff_tot_tmp*1e3, eta_diff_tot_y_tmp*1e3




                if dz == '30':
                    # wall viscosity as fn of z
                    eta_wall_z = []
                    z_eta_wall = []
                    for i in range(2,10):
                        print 'Reading in C_vv_y_{}_1_1_10000_z{}_{}.dat'.format(f,i,i+1)
                        C_vv_array = np.loadtxt("Water_Graphene/C_vv_y_{}_1_1_10000_z{}_{}.dat".format(f,i,i+1))
                        times =  C_vv_array[:,0]
                        C_vv_ave = C_vv_array[:,1]
                        diff_wall_z_tmp = diffusion(C_vv_ave, dt, time_conv, space_conv, 'C_vv_y')
                        eta_wall_z_tmp = eta_diff(diff_wall_z_tmp)
                        eta_wall_z.append(eta_wall_z_tmp*1e3)
                        z_eta_wall.append(i+0.5)


            except IOError:
                print 'File C_vv_{}.dat or similar does not exist.'.format(f)


            # viscosity profile from GK
            
            try:

                

                # effective shear viscosity (GK)
                print 'Reading in visc.{}'.format(f,int(dz)-1)
                eta_gk_tmp = visc_gk(f, height, 'etas','eff')
                eta_tot_gk.append(eta_gk_tmp*1e3)

                # parallel shear viscosity (GK)
                print 'Reading in visc.{}'.format(f,int(dz)-1)
                eta_gk_xy_tmp = visc_gk(f, height, 'etas', 'xy')
                eta_xy_gk.append(eta_gk_xy_tmp*1e3)

                # bulk shear viscosity (GK, acf)
                #print 'Reading in acfsv.{}'.format(f,int(dz)-1)
                #eta_gk_acf_tmp = eta_gk_acf(f, height)
                #eta_tot_gk_acf.append(eta_gk_acf_tmp*1e3)

                print 'The viscosities from GK are:', eta_gk_tmp*1e3, eta_gk_xy_tmp*1e3 #, eta_gk_acf_tmp*1e3

                # bulk viscosity
                print 'Reading in visc.{}'.format(f,int(dz)-1)
                kappa_gk_tmp = visc_gk(f, height, 'etab','nvt')
                kappa_tot_gk.append(kappa_gk_tmp*1e3)
                kappa_eta_tot.append(kappa_gk_tmp/eta_gk_xy_tmp)




            except IOError:
                print 'File visc.{} or similar does not exist.'.format(f)





            # viscosity profile

            # try:
            #     print 'Reading in stresseta.{}'.format(f)
            #     Z_P, P, Ptot, Pxy, Pxz, Pyz, Pxx, Pyy, Pzz, delz = stress_prof(f, 'eta')
            #     MID_P, LEFT_P, RIGHT_P = mid_point(Z_P,P)
            #     Z_P = Z_P-Z_P[MID_P]
            #     #DZ_eta = Z_P[RIGHT_P]-Z_P[LEFT]
            #     DZ_eta = height
            #     Eta_xy = eta(f,np.mean(Pxy,axis=1), DZ_eta)
            #     Eta_yz = eta(f,np.mean(Pyz,axis=1), DZ_eta)
            #     Eta_xz = eta(f,np.mean(Pxz,axis=1), DZ_eta)

            #     print 'Viscosities: ', Eta_xy, Eta_yz, Eta_xz

                # normal corrlen is too long. The length we have is only 422 data points
                '''
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
                '''
            # except IOError:
            #     print 'File stresseta.{}_1 does not exist.'.format(f)


    # Diffusion related properties
    print '\n#-----------------------Analysing diffusion coefficient------------------\n#'
    count_p = 0
    #diff_bulk = 2.86624087647e-09 # needs to be different for 1 and 2 d!!
    pre_name = ''
    for p in pre:
        pre_name += LABEL[p].split('$')[1]+'_'
        z_list = []
        d_list, d_list_tmp, d_list_err = [], [], []
        rho_list,rho_list_tmp, rho_list_err = [], [], []
        zr_diff = []
    

        z_name = 'z_'
        for dz in dZ:
            diff_list, diff_list_err = [], []
            f = 'spce_T298_z{}_eps1.0'.format(dz)

            '\n#-----------------------dz = {}------------------\n#'.format(dz)

            z_name += dz

            for r in region:
                diff_r_tmp = []
                for c in configs:
                    zr_diff.append((dz,r))

                    # diffusion
                    try:
                        print 'Reading in {}_{}_{}_{}.dat'.format(p,f,c,r)
                        C_vv_array = np.loadtxt("Water_Graphene/{}_{}_{}_{}.dat".format(p,f,c,r))
                        times =  C_vv_array[:,0]
                        C_vv_ave = C_vv_array[:,1]
                        diff_tmp = diffusion(C_vv_ave, dt, time_conv, space_conv, p)
                        diff_bulk = 2.86624087647e-09
                        #diff_list.append(diff_tmp)#/diff_bulk)
                        diff_r_tmp.append(diff_tmp/diff_bulk)
                        if r == 'wall' and p == 'C_vv_y':
                            print dz
                            RHOdat = RHOdict[dz]
                            RHO_max = rho_wall_max(RHOdat[0], RHOdat[1])
                            z_list.append(dz)
                            d_list_tmp.append(diff_tmp/diff_bulk)
                            rho_list_tmp.append(RHO_max)
                            print diff_tmp/diff_bulk, RHO_max



                    except IOError:
                        print 'File {}_{}_{}_{}.dat does not exist.'.format(p,f,c,r)
                diff_list.append(np.mean(diff_r_tmp))
                diff_list_err.append(np.std(diff_r_tmp)/len(diff_r_tmp))
                d_list.append(np.mean(d_list_tmp))
                d_list_err.append(np.std(d_list_tmp)/len(d_list_tmp))
                rho_list.append(np.mean(rho_list_tmp))
                rho_list_err.append(np.std(rho_list_tmp)/len(rho_list_tmp))
                count_list = range(1,len(diff_list)+1)

            ax1.errorbar(count_list, diff_list, yerr=diff_list_err, marker=markers[count_p], label=LABEL[p], markersize=12, linestyle='dashed')
            count_p += 1


    # average over input configurations
    if len(configs)>1:
        Z, Hlist, Hlist_err = averaging(Z_init,Hlist)
        Z, eta_diff_wall, eta_diff_wall_err = averaging(Z_init,eta_diff_wall)
        Z, eta_diff_tot, eta_diff_tot_err = averaging(Z_init,eta_diff_tot)
        Z, eta_diff_tot_y, eta_diff_tot_y_err = averaging(Z_init,eta_diff_tot_y)
        Z, eta_tot_gk, eta_tot_gk_err = averaging(Z_init,eta_tot_gk)
        Z, eta_xy_gk, eta_xy_gk_err = averaging(Z_init,eta_xy_gk)
        Z, diff_tot_y, diff_tot_y_err = averaging(Z_init,diff_tot_y)
        Z, diff_tot_y_msd, diff_tot_y_msd_err = averaging(Z_init,diff_tot_y_msd)
        Z, kappa_tot_gk, kappa_tot_gk_err = averaging(Z_init,kappa_tot_gk)
        Z, kappa_eta_tot, kappa_eta_tot_err = averaging(Z_init,kappa_eta_tot)

    else:
        Hlist_err = [0]*len(Hlist)
        eta_diff_wall_err = [0]*len(eta_diff_wall)
        eta_diff_tot_err = [0]*len(eta_diff_tot)
        eta_diff_tot_y_err = [0]*len(eta_diff_tot_y)
        eta_tot_gk_err = [0]*len(eta_tot_gk)
        eta_xy_gk_err = [0]*len(eta_xy_gk)
        diff_tot_y_err = [0]*len(diff_tot_y)
        kappa_tot_gk_err = [0]*len(kappa_tot_gk)
        kappa_eta_tot_err = [0]*len(kappa_eta_tot)

    print diff_tot_y_msd_err
    # print to csv for latex use
    with open('DATA/heights_{}.csv'.format(z_name), 'wb') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(("dz", "H"))
        writer.writerows(imap(lambda x,y: (y, round(x, 2)), Hlist, dZ))#izip(dZ,round(Hlist,2)))


    print '\n#-----------------------Final plotting------------------\n#'

    #ax1.set_xlabel('region')
    ax1.set_ylabel('$D_s/D_{\mathrm{iso}}$')
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

    ax3b = ax3.twinx()
    #if len(dZ) > 1:
    #    Hfit, diff_fit = exp_fit(Hlist, d_list, 1, 40, [-100, 10, 1.5])
    #    ax3.plot(Hfit, diff_fit, c=colours[0], linestyle='dashed', marker='None')
    #    Hfit, rho_fit = exp_fit(Hlist, rho_list, 1, 40)
    #    ax3b.plot(Hfit, rho_fit, c=colours[6], linestyle='dashed', marker='None')
    ax3.plot(Hlist, d_list, c=colours[0], linestyle='dashed', marker='D')
    ax3.set_ylabel('$D_{||}/D_{\mathrm{iso}}$')
    #ax3.set_ylim(0,1.3)
    ax3.yaxis.label.set_color(colours[0])
    for tl in ax3.get_yticklabels():   
        tl.set_color(colours[0])
    ax3b.plot(Hlist, rho_list, c=colours[6], linestyle='dashed', marker='o')
    #ax3b.set_ylabel('$\\rho_{\mathrm{wall}}$ (g/cm$^3$)', rotation=270)
    #ax3b.yaxis.label.set_color(colours[6])
    fig3.text(0.975, 0.5, '$\\rho_{\mathrm{wall}}$ (g/cm$^3$)', color=colours[6], ha='center', va='center', fontsize=24,rotation=270)
    for tl in ax3b.get_yticklabels():   
        tl.set_color(colours[6])
    ax3.set_xlabel('H (\AA)')
    #ax3.set_xlabel('$z-z_{\mathrm{mid}}$')
    #ax3.legend()

    ax4.set_ylabel('P (MPa)')
    ax4.set_xlabel('$z-z_{\mathrm{mid}}$')
    #ax2.set_xlim(Z_left, Z_right)
    ax4.set_xlim(-18, 18)
    ax4.legend()
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

    # ETA

    ax9.errorbar(Hlist, eta_diff_wall, yerr=eta_diff_wall_err, marker='o')
    ax9.plot(Hlist, len(Hlist)*[0.67], linestyle='dashed', label='Bulk')
    ax9.set_xlabel('H (\AA)')
    ax9.set_ylabel('$\eta_{\mathrm{wall}}$ (mPas)')
    ax9.legend()

    ax9a.set_yscale('log')
    ax9a.errorbar(Hlist, eta_diff_wall, yerr=eta_diff_wall_err, marker='o')
    ax9a.plot(Hlist, len(Hlist)*[0.67], linestyle='dashed', label='Bulk')
    ax9a.set_xlabel('H (\AA)')
    ax9a.set_ylabel('$\eta_{\mathrm{wall}}$ (mPas)')
    ax9a.legend()

    ax10.errorbar(Hlist, eta_diff_tot, yerr=eta_diff_tot_err, marker='o', label='isotropic')
    ax10.errorbar(Hlist, eta_diff_tot_y, yerr=eta_diff_tot_y_err, marker='o', label='xy')
    ax10.plot(Hlist, len(Hlist)*[0.67], linestyle='dashed', label='Bulk')
    ax10.set_xlabel('H (\AA)')
    ax10.set_ylabel('$\eta_{\mathrm{channel}}$ (mPas)')
    ax10.legend()

    ax10a.set_yscale('log')
    ax10a.errorbar(Hlist, eta_diff_tot, yerr=eta_diff_tot_err, marker='o', label='isotropic')
    ax10a.errorbar(Hlist, eta_diff_tot_y, yerr=eta_diff_tot_y_err, marker='o', label='xy')
    ax10a.plot(Hlist, len(Hlist)*[0.67], linestyle='dashed', label='Bulk')
    ax10a.set_xlabel('H (\AA)')
    ax10a.set_ylabel('$\eta_{\mathrm{channel}}$ (mPas)')
    ax10a.legend()

    if '30' in dZ:
        ax11.plot(z_eta_wall, eta_wall_z, marker='o')
        ax11.plot(z_eta_wall, len(z_eta_wall)*[0.67], linestyle='dashed', label='Bulk')
        ax11.set_xlabel('z (\AA)')
        ax11.set_ylabel('$\eta_{\mathrm{xy}}$ (mPas)')
        ax11.legend()
        fig11.savefig('PLOTS_C/eta_diff_wall_z.pdf',bbox_inches='tight')

    ax12.errorbar(Hlist, eta_tot_gk, yerr=eta_tot_gk_err, marker='o', label='$\eta_{\mathrm{eff}}$')
    ax12.errorbar(Hlist, eta_xy_gk, yerr=eta_xy_gk_err, marker='o', label='$\eta_{\mathrm{||}}$')
    #ax12.plot(Hlist, eta_tot_gk_acf, marker='o', label='ACF')
    ax12.plot(Hlist, len(Hlist)*[0.67], linestyle='dashed', label='Bulk')
    ax12.set_xlabel('H (\AA)')
    ax12.set_ylabel('$\eta_{\mathrm{GK}}$ (mPas)')
    ax12.legend()

    ax12a.set_yscale('log')
    ax12a.errorbar(Hlist, eta_tot_gk, yerr=eta_tot_gk_err, marker='o', label='$\eta_{\mathrm{eff}}$')
    ax12a.errorbar(Hlist, eta_xy_gk, yerr=eta_xy_gk_err, marker='o', label='$\eta_{\mathrm{||}}$')
    ax12a.plot(Hlist, len(Hlist)*[0.67], linestyle='dashed', label='Bulk')
    #ax12a.semilogy(Hlist, eta_tot_gk, marker='o', label='$\eta_{\mathrm{eff}}$')
    #ax12a.semilogy(Hlist, eta_xy_gk, marker='o', label='$\eta_{\mathrm{||}}$')
    #ax12a.semilogy(Hlist, len(Hlist)*[0.67], linestyle='dashed', label='Bulk')
    ax12a.set_xlabel('H (\AA)')
    ax12a.set_ylabel('$\eta_{\mathrm{GK}}$ (mPas)')
    ax12a.legend()

    #ax13.plot(Hlist, eta_tot_gk, marker='o', label='$\eta_{\mathrm{eff}}$')
    ax13.plot(1/np.array(diff_tot_y), eta_diff_tot_y, ls='None', marker='s', label='Stokes-Einstein')
    ax13.plot(1/np.array(diff_tot_y), eta_xy_gk, ls='None', marker='o', label='Green-Kubo')
    ax13.set_xlabel('$1/D_s$ ($10^{9}$s/m$^2$)')
    ax13.set_ylabel('$\eta$ (mPas)')
    ax13.set_xlim(0,10)
    ax13.set_ylim(0,25)
    #ax13.set_xlim(0.2,0.8)
    #ax13.set_ylim(0.4,1)
    ax13.legend()

    ax14.errorbar(Hlist, kappa_tot_gk, yerr=kappa_tot_gk_err, marker='o')
    ax14.plot(Hlist, len(Hlist)*[1.59], linestyle='dashed', label='Bulk')
    ax14.set_xlabel('H (\AA)')
    ax14.set_ylabel('$\kappa_{\mathrm{GK}}$ (mPas)')

    ax15.errorbar(Hlist, kappa_eta_tot, yerr=kappa_eta_tot_err, marker='o')
    ax15.plot(Hlist, len(Hlist)*[2.32], linestyle='dashed', label='Bulk')
    ax15.set_xlabel('H (\AA)')
    ax15.set_ylabel('$\kappa/\eta_{||}$ ')

    ax16.errorbar(Hlist, diff_tot_y, yerr=diff_tot_y_err, marker='o', label='VACF')
    ax16.errorbar(Hlist, diff_tot_y_msd, yerr=diff_tot_y_msd_err, marker='o', label='MSD')
    ax16.plot(Hlist, len(Hlist)*[2.866], linestyle='dashed', label='Bulk')
    ax16.set_xlabel('H (\AA)')
    ax16.set_ylabel('$D_s$ (10$^{-9}$ m$^2$/s)')
    ax16.legend()

    dlayer = 2.75
    ax17.plot(Hlist, np.array(Hlist)/dlayer, marker='o', linestyle='dashed')
    ax17.plot(Hlist,len(Hlist)*[1], linestyle='dashed')
    ax17.plot(Hlist,len(Hlist)*[2], linestyle='dashed')
    ax17.plot(Hlist,len(Hlist)*[3], linestyle='dashed')
    ax17.plot(Hlist,len(Hlist)*[4], linestyle='dashed')
    ax17.plot(Hlist,len(Hlist)*[5], linestyle='dashed')
    ax17.plot(Hlist,len(Hlist)*[6], linestyle='dashed')
    ax17.plot(Hlist,len(Hlist)*[7], linestyle='dashed')
    ax17.set_xlabel('H (\AA)')
    ax17.set_ylabel('H/d$_{\mathrm{layer}}$')
    ax17.legend()

    Hlayer_sample = np.linspace(1,np.max(Hlist)/dlayer,100)
    ax17a.plot(np.array(Hlist)/dlayer, map(layers.get, dZ), marker='o', linestyle='dashed')
    ax17a.plot(Hlayer_sample, f1(Hlayer_sample, 1, 0), linestyle='dotted')
    ax17a.set_ylabel('No. of layers')
    ax17a.set_xlabel('H/d$_{\mathrm{layer}}$')
    ax17a.legend()

    

    fig2.savefig('PLOTS_C/densprof_{}.pdf'.format(z_name),bbox_inches='tight')
    fig4.savefig('PLOTS_C/stressprof_{}.pdf'.format(z_name),bbox_inches='tight')
    fig3.savefig('PLOTS_C/diffdens_wall_{}.pdf'.format(z_name),bbox_inches='tight')
    fig9.savefig('PLOTS_C/eta_diff_wall_{}.pdf'.format(z_name),bbox_inches='tight')
    fig9a.savefig('PLOTS_C/eta_diff_wall_log_{}.pdf'.format(z_name),bbox_inches='tight')
    fig10.savefig('PLOTS_C/eta_diff_channel_{}.pdf'.format(z_name),bbox_inches='tight')
    fig10a.savefig('PLOTS_C/eta_diff_channel_log_{}.pdf'.format(z_name),bbox_inches='tight')
    fig12.savefig('PLOTS_C/eta_gk_channel_{}.pdf'.format(z_name),bbox_inches='tight')
    fig12a.savefig('PLOTS_C/eta_gk_channel_log_{}.pdf'.format(z_name),bbox_inches='tight')
    fig13.savefig('PLOTS_C/eta_vs_diff_{}.pdf'.format(z_name),bbox_inches='tight')
    fig14.savefig('PLOTS_C/kappa_gk_{}.pdf'.format(z_name),bbox_inches='tight')
    fig15.savefig('PLOTS_C/kappa_eta_{}.pdf'.format(z_name),bbox_inches='tight')
    fig16.savefig('PLOTS_C/diff_msd_vacf_{}.pdf'.format(z_name),bbox_inches='tight')
    fig17.savefig('PLOTS_C/height_layer_{}.pdf'.format(z_name),bbox_inches='tight')
    fig17a.savefig('PLOTS_C/height_layer2_{}.pdf'.format(z_name),bbox_inches='tight')
    fig1.savefig('PLOTS_C/diffusion_{}{}.pdf'.format(pre_name, z_name))
    
    print '\n#-----------------------Done------------------\n#'
    return

if __name__ == "__main__":
    sys.path.insert(0,'/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()