#! /usr/bin/env python
# Plots for confined water data

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
from matplotlib import ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def GetArgs():
    ## Parse command line arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', '--width', nargs='+', required=False, default='liquid', action='store',
                       help='Width of channel')
    parser.add_argument('-r', '--region', nargs='+', required=False, default='liquid', action='store',
                       help='Which region? wall, intermediate, bulk')
    parser.add_argument('-u', '--units', required=False, type=str, action='store', default=[0,100],
                        help='metal or real')
    parser.add_argument('-p', '--pre', nargs='+', required=False, type=str, action='store', default=['C_vv'],
                        help='C_vv or C_vv_y')
    parser.add_argument('-n', '--name', required=False, default='None', action='store',
                       help='Plot name')
    parser.add_argument('-s', '--start', nargs='+',required=False, default='None', action='store',
                       help='Starting configurations')
    parser.add_argument('-rho', '--rho', required=False, default='RHO1', action='store',
                       help='Average channel density')
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
    rho = args.rho

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
    fig_size_sq2 = (6,9)

    markers = ['o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P','X', '8', 'p',\
    'o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P','X', '8', 'p',\
    'o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P','X', '8', 'p']
    colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026',\
    '#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026',\
    '#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']
    LABEL = {'C_vv': '$xyz$', 'C_vv_y': '$xy$'}
    


    ls = ['--','-', '--', '-.', ':','-', '--', '-.', ':']

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
    # density profile errorbar
    fig2b = plt.figure(figsize=fig_size_sq)

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
    fig15 = plt.figure(figsize=fig_size_long)
    ax15  = fig15.add_axes([0.15,0.15,0.75,0.75])

    # diffusion
    fig16 = plt.figure(figsize=fig_size_sq)
    ax16  = fig16.add_axes([0.15,0.15,0.75,0.75])

    fig16a = plt.figure(figsize=fig_size_sq)
    ax16a  = fig16a.add_axes([0.15,0.15,0.75,0.75])

    # height layer
    fig17 = plt.figure(figsize=fig_size_sq)
    ax17  = fig17.add_axes([0.15,0.15,0.75,0.75])

    fig17a = plt.figure(figsize=fig_size_sq)
    ax17a  = fig17a.add_axes([0.15,0.15,0.75,0.75])

    # friction coeff
    fig18 = plt.figure(figsize=fig_size_sq)
    ax18  = fig18.add_axes([0.15,0.15,0.75,0.75])

    # fric and ls
    fig18a = plt.figure(figsize=fig_size_long)
    ax18a  = fig18a.add_axes([0.15,0.15,0.75,0.75])

    # work of adhesion
    fig19 = plt.figure(figsize=fig_size_sq)
    ax19  = fig19.add_axes([0.15,0.15,0.75,0.75])

    fig19a = plt.figure(figsize=fig_size_sq)
    ax19a  = fig19a.add_axes([0.15,0.15,0.75,0.75])

    # slip length
    fig20 = plt.figure(figsize=fig_size_sq)
    ax20  = fig20.add_axes([0.15,0.15,0.75,0.75])

    # slip length by aspect ratio
    fig21 = plt.figure(figsize=fig_size_sq)
    ax21  = fig21.add_axes([0.15,0.15,0.75,0.75])

    # enhancement
    fig22 = plt.figure(figsize=fig_size_sq)
    ax22  = fig22.add_axes([0.15,0.15,0.75,0.75])

    # scaled enhancement
    fig22a = plt.figure(figsize=fig_size_sq)
    ax22a  = fig22a.add_axes([0.15,0.15,0.75,0.75])

    # Ds/Wa channel length enhancement
    fig22b = plt.figure(figsize=fig_size_sq)
    ax22b  = fig22b.add_axes([0.15,0.15,0.75,0.75])

    # Ds and eta comb
    fig23 = plt.figure(figsize=fig_size_sq)
    #ax23  = fig23.add_axes([0.15,0.15,0.75,0.75])

    # Ds, eta, kappa/eta comb
    fig24 = plt.figure(figsize=(9,9))

    # VACF wall bulk
    fig25 = plt.figure(figsize=fig_size_sq)


    # VACF comp
    fig26 = plt.figure(figsize=fig_size_sq)
    ax26  = fig26.add_axes([0.15,0.15,0.75,0.75])

    # VACF comp msd
    fig26a = plt.figure(figsize=fig_size_sq)
    ax26a  = fig26a.add_axes([0.15,0.15,0.75,0.75])
    ax26a_i  = fig26a.add_axes([0.35,0.63,0.25,0.25])

    # Pressure (H)
    fig27 = plt.figure(figsize=fig_size_sq)
    ax27  = fig27.add_axes([0.15,0.15,0.75,0.75])

    # Pressure T max (H)
    fig27a = plt.figure(figsize=fig_size_sq)
    ax27a  = fig27a.add_axes([0.15,0.15,0.75,0.75])

    # Rho_ave (H)
    fig28 = plt.figure(figsize=fig_size_sq)
    ax28  = fig28.add_axes([0.15,0.15,0.75,0.75])

    # Pave vs Kappa (H)
    fig29 = plt.figure(figsize=fig_size_sq)
    ax29  = fig29.add_axes([0.15,0.15,0.75,0.75])

    # rho configs
    fig30 = plt.figure(figsize=fig_size_sq)
    ax30  = fig30.add_axes([0.15,0.15,0.75,0.75])

    # area density
    fig31 = plt.figure(figsize=fig_size_sq)
    ax31  = fig31.add_axes([0.15,0.15,0.75,0.75])


    # 2d dens
    #fig8 = plt.figure(figsize=fig_size_sq)

    if rho == 'RHO1':
        steps = 3000
        DIR = 'RHO1'
        layers =  {'6':1, '6.2': 1, '6.8':1, '7':1, '7.2':1, '7.5': 2, '7.8':2, \
            '8':2, '8.2':2, '8.5':2 ,'8.8':2, '9':2, '9.2':2, '9.5':2, '9.8':2,'10':2, \
            '10.2':2, '10.5':2, '10.8':2, '11':2, '11.5':3, '12':3, '12.5':3,'13':3, \
            '13.4':3,'14':4, '14.5':4,\
            '15':4, '16':4, '17':4, '20':5, '22':6, '25':7,\
            '27':8, '30':9, '33':10, '35':10, '40':12}
    else:
        steps = 10000
        DIR = 'RHO1.1'
        layers =  {'6':2, '6.2': 2, '6.5':2, '7':2, '7.2':2, '7.8':2, \
            '8':2, '8.2':2, '8.5':2 ,'8.8':2, '9':2, '9.2':2, '9.5':3, '9.8':3,'10':3, \
            '10.2':3, '10.8':3, '11':3, '11.5':3, '12.1':3, '12.5':4,'13.1':4, '13.5':4,'14.1':4, '14.5':4,\
            '15':4, '16':5, '17':5, '20':6, '22':6, '25':7,\
            '27':8, '30':9, '33':10, '35':10, '40':12}
    


    # Properties independent of dimensionality
    
    print '\n#-----------------------Analysing density, pressure and viscosities------------------\n#'
    RHOdict = {}
    Hlist, Hlist_vacf = [], []
    Z_init, Z_init_vacf = [], []
    oh_flag = 0
    eta0 = 0.67 #mPas
    k0 = 2.7 #Angstrom
    Pave_tot, PTmax, rhoave, area_rho, rhoeff = [], [], [], [], []
    eta_diff_wall, eta_diff_tot, eta_diff_tot_y, eta_tot_gk, eta_tot_gk_acf, eta_xy_gk = [], [], [], [], [], []
    kappa_tot_gk, kappa_eta_tot, fric, wa, Ls_fric, Ls_wa = [], [], [], [], [], []
    diff_tot_y, diff_tot, diff_tot_y_msd, diff_tot_y_msd_err = [], [], [], []
    count_tmp = 0
    for dz in dZ:
        for c in configs:
            Z_init.append(float(dz))
            count = 1
            count_list, diff_list = [], []
            f = 'spce_T298_z{}_eps1.0_{}'.format(dz,c)
            print '\n#-----------------------dz = {}, DEN = {}------------------\n#'.format(dz,c)

            # area density
            try:
                area_xy = geometry(DIR,f)[0]

                # find number of atoms
                file = open('Water_Graphene/{}/log.{}'.format(DIR,f),'r')

                data = file.read()
                data_lines = data.split('\n')

                DATA = []
                flag=0
                for i in range(len(data_lines)-1):
                    if data_lines[i].split() != [] and data_lines[i].split()[0] == 'group' and data_lines[i].split()[1] == 'oxygen':
                        numwater = int(data_lines[i+1].split()[0])
                        break
                area_dens = numwater/(area_xy*1e-2)
                area_rho.append(area_dens)

                RHOeff = avg_density('Water_Graphene', f, DIR)
                if DIR == 'RHO1.1':
                        RHOeff *= (float(dz)-0.5)/(float(dz)-3.19)
                rhoeff.append(RHOeff)
                

            except IOError:
                print 'File log.{} does not exist.'.format(f)


            # density profile
            try:
                print 'Reading in densprof.{}'.format(f)
                Z, RHO, RHOerr = read_densprof(DIR,f)
                MID, LEFT, RIGHT = mid_point(Z,RHO)
                # store data for use in diffusion plots
                RHOdict[dz] = (RHO, LEFT)
                RHOave = np.mean(RHO[LEFT:RIGHT])
                rhoave.append(RHOave)
                Z = Z-Z[MID]
                Z_left = Z[LEFT]
                Z_right = Z[RIGHT]
                Hlist.append(float(dz)-3.19)#Z_right-Z_left)
                height = Z_right-Z_left
                print 'Channel height:', height
                label = '$\Delta z = {}$ \AA'.format(dz)
                #RHOfit, roots = interpolate_derivative(Z, RHO)
                #roots = roots[(roots>Z_left) & (roots < Z_right)]
                #print 'Roots for {} are:'.format(dz), roots


                ax2.plot(Z, RHO, label=label)
                ax30.plot(Z, RHO, label=c)

                fig2a.clear()
                matplotlib.rcParams.update({'font.size': 30})
                ax2a  = fig2a.add_axes([0.15,0.15,0.75,0.75])
                ax2a.plot(RHO,Z)
                #ax2a.plot(RHOfit,Z, linestyle='dashed')
                ax2a.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax2a.set_ylabel('$z-z_{\mathrm{mid}}$ (\AA)')
                #fig2a.text(0.01, 0.54, '$z-z_{\mathrm{mid}}$ (\AA)', ha='center', va='center', rotation=90) #, fontsize=26)
                ax2a.set_ylim(Z_left-0.5, Z_right+0.5)#limits[dz][0], limits[dz][1])
                #ax2a.set_ylim(-8, 8)
                ax2a.set_xlim(0, 3)
                fig2a.savefig('PLOTS_C/{}/densprof_sideways_{}_{}.pdf'.format(DIR,dz,c),
                              bbox_inches='tight',
                              pad_inches=0.15)

                fig2b.clear()
                matplotlib.rcParams.update({'font.size': 30})
                ax2b  = fig2b.add_axes([0.15,0.15,0.75,0.75])
                ax2b.errorbar(Z,RHO, yerr=RHOerr)
                #ax2a.plot(RHOfit,Z, linestyle='dashed')
                ax2b.set_ylabel('$\\rho$ (g/cm$^3$)')
                ax2b.set_xlabel('$z-z_{\mathrm{mid}}$')
                ax2b.set_xlim(Z_left-0.5, Z_right+0.5)#limits[dz][0], limits[dz][1])
                #ax2a.set_ylim(-8, 8)
                #ax2b.set_xlim(0, 7)
                fig2b.savefig('PLOTS_C/{}/densprof_err_{}_{}.pdf'.format(DIR,dz,c),bbox_inches='tight')


            except IOError:
                print 'File densprof.{} does not exist.'.format(f)

            # oxygen and hydrogen densities
            try:
                if oh_flag == 1:
                    print 'Reading in denso.{}'.format(f)
                    RHO_O, CY, CX = read_densprof_2d(DIR,f,'o')

                    # compute area density
                    CYmin, CYmax = np.min(CY), np.max(CY)
                    CXmin, CXmax = np.min(CX), np.max(CX)

                    den_max = np.max(RHO_O)
                    den_min = np.min(RHO_O)
                    fig8 = plt.figure(figsize=(fig_size_sq)) 
                    ax8  = fig8.add_axes([0.1,0.15,0.8,0.75])
                    #fig1.text(0.44, 0.025, '$y$ (\AA)', ha='center', va='center', fontsize=26)
                    #ax8.set_xlim(ylo,yhi)
                    #ax8.set_ylim(zlo,zhi)
                    ax8.set_xlabel('$y$ (\AA)')
                    ax8.set_ylabel('$x$ (\AA)')
                    #plt.tick_params(pad=7)
                    ctest=ax8.contourf(CY, CX, np.transpose(RHO_O), cmap=cm.RdBu, levels=np.linspace(den_min,den_max,500))
                    cbar = fig8.colorbar(ctest)
                    tick_locator = ticker.MaxNLocator(nbins=5)
                    cbar.locator = tick_locator
                    cbar.update_ticks()
                    fig8.savefig('PLOTS_C/{}/denso_2d_z{}_{}.png'.format(DIR,dz,c),dpi=500,bbox_inches='tight')
                    fig8.clear()

                    print 'Reading in densh.{}'.format(f)
                    RHO_H, CY, CX = read_densprof_2d(DIR,f,'h')

                    den_max = np.max(RHO_H)
                    den_min = np.min(RHO_H)
                    fig8 = plt.figure(figsize=(fig_size_sq)) 
                    ax8  = fig8.add_axes([0.1,0.15,0.8,0.75])
                    #fig1.text(0.44, 0.025, '$y$ (\AA)', ha='center', va='center', fontsize=26)
                    #ax8.set_xlim(ylo,yhi)
                    #ax8.set_ylim(zlo,zhi)
                    ax8.set_xlabel('$y$ (\AA)')
                    ax8.set_ylabel('$x$ (\AA)')
                    #plt.tick_params(pad=7)
                    ctest=ax8.contourf(CY, CX, np.transpose(RHO_H), cmap=cm.RdBu, levels=np.linspace(den_min,den_max,500))
                    cbar = fig8.colorbar(ctest)
                    fig8.savefig('PLOTS_C/{}/densh_2d_z{}_{}.png'.format(DIR,dz,c),bbox_inches='tight')
                    fig8.clear()


            except IOError:
                print 'File denso.{} does not exist.'.format(f)

            # pressure profile
            try:
                print 'Reading in stress.{}'.format(f)
                Z_P, P, Ptot, Pxy, Pxz, Pyz, Pxx, Pyy, Pzz, delz = stress_prof(DIR,f, 'None')
                MID_P, LEFT_P, RIGHT_P = mid_point(Z_P,P)
                Z_P = Z_P-Z_P[MID_P]
                P_T_tot = ((np.mean(Pxx[LEFT_P:RIGHT_P])+
                    np.mean(Pyy[LEFT_P:RIGHT_P]))/2)
                PTmax.append((np.max(Pxx)+np.max(Pyy))/(2*1000))
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
                fig5.savefig('PLOTS_C/{}/stressdens_{}_{}.pdf'.format(DIR,dz,c))
            except:
                pass

            # pmf
            try:
                print 'Reading in pmf.{}'.format(f)
                Z_PMF, PMF = pmf_prof(DIR,f)
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
                fig7.savefig('PLOTS_C/{}/pmf_{}.pdf'.format(DIR,dz),bbox_inches='tight')

            except IOError:
                print 'File pmf.{}_1 does not exist.'.format(f)


            # Diffusion (MSD)
            
            try:
                
                diff_bulk = 2.86624087647e-09
                # diff from msd
                diff_tot_y_msd_tmp, diff_tot_y_msd_err_tmp, t_msd, msd_ave = diff_msd(DIR,f)
                diff_tot_y_msd.append(diff_tot_y_msd_tmp*1e9)
                diff_tot_y_msd_err.append(diff_tot_y_msd_err_tmp*1e9)


            except IOError:
                print 'MSD not working'

            # viscosity profile from Diffusion (VACF)
            
            try:

                
                # wall viscosity
                print 'Reading in C_vv_y_{}_wall.dat'.format(f)
                C_vv_array = np.loadtxt("Water_Graphene/{}/C_vv_y_{}_wall.dat".format(DIR,f))
                times_wall =  C_vv_array[:,0]
                C_vv_ave_wall = C_vv_array[:,1]
                diff_wall_tmp = diffusion(C_vv_ave_wall, dt, time_conv, space_conv, 'C_vv_y')
                eta_diff_wall_tmp = eta_diff(diff_wall_tmp)
                eta_diff_wall.append(eta_diff_wall_tmp*1e3)

                # bulk shear viscosity
                print 'Reading in C_vv_{}_1_{}_z0_{}.dat'.format(f,steps,int(float(dz))-1)
                C_vv_array = np.loadtxt("Water_Graphene/{}/C_vv_{}_1_{}_z0_{}.dat".format(DIR,f, steps, int(float(dz))-1))
                times =  C_vv_array[:,0]
                C_vv_ave = C_vv_array[:,1]
                diff_tot_tmp = diffusion(C_vv_ave, dt, time_conv, space_conv, 'C_vv')
                eta_diff_tot_tmp = eta_diff(diff_tot_tmp)
                eta_diff_tot.append(eta_diff_tot_tmp*1e3)
                diff_tot.append(diff_tot_tmp*1e9)

                # bulk shear viscosity surface
                # diff from vacf
                diff_bulk = 2.86624087647e-09
                print 'Reading in C_vv_y_{}_1_{}_z0_{}.dat'.format(f, steps, int(float(dz))-1)
                C_vv_array = np.loadtxt("Water_Graphene/{}/C_vv_y_{}_1_{}_z0_{}.dat".format(DIR,f, steps, int(float(dz))-1))
                times =  C_vv_array[:,0]
                C_vv_ave = C_vv_array[:,1]
                diff_tot_y_tmp = diffusion(C_vv_ave, dt, time_conv, space_conv, 'C_vv_y')
                # eta from diff vacf
                eta_diff_tot_y_tmp = eta_diff(diff_tot_y_tmp)
                eta_diff_tot_y.append(eta_diff_tot_y_tmp*1e3)
                diff_tot_y.append(diff_tot_y_tmp*1e9)

                Z_init_vacf.append(float(dz))
                Hlist_vacf.append(float(dz)-3.19)

                print 'The viscosities from diffusion are:', eta_diff_wall_tmp*1e3, eta_diff_tot_tmp*1e3, eta_diff_tot_y_tmp*1e3

                matplotlib.rc('font', size=24)
                fig25.clear()
                ax25  = fig25.add_axes([0.15,0.15,0.75,0.75])
                ax25.plot(times_wall, C_vv_ave_wall/C_vv_ave_wall[0], linestyle='dashed',c=colours[0], label='wall')
                ax25.plot(times, C_vv_ave/C_vv_ave[0], c=colours[2], label='channel')
                ax25.set_ylabel('$\Psi_{\mathrm{u,u}}$ (t)', fontsize=28)
                ax25.set_xlabel('t (ps)', fontsize=28)
                ax25.set_xlim(0,0.6)
                ax25.set_ylim(-0.3,1)
                ax25.legend()
                fig25.savefig('PLOTS_C/{}/VACF_wall_channel_{}_{}.pdf'.format(DIR,dz,c),bbox_inches='tight')
                


                print dz, len(dZ)
                if len(dZ) == 2:
                    if count_tmp == 0:
                        ccount = 0
                    else:
                        ccount = 2*count+4
                    if dz=='6':
                        Hlabel = 'liquid'
                    elif dz=='7.5':
                        Hlabel='frozen'
                    else:
                        Hlabel='H = {} \AA'.format(dz)
                    ax26.plot(times, C_vv_ave/C_vv_ave[0], linestyle = ls[count_tmp], c=colours[ccount],label=Hlabel)

                    ax26a.plot(times, C_vv_ave/C_vv_ave[0], linestyle = ls[count_tmp], c=colours[ccount],label=Hlabel)
                    ax26a_i.plot(np.array(t_msd)/1000, msd_ave, linestyle = ls[count_tmp], c=colours[ccount])
                    ax26a_i.set_xlabel('t (ns)')
                    ax26a_i.set_ylabel('MSD (\AA$^2$)')
                else:
                    ax26.plot(times, C_vv_ave/C_vv_ave[0], label='H = {} \AA'.format(dz))
                
                matplotlib.rc('font', size=22)
                if dz == '30' and rho != 'RHO1':
                    # wall viscosity as fn of z
                    eta_wall_z = []
                    z_eta_wall = []
                    for i in range(2,10):
                        print 'Reading in C_vv_y_{}_1_1_{}_z{}_{}.dat'.format(f,steps,i,i+1)
                        C_vv_array = np.loadtxt("Water_Graphene/{}/C_vv_y_{}_1_1_{}_z{}_{}.dat".format(DIR,f,steps,i,i+1))
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
                print 'Reading in visc.{}'.format(f)
                eta_gk_tmp = visc_gk(DIR,f, height, 'etas','eff')
                eta_tot_gk.append(eta_gk_tmp*1e3)

                # parallel shear viscosity (GK)
                print 'Reading in visc.{}'.format(f)
                eta_gk_xy_tmp = visc_gk(DIR,f, height, 'etas', 'xy')
                if rho=='RHO1.1':
                    print 'Scaling eta'
                    eta_gk_xy_tmp *= (float(dz)-3.19)/(float(dz)-0.5)
                eta_xy_gk.append(eta_gk_xy_tmp*1e3)

                # bulk shear viscosity (GK, acf)
                print 'Reading in acfsv.{}'.format(f,int(float(dz))-1)
                eta_gk_acf_tmp = eta_gk_acf(DIR,f, height)
                #eta_tot_gk_acf.append(eta_gk_acf_tmp*1e3)

                print 'The viscosities from GK are:', eta_gk_tmp*1e3, eta_gk_xy_tmp*1e3, eta_gk_acf_tmp*1e3


                # bulk viscosity
                print 'Reading in visc.{}'.format(f)
                if rho == 'RHO1':
                    kappa_gk_tmp = visc_gk(DIR,f, height, 'etab','xy')
                else:
                    kappa_gk_tmp = visc_gk(DIR,f, height, 'etab','nvt')
                kappa_tot_gk.append(kappa_gk_tmp*1e3)
                kappa_eta_tot.append(kappa_gk_tmp/eta_gk_xy_tmp)




            except IOError:
                print 'File visc.{} or similar does not exist.'.format(f)


            # fluid wall 
            
            try:

                

                # reading in fric x and fric y
                print 'Reading in fric.{}'.format(f)
                fric_x_tmp = fric_gk(DIR,f, height, 'x')
                fric_y_tmp = fric_gk(DIR,f, height, 'y')
                fric_ave_tmp = 0.5*(fric_x_tmp+fric_y_tmp)
                fric.append(fric_ave_tmp*1e-4)
                Ls_fric_tmp = eta_gk_xy_tmp/fric_ave_tmp
                Ls_fric.append(Ls_fric_tmp*1e9)



                # work of adhesion
                wa_tmp = wa_log(DIR,f)
                wa.append(wa_tmp*1e3)
                print 'WA:', wa_tmp*1e3
                print 'Ds:', diff_tot_y_msd_tmp*1e9
                print 'eta:', eta_gk_xy_tmp*1e3
                Ls_wa_tmp = ((eta_gk_xy_tmp*32.8)/(height/2))*(diff_tot_y_msd_tmp/wa_tmp)
                Ls_wa.append(Ls_wa_tmp*1e9)

                print 'The slip lengths are:', Ls_fric_tmp*1e9, Ls_wa_tmp*1e9, 'nm.'

            except IOError:
                print 'File fric.{} or similar does not exist.'.format(f)



            # Calculate the overall system pressure from log file
            try:
                thermdata = read_log(DIR,f)
                Tarr, _p, Pcum, _p1, _dp, _v = props(thermdata)
                Pave = Pcum[-1]/10 # in MPa
                Pave_tot.append(Pave)
            except IOError:
                print 'File log.{} does not exist.'.format(f)



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

                # # normal corrlen is too long. The length we have is only 422 data points
                
                # # Calculate viscosity
                # Pxy_list, Pxz_list, Pyz_list = [], [], []
                # print Pxy.shape
                # for i in range(Pxy.shape[1]):
                #     Pxy_list.append(eta(f,Pxy[:-1,i],delz))
                #     Pxz_list.append(eta(f,Pxz[:-1,i],delz))
                #     Pyz_list.append(eta(f,Pyz[:-1,i],delz))
                # Paniso = np.mean(np.array([Pxy_list, Pxz_list, Pyz_list]), axis=0)
                # #ax6.plot(Z_P, Paniso)


                # fig6.clear()
                # ax6  = fig6.add_axes([0.15,0.15,0.75,0.75])
                
                # ax6.plot(Z_P, Pxz_list, linestyle='--', marker=markers[1], label='$xz$')
                # ax6.plot(Z_P, Pyz_list, linestyle=':', marker=markers[2], label='$yz$')
                # ax6.plot(Z_P, Pxy_list, linestyle='-', marker=markers[0], label='$xy$')
                # ax6.set_ylabel('$\eta$')
                # ax6.set_xlabel('$z-z_{\mathrm{mid}}$')
                # #ax2.set_xlim(Z_left, Z_right)
                # #ax6.set_xlim(-17.4, 17.4)
                # ax6.legend()
                # #ax4.set_ylim(0, 4)
                # fig6.savefig('PLOTS_C/SV_z{}.pdf'.format(dz))
                
            # except IOError:
            #     print 'File stresseta.{}_1 does not exist.'.format(f)
            count_tmp += 1
        ax30.set_xlim(Z_left-0.5, Z_right+0.5)
        ax30.set_ylabel('$\\rho$ (g/cm$^3$)')
        ax30.set_xlabel('$z-z_{\mathrm{mid}}$')
        ax30.legend()
        fig30.savefig('PLOTS_C/{}/dens_z{}.pdf'.format(DIR,dz),bbox_inches='tight')
        fig30.clear()
        ax30  = fig30.add_axes([0.15,0.15,0.75,0.75])

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
                        C_vv_array = np.loadtxt("Water_Graphene/{}/{}_{}_{}_{}.dat".format(DIR,p,f,c,r))
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
        Z_final, Hlist, Hlist_err = averaging(Z_init,Hlist)
        Z_final_vacf, Hlist_vacf, Hlist_vacf_err = averaging(Z_init_vacf,Hlist_vacf)
        Z_final_vacf, eta_diff_wall, eta_diff_wall_err = averaging(Z_init_vacf,eta_diff_wall)
        Z_final_vacf, eta_diff_tot, eta_diff_tot_err = averaging(Z_init_vacf,eta_diff_tot)
        Z_final_vacf, eta_diff_tot_y, eta_diff_tot_y_err = averaging(Z_init_vacf,eta_diff_tot_y)
        Z_final, eta_tot_gk, eta_tot_gk_err = averaging(Z_init,eta_tot_gk)
        Z_final, eta_xy_gk, eta_xy_gk_err = averaging(Z_init,eta_xy_gk)
        Z_final_vacf, diff_tot_y, diff_tot_y_err = averaging(Z_init_vacf,diff_tot_y)
        Z_final, diff_tot_y_msd, diff_tot_y_msd_err = averaging(Z_init,diff_tot_y_msd)
        Z_final, kappa_tot_gk, kappa_tot_gk_err = averaging(Z_init,kappa_tot_gk)
        Z_final, kappa_eta_tot, kappa_eta_tot_err = averaging(Z_init,kappa_eta_tot)
        Z_final, fric, fric_err = averaging(Z_init,fric)
        Z_final, wa, wa_err = averaging(Z_init,wa)
        Z_final, Ls_fric, Ls_fric_err = averaging(Z_init,Ls_fric)
        Z_final, Ls_wa, Ls_wa_err = averaging(Z_init,Ls_wa)
        Z_final, Pave_tot, Pave_err = averaging(Z_init,Pave_tot)
        Z_final, PTmax, PTmax_err = averaging(Z_init,PTmax)
        Z_final, rhoave, rhoave_err = averaging(Z_init,rhoave)
        Z_final, area_rho, area_rho_err = averaging(Z_init, area_rho)
        Z_final, rhoeff, rhoeff_err = averaging(Z_init, rhoeff)

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
        fric_err = [0]*len(fric)
        wa_err = [0]*len(wa)
        Ls_fric_err = [0]*len(Ls_fric)
        Ls_wa_err = [0]*len(Ls_wa)
        Pave_err = [0]*len(Pave_tot)
        PTmax_err = [0]*len(PTmax)
        rhoave_err = [0]*len(rhoave)
        area_rho_err = [0]*len(area_rho)
        Z_final = Z_init
        Z_final_vacf = Z_init_vacf

    print '\n#-----------------------Average values for fluid wall---------------#'
    print 'The average WA is:', np.mean(wa), '+/-', stats.sem(wa)
    print 'The average Ds is:', np.mean(diff_tot_y_msd), '+/-', stats.sem(diff_tot_y_msd)
    print 'The average Eta is:', np.mean(eta_xy_gk), '+/-', stats.sem(eta_xy_gk)
    print 'The average fric is:', np.mean(fric), '+/-', stats.sem(fric) 
    print 'The average Ls fric is:', np.mean(Ls_fric), '+/-', stats.sem(Ls_fric)

    print Z_final, rhoeff
    print eta_xy_gk, eta_xy_gk_err
    print kappa_tot_gk, kappa_tot_gk_err
    print kappa_eta_tot, kappa_eta_tot_err

    # in SI units
    WA_mean = np.mean(wa)*1e-3
    Ds_mean = np.mean(diff_tot_y_msd)*1e-9
    Eta_mean = np.mean(eta_xy_gk)*1e-3
    Ls_fric_mean = np.mean(Ls_fric)

    # Aspect ratio array
    LD = np.linspace(0,1000, 100)
    Ls_LD = (Eta_mean*LD)*(Ds_mean/WA_mean)*1e9

    # enhancement
    # fixed L length of 1000 A
    D = np.linspace(0,1000,1000)*1e-10
    L_tmp = 10000*1e-10
    enh_fric = 1+3*(Ls_fric_mean*1e-9/D)
    enh_LD = 1+(3/D**2)*(Eta_mean*L_tmp)*(Ds_mean/WA_mean)
    enh_LD_100 = 1+(3/D**2)*(Eta_mean*100*1e-10)*(Ds_mean/WA_mean)
    enh_LD_1000 = 1+(3/D**2)*(Eta_mean*1000*1e-10)*(Ds_mean/WA_mean)
    enh_LD_10000 = 1+(3/D**2)*(Eta_mean*10000*1e-10)*(Ds_mean/WA_mean)
    scale = WA_mean/(L_tmp*Ds_mean)*1e-12
    enh_LD_scaled = enh_LD*scale


    ###################################
    Z_pad_min = np.min(Z_final)-0.5
    Z_pad_max = np.max(Z_final)+0.5
    Z_final_pad = np.linspace(Z_pad_min, Z_pad_max, 10)


    # print to csv for latex use
    with open('DATA/{}/heights_{}.csv'.format(DIR, z_name), 'wb') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(("dz", "H"))
        writer.writerows(imap(lambda x, y: (y, round(x, 2)), Hlist, dZ))#izip(dZ,round(Hlist,2)))

    # print to csv for latex use
    with open('DATA/{}/rho_H_{}.csv'.format(DIR, z_name), 'wb') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(("rho", "H"))
        writer.writerows(imap(lambda x, y: (y, x), Z_final, rhoeff))#izip(dZ,round(Hlist,2)))


    with open('DATA/{}/viscosity_{}.csv'.format(DIR, z_name), 'wb') as csvfile:
        writer = csv.writer(csvfile)
        if DIR == 'RHO1':
            writer.writerow(("height", "rhoave", "diff", "shear", "bulk", "ratio"))
            writer.writerows(imap(
                             lambda a, b, c, d, e, f, g, h, i, j:
                             (round(a, 1),
                              "{}".format(round(j, 3)),
                              "{} \pm {}".format(round(b, 3), round(c, 3)),
                              "{} \pm {}".format(round(d, 2), round(e, 2)),
                              "{} \pm {}".format(round(f, 2), round(g, 2)),
                              "{} \pm {}".format(round(h, 1), round(i, 1))),
                             Z_final, diff_tot_y_msd, diff_tot_y_msd_err,
                             eta_xy_gk, eta_xy_gk_err,
                             kappa_tot_gk, kappa_tot_gk_err,
                             kappa_eta_tot, kappa_eta_tot_err,
                             rhoeff))
        elif DIR == 'RHO1.1':

            writer.writerow(("height", "rhoave", "diff", "shear"))
            writer.writerows(imap(
                             lambda a, b, c, d, e, f:
                             (round(a, 1),
                              "{}".format(round(f, 3)),
                              "{:6f} \pm {:6f}".format(round(b, 6), round(c, 6)),
                              "{} \pm {}".format(round(d, 2), round(e, 2))),
                             Z_final, diff_tot_y_msd, diff_tot_y_msd_err,
                             eta_xy_gk, eta_xy_gk_err,
                             rhoeff))

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
    #ax3.errorbar(Z_final, d_list, yerr=d_list_err, c=colours[0], linestyle='dashed', marker='D')
    ax3.set_ylabel('$D_{||}/D_{\mathrm{iso}}$')
    #ax3.set_ylim(0,1.3)
    ax3.yaxis.label.set_color(colours[0])
    for tl in ax3.get_yticklabels():   
        tl.set_color(colours[0])
    ax3b.errorbar(Z_final, rho_list, yerr=rho_list_err, c=colours[6], linestyle='dashed', marker='o')
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

    # fit curve for squeeze film equation
    #eta_squeeze = squeeze_visc(Z_final, k0)
    #eta_squeeze_fit = squeeze_visc(Z_final, 1.7)

    ax9.errorbar(Z_final_vacf, eta_diff_wall, yerr=eta_diff_wall_err, marker='o')
    ax9.plot(Z_final_vacf, len(Z_final_vacf)*[0.67], linestyle='dashed', label='Bulk')
    ax9.set_xlabel('H (\AA)')
    ax9.set_ylabel('$\eta_{\mathrm{wall}}$ (mPas)')
    ax9.legend()

    ax9a.set_yscale('log')
    ax9a.errorbar(Z_final_vacf, eta_diff_wall, yerr=eta_diff_wall_err, marker='o')
    ax9a.plot(Z_final_vacf, len(Z_final_vacf)*[0.67], linestyle='dashed', label='Bulk')
    ax9a.set_xlabel('H (\AA)')
    ax9a.set_ylabel('$\eta_{\mathrm{wall}}$ (mPas)')
    ax9a.legend()

    ax10.errorbar(Z_final_vacf, eta_diff_tot, yerr=eta_diff_tot_err, marker='o', label='isotropic')
    ax10.errorbar(Z_final_vacf, eta_diff_tot_y, yerr=eta_diff_tot_y_err, marker='o', label='xy')
    ax10.plot(Z_final_vacf, len(Z_final)*[0.67], linestyle='dashed', label='Bulk')
    ax10.set_xlabel('H (\AA)')
    ax10.set_ylabel('$\eta_{\mathrm{channel}}$ (mPas)')
    ax10.legend()

    ax10a.set_yscale('log')
    ax10a.errorbar(Z_final_vacf, eta_diff_tot, yerr=eta_diff_tot_err, marker='o', label='isotropic')
    ax10a.errorbar(Z_final_vacf, eta_diff_tot_y, yerr=eta_diff_tot_y_err, marker='o', label='xy')
    ax10a.plot(Z_final_vacf, len(Z_final_vacf)*[0.67], linestyle='dashed', label='Bulk')
    ax10a.set_xlabel('H (\AA)')
    ax10a.set_ylabel('$\eta_{\mathrm{channel}}$ (mPas)')
    ax10a.legend()

    #if '30' in dZ:
    #    ax11.plot(z_eta_wall, eta_wall_z, marker='o')
    #    ax11.plot(z_eta_wall, len(z_eta_wall)*[0.67], linestyle='dashed', label='Bulk')
    #    ax11.set_xlabel('z (\AA)')
    #    ax11.set_ylabel('$\eta_{\mathrm{xy}}$ (mPas)')
    #    ax11.legend()
    #    fig11.savefig('PLOTS_C/eta_diff_wall_z.pdf',bbox_inches='tight')

    ax12.errorbar(Z_final, eta_tot_gk, yerr=eta_tot_gk_err, marker='o', label='$\eta_{\mathrm{eff}}$')
    ax12.errorbar(Z_final, eta_xy_gk, yerr=eta_xy_gk_err, marker='o', label='$\eta_{\mathrm{||}}$')
    #ax12.plot(Hlist, eta_tot_gk_acf, marker='o', label='ACF')
    ax12.plot(Z_final, len(Z_final)*[0.67], linestyle='dashed', label='Bulk')
    ax12.set_xlabel('H (\AA)')
    ax12.set_ylabel('$\eta_{\mathrm{GK}}$ (mPas)')
    ax12.legend()

    ax12a.set_yscale('log')
    #ax12a.errorbar(Hlist, eta_tot_gk, yerr=eta_tot_gk_err, marker='o', label='$\eta_{\mathrm{eff}}$')
    ax12a.errorbar(Z_final, eta_xy_gk, yerr=eta_xy_gk_err, marker='o', c=colours[0])
    ax12a.plot(Z_final, len(Z_final)*[0.67], linestyle='dashed', c=colours[2], label='Bulk')
    #ax12a.plot(Hlist, eta_squeeze, linestyle='dotted', c=colours[4],label='Squeeze film eq ($k_0=2.7 \AA$)')
    #ax12a.plot(Hlist, eta_squeeze_fit, c=colours[6],linestyle='dotted', label='Squeeze film eq ($k_0=1.7 \AA$)')
    #ax12a.semilogy(Hlist, eta_tot_gk, marker='o', label='$\eta_{\mathrm{eff}}$')
    #ax12a.semilogy(Hlist, eta_xy_gk, marker='o', label='$\eta_{\mathrm{||}}$')
    #ax12a.semilogy(Hlist, len(Hlist)*[0.67], linestyle='dashed', label='Bulk')
    ax12a.set_xlabel('H (\AA)')
    ax12a.set_ylabel('$\eta$ (mPas)')
    ax12a.legend()

    #ax13.plot(Hlist, eta_tot_gk, marker='o', label='$\eta_{\mathrm{eff}}$')
    ax13.plot(1/np.array(diff_tot_y), eta_diff_tot_y, ls='None', marker='s', label='Stokes-Einstein')
    ax13.plot(1/np.array(diff_tot_y), eta_xy_gk, ls='None', marker='o', label='Green-Kubo')
    ax13.set_xlabel('$1/D_s$ ($10^{9}$s/m$^2$)')
    ax13.set_ylabel('$\eta$ (mPas)')
    #ax13.set_xlim(0,10)
    #ax13.set_ylim(0,25)
    ax13.set_xlim(0.2,2.0)
    ax13.set_ylim(0,5)
    ax13.legend()

    ax14.set_yscale('log')
    ax14.errorbar(Z_final, kappa_tot_gk, yerr=kappa_tot_gk_err, marker='o')
    ax14.plot(Z_final, len(Z_final)*[1.59], linestyle='dashed', label='Bulk')
    ax14.set_xlabel('H (\AA)')
    ax14.set_ylabel('$\kappa$ (mPas)')

    ax15.errorbar(Z_final, kappa_eta_tot, yerr=kappa_eta_tot_err, marker='o')
    ax15.plot(Z_final_pad, len(Z_final_pad)*[2.32], linestyle='dashed', label='Bulk')
    ax15.set_xlim(Z_pad_min, Z_pad_max)
    ax15.set_xlabel('H (\AA)')
    ax15.set_ylabel('$\kappa/\eta$ ')

    ax16.errorbar(Z_final_vacf, diff_tot_y, yerr=diff_tot_y_err, marker='o', label='VACF')
    ax16.errorbar(Z_final, diff_tot_y_msd, yerr=diff_tot_y_msd_err, marker='o', label='MSD')
    ax16.plot(Z_final, len(Z_final)*[2.866], linestyle='dashed', label='Bulk')
    ax16.set_xlabel('H (\AA)')
    ax16.set_ylabel('$D_s$ (10$^{-9}$ m$^2$/s)')
    ax16.legend()

    ax16a.set_yscale('log')
    ax16a.errorbar(Z_final_vacf, diff_tot_y, yerr=diff_tot_y_err, marker='o', label='VACF')
    ax16a.errorbar(Z_final, diff_tot_y_msd, yerr=diff_tot_y_msd_err, marker='o', label='MSD')
    ax16a.plot(Z_final, len(Z_final)*[2.866], linestyle='dashed', label='Bulk')
    ax16a.set_xlabel('H (\AA)')
    ax16a.set_ylabel('$D_s$ (10$^{-9}$ m$^2$/s)')
    ax16a.legend()

    dlayer = 2.75
    ax17.plot(Z_final, np.array(Z_final)/dlayer, marker='o', linestyle='dashed')
    ax17.plot(Z_final,len(Z_final)*[1], linestyle='dashed')
    ax17.plot(Z_final,len(Z_final)*[2], linestyle='dashed')
    ax17.plot(Z_final,len(Z_final)*[3], linestyle='dashed')
    ax17.plot(Z_final,len(Z_final)*[4], linestyle='dashed')
    ax17.plot(Z_final,len(Z_final)*[5], linestyle='dashed')
    ax17.plot(Z_final,len(Z_final)*[6], linestyle='dashed')
    ax17.plot(Z_final,len(Z_final)*[7], linestyle='dashed')
    ax17.set_xlabel('H (\AA)')
    ax17.set_ylabel('H/d$_{\mathrm{layer}}$')
    ax17.legend()

    Hlayer_sample = np.linspace(1,np.max(Z_final)/dlayer,100)
    ax17a.plot(np.array(Z_final)/dlayer, map(layers.get, dZ), marker='o', linestyle='dashed')
    ax17a.plot(Hlayer_sample, f1(Hlayer_sample, 1, 0), linestyle='dotted')
    ax17a.set_ylabel('No. of layers')
    ax17a.set_xlabel('H/d$_{\mathrm{layer}}$')
    ax17a.legend()

    ax18.errorbar(Z_final, fric, yerr=fric_err, marker='o', linestyle='dashed')
    ax18.set_xlabel('H (\AA)')
    ax18.set_ylabel('$\zeta$ ($10^4$ Ns/m$^3$)')
    ax18.legend()

    ax18b = ax18a.twinx()
    ax18a.set_xlabel('H (\AA)')
    ax18a.set_ylabel('$L_{s,\zeta}$ (nm)') 
    ax18a.yaxis.label.set_color(colours[0])
    for tl in ax18a.get_yticklabels():   
        tl.set_color(colours[0])
    fig18a.text(0.96, 0.5, '$\zeta$ ($10^4$ Ns/m$^3$)', color=colours[6], ha='center', va='center', fontsize=22,rotation=270)
    for tl in ax18b.get_yticklabels():   
        tl.set_color(colours[6])
    ax18b.errorbar(Z_final, fric, yerr=fric_err, marker='o', c=colours[6],linestyle='dashed',zorder=0)
    ax18a.errorbar(Z_final, Ls_fric, yerr=Ls_fric_err, marker='D', c=colours[0],linestyle='dashed',zorder=1)
    ax18a.set_zorder(ax18b.get_zorder()+1) # put ax in front of ax2 
    ax18a.patch.set_visible(False) # hide the 'canvas' 
    #ax18a.legend()

    ax19.errorbar(Z_final, wa, yerr=wa_err, marker='o', linestyle='dashed')
    ax19.set_xlabel('H (\AA)')
    ax19.set_ylabel('$W_A$ ($10^3$ N/m)')
    ax19.legend()

    ax19a.errorbar(Z_final, wa/(WA_mean*1e3), yerr=wa_err/(WA_mean*1e3), marker='o', c=colours[0], linestyle='dashed', label='$W_{A,\mathrm{ave}}$ = %.1f $\\times 10^3$N/m'%(WA_mean*1e3))
    ax19a.errorbar(Z_final, diff_tot_y_msd/(Ds_mean*1e9), yerr=diff_tot_y_msd_err/(Ds_mean*1e9), marker='D', c=colours[4], linestyle='dashed', label='$D_{s,\mathrm{ave}}$ = %.1f $\\times 10^{-9}$m$^2$/s'%(Ds_mean*1e9))
    ax19a.errorbar(Z_final, eta_xy_gk/(Eta_mean*1e3), yerr=eta_xy_gk_err/(Eta_mean*1e3), marker='s', c=colours[6], linestyle='dashed', label='$\eta_{\mathrm{ave}}$ = %.1f mPas'%(Eta_mean*1e3))
    ax19a.set_xlabel('H (\AA)')
    ax19a.set_ylabel('$X/X_{\mathrm{ave}}$')
    ax19a.set_ylim(0.6,1.4)
    ax19a.legend()

    ax20.errorbar(Z_final, Ls_fric, yerr=Ls_fric_err, marker='o', linestyle='dashed', label='$\zeta$')
    #ax20.errorbar(Hlist, Ls_wa, yerr=Ls_wa_err, marker='o', linestyle='dashed', label='$D_s/W_A$')
    ax20.set_xlabel('H (\AA)')
    ax20.set_ylabel('$L_s$ (nm)')
    ax20.legend()


    ax21.plot(LD, Ls_LD, marker='o', linestyle='dashed', label='$D_s/W_A$')
    ax21.plot(LD, len(LD)*[Ls_fric_mean],  marker='None', linestyle='dotted', label='$\zeta$')
    ax21.set_xlabel('L/D')
    ax21.set_ylabel('$L_s$ (nm)')
    ax21.legend()

    ax22.set_yscale('log')
    ax22.plot(D*1e10, enh_fric, marker='o', linestyle='dashed', label='$\zeta$')
    ax22.plot(D*1e10, enh_LD, marker='D', linestyle='dashed', label='$D_s/W_A$')
    ax22.set_xlabel('D (\AA)')
    ax22.set_ylabel('$\\varepsilon$')
    ax22.set_xlim(0,40)
    ax22.legend()


    ax22a.set_yscale('log')
    ax22a.set_xscale('log')
    #ax22.plot(D, enh_fric, marker='o', linestyle='dashed', label='$\zeta$')
    ax22a.plot(D*1e10, enh_LD_scaled, marker='o', linestyle='dashed', label='$D_s/W_A$')
    ax22a.set_xlabel('D (\AA)')
    ax22a.set_ylabel('$\\varepsilon \\times W_A/(LD_s) ~(10^{12}$Nsm$^{-4})$')
    ax22a.legend()

    ax22b.set_yscale('log')
    ax22b.plot(D*1e10, enh_LD_100, marker='o', linestyle='dashed', label='L = 10 nm')
    ax22b.plot(D*1e10, enh_LD_1000, marker='D', linestyle='dashed', label='L = 100 nm')
    ax22b.plot(D*1e10, enh_LD_10000, marker='s', linestyle='dashed', label='L = 1000 nm')
    ax22b.set_xlabel('D (\AA)')
    ax22b.set_ylabel('$\\varepsilon$')
    ax22b.set_xlim(0,40)
    ax22b.legend()

    ax23a = fig23.add_subplot(211)
    ax23a.set_yscale('log')
    ax23a.errorbar(Z_final, diff_tot_y_msd, yerr=diff_tot_y_msd_err, marker='o')
    ax23a.plot(Z_final_pad, len(Z_final_pad)*[2.866], linestyle='dashed')
    #ax23a.axes.get_xaxis().set_ticks([])
    ax23a.tick_params(labelbottom='off')   
    ax23a.set_ylabel('$D_s$ (10$^{-9}$ m$^2$/s)')
    ax23a.set_xlim(Z_pad_min, Z_pad_max)
    ax23a.legend()

    ax23b = fig23.add_subplot(212)
    ax23b.set_yscale('log')
    ax23b.errorbar(Z_final, eta_xy_gk, yerr=eta_xy_gk_err, marker='o')
    ax23b.plot(Z_final_pad, len(Z_final_pad)*[0.67], linestyle='dashed')
    ax23b.set_xlabel('H (\AA)')
    ax23b.set_ylabel('$\eta$ (mPas)', labelpad=20)
    ax23b.set_xlim(Z_pad_min, Z_pad_max)
    ax23b.legend()

    ax24a = fig24.add_subplot(311)
    ax24a.set_yscale('log')
    ax24a.errorbar(Z_final, diff_tot_y_msd, yerr=diff_tot_y_msd_err, marker='o')
    ax24a.plot(Z_final_pad, len(Z_final_pad)*[2.866], linestyle='dashed')
    #ax23a.axes.get_xaxis().set_ticks([])
    ax24a.tick_params(labelbottom='off')   
    ax24a.set_ylabel('$D_s$ (10$^{-9}$ m$^2$/s)')
    ax24a.set_xlim(Z_pad_min, Z_pad_max)
    ax24a.legend()

    ax24b = fig24.add_subplot(312)
    ax24b.set_yscale('log')
    ax24b.errorbar(Z_final, eta_xy_gk, yerr=eta_xy_gk_err, marker='o')
    ax24b.plot(Z_final_pad, len(Z_final_pad)*[0.67], linestyle='dashed')
    ax24b.set_ylabel('$\eta$ (mPas)', labelpad=20)
    ax24b.tick_params(labelbottom='off')   
    ax24b.set_xlim(Z_pad_min, Z_pad_max)
    ax24b.legend()

    ax24c = fig24.add_subplot(313)
    ax24c.set_yscale('log')
    ax24c.errorbar(Z_final, kappa_tot_gk, yerr=kappa_tot_gk_err, marker='o')
    ax24c.plot(Z_final_pad, len(Z_final_pad)*[1.59], linestyle='dashed')
    ax24c.set_xlabel('H (\AA)')
    ax24c.set_ylabel('$\kappa$ (mPas)', labelpad=20)
    ax24c.set_xlim(Z_pad_min, Z_pad_max)
    ax24c.legend()

    matplotlib.rc('font', size=28)
    ax26.set_xlabel('t (ps)', fontsize=28)
    ax26.set_ylabel('$\Psi_{\mathrm{u,u}}$ (t)', fontsize=28 )
    ax26.set_xlim(0,0.6)
    ax26.set_ylim(-0.3,1)
    ax26.legend()

    matplotlib.rc('font', size=24)
    ax26a.set_xlabel('t (ps)', fontsize=28)
    ax26a.set_ylabel('$C_{\mathrm{u,u}}$ (t)', fontsize=28 )
    ax26a.set_xlim(0,0.6)
    ax26a.set_ylim(-0.3,1)
    ax26a.legend()

    ax27.errorbar(Z_final, Pave_tot, yerr=Pave_err, marker='o')
    ax27.set_xlabel('H (\AA)')
    ax27.set_ylabel('P (MPa)')
    ax27.set_xlim(Z_pad_min, Z_pad_max)
    ax27.legend()

    ax27a.errorbar(Z_final, PTmax, yerr=PTmax_err, marker='o')
    ax27a.set_xlabel('H (\AA)')
    ax27a.set_ylabel('$P_{\mathrm{T, max}}$ (GPa)')
    ax27a.set_xlim(Z_pad_min, Z_pad_max)
    ax27a.legend()

    ax28.errorbar(Z_final, rhoave, yerr=rhoave_err, marker='o')
    ax28.set_xlabel('H (\AA)')
    ax28.set_ylabel('$\\rho$ (g/cm$^3$)')
    ax28.set_xlim(Z_pad_min, Z_pad_max)
    ax28.legend()

    ax29.set_yscale('log')
    ax29.errorbar(Pave_tot, kappa_tot_gk, yerr=kappa_tot_gk_err, xerr=Pave_err, marker='o', linestyle='None')
    ax29.set_ylabel('$\kappa$ (mPas)')
    ax29.set_xlabel('P (MPa)')
    #ax29.set_xlim(Z_pad_min, Z_pad_max)
    ax29.legend()

    ax31.errorbar(Z_final, area_rho, yerr=area_rho_err, marker='o')
    ax31.set_xlabel('H (\AA)')
    ax31.set_ylabel('$\\rho_{\mathrm{area}}$ (nm$^{-2}$)')
    ax31.set_xlim(Z_pad_min, Z_pad_max)
    ax31.legend()
    
    fig1.savefig('PLOTS_C/{}/diffusion_{}{}.pdf'.format(DIR,pre_name, z_name))
    fig2.savefig('PLOTS_C/{}/densprof_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig4.savefig('PLOTS_C/{}/stressprof_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig3.savefig('PLOTS_C/{}/diffdens_wall_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig9.savefig('PLOTS_C/{}/eta_diff_wall_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig9a.savefig('PLOTS_C/{}/eta_diff_wall_log_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig10.savefig('PLOTS_C/{}/eta_diff_channel_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig10a.savefig('PLOTS_C/{}/eta_diff_channel_log_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig12.savefig('PLOTS_C/{}/eta_gk_channel_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig12a.savefig('PLOTS_C/{}/eta_gk_channel_log_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig13.savefig('PLOTS_C/{}/eta_vs_diff_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig14.savefig('PLOTS_C/{}/kappa_gk_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig15.savefig('PLOTS_C/{}/kappa_eta_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig16.savefig('PLOTS_C/{}/diff_msd_vacf_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig16a.savefig('PLOTS_C/{}/diff_msd_vacf_log_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig17.savefig('PLOTS_C/{}/height_layer_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig17a.savefig('PLOTS_C/{}/height_layer2_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig18.savefig('PLOTS_C/{}/fric_gk_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig18a.savefig('PLOTS_C/{}/fric_comb_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig19.savefig('PLOTS_C/{}/wa_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig19a.savefig('PLOTS_C/{}/WADseta_ave_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig20.savefig('PLOTS_C/{}/Ls_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig21.savefig('PLOTS_C/{}/Ls_aspect_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig22.savefig('PLOTS_C/{}/enh_fixedL{}_{}.pdf'.format(DIR,L_tmp*1e10,z_name),bbox_inches='tight')
    fig22a.savefig('PLOTS_C/{}/enh_fixedL_scaled_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig22b.savefig('PLOTS_C/{}/enh_DsWA_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig23.savefig('PLOTS_C/{}/eta_diff_log_sub_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig24.savefig('PLOTS_C/{}/eta_diff_kappa_log_sub_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig26.savefig('PLOTS_C/{}/VACF_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig26a.savefig('PLOTS_C/{}/VACF_MSD_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig27.savefig('PLOTS_C/{}/Pave_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig27a.savefig('PLOTS_C/{}/PTmax_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig28.savefig('PLOTS_C/{}/rhoave_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig29.savefig('PLOTS_C/{}/Pave_vs_kappa_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    fig31.savefig('PLOTS_C/{}/area_dens_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    print '\n#-----------------------Done------------------\n#'

    return

if __name__ == "__main__":
    sys.path.insert(0,'/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()