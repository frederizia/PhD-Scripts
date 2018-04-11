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
from class_bulk_props import bulk_properties

def averaging3(k,v,e):
    averages = {}
    errors = {}
    counts = {}
    for name, value, error in zip(k, v, e):
        if name in averages:
            averages[name].append(value)
            errors[name].append(error)
            counts[name] += 1
        else:
            averages[name] = [value]
            errors[name] = [error]
            counts[name] = 1
    for name in averages:
        ave_array = np.array(averages[name])
        ave_val = np.mean(ave_array)
        err_array = np.array(errors[name])
        ave_error = np.sqrt(np.sum(np.square(err_array)))/counts[name]
        averages[name] = ave_val
        errors[name] = ave_error 

    k = averages.keys()
    v = averages.values()
    e = errors.values()
    return k,v,e

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
    parser.add_argument('-o', '--opt', required=False, type=str, action='store', default=[0,100],
                        help='msd or vacf')
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
    opt = args.opt

    rho = 'RHO1'
    avg = 0


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

    if rho == 'RHO1':
        steps = 3000
        DIR = 'RHO1'
        layers =  {'6':1, '6.2': 1, '6.8':1, '7':1, '7.2':2, '7.5': 2, '7.8':2, \
            '8':2, '8.2':2, '8.5':2 ,'8.8':2, '9':2, '9.2':2, '9.5':2, '9.8':2,'10':2, \
            '10.2':2, '10.5':2, '10.8':2, '11':2, '11.5':3, '12':3, '12.5':4,'13':3, \
            '13.4':3,'14':4, '14.5':4,\
            '15':4, '16':4, '17':4, '20':5, '22':6, '25':7,\
            '27':8, '30':9, '33':10, '35':10, '40':12}

    # ------------------- Initialise figures -----------------------

    # eta GK x Ds confined
    fig1 = plt.figure(figsize=fig_size_sq)
    ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])

    # eta GK vs Ds
    fig13 = plt.figure(figsize=fig_size_sq)
    ax13  = fig13.add_axes([0.15,0.15,0.75,0.75])


    # Properties independent of dimensionality
    
    print '\n#-----------------------Analysing density, pressure and viscosities------------------\n#'
    RHOdict = {}
    Hlist, Hlist_vacf = [], []
    Z_init, Z_init_vacf = [], []
    z_name = 'z_'
    eta_diff_wall, eta_diff_tot, eta_diff_tot_y, eta_diff_tot_y_msd, eta_tot_gk, eta_tot_gk_acf, eta_xy_gk = [], [], [], [], [], [], []
    kappa_tot_gk, kappa_eta_tot = [], []
    diff_tot_y, diff_tot, diff_tot_y_msd, diff_tot_y_msd_err = [], [], [], []
    for dz in dZ:
        z_name += dz
        for c in configs:
            Z_init.append(float(dz))
            count = 1
            count_list, diff_list = [], []
            f = 'spce_T298_z{}_eps1.0_{}'.format(dz,c)
            print '\n#-----------------------dz = {}, DEN = {}------------------\n#'.format(dz,c)

            # density profile
            try:
                print 'Reading in densprof.{}'.format(f)
                Z, RHO = read_densprof(DIR,f)
                MID, LEFT, RIGHT = mid_point(Z,RHO)
                # store data for use in diffusion plots
                RHOdict[dz] = (RHO, LEFT)
                Z = Z-Z[MID]
                Z_left = Z[LEFT]
                Z_right = Z[RIGHT]
                Hlist.append(Z_right-Z_left)
                height = Z_right-Z_left
                print 'Channel height:', height


            except IOError:
                print 'File densprof.{} does not exist.'.format(f)

 
            # viscosity profile from Diffusion (VACF)
            
            try:

                
                # wall viscosity
                print 'Reading in C_vv_y_{}_wall.dat'.format(f)
                C_vv_array = np.loadtxt("Water_Graphene/{}/C_vv_y_{}_wall.dat".format(DIR,f))
                times =  C_vv_array[:,0]
                C_vv_ave = C_vv_array[:,1]
                diff_wall_tmp = diffusion(C_vv_ave, dt, time_conv, space_conv, 'C_vv_y')
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

            # Diffusion (MSD)
            
            try:
                print 'Calculating MSD for {}'.format(f)
                diff_bulk = 2.86624087647e-09
                # diff from msd
                diff_tot_y_msd_tmp, diff_tot_y_msd_err_tmp = diff_msd(DIR,f)
                diff_tot_y_msd.append(diff_tot_y_msd_tmp*1e9)
                diff_tot_y_msd_err.append(diff_tot_y_msd_err_tmp*1e9)

                eta_diff_tot_y_msd_tmp = eta_diff(diff_tot_y_msd_tmp)
                eta_diff_tot_y_msd.append(eta_diff_tot_y_msd_tmp*1e3)

            except IOError:
                print 'MSD not working'
            # viscosity profile from GK
            
            try:
                # effective shear viscosity (GK)
                print 'Reading in visc.{}'.format(f)
                eta_gk_tmp = visc_gk(DIR,f, height, 'etas','eff')
                eta_tot_gk.append(eta_gk_tmp*1e3)

                # parallel shear viscosity (GK)
                print 'Reading in visc.{}'.format(f)
                eta_gk_xy_tmp = visc_gk(DIR,f, height, 'etas', 'xy')
                eta_xy_gk.append(eta_gk_xy_tmp*1e3)


                print 'The viscosities from GK are:', eta_gk_tmp*1e3, eta_gk_xy_tmp*1e3 #, eta_gk_acf_tmp*1e3

            except IOError:
                print 'File visc.{} or similar does not exist.'.format(f)



    # average over input configurations
    if avg == 1:
        Z_final_vacf, eta_diff_wall, eta_diff_wall_err = averaging(Z_init_vacf,eta_diff_wall)
        Z_final_vacf, eta_diff_tot, eta_diff_tot_err = averaging(Z_init_vacf,eta_diff_tot)
        Z_final_vacf, eta_diff_tot_y, eta_diff_tot_y_err = averaging(Z_init_vacf,eta_diff_tot_y)
        Z_final_vacf, eta_diff_tot_y_msd, eta_diff_tot_y_msd_err = averaging(Z_init_vacf,eta_diff_tot_y_msd)
        Z_final, eta_tot_gk, eta_tot_gk_err = averaging(Z_init,eta_tot_gk)
        Z_final, eta_xy_gk, eta_xy_gk_err = averaging(Z_init,eta_xy_gk)
        Z_final_vacf, diff_tot_y, diff_tot_y_err = averaging(Z_init_vacf,diff_tot_y)
        Z_final, diff_tot_y_msd, diff_tot_y_msd_err = averaging(Z_init,diff_tot_y_msd)

    else:
        eta_diff_wall_err = [0]*len(eta_diff_wall)
        eta_diff_tot_err = [0]*len(eta_diff_tot)
        eta_diff_tot_y_err = [0]*len(eta_diff_tot_y)
        eta_tot_gk_err = [0]*len(eta_tot_gk)
        eta_xy_gk_err = [0]*len(eta_xy_gk)
        diff_tot_y_err = [0]*len(diff_tot_y)
        Z_final = Z_init



    m = 'spce'
    fluid = 'Water'
    temp = '298'
    den = ['5','6']
    press = [1, 101, 201, 301, 401, 501, 601, 701, 801]
    print '\n\n++++++++++++++++++++++++ %s +++++++++++++++++++++\n' %m

    shear_dat = []
    shear_msd_dat = []
    rho_dat_init = []
    shear_vacf_dat = []
    shear_err = []
    diff_dat = []
    diff_vacf_dat = []

    for p in press:
        for i in den:
            print 'Evaluating %s, DEN=%s and P=%s.'%(m,i,p)
            try:
                results = bulk_properties(fluid,m,temp,p,i)
                rho_dat_init.append(results.rho())
                shear_dat.append(results.shear())
                shear_msd_dat.append(results.shear_diff('msd', 'None'))
                shear_vacf_dat.append(results.shear_diff('vacf', steps))
                diff_dat.append(results.diff()[0]*1e9)
                diff_vacf_dat.append(results.diff2(steps=steps)*1e9)
                #print results.rho(), results.bulk2()[0], results.bulk2()[1]
                print 'Data collection successful.\n'
            except:
                print 'ERROR: Maybe no data for %s, DEN=%s and P=%s exists.\n'%(m,c,p)


    if avg == 1:
        rho_dat, shear_msd_dat, shear_msd_err   = averaging(rho_dat_init, shear_msd_dat)
        rho_dat, shear_dat, shear_err   = averaging(rho_dat_init, shear_dat)
        rho_dat, shear_vacf_dat, shear_vacf_err   = averaging(rho_dat_init, shear_vacf_dat)#, shear_err)
        rho_dat, diff_vacf_dat, diff_vacf_err   = averaging(rho_dat_init, diff_vacf_dat)
        rho_dat, diff_dat, diff_err   = averaging(rho_dat_init, diff_dat)
    else:
        rho_dat = rho_dat_init


    # add the two sets

    if opt=='msd':
        diff_confined = np.array(diff_tot_y_msd)
        diff_final = list(diff_tot_y_msd) + list(diff_dat)
        shear_se_final = list(eta_diff_tot_y_msd) + list(shear_msd_dat)

    elif opt == 'vacf':
        diff_confined = np.array(diff_tot_y)
        diff_final = list(diff_tot_y) + list(diff_vacf_dat)
        shear_se_final = list(eta_diff_tot_y) + list(shear_vacf_dat)

    print len(diff_final), len(shear_se_final)



    shear_gk_final = list(eta_xy_gk) + list(shear_dat)
    #shear_gk_err_final = eta_xy_gk_err + list(shear_err)

    # determine the product
    shear_gk_confined = np.array(eta_xy_gk)
    shear_diff_confined = shear_gk_confined*diff_confined



    # average
    Z_confined, shear_diff_confined, shear_diff_err = averaging(Z_init, shear_diff_confined)
    shear_diff_conf_avg = np.mean(shear_diff_confined[np.where(Z_confined>=9)])
    shear_diff_conf_err = stats.sem(shear_diff_confined[np.where(Z_confined>=9)])
    shear_diff_dev = abs(shear_diff_confined-shear_diff_conf_avg)
    print 'The average product eta x Ds is:', shear_diff_conf_avg, '+/-', shear_diff_conf_err
    


    xfit_gk, shear_fit_gk, slope_gk, slope_gk_err = straight_fit(1e9/np.array(diff_final), 1e-3*np.array(shear_gk_final), 1e9/np.max(diff_final), 1e9/np.min(diff_final))
    xfit_se, shear_fit_se, slope_se, slope_se_err = straight_fit(1e9/np.array(diff_final), 1e-3*np.array(shear_se_final), 1e9/np.max(diff_final), 1e9/np.min(diff_final))

    # effect diameter
    kB = 1.38*1e-23
    T=int(temp)
    convert = 1e10
    factor = (kB*T)/(3*np.pi)

    alpha_se = convert*factor/slope_se # should be 1.7
    alpha_se_err = convert*factor*(-1)*(slope_se_err/slope_se**2)
    alpha_gk = convert*factor/slope_gk
    alpha_gk_err = convert*factor*(-1)*(slope_gk_err/slope_gk**2)
    xfit_gk, shear_fit_gk = 1e-9*np.array(xfit_gk), 1e3*np.array(shear_fit_gk)
    xfit_se, shear_fit_se = 1e-9*np.array(xfit_se), 1e3*np.array(shear_fit_se)
    print 'The effective diameters are:', alpha_se, '+/-',alpha_se_err, 'A (Diff) and', alpha_gk, '+/-',alpha_gk_err, 'A (GK)'

    print '\n#-----------------------Final plotting------------------\n#'

    # eta x D product
    ax1.set_yscale('log')
    ax1.plot([np.min(Z_confined)-1,np.max(Z_confined)+1],[shear_diff_conf_avg,shear_diff_conf_avg], linestyle='dashed', marker='None')
    cm = plt.cm.get_cmap('RdYlBu_r')
    sc = ax1.scatter(Z_confined, shear_diff_confined,c=shear_diff_dev, vmin=0, vmax=1.5, s=70, cmap=cm)
    #create colorbar according to the scatter plot
    clb = plt.colorbar(sc)
    #create errorbar plot and return the outputs to a,b,c
    a,b,c = ax1.errorbar(Z_confined, shear_diff_confined, yerr=shear_diff_err, marker='None', linestyle='None', capsize=0,zorder=0)
    #convert time to a color tuple using the colormap used for scatter
    err_color = clb.to_rgba(shear_diff_dev)
    #adjust the color of c[0], which is a LineCollection, to the colormap
    c[0].set_color(err_color)

    ax1.set_xlabel('$\Delta z$ (\AA)')
    ax1.set_xlim(np.min(Z_confined)-1,np.max(Z_confined)+1)
    ax1.set_ylabel('$\eta\\times D_s$ ($10^6$ m$^2$Pa)')

    #ax13.plot(Hlist, eta_tot_gk, marker='o', label='$\eta_{\mathrm{eff}}$')
    #ax13.plot(1/np.array(diff_final), shear_se_final, ls='None', marker='s', c=colours[0], label='Stokes-Einstein')
    ax13.plot(1/np.array(diff_final), shear_gk_final, ls='None', marker='o', c=colours[0], label='Green-Kubo')
    ax13.plot(xfit_gk, shear_fit_gk, linestyle='dashed',c=colours[0],label='fit')
    ax13.plot(xfit_se, shear_fit_se, linestyle='dotted',c=colours[7], label='$\\alpha = 1.7$ \AA')
    ax13.set_xlabel('$1/D_s$ ($10^{9}$s/m$^2$)')
    ax13.set_ylabel('$\eta$ (mPas)')
    #ax13.set_xlim(0,10)
    #ax13.set_ylim(0,25)
    #ax13.set_xlim(0.2,1)
    #ax13.set_ylim(0.5,2.5)
    ax13.legend()

    
    fig1.savefig('PLOTS_C/{}/stokes_einstein_product_{}_{}.pdf'.format(DIR,opt,z_name),bbox_inches='tight')
    fig13.savefig('PLOTS_C/{}/stokes_einstein_{}_{}.pdf'.format(DIR,opt,z_name),bbox_inches='tight')
    
    print '\n#-----------------------Done------------------\n#'
    return

if __name__ == "__main__":
    sys.path.insert(0,'/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()