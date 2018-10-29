#! /usr/bin/env python
# Plots for confined water data

import argparse
import matplotlib.pyplot as plt
import matplotlib
import sys
import csv
import numpy as np
from scipy.integrate import simps
from itertools import izip, imap
from confined_tools import *
from matplotlib import cm
from matplotlib import ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def GetArgs():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', '--width', nargs=1, required=False,
                        default=6.8, action='store',
                        help='Width of channel')
    parser.add_argument('-r', '--rhoflag', nargs='+', required=False,
                        default='liquid', action='store',
                        help='Density flag')
    parser.add_argument('-s', '--start', nargs='+', required=False,
                        default=1, action='store',
                        help='Starting configurations')
    parser.add_argument('-rho', '--rho', required=False, default='RHO1',
                        action='store',
                        help='Average channel density')
    args = parser.parse_args()
    return args


def main():

    args = GetArgs()

    dz = args.width[0]
    configs = args.start
    rho = args.rho
    rhoflag = args.rhoflag


    # Definitions
    fig_size_long = (14, 5)
    fig_size_sq = (9, 7)
    fig_size_sq2 = (7, 9)

    markers = ['o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P', 'X', '8', 'p',
               'o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P', 'X', '8', 'p',
               'o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P', 'X', '8', 'p']
    colours = ['#313695', '#4575b4', '#74add1',
               '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026',
               '#313695', '#4575b4', '#74add1',
               '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026',
               '#313695', '#4575b4', '#74add1',
               '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']
    LABEL = {'C_vv': '$xyz$', 'C_vv_y': '$xy$'}

    ls = ['--', '-', '--', '-.', ':', '-', '--', '-.', ':']

    # ------------------- Initialise figures -----------------------

    # rho (z)
    fig1 = plt.figure(figsize=fig_size_sq)
    ax1 = fig1.add_axes([0.15, 0.15, 0.75, 0.75])

    # P(z)
    fig2 = plt.figure(figsize=fig_size_long)
    ax2 = fig2.add_axes([0.15, 0.15, 0.75, 0.75])

    # rho vs PT
    fig3a = plt.figure(figsize=fig_size_sq)
    ax3a = fig3a.add_axes([0.15, 0.15, 0.75, 0.75])

    # rho vs P
    fig3b = plt.figure(figsize=fig_size_sq)
    ax3b = fig3b.add_axes([0.15, 0.15, 0.75, 0.75])

    # rho vs area density
    fig4 = plt.figure(figsize=fig_size_sq)
    ax4 = fig4.add_axes([0.15, 0.15, 0.75, 0.75])

    # Ds, eta, kappa comb
    fig5 = plt.figure(figsize=(9, 9))

    # kappa vs P
    fig6 = plt.figure(figsize=fig_size_sq)
    ax6 = fig6.add_axes([0.15, 0.15, 0.75, 0.75])

    # kappa/eta
    fig9 = plt.figure(figsize=fig_size_sq)
    ax9 = fig9.add_axes([0.15, 0.15, 0.75, 0.75])

    # R vs rho
    fig10 = plt.figure(figsize=fig_size_sq)
    ax10 = fig10.add_axes([0.15, 0.15, 0.75, 0.75])

    steps = 3000
    DIR = 'RHOtest'

    # Properties independent of dimensionality

    print '\n#-----------------------Analysing density, pressure and viscosities------------------\n#'
    RHOdict = {}
    R_init, R_init_vacf = [], []
    oh_flag = 0
    eta0 = 0.67  # mPas
    Pave_tot, PTmax, rhoave, area_rho = [], [], [], []
    RHOeff_list = []
    eta_tot_gk, eta_tot_gk_acf, eta_xy_gk = [], [], []
    kappa_tot_gk, kappa_eta_tot = [], []
    diff_tot_y_msd, diff_tot_y_msd_err = [], []
    count_tmp = 0
    r_name = 'r'
    for r in rhoflag:
        r_name += r
        for c in configs:
            R_init.append(float(r))
            count = 1
            count_list, diff_list = [], []
            f = 'spce_T298_z{}_r{}_eps1.0_{}'.format(dz, r, c)
            print '\n#-----------------------dz = {}, r = {}, DEN = {}------------------\n#'.format(dz, r, c)

            # area density
            try:
                area_xy = geometry(DIR, f)[0]

                # find number of atoms
                file = open('Water_Graphene/{}/log.{}'.format(DIR, f), 'r')

                data = file.read()
                data_lines = data.split('\n')

                DATA = []
                flag = 0
                for i in range(len(data_lines)-1):
                    if (data_lines[i].split() != []
                            and data_lines[i].split()[0] == 'group'
                            and data_lines[i].split()[1] == 'oxygen'):
                        numwater = int(data_lines[i+1].split()[0])
                        break
                area_dens = numwater/(area_xy*1e-2)
                area_rho.append(area_dens)

            except IOError:
                print 'File log.{} does not exist.'.format(f)

            # average effective density
            try:
                # open density file
                file = open('Water_Graphene/{}/dens.{}'.format(DIR, f), 'r')
                _, __, data = (file.readline(), file.readline(),
                               file.readline())

                RHOeff = float(data.split()[1])
                RHOeff_list.append(RHOeff)
                print 'Effective density:', RHOeff
            except IOError:
                print 'File dens.{} does not exist.'.format(f)

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
                height = Z_right-Z_left
                label = '$\rrho = {}$ g/cm$^3$'.format(RHOeff)

                fig1.clear()
                ax1 = fig1.add_axes([0.15, 0.15, 0.75, 0.75])
                ax1.plot(Z, RHO, label=label)
                ax1.set_ylabel('$\\rho$ (g/cm$^3$)')
                ax1.set_xlabel('$z-z_{\mathrm{mid}}$')
                ax1.set_xlim(Z_left-0.5, Z_right+0.5)
                fig1.savefig('PLOTS_C/{}/densprof_z{}_r{}_{}.pdf'.format(DIR, dz, r, c),
                             bbox_inches='tight')

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
                    fig7 = plt.figure(figsize=(fig_size_sq))
                    ax7 = fig7.add_axes([0.1, 0.15, 0.8, 0.75])
                    ax7.set_xlabel('$y$ (\AA)')
                    ax7.set_ylabel('$x$ (\AA)')
                    ctest = ax7.contourf(CY, CX, np.transpose(RHO_O),
                                         cmap=cm.RdBu,
                                         levels=np.linspace(den_min,
                                                            den_max,
                                                            500))
                    cbar = fig7.colorbar(ctest)
                    tick_locator = ticker.MaxNLocator(nbins=5)
                    cbar.locator = tick_locator
                    cbar.update_ticks()
                    fig7.savefig('PLOTS_C/{}/denso_2d_z{}_{}.png'.format(DIR,
                                                                         dz,
                                                                         c),
                                 dpi=500, bbox_inches='tight')
                    fig7.clear()

                    print 'Reading in densh.{}'.format(f)
                    RHO_H, CY, CX = read_densprof_2d(DIR, f, 'h')

                    den_max = np.max(RHO_H)
                    den_min = np.min(RHO_H)
                    fig8 = plt.figure(figsize=(fig_size_sq))
                    ax8 = fig8.add_axes([0.1, 0.15, 0.8, 0.75])

                    ax8.set_xlabel('$y$ (\AA)')
                    ax8.set_ylabel('$x$ (\AA)')
                    ctest = ax8.contourf(CY, CX, np.transpose(RHO_H),
                                         cmap=cm.RdBu,
                                         levels=np.linspace(den_min,
                                                            den_max,
                                                            500))
                    cbar = fig8.colorbar(ctest)
                    fig8.savefig('PLOTS_C/{}/densh_2d_z{}_{}.png'.format(DIR,
                                                                         dz,
                                                                         c),
                                 bbox_inches='tight')
                    fig8.clear()

            except IOError:
                print 'File denso.{} does not exist.'.format(f)

            # pressure profile
            try:
                print 'Reading in stress.{}'.format(f)
                (Z_P, P, Ptot, Pxy, Pxz,
                 Pyz, Pxx, Pyy, Pzz, delz) = stress_prof(DIR, f, 'None')
                MID_P, LEFT_P, RIGHT_P = mid_point(Z_P, P)
                Z_P = Z_P-Z_P[MID_P]
                P_T_ave = ((Pxx[LEFT_P:RIGHT_P] +
                           Pyy[LEFT_P:RIGHT_P])/2)
                P_T_tot = ((np.mean(Pxx[LEFT_P:RIGHT_P]) +
                           np.mean(Pyy[LEFT_P:RIGHT_P]))/2)
                PTmax.append((np.max(Pxx)+np.max(Pyy))/(2*1000))
                label = '$\rrho = {}$ g/cm$^3$'.format(RHOeff)

                fig2.clear()
                ax2 = fig2.add_axes([0.15, 0.15, 0.75, 0.75])
                ax2.plot(Z_P[LEFT_P:RIGHT_P], P_T_ave, label=label)
                ax2.set_ylabel('$P_T$ (MPa)')
                ax2.set_xlabel('$z-z_{\mathrm{mid}}$')
                ax2.set_xlim(Z_P[LEFT_P]-0.5, Z_P[RIGHT_P]+0.5)
                fig2.savefig('PLOTS_C/{}/pt_z{}_r{}_{}.pdf'.format(DIR,
                                                                   dz,
                                                                   r,
                                                                   c),
                             bbox_inches='tight')

            except IOError:
                print 'File stress.{} does not exist.'.format(f)

            # Diffusion (MSD)
            try:
                diff_bulk = 2.86624087647e-09
                # diff from msd
                (diff_tot_y_msd_tmp, diff_tot_y_msd_err_tmp,
                    t_msd, msd_ave) = diff_msd(DIR, f)
                diff_tot_y_msd.append(diff_tot_y_msd_tmp*1e9)
                diff_tot_y_msd_err.append(diff_tot_y_msd_err_tmp*1e9)

            except IOError:
                print 'MSD not working'

            # viscosity profile from GK
            try:

                # effective shear viscosity (GK)
                print 'Reading in visc.{}'.format(f)
                eta_gk_tmp = visc_gk(DIR, f, height, 'etas', 'eff')
                eta_tot_gk.append(eta_gk_tmp*1e3)

                # parallel shear viscosity (GK)
                print 'Reading in visc.{}'.format(f)
                eta_gk_xy_tmp = visc_gk(DIR, f, height, 'etas', 'xy')
                eta_xy_gk.append(eta_gk_xy_tmp*1e3)

                # bulk shear viscosity (GK, acf)
                print 'Reading in acfsv.{}'.format(f, int(float(dz))-1)
                eta_gk_acf_tmp = eta_gk_acf(DIR, f, height)
                # eta_tot_gk_acf.append(eta_gk_acf_tmp*1e3)

                print 'The viscosities from GK are:', (eta_gk_tmp*1e3,
                                                       eta_gk_xy_tmp*1e3,
                                                       eta_gk_acf_tmp*1e3)

                # bulk viscosity
                print 'Reading in visc.{}'.format(f)
                kappa_gk_tmp = visc_gk(DIR, f, height, 'etab', 'xy')
                kappa_tot_gk.append(kappa_gk_tmp*1e3)
                kappa_eta_tot.append(kappa_gk_tmp/eta_gk_xy_tmp)

            except IOError:
                print 'File visc.{} or similar does not exist.'.format(f)

            # Calculate the overall system pressure from log file
            try:
                thermdata = read_log(DIR, f)
                Tarr, _p, Pcum, _p1, _dp, _v = props(thermdata)
                Pave = Pcum[-1]/10  # in MPa
                Pave_tot.append(Pave)
            except IOError:
                print 'File log.{} does not exist.'.format(f)

            count_tmp += 1

    # average over input configurations
    if len(configs) > 1:
        rhoeff, eta_tot_gk, eta_tot_gk_err = averaging(RHOeff_list, eta_tot_gk)
        rhoeff, eta_xy_gk, eta_xy_gk_err = averaging(RHOeff_list, eta_xy_gk)
        rhoeff, diff_tot_y_msd, diff_tot_y_msd_err = averaging(RHOeff_list, diff_tot_y_msd)
        rhoeff, kappa_tot_gk, kappa_tot_gk_err = averaging(RHOeff_list, kappa_tot_gk)
        rhoeff, kappa_eta_tot, kappa_eta_tot_err = averaging(RHOeff_list, kappa_eta_tot)
        rhoeff, Pave_tot, Pave_err = averaging(RHOeff_list, Pave_tot)
        rhoeff, PTmax, PTmax_err = averaging(RHOeff_list, PTmax)
        rhoeff, rhoave, rhoave_err = averaging(RHOeff_list, rhoave)
        rhoeff, R_final, R_final_err = averaging(RHOeff_list, R_init)
        rhoeff, area_rho, area_rho_err = averaging(RHOeff_list, area_rho)

    else:
        eta_tot_gk_err = [0]*len(eta_tot_gk)
        eta_xy_gk_err = [0]*len(eta_xy_gk)
        diff_tot_y_msd_err = [0]*len(diff_tot_y)
        kappa_tot_gk_err = [0]*len(kappa_tot_gk)
        kappa_eta_tot_err = [0]*len(kappa_eta_tot)
        Pave_err = [0]*len(Pave_tot)
        PTmax_err = [0]*len(PTmax)
        rhoave_err = [0]*len(rhoave)
        area_rho_err = [0]*len(area_rho)
        rhoeff = RHOeff_list
        R_final = R_init

    ###################################
    R_pad_min = np.min(rhoeff)-0.05
    R_pad_max = np.max(rhoeff)+0.05
    R_final_pad = np.linspace(R_pad_min, R_pad_max, 10)

    # print to csv for latex use
    with open('DATA/{}/rho_H_{}.csv'.format(DIR, r_name), 'wb') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(("rho", "H"))
        writer.writerows(imap(lambda x, y: (y, x), len(rhoeff)*[6.8], rhoeff))#izip(dZ,round(Hlist,2)))

    with open('DATA/{}/viscosity_{}.csv'.format(DIR, r_name), 'wb') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(("rhoave", "diff", "shear", "bulk", "ratio"))
        writer.writerows(imap(
                         lambda b, c, d, e, f, g, h, i, j:
                         ("{}".format(round(j, 3)),
                          "{} \pm {}".format(round(b, 3), round(c, 3)),
                          "{} \pm {}".format(round(d, 2), round(e, 2)),
                          "{} \pm {}".format(round(f, 2), round(g, 2)),
                          "{} \pm {}".format(round(h, 1), round(i, 1))),
                         diff_tot_y_msd, diff_tot_y_msd_err,
                         eta_xy_gk, eta_xy_gk_err,
                         kappa_tot_gk, kappa_tot_gk_err,
                         kappa_eta_tot, kappa_eta_tot_err,
                         rhoeff))

    print '\n#-----------------------Final plotting------------------\n#'

    ax3a.errorbar(rhoeff, PTmax, yerr=PTmax_err, marker='o')
    ax3a.set_xlabel('$\\rho$ (g/cm$^3$)')
    ax3a.set_ylabel('$P_{\mathrm{T, max}}$ (GPa)')
    ax3a.set_xlim(R_pad_min, R_pad_max)
    ax3a.legend()

    ax3b.errorbar(rhoeff, Pave_tot, yerr=Pave_err, marker='o')
    ax3b.set_xlabel('$\\rho$ (g/cm$^3$)')
    ax3b.set_ylabel('P (MPa)')
    ax3b.set_xlim(R_pad_min, R_pad_max)
    ax3b.legend()

    ax4.errorbar(rhoeff, area_rho, yerr=area_rho_err, marker='o')
    ax4.set_xlabel('$\\rho$ (g/cm$^3$)')
    ax4.set_ylabel('$\\rho_{\mathrm{area}}$ (nm$^{-2}$)')
    ax4.set_xlim(R_pad_min, R_pad_max)
    ax4.legend()

    ax5a = fig5.add_subplot(311)
    ax5a.set_yscale('log')
    ax5a.errorbar(rhoeff, diff_tot_y_msd, yerr=diff_tot_y_msd_err, marker='o')
    ax5a.plot(R_final_pad, len(R_final_pad)*[2.866], linestyle='dashed')
    # ax23a.axes.get_xaxis().set_ticks([])
    ax5a.tick_params(labelbottom='off')
    ax5a.set_ylabel('$D_s$ (10$^{-9}$ m$^2$/s)')
    ax5a.set_xlim(R_pad_min, R_pad_max)
    ax5a.legend()

    ax5b = fig5.add_subplot(312)
    ax5b.set_yscale('log')
    ax5b.errorbar(rhoeff, eta_xy_gk, yerr=eta_xy_gk_err, marker='o')
    ax5b.plot(R_final_pad, len(R_final_pad)*[0.67], linestyle='dashed')
    ax5b.set_ylabel('$\eta$ (mPas)', labelpad=20)
    ax5b.tick_params(labelbottom='off')
    ax5b.set_xlim(R_pad_min, R_pad_max)
    ax5b.legend()

    ax5c = fig5.add_subplot(313)
    ax5c.set_yscale('log')
    ax5c.errorbar(rhoeff, kappa_tot_gk, yerr=kappa_tot_gk_err, marker='o')
    ax5c.plot(R_final_pad, len(R_final_pad)*[1.59], linestyle='dashed')
    ax5c.set_xlabel('$\\rho$ (g/cm$^3$)')
    ax5c.set_ylabel('$\kappa$ (mPas)', labelpad=20)
    ax5c.set_xlim(R_pad_min, R_pad_max)
    ax5c.legend()

    ax6.errorbar(Pave_tot, kappa_tot_gk, yerr=kappa_tot_gk_err, marker='o', linestyle=None)
    ax6.set_xlabel('P (MPa)')
    ax6.set_ylabel('$\kappa$ (mPas)')

    ax9.errorbar(rhoeff, kappa_eta_tot, yerr=kappa_eta_tot_err, marker='o')
    ax9.plot(R_final_pad, len(R_final_pad)*[2.32], linestyle='dashed', label='Bulk')
    ax9.set_xlim(R_pad_min, R_pad_max)
    ax9.set_xlabel('$\\rho$ (g/cm$^3$)')
    ax9.set_ylabel('$\kappa/\eta$ ')

    ax10.plot(rhoeff, R_final, linestyle='None', marker='o')
    ax10.set_xlim(R_pad_min, R_pad_max)
    ax10.set_xlabel('$\\rho$ (g/cm$^3$)')
    ax10.set_ylabel('r')

    fig3a.savefig('PLOTS_C/{}/rho_vs_pt_{}.pdf'.format(DIR, r_name),
                  box_inches='tight')
    fig3b.savefig('PLOTS_C/{}/rho_vs_p_{}.pdf'.format(DIR, r_name),
                  box_inches='tight')
    fig4.savefig('PLOTS_C/{}/area_dens_{}.pdf'.format(DIR, r_name),
                 bbox_inches='tight')
    fig5.savefig('PLOTS_C/{}/eta_diff_kappa_log_{}.pdf'.format(DIR, r_name),
                 bbox_inches='tight')
    fig6.savefig('PLOTS_C/{}/Pave_vs_kappa_{}.pdf'.format(DIR, r_name),
                 bbox_inches='tight')
    fig9.savefig('PLOTS_C/{}/kappa_eta_vs_rho_{}.pdf'.format(DIR, r_name),
                 bbox_inches='tight')
    fig10.savefig('PLOTS_C/{}/r_vs_rho_{}.pdf'.format(DIR, r_name),
                 bbox_inches='tight')
    print '\n#-----------------------Done------------------\n#'

    return


if __name__ == "__main__":
    sys.path.insert(0, '/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()
