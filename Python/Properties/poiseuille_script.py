#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
from poiseuille_tools import *
import sys


def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sep', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-fl', '--fluid', nargs=1, required=False, default='Water', action='store')
    parser.add_argument('-w', '--wall', nargs=1, required=False, default='Graphene', action='store')
    parser.add_argument('-e', '--eps', nargs=1, type=float,required=False, default='[1.0]', action='store')
    parser.add_argument('-m', '--model', nargs=1, required=False, default='spce', action='store')
    parser.add_argument('-f', '--force', required=False, nargs='+', default='1.0', action='store')
    parser.add_argument('-png', '--png', nargs=1, type=str, required=False, default='n', action='store')
    parser.add_argument("-r", type=str, nargs=1, help="Rerun", default='n')
    args = parser.parse_args()
    return args

def main():
    args        = GetArgs()
    sep         = args.sep
    fluid       = args.fluid[0]
    wall        = args.wall[0]
    model       = args.model[0]
    eps         = args.eps[0]
    force       = args.force
    png         = args.png[0]
    rerun       = args.r[0]

    # Definitions
    fig_size_long = (14,5)
    fig_size_sq   = (9,7)
    if png == 'n':
        EXT = 'PDF'
        ext = 'pdf'
    else:
        EXT = 'PNG'
        ext = 'png'


    eta = 0.7*1e-3
    

    markers = ['o', 'D', 's', 'v', '^', 'd', '*']
    #colours = ['#3a42d4', '#d43a3a', '#60c1dc',\
    #'#60dc96', '#dc6060', '#addc60', '#dc6079', '#60dcd4']
    colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']

    # Enh vs L/D from EMD
    fig1 = plt.figure(figsize=fig_size_sq) 
    ax1  = fig1.add_axes([0.1,0.15,0.8,0.75])
    # enh WADs only
    fig1a = plt.figure(figsize=fig_size_sq) 
    ax1a  = fig1a.add_axes([0.1,0.15,0.8,0.75])
    # enh ALL
    fig1b = plt.figure(figsize=fig_size_sq) 
    ax1b  = fig1b.add_axes([0.1,0.15,0.8,0.75])
    # Ls vs D from EMD
    fig2 = plt.figure(figsize=fig_size_sq) 
    ax2  = fig2.add_axes([0.1,0.15,0.8,0.75])
    # WA/Ds vs fric
    fig3 = plt.figure(figsize=fig_size_sq) 
    ax3  = fig3.add_axes([0.1,0.15,0.8,0.75])
    # Eta vs fric
    fig5 = plt.figure(figsize=fig_size_sq) 
    ax5  = fig5.add_axes([0.1,0.15,0.8,0.75])

    Ls_WADs = []
    Ls_fric = []
    Ls_enh = []
    Ls_enh_avg = []
    Ls_enh_err = []
    Ls_enh_mc = []
    Ls_enh_mc_avg = []
    Ls_enh_mc_err = []
    enh_WADs = []
    enh_fric = []
    enh_force = []
    enh_force_mc = []
    enh_force_err = []
    enh_force_mc_err = []
    Eta_force_avg = []
    Eta_force_avg_err = []
    Eta_EMD = []

    LD, LH = [], []
    D = []
    H = []
    Dnemd = []
    count_s = 0
    name_s = 'eps{}'.format(eps)
    for s in sep:
        name_s += s
        enh_force_tmp = []
        enh_force_mc_tmp = []
        Ls_enh_tmp = []
        Ls_enh_mc_tmp = []
        Eta_force = []
        # Equilibrium properties
        print '--------------------- s = {} ----------------------'.format(s)
        print 'Calculating equilibrium properties for s=',s,'A...'

        # channel dimensions
        Lval, Hval, Ylo, Yhi, Zlo, Zhi, Xwidth, RHOave = dimensions(fluid,wall,model,eps,s,rerun)
        Lval *=1e-10
        Dval = (Hval/2)*1e-10
        D.append(Dval*1e10)
        H.append(Hval)
        LD.append(Lval/Dval)
        LH.append(Lval/(Hval*1e-10))

        print 'L/D =', Lval/Dval
        print 'Eta =', eta



        # WA in J/m^2
        WA_val, WA_err = wa(fluid,wall,model,eps,s)

        # Ds in m^2/s
        Ds_val = diff(fluid,wall,model,eps,s) 
        print 'D_s/W_A =', Ds_val/WA_val

        #Eta
        alpha = 1.7*1e-10
        T = 298
        kB = 1.38*1e-23
        Eta_EMD_val = (kB*T)/(3*np.pi*alpha*Ds_val)
        Eta_EMD.append(Eta_EMD_val*1e3)
        print 'Eta_wall =', Eta_EMD_val*1e3, 'mPas'


        # Ls from WA/Ds
        Ls_WADs_val = -(eta/2)*(Lval/Dval)*(Ds_val/WA_val)
        Ls_WADs.append(abs(Ls_WADs_val)*1e9)
        print 'Ls_WADs = ', Ls_WADs_val*1e9, 'nm'

        # enhancement
        enh_WADs_val = enh_WADs_fn(Lval, Dval, WA_val, Ds_val, eta)
        enh_WADs.append(enh_WADs_val)
        print 'Predicted enhancement WADs = ', enh_WADs_val
        # error
        # Add to plot

        # Friction
        fric_val, fric_err = fric(fluid,wall,model,eps,s)

        # Ls from friction
        Ls_fric_val = eta/fric_val
        # Ls_err
        # Add to plot
        Ls_fric.append(abs(Ls_fric_val)*1e9)
        enh_fric_val = 1+3*(Ls_fric_val/Dval)
        enh_fric.append(enh_fric_val)
        print 'Ls_fric = ', Ls_fric_val*1e9, 'nm'
        print 'Predicted enhancement fric = ', enh_fric_val


        # WA/Ds vs fric

        # WA/Ds vs fric
        fig4 = plt.figure(figsize=fig_size_sq) 
        ax4  = fig4.add_axes([0.1,0.15,0.8,0.75])
        dP_list, dP_err_list = [], []
        Q_list, Qmc_list, QP_list = [], [], []
           

        # Non-equilibirum properties
        print '----------------------------------------'
        print 'Calculating non-equilibrium properties for s=',s,'A...'
        for f in force:
            try:
                print '--------------------- f = {} ----------------------'.format(f)
                
                Dnemd.append(Dval*1e10)
                # pressure drop
                dP, dP_err = delP(fluid,wall,model,eps,s,f)
                dP_list.append(abs(dP))
                dP_err_list.append(dP_err)

                # Q and velocity plot
                Q, Z, VEL = flux(fluid,wall,model,eps,s,f,Ylo,Yhi,Zlo, Zhi,rerun)
                Q_list.append(Q)

                # Q mol count
                Qmc = Qmolcount(fluid,wall,model,eps,s,f,RHOave,Xwidth)
                Qmc_list.append(Qmc)

                # QP
                QP = vol_flow_P(dP, eta, Lval, Dval)
                QP_list.append(QP)

                # enhancement
                enh_force_val = Q/QP
                print 'Flow enhancement = ', enh_force_val
                enh_force_tmp.append(enh_force_val)

                # enhancement mc
                enh_force_mc_val = Qmc/QP
                print 'Flow enhancement MC = ', enh_force_mc_val
                enh_force_mc_tmp.append(enh_force_mc_val)

                # enhancement
                Ls_enh_val = (Dval/3)*(enh_force_val-1)
                Ls_enh.append(Ls_enh_val*1e9)
                Ls_enh_tmp.append(Ls_enh_val*1e9)
                print 'Ls_enh = ', Ls_enh_val*1e9, 'nm'

                # enhancement
                Ls_enh_mc_val = (Dval/3)*(enh_force_mc_val-1)
                Ls_enh_mc.append(Ls_enh_mc_val*1e9)
                Ls_enh_mc_tmp.append(Ls_enh_mc_val*1e9)
                print 'Ls_enh_mc = ', Ls_enh_mc_val*1e9, 'nm'


                EtaVal = viscosity(fluid,wall,model,eps,s,f,Z,VEL,dP,Lval,Dval)
                Eta_force.append(EtaVal*1e3)
            except:
                print 'No data for s={}, f={}.'.format(s,f)
        enh_force_avg_val = np.mean(np.array(enh_force_tmp))
        enh_force_err_val = stats.sem(enh_force_tmp)
        enh_force_mc_avg_val = np.mean(np.array(enh_force_mc_tmp))
        enh_force_mc_err_val = stats.sem(enh_force_mc_tmp)
        Ls_enh_avg_val = np.mean(np.array(Ls_enh_tmp))
        Ls_enh_err_val = stats.sem(Ls_enh_tmp)
        Ls_enh_mc_avg_val = np.mean(np.array(Ls_enh_mc_tmp))
        Ls_enh_mc_err_val = stats.sem(Ls_enh_mc_tmp)
        enh_force.append(enh_force_avg_val)
        enh_force_mc.append(enh_force_mc_avg_val)
        enh_force_err.append(enh_force_err_val)
        enh_force_mc_err.append(enh_force_err_val)
        Ls_enh_avg.append(Ls_enh_avg_val)
        Ls_enh_err.append(Ls_enh_err_val)
        Ls_enh_mc_avg.append(Ls_enh_mc_avg_val)
        Ls_enh_mc_err.append(Ls_enh_mc_err_val)
        Eta_force_avg_val = np.mean(np.array(Eta_force))
        Eta_force_err_val = stats.sem(Eta_force)
        Eta_force_avg.append(Eta_force_avg_val)
        Eta_force_avg_err.append(Eta_force_err_val)
        print '------------------------------------'
        print 'Average enh_force = ', enh_force_avg_val, '+/-', enh_force_err_val
        print 'Average Ls_enh = ', Ls_enh_avg_val, '+/-', Ls_enh_err_val, 'nm.'
        print 'Average enh_force_mc = ', enh_force_mc_avg_val, '+/-', enh_force_mc_err_val
        print 'Average Ls_enh_mc = ', Ls_enh_mc_avg_val, '+/-', Ls_enh_mc_err_val, 'nm.'
        print 'Average Eta_force = ', Eta_force_avg_val, '+/-', Eta_force_err_val, 'mPas.'
        print '------------------------------------'

        ax4.errorbar(dP_list, Q_list, xerr=dP_err_list, label = '$Q_{MD}$', marker=markers[0], c=colours[0], ls='dashed')
        ax4.errorbar(dP_list, Qmc_list, xerr=dP_err_list, label = '$Q_{count}$', marker=markers[4], c=colours[4], ls='dashed')
        ax4.errorbar(dP_list, QP_list, xerr=dP_err_list, label = '$Q_{P}$', marker=markers[2], c=colours[2], ls='dashed')
        ax4.set_xlabel('$\Delta$P (Pa)')
        ax4.set_ylabel('Q (m$^2$/s)')
        ax4.legend()
        fig4.savefig('PLOTS/{}/Q_QP_vs_dP_s{}.{}'.format(EXT,s,ext))


        # Calculate Ls from 
        count_s +=1

    enh_WADs_mod, LH_mod = [], []
    for i in range(1, 30):
        enh_WADs_mod.append(enh_WADs_fn(66.6*1e-10, i*1e-10, -10.6, 2.4e-9, eta))
        LH_mod.append(66.6/i)

    # Calculate average of Ls from fit to Q/QP vs delP plot

    # Calculate Ls from Q/QP vs D plot

    # PLOTTING
    # enhancement EMD
    ax1.plot(LH, enh_fric, marker=markers[0], c=colours[0], ls='dashed', label='$\lambda$')
    ax1.plot(LH, enh_WADs, marker=markers[1], c=colours[2], ls='dashed', label='$D_s/W_A$')
    ax1.set_xlabel('L/H')
    ax1.set_ylabel('$\epsilon$')
    ax1.legend()

    # enh WADs only
    ax1a.plot(LH, enh_WADs, marker=markers[1], c=colours[2], ls='dashed', label='$D_s/W_A$')
    ax1a.plot(LH_mod, enh_WADs_mod, c=colours[2], label='model')
    ax1a.set_xlabel('L/H')
    ax1a.set_ylabel('$\epsilon$')
    ax1a.legend()

    # enh All
    ax1b.plot(LH, enh_fric, marker=markers[0], c=colours[0], ls='dashed', label='$\lambda$')
    ax1b.plot(LH, enh_WADs, marker=markers[1], c=colours[2], ls='dashed', label='$D_s/W_A$')
    ax1b.errorbar(LH, enh_force, yerr=enh_force_err, marker=markers[2], c=colours[4], ls='dashed', label='$Q_{MD}$')
    ax1b.errorbar(LH, enh_force_mc, yerr=enh_force_mc_err, marker=markers[3], c=colours[6], ls='dashed', label='$Q_{count}$')
    ax1b.set_xlabel('L/H')
    ax1b.set_ylabel('$\epsilon$')
    ax1b.legend()

    # Ls EMD
    ax2.plot(H, Ls_fric, marker=markers[0], c=colours[0], ls='dashed', label='$\lambda$')
    ax2.plot(H, Ls_WADs, marker=markers[1], c=colours[2], ls='dashed', label='$D_s/W_A$')
    ax2.set_xlabel('H (\AA)')
    ax2.set_ylabel('$L_s (\mathrm{nm})$')
    ax2.legend()

    # Ls all
    ax3.plot(H, Ls_fric, marker=markers[0], c=colours[0], ls='dashed', label='$\lambda$')
    ax3.plot(H, Ls_WADs, marker=markers[1], c=colours[2], ls='dashed', label='$D_s/W_A$')
    ax3.errorbar(H, Ls_enh_avg, yerr=Ls_enh_err, marker=markers[2], c=colours[4], ls='dashed', label='$\epsilon$')
    ax3.errorbar(H, Ls_enh_mc_avg, yerr=Ls_enh_mc_err, marker=markers[3], c=colours[6], ls='dashed', label='$\epsilon$')
    #ax3.plot(Dnemd, Ls_enh, marker=markers[2], c=colours[4], ls='dashed', label='$\epsilon$')
    ax3.set_xlabel('H (\AA)')
    ax3.set_ylabel('$L_s (\mathrm{nm})$')
    ax3.legend()

    # Eta
    if eps == 5.0:
        ax5.errorbar(H, Eta_force_avg, yerr=Eta_force_avg_err, marker=markers[0], c=colours[0], label='NEMD')
    ax5.plot(H, Eta_EMD, marker=markers[1], c=colours[1], label='EMD, surface')
    ax5.set_xlabel('H (\AA)')
    ax5.set_ylabel('$\eta (\mathrm{mPas})$')
    ax5.legend()

    # SAVEFIG
    fig1.savefig('PLOTS/{}/enh_EMD_vs_D_s{}.{}'.format(EXT,name_s,ext),bbox_inches='tight')
    fig1a.savefig('PLOTS/{}/enh_WADs_vs_D_s{}.{}'.format(EXT,name_s,ext),bbox_inches='tight')
    fig1b.savefig('PLOTS/{}/enh_all_vs_D_s{}.{}'.format(EXT,name_s,ext),bbox_inches='tight')
    fig2.savefig('PLOTS/{}/Ls_EMD_vs_H_s{}.{}'.format(EXT,name_s,ext),bbox_inches='tight')
    fig3.savefig('PLOTS/{}/Ls_all_vs_H_s{}.{}'.format(EXT,name_s,ext),bbox_inches='tight')
    fig5.savefig('PLOTS/{}/Eta_vs_H_s{}.{}'.format(EXT,name_s,ext),bbox_inches='tight')

    return

if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()