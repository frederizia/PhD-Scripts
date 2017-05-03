#!/usr/bin/env python
'''Script to analyse flow rates compared to HP'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import sys
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
from tools import *
import scipy.integrate as si
import re

def GetArgs():
    # Parse arguments form command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-var", type=str, nargs=1, \
                        help="Variable to be considered ",required=True)
    parser.add_argument("-H", type=str, nargs='+', \
                        help="Height ",required=True)
    parser.add_argument("-f", type=str, nargs='+', \
                        help="force",required=True)
    parser.add_argument("-e", type=str, nargs='+', \
                        help="epsilon fluid-solid",required=True)
    parser.add_argument("-r", type=str, nargs=1, \
                        help="liquid density",required=True)
    parser.add_argument("-x", type=str, nargs='+', \
                        help="Xtra",required=True)
    args = parser.parse_args()
    return args

def shear_17(rho):
    #rho = np.array(map(float,rho))
    return 0.167*np.exp(1.715*rho)+0.008242*np.exp(6.292*rho)


def main():
    args = GetArgs()
    var = args.var[0]

    if var == 'eps':
        H = args.H[0]
        F = args.f[0]
        VAR = args.e
        RHOF = args.r[0]
        Xtra = args.x[0]
        name_nemd = 'H%s_eps%s_f%s_rhof%s'%(H,'%s',F,RHOF)
        name_plot = 'eps_H%s_f%s_rhof%s'%(H,F,RHOF)
        var_lgd = '\epsilon'
        logfile = 'H%s/eps%s/log.LJ_nemd_H%s_rhof%s_f%s'%(H,'%s',H,RHOF,F)

    elif var == 'H':
        VAR = args.H
        F = args.f[0]
        EPS = args.e[0]
        RHOF = args.r[0]
        Xtra = args.x[0]
        name_nemd = 'H%s_eps%s_f%s_rhof%s'%('%s',EPS,F,RHOF)
        name_plot = 'H_eps%s_f%s_rhof%s'%(EPS,F,RHOF)
        var_lgd = var
        logfile = 'H%s/eps%s/log.LJ_nemd_H%s_rhof%s_f%s'%('%s',EPS,'%s',RHOF,F)

    elif var == 'F':
        VAR = args.f
        H = args.H[0]
        EPS = args.e[0]
        RHOF = args.r[0]
        Xtra = args.x[0]
        name_nemd = 'H%s_eps%s_f%s_rhof%s'%(H,EPS,'%s',RHOF)
        name_plot = 'f_H%s_eps%s_rhof%s'%(H,EPS,RHOF)
        var_lgd = 'f'
        logfile = 'H%s/eps%s/log.LJ_nemd_H%s_rhof%s_f%s'%(H,EPS,H,RHOF,'%s')

    elif var == 'rhos':
        VAR = args.x
        H = args.H[0]
        EPS = args.e[0]
        RHOF = args.r[0]
        F = args.f[0]
        Xtra = 'None'
        name_nemd = 'H%s_eps%s_f%s_rhof%s_%s'%(H,EPS,F,RHOF,'%s')
        name_plot = 'rhos_H%s_eps%s_f%s_rhof%s'%(H,EPS,F,RHOF)
        var_lgd = '\\rho_s'
        logfile = 'H%s/eps%s/%s/log.LJ_nemd_H%s_rhof%s_f%s'%(H,EPS,'%s',H,RHOF,F)

    elif var == 'rhof':
        VAR = args.x
        H = args.H[0]
        EPS = args.e[0]
        RHOF = '%s'
        F = args.f[0]
        Xtra = 'None'
        name_nemd = 'H%s_eps%s_f%s_%s_%s'%(H,EPS,F,RHOF,'%s')
        name_plot = 'rhof_H%s_eps%s_f%s'%(H,EPS,F)
        var_lgd = '\\rho_f'
        logfile = 'H%s/eps%s/%s/log.LJ_nemd_H%s_%s_f%s'%(H,EPS,'%s',H,'%s',F)

    # If we are in the rhos folder
    if Xtra != 'None' and var != 'rhos' and var != 'rhof':
        name_nemd = '%s_%s'%(name_nemd,Xtra)
        name_plot = '%s_%s'%(name_plot,Xtra)
    if var == 'rhos':
        VAR_ticks = []
        for v in VAR:
            tick_tmp = re.search('rhos(.*)',v).group(1)
            VAR_ticks.append(tick_tmp)
    elif var == 'rhof':
        VAR_ticks = []
        for v in VAR:
            tick_tmp = re.search('rhof(.*)',v).group(1)
            VAR_ticks.append(tick_tmp)
    else:
        VAR_ticks = VAR

    rhovals = []
    uvals = []
    umax = [] 
    VARs = []
    massflow = []
    deltaP = []
    colours = ['b','g','r','k','c']
    i_tmp = 0

    for i in range(len(VAR)):
        try:
            V = VAR[i]
            if V[0:4] == 'rhof':	
                name_tmp = name_nemd %(V,V)
            else:
                name_tmp = name_nemd %V

            print 'Calculating flow rates for:'
            print var_lgd, ':', V

            rho_tmp = np.loadtxt('DATA/rho_%s.dat'%name_tmp)
            u_tmp = np.loadtxt('DATA/u_%s.dat'%name_tmp)
            y_tmp = np.loadtxt('DATA/y_%s.dat'%name_tmp)
            x_tmp = np.loadtxt('DATA/x_%s.dat'%name_tmp)		


            if var[0:4] == 'rhof' or var=='H':
                logfile_tmp = logfile % (V,V)
            else:
                logfile_tmp = logfile % (V)


            umax.append(np.max(np.average(u_tmp,axis=1)))

            #try finding midpoint
            midp_tmp, left, right = mid_point(y_tmp, y_tmp, rho_tmp)

            buff = int(left/3)
            LEFT_minus = left - buff
            RIGHT_plus = right + buff
            midp_new = midp_tmp-(LEFT_minus)

            y_tmp = y_tmp[LEFT_minus:RIGHT_plus]
            rho_tmp = rho_tmp[LEFT_minus:RIGHT_plus,:]
            u_tmp = u_tmp[LEFT_minus:RIGHT_plus,:]		

            y_mid_tmp = y_tmp[midp_new]
            for j in range(len(y_tmp)):
                y_tmp[j] -= y_mid_tmp


            rhovals.append(rho_tmp)
            uvals.append(u_tmp)
            VARs.append(V)
            yvals = y_tmp

            #--------------------Define values------------------

            
            # Calculation of quantities needed for HP flow rate
            # Find pressure in left and right bulk
            pl = np.array(map(float,read_log(logfile_tmp, 'Pleft')[1]))
            pr = np.array(map(float,read_log(logfile_tmp, 'Pright')[1]))
            Pleft_tmp = np.mean(pl)
            Pright_tmp = np.mean(pr)
            deltaP_tmp = Pright_tmp-Pleft_tmp
            deltaP.append(deltaP_tmp) # pressure drop

            
            #Find average density
            rho_avg_tmp = np.mean(rho_tmp)

            # Take viscosity from relations
            eta_tmp = shear_17(rho_avg_tmp)
            # Length of channel
            L_tmp = x_tmp[-1]
            # Height/2
            wallpos = wall_pos(rho_tmp,y_tmp)
            D_tmp = yvals[len(yvals)-wallpos[2]]
            #print D_tmp, yvals[-1]

            # Diffusivity
            #Ds_dict  = {'1.0': 0.0018, '0.5': 0.0019, '2.0': 0.0017}
            Ds_dict  = {'1.0': 0.2, '0.5': 0.4, '2.0': 0.1} # values currently made up
            # Work of adhesion
            WA_dict = {'1.0': -3.4, '0.5': -1.2, '2.0': -10.9} # values currently made up
            if var == 'eps':
                WA = WA_dict[V]
                D_s = Ds_dict[V]
            else:
                WA = WA_dict[EPS]
                D_s = Ds_dict[EPS]

            print ''
            print 'Parameters:'
            print 'Pressure drop:', deltaP_tmp
            print 'Average density:', rho_avg_tmp
            print 'Average viscosity:', eta_tmp
            print 'Channel length:', L_tmp
            print 'Channel height:', D_tmp*2
            print 'Diffusivity:', D_s
            print 'Work of Adhesion:', WA
            #-------------------------Flow rates-------------------

            dx_tmp = yvals[1]-yvals[0]


            # VOLUMETRIC FLOW RATE

            # MD
            Qmd  = vol_flow_uavg(u_tmp, dx_tmp)

            # HP
            u_hp = u_HP(y_tmp,deltaP_tmp,eta_tmp,L_tmp,D_tmp)
            #Qhp = vol_flow_uavg_HP(u_hp, dx_tmp)
            Qhp = vol_flow_HP(deltaP_tmp, eta_tmp, L_tmp, D_tmp)

            # model
            Qmodel = vol_flow_model(deltaP_tmp, eta_tmp, L_tmp, D_tmp, D_s, WA)

            print 'The calculated volumetric flow rates are:', Qmd, '(MD) and', Qhp, '(HP)'

            Qratio = Qmd/Qhp

            print 'The enhancement is:', Qratio

            # Theortical enhancement
            #Qr_theo = 1 - 3*(D_s/WA)*(eta_tmp*L_tmp/D_tmp**2)
            print 'The predicted enhancement is', Qmodel/Qhp #Qr_theo
            print ' '


            # MASS FLOW RATE

            # MD
            mdot_avg = mass_flow_avg(u_tmp, rho_tmp, dx_tmp)

            # HP
            mdot_hp = rho_avg_tmp*Qhp

            print 'The average mass flow rate is', mdot_avg, '(MD) and', mdot_hp,'(HP)'

            mdot_ratio = mdot_avg/mdot_hp

            print 'The enhancement is', mdot_ratio

            print ' '

            
            i_tmp +=1

        except IOError:
            print name_tmp
            continue

    '''
    # massflow and pressure drop

    ind = np.arange(1,len(VAR)+1)
    width = 0.2
    mdot_bar = ax3.bar(ind, massflow, width, color = 'b')
    ax4 = ax3.twinx()
    dP_bar = ax4.bar(ind+width, deltaP, width, color = 'g')
    ax3.set_ylabel('$\dot{m}/\\tau^{-1}$', color='b')
    for tl in ax3.get_yticklabels():
        tl.set_color('b')
    ax4.set_ylabel('$\Delta P/(\sigma^3/\epsilon)$',color='g')
    for tl in ax4.get_yticklabels():
        tl.set_color('g')
    ax3.set_xlabel('$%s$'%var_lgd)
    ax3.set_xticks(ind+width)
    ax3.set_xticklabels(VAR_ticks)

    m,b = np.polyfit(deltaP, massflow,1)
    mf_fit = []
    for dp in deltaP:
    	mf_fit.append(m*dp+b)
    ax5.plot(deltaP, massflow, marker = 'D',c='b', linestyle = 'None')
    ax5.plot(deltaP, mf_fit,c='b',lw = 2, linestyle = 'dashed')
    ax5.set_xlabel('$\Delta P/(\sigma^3/\epsilon)$')
    #plt.xlim(l_lim-0.3,r_lim+0.3)
    ax5.set_ylabel('$\dot{m}/\\tau^{-1}$')


    ax6.set_ylim(0, np.max(umax)+0.05)


    fig1.savefig('PLOTS/comb_rhovals_avg_Y_%s.pdf'%name_plot)
    fig2.savefig('PLOTS/comb_uvals_avg_Y_%s.pdf'%name_plot)
    fig3.savefig('PLOTS/comb_massflow_%s.pdf'%name_plot)
    fig5.savefig('PLOTS/comb_massflow_straight_%s.pdf'%name_plot)
    fig6.savefig('PLOTS/comb_uvals_avg_Y_HP_%s.pdf'%name_plot)
    '''
    return


if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS')
    main()

