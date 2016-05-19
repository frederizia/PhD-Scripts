#!/usr/bin/env python
'''Script to compare different MD simulations based on certain parameters'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
from tools import *
import scipy.integrate as si
import re


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
        logfile = 'H%/eps%/log.LJ_nemd_H%s_rhof%s_f%s'%(H,'%s',H,RHOF,F)

elif var == 'H':
	VAR = args.H
	F = args.f[0]
	EPS = args.e[0]
	RHOF = args.r[0]
        Xtra = args.x[0]
	name_nemd = 'H%s_eps%s_f%s_rhof%s'%('%s',EPS,F,RHOF)
	name_plot = 'H_eps%s_f%s_rhof%s'%(EPS,F,RHOF)
	var_lgd = var
	logfile = 'H%/eps%/log.LJ_nemd_H%s_rhof%s_f%s'%('%s',EPS,'%s',RHOF,F)

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
        logfile = 'H%/eps%/%s/log.LJ_nemd_H%s_rhof%s_f%s'%(H,EPS,'%s',H,RHOF,F)

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
        logfile = 'H%/eps%/%s/log.LJ_nemd_H%s_%s_f%s'%(H,EPS,'%s',H,'%s',F)

# If we are in the rhos folder
if Xtra != 'None' and var != 'rhos' and var != 'rhof':
        name_nemd = '%s_%s'%(name_nemd,Xtra)
        name_plot = '%s_%s'%(name_plot,Xtra)
if var == 'rhos':
	VAR_ticks = []
	for v in VAR:
		tick_tmp = re.search('rhos(.*)',v).group(1)
		VAR_ticks.append(tick_tmp)
if var == 'rhof':
        VAR_ticks = []
        for v in VAR:
                tick_tmp = re.search('rhof(.*)',v).group(1)
                VAR_ticks.append(tick_tmp)
else:
	VAR_ticks = VAR


rhovals = []
uvals = []
VARs = []
massflow = []
deltaP = []
colours = ['b','g','r','k','c']
i_tmp = 0

def P_17(rho):
        #rho = np.array(map(float,rho))
        P = 0.163*np.exp(-1.078*rho)+0.02561*np.exp(6.305*rho)
        return P

def shear_17(rho):
        #rho = np.array(map(float,rho))
        return 0.167*np.exp(1.715*rho)+0.008242*np.exp(6.292*rho)

def HP(yvals,delta_P,eta,L,D):
        u = ((delta_P)/(2*eta*L))*(D**2-yvals**2)
        return u


matplotlib.rcParams.update({'font.size': 19})
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
fig5, ax5 = plt.subplots()
fig6, ax6 = plt.subplots()
#fig3 = plt.figure()
#ax3 = fig3.add_subplot(111)


for i in range(len(VAR)):
	try:
		V = VAR[i]
		if V[0:4] == 'rhof':	
			name_tmp = name_nemd %(V,V)
		else:
			name_tmp = name_nemd %V
		rho_tmp = np.loadtxt('DATA/rho_%s.dat'%name_tmp)
		u_tmp = np.loadtxt('DATA/u_%s.dat'%name_tmp)
		y_tmp = np.loadtxt('DATA/y_%s.dat'%name_tmp)
		x_tmp = np.loadtxt('DATA/x_%s.dat'%name_tmp)		

		if var[0:4] == 'rhof' or var=='H':
			logfile_tmp = logfile % (V,V)
		else:
			logfile_tmp = logfile % (V)

		
			
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

		# average mass flow
		mdot_tmp = mass_flow(u_tmp,rho_tmp)
		mdot_int = si.simps(mdot_tmp)
		massflow.append(mdot_int)
		#massflow.append(np.average(mdot_tmp))

                # Calculation of quantities needed for HP
                eta_tmp = shear_17(np.mean(rho_tmp[:,:])) # eta
                deltaP_tmp = rho_init(rho_tmp)[1]
                deltaP.append(deltaP_tmp) # pressure drop
                L_tmp = x_tmp[-1]
                D_tmp = yvals[-1]

	
		# rho average over L
		ax1.plot(y_tmp, np.average(rho_tmp[:,50:70],axis=1), c=colours[i_tmp], linestyle = 'dashed', label = '$%s = \mathrm{%s}$'%(var_lgd,VAR_ticks[i]))
		ax1.set_xlabel('$Y$ / $\sigma$')
		#plt.xlim(l_lim-0.3,r_lim+0.3)
		ax1.set_ylabel('$\\rho$ / $\sigma^{-3}$')
		ax1.legend()

		#plt.show()

		# u average over L

		ax2.plot(y_tmp, np.average(u_tmp[:,50:150],axis=1), c=colours[i_tmp], linestyle = 'dashed', label = '$%s = \mathrm{%s}$'%(var_lgd,VAR_ticks[i]))
		ax2.set_xlabel('$Y$ / $\sigma$')
		#plt.xlim(l_lim-0.3,r_lim+0.3)
		ax2.set_ylabel('$u_x$ / $(\epsilon/m)^{(1/2)}$')
		ax2.legend()


                # u average over L - comparison to HP

                ax6.plot(y_tmp, np.average(u_tmp[:,50:150],axis=1), c=colours[i_tmp], linestyle = 'dashed', label = '$%s = \mathrm{%s}$'%(var_lgd,VAR_ticks[i]))
                ax6.plot(y_tmp, HP(y_tmp,deltaP_tmp,eta_tmp,L_tmp,D_tmp), c=colours[i_tmp]) 
		ax6.set_xlabel('$Y$ / $\sigma$')
                #plt.xlim(l_lim-0.3,r_lim+0.3)
                ax6.set_ylabel('$u_x$ / $(\epsilon/m)^{(1/2)}$')
                ax6.legend()
		
		#plt.show()
		i_tmp +=1
		
	except IOError:
		print name_tmp
		continue


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


fig1.savefig('PLOTS/comb_rhovals_avg_Y_%s.pdf'%name_plot)
fig2.savefig('PLOTS/comb_uvals_avg_Y_%s.pdf'%name_plot)
fig3.savefig('PLOTS/comb_massflow_%s.pdf'%name_plot)
fig5.savefig('PLOTS/comb_massflow_straight_%s.pdf'%name_plot)
fig6.savefig('PLOTS/comb_uvals_avg_Y_HP_%s.pdf'%name_plot)
