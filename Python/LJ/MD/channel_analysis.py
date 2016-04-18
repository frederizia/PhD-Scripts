#!/usr/bin/env python

'''Reads in LAMMPS data files and organises it. Then plots both contour, 3d and simple data.'''
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


# Parse arguments form command line
parser = argparse.ArgumentParser()

parser.add_argument("-H", type=str, nargs=1, \
                    help="Channel width",required=True)
parser.add_argument("-f", type=str, nargs=1, \
                    help="force",required=True)
parser.add_argument("-e", type=str, nargs=1, \
                    help="epsilon fluid-solid",required=True)
parser.add_argument("-r", type=str, nargs=1, \
                    help="liquid density",required=True)
parser.add_argument("-x", type=str, nargs=1, \
                    help="Extra info",required=True, default='None')



args = parser.parse_args()
H = args.H[0]
F = args.f[0]
EPS = args.e[0]
RHOF = args.r[0]
Xtra = args.x[0]

#----------------DATA--------------------

name_nemd = 'H%s_eps%s_f%s_rhof%s'%(H,EPS,F,RHOF)
name_equ = 'H%s_eps%s_rhof%s'%(H,EPS,RHOF)

if Xtra != 'None':
	name_nemd = 'H%s_eps%s_f%s_rhof%s_%s'%(H,EPS,F,RHOF,Xtra)
	name_equ = 'H%s_eps%s_rhof%s_%s'%(H,EPS,RHOF,Xtra)

# NEMD data
try:
	rhovals = np.loadtxt('DATA/rho_%s.dat'%name_nemd)
	uvals = np.loadtxt('DATA/u_%s.dat'%name_nemd)
	xvals = np.loadtxt('DATA/x_%s.dat'%name_nemd)
	yvals = np.loadtxt('DATA/y_%s.dat'%name_nemd)
except IOError:
	xvals, yvals, rhovals = read_data(H,F, EPS, RHOF,'nemd','rho', Xtra)
	xvals, yvals, uvals = read_data(H,F, EPS, RHOF, 'nemd', 'u', Xtra)
	rhovals = np.transpose(rhovals)
	uvals = np.transpose(uvals)

	xvals, yvals = xvals-xvals[0], yvals-yvals[0]

	np.savetxt('DATA/rho_%s.dat'%name_nemd, rhovals, fmt="%s", header = "%i %i %.3f %.3f" % (rhovals.shape[0], rhovals.shape[1],yvals[-1],xvals[-1]))
	np.savetxt('DATA/u_%s.dat'%name_nemd, uvals, fmt="%s", header = "%i %i %.3f %.3f" % (uvals.shape[0], uvals.shape[1],yvals[-1],xvals[-1]))
	np.savetxt('DATA/x_%s.dat'%name_nemd, xvals, fmt="%s")
	np.savetxt('DATA/y_%s.dat'%name_nemd, yvals, fmt="%s")


# manipulate data

yvals_new, rho_new, l_lim, r_lim = symmetry(xvals,yvals,rhovals)
yvals_new_u, u_new, l_lim_u, r_lim_u = symmetry(xvals,yvals,uvals)
MIDP = mid_point(xvals,yvals,rhovals)[0]

np.savetxt('DATA/rho_sym_%s.dat'%name_nemd, rho_new, fmt="%s", header = "%i %i %.3f %.3f" % (rho_new.shape[1], rho_new.shape[0],xvals[-1],yvals_new[-1]))
np.savetxt('DATA/u_sym_%s.dat'%name_nemd, rho_new, fmt="%s", header = "%i %i %.3f %.3f" % (u_new.shape[1], u_new.shape[0],xvals[-1],yvals_new[-1]))

rho_simple_sym = rho_init(rho_new)
# try transposing it
np.savetxt('DATA/rho_simple_sym_%s.dat'%name_nemd, rho_simple_sym, fmt="%s", header = "%i %i %.3f %.3f" % (rho_new.shape[1], rho_new.shape[0],xvals[-1],yvals_new[-1]))

rhovals_avg_Y = np.average(rhovals, axis=1)
rhovals_new_avg_Y = np.average(rho_new, axis=1)
uvals_avg_Y = np.average(uvals, axis=1)


# EMD data
try:
	rho_equ = np.loadtxt('DATA/rho_equ_%s.dat'%name_equ)
	xvals_equ = np.loadtxt('DATA/x_equ_%s.dat'%name_equ)
	yvals_equ = np.loadtxt('DATA/y_equ_%s.dat'%name_equ)
except IOError:
	xvals_equ, yvals_equ, rho_equ = read_data(H,F, EPS, RHOF,'equ','rho', Xtra)
	rho_equ = np.transpose(rho_equ)

	xvals_equ, yvals_equ = xvals_equ-xvals_equ[0], yvals_equ-yvals_equ[0]
	np.savetxt('DATA/rho_equ_%s.dat'%name_equ, rho_equ, fmt="%s", header = "%i %i %.3f %.3f" % (rho_equ.shape[0], rho_equ.shape[1],yvals_equ[-1],xvals_equ[-1]))
	np.savetxt('DATA/x_equ_%s.dat'%name_equ, xvals_equ, fmt="%s")
	np.savetxt('DATA/y_equ_%s.dat'%name_equ, yvals_equ, fmt="%s")

rhovals_avg_equ_Y = np.average(rho_equ, axis=1)
MIDP_EQU = mid_point(xvals_equ,yvals_equ,rho_equ)[0]

#-------------------------------PLOTTING--------------------

matplotlib.rcParams.update({'font.size': 19})
matplotlib.rc('text', usetex=True)


# DENSITY 2D

den_max = np.max(rhovals)
den_min = np.min(rhovals)
fig = plt.figure()
fig.text(0.44, 0.025, '$x$/ $\sigma$', ha='center', va='center', fontsize=26)
plt.ylim(l_lim-0.3,r_lim+0.3)
plt.ylabel('$y$/ $\sigma$', fontsize=26)
plt.tick_params(pad=7)
ctest=plt.contourf(xvals, yvals, rhovals, cmap=cm.RdBu, levels=np.linspace(den_min,den_max,500))
plt.colorbar()
plt.savefig('PLOTS/den_2d_%s.png'%name_nemd)


# test plot
den_max = np.max(rho_simple_sym)
den_min = np.min(rho_simple_sym)
fig = plt.figure()
fig.text(0.44, 0.025, '$x$/ $\sigma$', ha='center', va='center', fontsize=26)
plt.ylim(l_lim-0.3,r_lim+0.3)
plt.ylabel('$y$/ $\sigma$', fontsize=26)
plt.tick_params(pad=7)
ctest=plt.contourf(xvals, yvals_new, rho_simple_sym, cmap=cm.RdBu, levels=np.linspace(den_min,den_max,500))
plt.colorbar()
plt.savefig('PLOTS/den_2d_simple_%s.png'%name_nemd)



#plt.show()

# DENSITY 3D

fig = plt.figure()
Xvals, Yvals = np.meshgrid(xvals, yvals)
ax = Axes3D(fig)
#ax.contourf3D(Xvals, Yvals, rhovals)
ax.plot_surface(Xvals, Yvals, rhovals, rstride=1, cstride=1, lw=0, alpha=1.0, cmap = cm.RdBu)
ax.set_xlabel('X / $\sigma$')
ax.set_ylabel('Y / $\sigma$')
ax.set_zlabel('$\\rho$ / $\sigma^{-3}$')
fig.add_axes(ax)
plt.savefig('PLOTS/den_3d_%s.png'%name_nemd)

#plt.xlim(35,65)
#plt.ylim(20,21.3+float(DELX[i]))

#plt.show()

# rho selection
fig = plt.figure()
ax = fig.gca(projection='3d')

for i in range(0,len(xvals),15):
	yplot = np.empty(len(yvals))
	yplot.fill(xvals[i])
	rhovals_tmp = np.average(rhovals[:,i:10+i], axis=1)
	ax.plot(yvals, yplot, zs = rhovals_tmp, c='b')#zs=rhovals[:,i])
for i in range(0,len(xvals_equ),15):
	yplot_equ = np.empty(len(yvals_equ))
	yplot_equ.fill(xvals_equ[i])
	rhovals_tmp_equ = np.average(rho_equ[:,i:10+i], axis=1)
	ax.plot(yvals_equ, yplot_equ, zs = rhovals_tmp_equ, c='r')#zs=rhovals[:,i])
ax.set_xlim3d(l_lim-0.3,r_lim+0.3)
ax.set_ylim3d(xvals[0],xvals[-1])
ax.set_xlabel('$y$ / $\sigma$')
ax.set_ylabel('$x$ / $\sigma$')
plt.savefig('PLOTS/rhovals_sel_%s.png'%name_nemd)	
#plt.show()




# rho average over L
fig = plt.figure()

plt.plot(yvals, rhovals_avg_Y)
plt.plot(yvals, np.average(rhovals[:,:10],axis=1), c='g', linestyle = 'dashed', label = '$\mathrm{Entrance}$')
plt.plot(yvals, np.average(rhovals[:,-10:],axis=1), c='r', linestyle = 'dashed', label = '$\mathrm{Exit}$')
plt.xlabel('$Y$ / $\sigma$')
plt.xlim(l_lim-0.3,r_lim+0.3)
plt.ylabel('$\\rho$ / $\sigma^{-3}$')
plt.legend()
plt.savefig('PLOTS/rhovals_avg_Y_%s.png'%name_nemd)
#plt.show()

# u average over L
fig = plt.figure()

plt.plot(yvals, uvals_avg_Y)
plt.plot(yvals, np.average(uvals[:,:10],axis=1), c='g', linestyle = 'dashed', label = '$\mathrm{Entrance}$')
plt.plot(yvals, np.average(uvals[:,-10:],axis=1), c='r', linestyle = 'dashed', label = '$\mathrm{Exit}$')
plt.xlabel('$Y$ / $\sigma$')
plt.xlim(l_lim-0.3,r_lim+0.3)
plt.ylabel('$u_x$ / $(\epsilon/m)^{(1/2)}$')
plt.legend()
plt.savefig('PLOTS/uvals_avg_Y_%s.png'%name_nemd)
#plt.show()

# rho average over L symmetric
fig = plt.figure()

plt.plot(yvals_new, rhovals_new_avg_Y)
plt.plot(yvals_new, np.average(rho_new[:,:10],axis=1), c='g', linestyle = 'dashed', label = '$\mathrm{Entrance}$')
plt.plot(yvals_new, np.average(rho_new[:,-10:],axis=1), c='r', linestyle = 'dashed', label = '$\mathrm{Exit}$')
plt.xlabel('$Y$ / $\sigma$')
plt.xlim(l_lim-0.3,r_lim+0.3)
plt.ylabel('$\\rho_{\mathrm{sym}}$ / $\sigma^{-3}$')
plt.legend()
plt.savefig('PLOTS/rhovals_sym_avg_Y_%s.png'%name_nemd)
#plt.show()


# rho average over L equilibrium
fig = plt.figure()

plt.plot(yvals_equ, rhovals_avg_equ_Y, label='$\\rho_{\mathrm{EMD}}$')
plt.plot(yvals, rhovals_avg_Y, label='$\\rho_{\mathrm{NEMD}}$')
plt.xlabel('$Y$ / $\sigma$')
plt.xlim(l_lim-0.3,r_lim+0.3)
plt.ylabel('$\\rho$ / $\sigma^{-3}$')
plt.legend()
plt.savefig('PLOTS/rhovals_avg_equ_Y_%s.png'%name_nemd)
#plt.show()




# rho centre
fig = plt.figure()
plt.plot(xvals, rhovals[MIDP,:], label='$\mathrm{NEMD}$')
plt.plot(xvals_equ, rho_equ[MIDP_EQU,:], label='$\mathrm{EMD}$')
plt.xlabel('$x$ / $\sigma$')
plt.ylabel('$\\rho_{\mathrm{centre}}$  / $\sigma^{-3}$')
plt.legend()
#plt.ylim(0.1,0.7)
plt.savefig('PLOTS/rhovals_centre_%s.png'%name_nemd)
#plt.show()

# u centre
fig = plt.figure()
plt.plot(xvals, uvals[MIDP,:], label='$\mathrm{NEMD}$')
plt.xlabel('$x$ / $\sigma$')
plt.ylabel('$u_x$ centre / $(\epsilon/m)^{(1/2)}$')
plt.legend()
#plt.ylim(0.1,0.7)
plt.savefig('PLOTS/uvals_centre_%s.png'%name_nemd)
#plt.show()

# save data for fitting

save_data = np.transpose(np.array([xvals, rhovals[len(yvals)/2.,:]]))
np.savetxt('DATA/x_vs_rho_%s.dat'%name_nemd, save_data,fmt="%s",)
print np.mean(rhovals[MIDP,:])



