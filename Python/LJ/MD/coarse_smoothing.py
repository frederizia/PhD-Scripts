#!/usr/bin/env python
'''Code that takes an MD profile and attempts to smooth it such that it is easier to use in a finite difference code with derivatives'''

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
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import LSQBivariateSpline

# Parse arguments form command line
parser = argparse.ArgumentParser()


parser.add_argument("-H", type=str, nargs=1, \
                    help="Height ",required=True)
parser.add_argument("-f", type=str, nargs=1, \
                    help="force",required=True)
parser.add_argument("-e", type=str, nargs=1, \
                    help="epsilon fluid-solid",required=True)
parser.add_argument("-r", type=str, nargs=1, \
                    help="liquid density",required=False, default = ['0.5'])
parser.add_argument("-x", type=str, nargs=1, \
                    help="Extra info",required=True, default='None')



args = parser.parse_args()

H = args.H[0]
F = args.f[0]
EPS = args.e[0]
RHOF = args.r[0]
Xtra = args.x[0]
if Xtra != 'None':
	name_nemd = 'H%s_eps%s_f%s_rhof%s_%s'%(H,EPS,F,RHOF,Xtra)
else:
	name_nemd = 'H%s_eps%s_f%s_rhof%s'%(H,EPS,F,RHOF)

# Read in data and define x,y
try:
	rhoinit = np.loadtxt('DATA/rho_sym_%s.dat'%(name_nemd))
except:
	raise Exception('Run channel_analysis.py first')

yvals = np.arange(rhoinit.shape[0])
xvals = np.arange(rhoinit.shape[1])
Xvals, Yvals = np.meshgrid(xvals, yvals)

fac = 10.0

# Define finer meshgrid
yvals_fine = np.linspace(0,rhoinit.shape[0], rhoinit.shape[0]*fac)
xvals_fine = np.linspace(0,rhoinit.shape[1], rhoinit.shape[1]*fac)
Xvals_fine, Yvals_fine = np.meshgrid(xvals_fine, yvals_fine)

'''x, y = np.mgrid[-1:1:135j, -1:1:109j]
z = (x+y) * np.exp(-6.0*(x*x+y*y))
print x.shape,z.shape

xnew, ynew = np.mgrid[-1:1:270j, -1:1:218j]
tck = interpolate.bisplrep(x, y, z, s=0)
znew = interpolate.bisplev(xnew[:,0], ynew[0,:], tck)

print xnew.shape,znew.shape'''

print rhoinit.shape, yvals.size, xvals.size

# Perform interpolation
#tck = interpolate.bisplrep(Xvals, Yvals, rhoinit)
#rhointer = interpolate.bisplev(Xvals_fine[:,0], Yvals_fine[0,:], tck)
inter = RectBivariateSpline(yvals, xvals, rhoinit)
rhointer = inter(yvals_fine, xvals_fine)

#rhointer = LSQBivariateSpline(Yvals.ravel(), Xvals.ravel(), rhoinit.ravel(),Yvals_fine.ravel(),Xvals_fine.ravel())
#rhointer = rhointer.reshape(yvals_fine.shape[0],yvals_fine.shape[1])


print rhointer

matplotlib.rcParams.update({'font.size': 19})
matplotlib.rc('text', usetex=True)

den_max = np.max(rhoinit)
den_min = np.min(rhoinit)
fig = plt.figure()
fig.text(0.44, 0.025, '$x$', ha='center', va='center', fontsize=26)
#plt.ylim(l_lim-0.3,r_lim+0.3)
plt.ylabel('$y$', fontsize=26)
plt.tick_params(pad=7)
ctest=plt.contourf(xvals, yvals, rhoinit, cmap=cm.RdBu, levels=np.linspace(den_min,den_max,500))
plt.colorbar()
plt.savefig('smoothing_orig_den_2d_%s_init.png'%name_nemd)
#plt.show()

den_max = np.max(rhointer)
den_min = np.min(rhointer)
fig = plt.figure()
fig.text(0.44, 0.025, '$x$', ha='center', va='center', fontsize=26)
#plt.ylim(l_lim-0.3,r_lim+0.3)
plt.ylabel('$y$', fontsize=26)
plt.tick_params(pad=7)
ctest=plt.contourf(xvals_fine, yvals_fine, rhointer, cmap=cm.RdBu, levels=np.linspace(den_min,den_max,500))
plt.colorbar()
plt.savefig('smoothing_den_2d_%s_init.png'%name_nemd)
plt.show()

#test plot
test_pos=100
plt.plot(yvals, rhoinit[:,test_pos], label='init')
plt.plot(yvals_fine, rhointer[:,test_pos*fac], label='inter')
plt.legend()
plt.show()
