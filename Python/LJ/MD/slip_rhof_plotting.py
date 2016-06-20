#!/usr/bin/env python

'''Code to plot the slip rate vs velocity/shear rate for different set ups'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import re

# Parse arguments form command line
parser = argparse.ArgumentParser()
parser.add_argument("-rhof", type=str, nargs='+', \
                    help="Fluid density",required=True)
parser.add_argument("-n", type=str, nargs=1, \
                    help="Epsilon, sigma  or type",required=True)
parser.add_argument("-t", type=str, nargs=1, \
                    help="thermostat",required=True)
args = parser.parse_args()
rhof = args.rhof
name = args.n[0]
thermo = args.t[0]

plt_name = '%s_%s'%(name,thermo) 

slip = []
error = []
RHO = []
for r in rhof:
	RHO_tmp = []
	slip_tmp = []
	data_tmp = np.loadtxt('DATA/slip_avg_%s_rhof%s_%s.dat'%(name,r,thermo))
	slip.append(data_tmp[1])
	error.append(data_tmp[2])
	RHO.append(data_tmp[0])


#---------------PLOTTING------------------

matplotlib.rcParams.update({'font.size': 15})
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

fig = plt.figure()
plt.errorbar(RHO, slip, yerr=error,linestyle = 'dashed',lw=2.0, marker = 'D')
plt.legend(loc='lower right')
plt.xlabel('$\\rho_f$')
plt.ylabel('$L_s$')
plt.savefig('PLOTS/comb_rhof_%s.pdf'%(plt_name))
plt.show()

