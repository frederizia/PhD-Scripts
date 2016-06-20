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
parser.add_argument("-v", type=str, nargs='+', \
                    help="Velocity",required=True)
parser.add_argument("-n", type=str, nargs='+', \
                    help="Epsilon, sigma  or type",required=True)
parser.add_argument("-t", type=str, nargs=1, \
                    help="thermostat",required=True)
args = parser.parse_args()
vel = args.v
name = args.n
thermo = args.t[0]

plt_name = '%s'%(thermo) 

slip = []
V = []
for n in name:
	V_tmp = []
	slip_tmp = []
	plt_name += '_'+n
	data_tmp = np.loadtxt('DATA/slip_%s_%s.dat'%(n,thermo))
	for i in range(len(vel)):
		print data_tmp[i][1:]
		if data_tmp[i][0:1] == float(vel[i]):
			slip_tmp.append(data_tmp[i][1:][0])
			V_tmp.append(data_tmp[i][0:1][0])
	slip.append(slip_tmp)
	V.append(V_tmp)


#---------------PLOTTING------------------

matplotlib.rcParams.update({'font.size': 15})
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

fig = plt.figure()
for i in range(len(name)):
	EPS = re.search('eps(.*)',name[i]).group(1)
	RHOS = re.search('rhos(.*)_',name[i]).group(1)
	#print vel, slip[i], len(V[i]), len(slip[i])
	plt.plot(V[i], slip[i], linestyle = 'dashed',lw=2.0, marker = 'D', label='$\\rho_s$ = %s, $\epsilon$ = %s'%(RHOS,EPS))
plt.legend(loc='lower right')
plt.savefig('PLOTS/comb_%s.pdf'%plt_name)
plt.show()

