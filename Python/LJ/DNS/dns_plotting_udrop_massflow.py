#!/usr/bin/python

"""Script written to plot the velocity change across channel and comparing the mass flow for different set ups with regard to the viscosities."""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import ast
from matplotlib import cm
import pylab as p
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
import argparse
import re
from matplotlib import rc





# Parse arguments form command line
parser = argparse.ArgumentParser()


parser.add_argument("-f", type=str, nargs='+', \
                    help="folder",required=True, default='XY')
parser.add_argument("-n", type=int, nargs=2, \
                    help="points for L, H",required=True)
parser.add_argument("-p", type=str, nargs=1, \
                    help="folder",required=False, default='n')
parser.add_argument("-l", type=str, nargs='+', \
                    help="labels",required=False, default='XY')
parser.add_argument("-t", type=str, nargs=1, \
                    help="t",required=True, default='n')
#parser.add_argument("-xy", type=int, nargs=2, \
 #                   help="X (length), y (height)",required=True)

args = parser.parse_args()
folder = args.f # folder
labels = args.l
ptype = args.t[0]
nptsX = args.n[0] #100
nptsY = args.n[1] #60 
plotting = args.p[0] # whether to plot or not
#domainLength = args.xy[0] #115
#domainHeight = args.xy[1] #1

cases = int(len(folder)/2)
folder = np.array(folder).reshape((cases,2))

# find geometry
domainLength = int(re.search('x(.*)y',folder[0][1]).group(1))
domainHeight = float(re.search('y(.*)',folder[0][1]).group(1))

print cases

delta_P = []#float(re.search('P(.*)x',folder[0][1]).group(1))-1

massflow = [1.08251621954, 0.9854034675, 0.679328986387]




var_list = ['u','v','rho']
uvals = []
vvals = []
rhovals = []

plotdir = '/home/fred/Dropbox/PROJECT/Presentations/APS Boston 2015/PLOTS'
col = ['b','g','r']
muvals = {'LJ':0.92695, 'water':1.0}



# read in data

for variable in var_list:
	for i in range(cases):
		if folder[i][0] == 't':
			f = 'LJ_Hollandetal/transient_261015'
		elif folder[i][0] == 'p':
			f = 'PIETRO/uvrho/scaled'
		filename = "{}/GRAD/{}/dns-{}Output-{}-{}.txt".format(f, folder[i][1], variable, nptsX, nptsY)
		f = open(filename,'r')
		data = f.read()

		# only from data split [2] to [2]
		temp_data = data.split('\n')[2:-2]

		for j in range(len(temp_data)):
			temp_data[j] = map(float, temp_data[j].split())
		#print temp_data[-1], len(temp_data)
		if variable == 'u':
			uvals.append(temp_data)
			delta_P.append(float(re.search('P(.*)x',folder[i][1]).group(1))-1)
		if variable == 'v':
			vvals.append(temp_data)
		if variable == 'rho':
			rhovals.append(temp_data)


'''
# initial data values

uvals_avg = np.average(uvals, axis=0)
rhovals_avg_L = np.average(rhovals, axis=0)
rhovals_avg_H = np.average(rhovals, axis=1)
uvals_centre = np.array(uvals)[:,nptsY-1]
rhovals_centre = np.array(rhovals)[:,nptsY-1]

'''


# define L and H grid


# think X,Y may be the wrong way around
Lvals = np.linspace(0, domainLength, nptsX)
Hvals = np.linspace(0, domainHeight, nptsY)[::-1]

H, L = np.meshgrid(Hvals, Lvals)


def HP(yvals,delta_P,mu,L,D):
	u = (delta_P)/(2*L*mu)*(D**2-yvals**2)
	return u
	


#----------------------------integrate rho*u over y--------------------------


def integrate(xval, rho, u):
	N = 2*nptsY-1
	h = domainLength/N
	rhoval = np.array(rho)[xval,:]
	rhoval = np.concatenate((rhoval, rhoval[::-1]),axis=1)
    	u = np.array(u)[xval,:]
    	u = np.concatenate((u, u[::-1]),axis=1)
	sum_rhou = 0
	for i in range(1, int(N/2)):
		sum_rhou += rhoval[2*i-2]*u[2*i-2]
        	sum_rhou += 4*rhoval[2*i-1]*u[2*i-1]
        	sum_rhou += rhoval[2*i]*u[2*i]
	final = (h/3)*sum_rhou
	return final



#print integrate(20, rhovals,uvals), integrate(40, rhovals,uvals)

print 'L', domainLength, 'H', domainHeight

#----------------------------------------------#


matplotlib.rcParams.update({'font.size': 19})
rc('text', usetex=True)


if plotting == 'y':
	# plotting

	# u drop
	delu = []
	ind = np.arange(1,cases+1)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	width = 0.2
	delu_lambda0 = np.array(uvals[0])[-1,nptsY-1]-np.array(uvals[0])[0,nptsY-1]
	for i in range(cases):
		delu.append((np.array(uvals[i])[-1,nptsY-1]-np.array(uvals[i])[0,nptsY-1])/delu_lambda0)
	u_bar = ax.bar(ind+width/2, delu, width, color = 'b')
	#P_bar = ax.bar(ind+width, delta_P, width, color='r')	
	#plt.xlabel('$H$')
	plt.ylabel('$\Delta u$/$\Delta u(\lambda_0)$')
	#plt.ylim(0.0007,0.00095)
	#ax.legend((u_bar[0],P_bar[0]),('$\Delta u$','$\Delta P$'))
	ax.set_xticks(ind+width)
	ax.set_xticklabels(labels)
	#plt.ylim(0,0.018)
	plt.legend()
	#plt.savefig('{}/delta_{}_{}_{}.pdf'.format(plotdir,nptsX,nptsY,ptype))
	plt.show()


	# mass flow
    	mdot = []
    	ind = np.arange(1,cases+1)
    	fig = plt.figure()
    	ax = fig.add_subplot(111)
    	width = 0.2
    	for i in range(cases):
        	mdot.append(integrate(55,rhovals[i],uvals[i]))
    	mdot_bar = ax.bar(ind+width/2, mdot, width, color = 'b')
	#P_bar = ax.bar(ind+width, delta_P, width, color='r')
	#plt.xlabel('$H$')
    	plt.ylabel('$\mathrm{Mass}$ $\mathrm{flow}$')
	#plt.ylim(0.0007,0.00095)
	#ax.legend((u_bar[0],P_bar[0]),('$\Delta u$','$\Delta P$'))
    	ax.set_xticks(ind+width)
    	ax.set_xticklabels(labels)
	#plt.ylim(0,0.018)
    	plt.legend()
#plt.savefig('{}/massflow_{}_{}_{}.pdf'.format(plotdir,nptsX,nptsY,ptype))
    	plt.show()




