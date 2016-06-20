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


# Plotting data for Pietro's relations, here just the initial distributions

# NOTE: should be run from '~/Dropbox/DNS\ Solver/PIETRO'


# Parse arguments form command line
parser = argparse.ArgumentParser()


parser.add_argument("-f", type=str, nargs=2, \
                    help="folder",required=True, default='XY')
parser.add_argument("-n", type=int, nargs=2, \
                    help="points for L, H",required=True)
#parser.add_argument("-xy", type=int, nargs=2, \
 #                   help="X (length), y (height)",required=True)

args = parser.parse_args()
folder = args.f[0] # folder
ftype = args.f[1] # type
nptsX = args.n[0] #100
nptsY = args.n[1] #60 
#domainLength = args.xy[0] #115
#domainHeight = args.xy[1] #1


# Determine the geometry of the system from the folder name
domainLength = int(re.search('x(.*)y',ftype).group(1))
domainHeight = float(re.search('y(.*)',ftype).group(1))



#var_list = ['u','rho','v']
var_list = ['u','rho']
ALL_DATA = []





# read in data

for variable in var_list:
	filename = "uvrho/scaled/{}/{}/dns-{}Output-init-{}-{}.txt".format(folder, ftype, variable, nptsX, nptsY)
	f = open(filename,'r')
	data = f.read()

	# only from data split [2] to [2]
	temp_data = data.split('\n')[2:-2]

	for i in range(len(temp_data)):
		temp_data[i] = map(float, temp_data[i].split())
	#print temp_data[-1], len(temp_data)
	ALL_DATA.append(temp_data)


uvals = ALL_DATA[0]
rhovals = ALL_DATA[1]
#vvals = ALL_DATA[2]







# initial data values

uvals_avg = np.average(uvals, axis=0)
rhovals_avg_L = np.average(rhovals, axis=0)
rhovals_avg_H = np.average(rhovals, axis=1)
uvals_centre = np.array(uvals)[:,nptsY-1]
rhovals_centre = np.array(rhovals)[:,nptsY-1]

#print uvals[0][59], uvals[1][59], np.array(uvals)[:,59], len(np.array(uvals)[:,59])


rhovals_start = rhovals[2]
rhovals_end = rhovals[nptsY-2]


# define L and H grid


# think X,Y may be the wrong way around
Lvals = np.linspace(0, domainLength, nptsX)
Hvals = np.linspace(0, domainHeight, nptsY)[::-1]

H, L = np.meshgrid(Hvals, Lvals)


#----------------------------------------------#

# plotting

# u
fig = plt.figure()
ctest=plt.contourf(L, H, uvals, cmap=cm.RdBu, levels=np.linspace(np.amin(uvals),np.amax(uvals),500))
plt.colorbar()
plt.xlabel('L')
plt.ylabel('H')
plt.savefig('uvrho/scaled/{}/{}/uvals_in_{}_{}_{}.png'.format(folder,ftype, nptsX,nptsY,ftype))
plt.show()

'''
# v
fig = plt.figure()
ctest=plt.contourf(L, H, vvals, cmap=cm.RdBu, levels=np.linspace(np.amin(vvals),np.amax(vvals),500))
plt.colorbar()
plt.xlabel('L')
plt.ylabel('H')
plt.savefig('vvals.pdf')
plt.show()
'''

# rho
fig = plt.figure()
ctest=plt.contourf(L, H, rhovals, cmap=cm.RdBu, levels=np.linspace(np.amin(rhovals),np.amax(rhovals),500))
plt.colorbar()
plt.xlabel('L')
plt.ylabel('H')
plt.savefig('uvrho/scaled/{}/{}/rhovals_in_{}_{}_{}.png'.format(folder,ftype,nptsX,nptsY,ftype))
plt.show()

# u centre
fig = plt.figure()
plt.plot(Lvals, uvals_centre, label='final')
plt.xlabel('L')
plt.ylabel('u centre')
plt.savefig('uvrho/scaled/{}/{}/uvals_centre_in_{}_{}_{}.png'.format(folder,ftype,nptsX,nptsY,ftype))
plt.show()

# rho centre
fig = plt.figure()
plt.plot(Lvals, rhovals_centre, label='final')
plt.xlabel('L')
plt.ylabel('$\\rho$ centre')
plt.savefig('uvrho/scaled/{}/{}/rhovals_centre_in_{}_{}_{}.png'.format(folder,ftype,nptsX,nptsY,ftype))
plt.show()


# u init
fig = plt.figure()
plt.plot(Hvals, np.array(uvals)[0,:], label='initial')
plt.plot(Hvals, np.array(uvals)[-1,:], label='final')
plt.plot(Hvals, np.array(uvals)[5,:], label='10')
plt.xlabel('H')
plt.ylabel('u init')
plt.legend()
plt.savefig('uvrho/scaled/{}/{}/uvals_init_in_{}_{}_{}.png'.format(folder,ftype,nptsX,nptsY,ftype))
plt.show()


# rho init
fig = plt.figure()
plt.plot(Hvals, np.array(rhovals)[0,:], label='initial')
#plt.plot(Hvals, np.array(uvals)[-1,:], label='final')
#plt.plot(Hvals, np.array(rhovals)[30,:], label='n=30')
#plt.plot(Hvals, np.array(rhovals)[59,:], label='n=59/final')
plt.legend()
plt.xlabel('H')
plt.ylabel('$\\rho$ init')
plt.savefig('uvrho/scaled/{}/{}/rhovals_init_in_{}_{}_{}.png'.format(folder,ftype,nptsX,nptsY,ftype))
plt.show()



'''
# u average
fig = plt.figure()

plt.plot(Hvals, uvals_avg, label='final')
plt.legend(loc='lower left')
plt.xlabel('H')
plt.ylabel('u')
plt.savefig('{}/uvals_avg_{}_{}.pdf'.format(folder,nptsX,nptsY,folder))
plt.show()

# rho average over L
fig = plt.figure()

plt.plot(Hvals, rhovals_avg_L, label='final')
plt.xlabel('H')
plt.ylabel('$\\rho$')
plt.savefig('{}/rhovals_avg_L_{}_{}_{}.pdf'.format(folder,nptsX,nptsY,folder))
plt.show()

# rho average over H
fig = plt.figure()
plt.plot(Lvals, rhovals_avg_H)
plt.xlabel('L')
plt.ylabel('$\\rho$')
plt.savefig('{}/rhovals_avg_H_{}_{}_{}.pdf'.format(folder,nptsX,nptsY,folder))
plt.show()

#rhovals start and end
fig = plt.figure()
plt.plot(Hvals, rhovals_start)
plt.plot(Hvals, rhovals_end)
plt.xlabel('H')
plt.ylabel('$\\rho$')
plt.savefig('rhovals_s_e.pdf')
plt.show()'''


