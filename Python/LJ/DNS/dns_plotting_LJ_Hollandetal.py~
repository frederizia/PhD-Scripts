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




# Parse arguments form command line
parser = argparse.ArgumentParser()

parser.add_argument("-n", type=int, nargs=2, \
                    help="points for L, H",required=True)
parser.add_argument("-u", type=str, nargs=1, \
                    help="units",required=True, default='SI')
parser.add_argument("-s", type=str, nargs=1, \
                    help="style",required=True, default='conv')
parser.add_argument("-bc", type=str, nargs=1, \
                    help="rho BC",required=True, default='soft')

args = parser.parse_args()
nptsL = args.n[0] #150
nptsH = args.n[1] #40
units = args.u[0] # SI or scaled
style = args.s[0] # conv or cont
rho_BC = args.bc[0] # noflux or soft


# define geometry of box

domainLength = 2.3*10**(-9)
domainHeight = 2*10**(-7)
length_scale = 2*10**(-9)




# read in data

if units == 'scaled':
	filename = "{}/dns-velocityOutput-{}-{}-1-100-1-{}-{}.txt".format(units,style, rho_BC, nptsL, nptsH)
elif units == 'HF':
	filename = "../LJ_{}/dns-velocityOutput-{}-{}-1-0-100-0-{}-{}.txt".format(units,style, rho_BC, nptsL, nptsH)
elif units == 'new' and rho_BC == 'soft':
	filename = "dns-velocityOutput-{}-{}-1-0-100--{}-{}.txt".format(style, rho_BC, nptsL, nptsH)
elif units == 'new' or units == 'inert' and rho_BC == 'noflux':
	filename = "dns-velocityOutput-{}-{}-{}-{}.txt".format(style, rho_BC, nptsL, nptsH)
else:
	filename = "dns-velocityOutput-{}-{}-1-0-100-0-{}-{}.txt".format(style, rho_BC, nptsL, nptsH)

#uvals = []
#vvals = []
#rhovals = []
    
f = open(filename,'r')
data = f.read()

if style == 'conv':    
	data_final = data.split('\n')[0]
	data_final = ast.literal_eval(data_final[:-1])
	uvals = np.array(data_final[1])
	vvals = np.array(data_final[2])
	rhovals = np.array(data_final[3])
elif style == 'cont':
	# final data values
	data_final = data.split('\n')[-3]
	data_final = ast.literal_eval(data_final[:-1])
	uvals = np.array(data_final[2])
	vvals = np.array(data_final[3])
	rhovals = np.array(data_final[4])

	# initial data values
	data_init = data.split('\n')[0]
	data_init = ast.literal_eval(data_init[:-1])
	uvals_init = np.array(data_init[2])
	vvals_init = np.array(data_init[3])
	rhovals_init = np.array(data_init[4])
	uvals_avg_init = np.average(uvals_init, axis=0)
	rhovals_avg_L_init = np.average(rhovals_init, axis=0)

	nptsL_init = uvals_init.shape[0]
	nptsH_init = uvals_init.shape[1]

elif style == 'HF':    
	data_final = data.split('\n')[0]
	data_final = ast.literal_eval(data_final[:-1])
	uvals = np.array(data_final[2])
	vvals = np.array(data_final[3])
	rhovals = np.array(data_final[4])
elif style == 'new' or style == 'simple':    
	# final data values
	data_final = data.split('\n')[-3]
	data_final = ast.literal_eval(data_final[:-1])
	uvals = np.array(data_final[2])
	vvals = np.array(data_final[3])
	rhovals = np.array(data_final[4])


	# initial data values
	data_init = data.split('\n')[0]
	data_init = ast.literal_eval(data_init[:-1])
	uvals_init = np.array(data_init[2])
	vvals_init = np.array(data_init[3])
	rhovals_init = np.array(data_init[4])
	uvals_avg_init = np.average(uvals_init, axis=0)
	rhovals_avg_L_init = np.average(rhovals_init, axis=0)


	nptsL_init = uvals_init.shape[0]
	nptsH_init = uvals_init.shape[1]

nptsL = uvals.shape[0]
nptsH = uvals.shape[1]


uvals_avg = np.average(uvals, axis=0)
rhovals_avg_L = np.average(rhovals, axis=0)
rhovals_avg_H = np.average(rhovals, axis=1)


rhovals_start = rhovals[2]
rhovals_end = rhovals[nptsL-2]


# define L and H grid
if units == 'scaled':
	if style == 'cont':
		domainLength = data_final[0]
		domainHeight = data_final[1]

		domainLength_init = data_init[0]
		domainHeight_init = data_init[1]
		print domainLength, domainHeight

		Lvals_init = np.linspace(0, domainLength_init, nptsL_init)
		Hvals_init = np.linspace(0, domainHeight_init, nptsH_init)[::-1]
		H_init, L_init = np.meshgrid(Hvals_init, Lvals_init)
	else:
		domainLength = domainLength/length_scale
		domainHeight = domainHeight/length_scale
		print domainLength, domainHeight		
	Lvals = np.linspace(0, domainLength, nptsL)[::-1]
	Hvals = np.linspace(0, domainHeight, nptsH)[::-1]


elif units == 'SI':

	if style == 'cont':
		domainLength = data_final[0]
		domainHeight = data_final[1]

		domainLength_init = data_init[0]
		domainHeight_init = data_init[1]
		print domainLength, domainHeight

		Lvals = np.linspace(0, domainLength, nptsL)
		Hvals = np.linspace(0, domainHeight, nptsH)[::-1]
		Lvals_init = np.linspace(0, domainLength_init, nptsL_init)
		Hvals_init = np.linspace(0, domainHeight_init, nptsH_init)[::-1]
		H_init, L_init = np.meshgrid(Hvals_init, Lvals_init)
		
	else:
		print domainLength, domainHeight

		Lvals = np.linspace(0, domainLength*10**(9), nptsL)
		Hvals = np.linspace(0, domainHeight*10**(9), nptsH)[::-1]

elif units == 'new' or units == 'simple' or units == 'inert':

	if style == 'new' or style == 'simple':
		domainLength = data_final[0]
		domainHeight = data_final[1]

		domainLength_init = data_init[0]
		domainHeight_init = data_init[1]
		print domainLength, domainHeight

		Lvals = np.linspace(0, domainLength, nptsL)
		Hvals = np.linspace(0, domainHeight, nptsH)[::-1]
		Lvals_init = np.linspace(0, domainLength_init, nptsL_init)
		Hvals_init = np.linspace(0, domainHeight_init, nptsH_init)[::-1]
		H_init, L_init = np.meshgrid(Hvals_init, Lvals_init)
		
	else:
		print domainLength, domainHeight

		Lvals = np.linspace(0, domainLength*10**(9), nptsL)
		Hvals = np.linspace(0, domainHeight*10**(9), nptsH)[::-1]

elif units == 'HF':

	if style == 'cont':
		domainLength = data_final[0]
		domainHeight = data_final[1]

		domainLength_init = data_init[0]
		domainHeight_init = data_init[1]
		print domainLength, domainHeight

		Lvals = np.linspace(0, domainLength, nptsL)
		Hvals = np.linspace(0, domainHeight, nptsH)[::-1]
		Lvals_init = np.linspace(0, domainLength_init, nptsL)
		Hvals_init = np.linspace(0, domainHeight_init, nptsH)[::-1]
		H_init, L_init = np.meshgrid(Hvals_init, Lvals_init)
	else:
		domainLength = data_final[0]
		domainHeight = data_final[1]

		print domainLength, domainHeight

		Lvals = np.linspace(0, domainLength, nptsL)
		Hvals = np.linspace(0, domainHeight, nptsH)[::-1]



H, L = np.meshgrid(Hvals, Lvals)


#----------------------------------------------#

# plotting

# u
fig = plt.figure()
ctest=plt.contourf(L, H, uvals, cmap=cm.RdBu, levels=np.linspace(np.amin(uvals),np.amax(uvals),500))
plt.colorbar()
plt.xlabel('L')
plt.ylabel('H')
plt.savefig('uvals_{}_{}_{}_{}_{}.png'.format(units,style, rho_BC, nptsL, nptsH))
plt.show()

'''# v
fig = plt.figure()
ctest=plt.contourf(L, H, vvals, cmap=cm.RdBu, levels=np.linspace(np.amin(vvals),np.amax(vvals),500))
plt.colorbar()
plt.xlabel('L')
plt.ylabel('H')
plt.savefig('vvals.pdf')
plt.show()'''

# rho
fig = plt.figure()
ctest=plt.contourf(L, H, rhovals, cmap=cm.RdBu, levels=np.linspace(np.amin(rhovals),np.amax(rhovals),500))
plt.colorbar()
plt.xlabel('L')
plt.ylabel('H')
plt.savefig('rhovals_{}_{}_{}_{}_{}.png'.format(units,style, rho_BC, nptsL, nptsH))
plt.show()

print len(Hvals_init), len(uvals_avg_init)

# u average
fig = plt.figure()
if style == 'cont' or style == 'new':
	plt.plot(Hvals_init, uvals_avg_init, label='initial')
	plt.plot(Hvals, uvals_avg, label='final')
	plt.legend(loc='lower left')
else:
	plt.plot(Hvals, uvals_avg)
plt.xlabel('H')
plt.ylabel('u')
plt.savefig('uvals_avg_{}_{}_{}_{}_{}.pdf'.format(units,style, rho_BC, nptsL, nptsH))
plt.show()

# rho average over L
fig = plt.figure()
if style == 'cont' or style == 'new':
	plt.plot(Hvals, rhovals_avg_L, label='final')
	plt.plot(Hvals_init, rhovals_avg_L_init, label='initial')
	plt.legend(loc='lower left')
else:
	plt.plot(Hvals, rhovals_avg_L)
plt.xlabel('H')
plt.ylabel('$\\rho$')
plt.savefig('rhovals_avg_L_{}_{}_{}_{}_{}.pdf'.format(units,style, rho_BC, nptsL, nptsH))
plt.show()

# rho average over H
fig = plt.figure()
plt.plot(Lvals, rhovals_avg_H)
plt.xlabel('L')
plt.ylabel('$\\rho$')
plt.savefig('rhovals_avg_H_{}_{}_{}_{}_{}.pdf'.format(units,style, rho_BC, nptsL, nptsH))
plt.show()

'''#rhovals start and end
fig = plt.figure()
plt.plot(Hvals, rhovals_start)
plt.plot(Hvals, rhovals_end)
plt.xlabel('H')
plt.ylabel('$\\rho$')
plt.savefig('rhovals_s_e.pdf')
plt.show()'''



