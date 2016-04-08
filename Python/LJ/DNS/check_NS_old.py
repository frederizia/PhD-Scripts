"""Determine whether the Navier-Stokes like relations are acutally satisfied in this version of the code.

Should be run from ~/Dropbox/PROJECT/DNS\ Solver/LJ_Hollandetal"""
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
elif units == 'new' and rho_BC == 'noflux':
	filename = "dns-velocityOutput-{}-{}-{}-{}.txt".format(style, rho_BC, nptsL, nptsH)
elif units == 'inert':
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
elif style == 'HF':    
	data_final = data.split('\n')[0]
	data_final = ast.literal_eval(data_final[:-1])
	uvals = np.array(data_final[2])
	vvals = np.array(data_final[3])
	rhovals = np.array(data_final[4])
elif style == 'new':    
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



print uvals.shape
nptsL = uvals.shape[0]
nptsH = uvals.shape[1]


# define L and H grid
if units == 'scaled':
	if style == 'cont':
		domainLength = data_final[0]
		domainHeight = data_final[1]

		domainLength_init = data_init[0]
		domainHeight_init = data_init[1]
		print domainLength, domainHeight

		Lvals_init = np.linspace(0, domainLength_init, nptsL)
		Hvals_init = np.linspace(0, domainHeight_init, nptsH)[::-1]
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
		Lvals_init = np.linspace(0, domainLength_init, nptsL)
		Hvals_init = np.linspace(0, domainHeight_init, nptsH)[::-1]
		H_init, L_init = np.meshgrid(Hvals_init, Lvals_init)
		
	else:
		print domainLength, domainHeight

		Lvals = np.linspace(0, domainLength*10**(9), nptsL)
		Hvals = np.linspace(0, domainHeight*10**(9), nptsH)[::-1]

elif units == 'new':

	if style == 'new':
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

#-------------------CALC P and MU--------------------
print uvals.shape

mu = np.zeros((rhovals.shape))
p = np.zeros((rhovals.shape))
rhou = np.zeros((rhovals.shape))
rhov = np.zeros((rhovals.shape))

for i in range(nptsL):
	for j in range(nptsH):
		mu[i][j]= 7.96*pow(10,-10)*pow(rhovals[i][j],2.0) - 1.774*pow(10,-6)*rhovals[i][j] + 0.001106 #from Holland et al 2015 in SI
		p[i][j]= 0.001559*pow(10,-6)*pow(rhovals[i][j],2.0) - 3.387*pow(10,-6)*rhovals[i][j] + 2020.6*pow(10,-6)
		rhou[i][j] = rhovals[i][j]*uvals[i][j]
		rhov[i][j] = rhovals[i][j]*vvals[i][j]

#-------------------DERIVATIVES----------------------

# Initialise arrays

dX = domainLength/nptsL
dY = domainHeight/nptsH

ux = np.zeros((nptsL-2,nptsH-2))
uy = np.zeros((nptsL-2,nptsH-2))
uxx = np.zeros((nptsL-2,nptsH-2))
uyy = np.zeros((nptsL-2,nptsH-2))
uxy = np.zeros((nptsL-2,nptsH-2))

vx = np.zeros((nptsL-2,nptsH-2))
vy = np.zeros((nptsL-2,nptsH-2))
vxx = np.zeros((nptsL-2,nptsH-2))
vyy = np.zeros((nptsL-2,nptsH-2))
vxy = np.zeros((nptsL-2,nptsH-2))

rhoux = np.zeros((nptsL-2,nptsH-2))
rhovy = np.zeros((nptsL-2,nptsH-2))

px = np.zeros((nptsL-2,nptsH-2))
py = np.zeros((nptsL-2,nptsH-2))

mux = np.zeros((nptsL-2,nptsH-2))
muy = np.zeros((nptsL-2,nptsH-2))



# Determine derivatives
for i in range(1,nptsL-1):
	for j in range(1,nptsH-1):
		ux[i-1][j-1] = (uvals[i+1][j]-uvals[i-1][j])/(2*dX)
		uy[i-1][j-1] = (uvals[i][j+1]-uvals[i][j-1])/(2*dY)
		uxx[i-1][j-1] = (uvals[i+1][j] - 2*uvals[i][j] +uvals[i-1][j])/dX**2
		uyy[i-1][j-1] = (uvals[i][j+1] - 2*uvals[i][j] +uvals[i][j-1])/dY**2
		uxy[i-1][j-1] = (uvals[i+1][j+1] - uvals[i+1][j-1] - uvals[i-1][j+1] + uvals[i-1][j-1])/(4*dX*dY)

		vx[i-1][j-1] = (vvals[i+1][j]-vvals[i-1][j])/(2*dX)
		vy[i-1][j-1] = (vvals[i][j+1]-vvals[i][j-1])/(2*dY)
		vxx[i-1][j-1] = (vvals[i+1][j] - 2*vvals[i][j] +vvals[i-1][j])/dX**2
		vyy[i-1][j-1] = (vvals[i][j+1] - 2*vvals[i][j] +vvals[i][j-1])/dY**2
		vxy[i-1][j-1] = (vvals[i+1][j+1] - vvals[i+1][j-1] - vvals[i-1][j+1] + vvals[i-1][j-1])/(4*dX*dY)

		rhoux[i-1][j-1] = (rhou[i+1][j]-rhou[i-1][j])/(2*dX)
		rhovy[i-1][j-1] = (rhov[i][j+1]-rhov[i][j-1])/(2*dY)

		px[i-1][j-1] = (p[i+1][j]-p[i-1][j])/(2*dX)
		py[i-1][j-1] = (p[i][j+1]-p[i][j-1])/(2*dY)

		mux[i-1][j-1] = (mu[i+1][j]-mu[i-1][j])/(2*dX)
		muy[i-1][j-1] = (mu[i][j+1]-mu[i][j-1])/(2*dY)



#----------------------------NS -----------------------------------

RHS_x = np.zeros((nptsL-2,nptsH-2))
RHS_y = np.zeros((nptsL-2,nptsH-2))

LHS_x = np.zeros((nptsL-2,nptsH-2))
LHS_y = np.zeros((nptsL-2,nptsH-2))

cont_eqn = np.zeros((nptsL-2,nptsH-2))

# Compute actual NS equations
for i in range(nptsL-2):
	for j in range(nptsH-2):
		RHS_x[i][j] = -px[i][j] + mu[i+1][j+1]*(uxx[i][j]+uyy[i][j]) + mux[i][j]*(ux[i][j] + vy[i][j])
		RHS_y[i][j] = -py[i][j] + mu[i+1][j+1]*(vxx[i][j]+vyy[i][j]) + muy[i][j]*(ux[i][j] + vy[i][j])

		LHS_x[i][j] = rhovals[i+1][j+1] * (uvals[i+1][j+1]*ux[i][j] + vvals[i+1][j+1]*uy[i][j])
		LHS_y[i][j] = rhovals[i+1][j+1] * (uvals[i+1][j+1]*vx[i][j] + vvals[i+1][j+1]*vy[i][j])
		cont_eqn[i][j] = rhoux[i][j] + rhovy[i][j]

		print 'x: ', RHS_x[i][j] - LHS_x[i][j],' y: ',  RHS_y[i][j] - LHS_y[i][j],' cont: ', cont_eqn[i][j]

RHS_x_sum = np.sum(RHS_x)
LHS_x_sum = np.sum(LHS_x)
RHS_y_sum = np.sum(RHS_y)
LHS_y_sum = np.sum(LHS_y)
cont_eqn_sum = np.sum(cont_eqn)


print 'Sums'
print 'x: ', RHS_x_sum - LHS_x_sum,' y: ',  RHS_y_sum - LHS_y_sum,' cont: ', cont_eqn_sum





