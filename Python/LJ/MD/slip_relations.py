#!/usr/bin/env python

'''Code to determine the slip length as a function of shear rate'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
from scipy.interpolate import splrep, splev

# Parse arguments form command line
parser = argparse.ArgumentParser()
parser.add_argument("-v", type=str, nargs='+', \
                    help="Velocity",required=True)
parser.add_argument("-n", type=str, nargs=1, \
                    help="Epsilon, sigma  or type",required=True)
parser.add_argument("-t", type=str, nargs=1, \
                    help="thermostat",required=True)
parser.add_argument("-rhos", type=str, nargs=1, \
                    help="Rhos",required=False, default = 'None')

args = parser.parse_args()
vel = args.v
name = args.n[0]
thermo = args.t[0]
rhos = args.rhos[0]

plt_name = '%s_%s'%(name, thermo) 
if name != 'Nature1997':
	plt_name = 'rhos%s_%s_%s'%(rhos,name, thermo)


def read_data(V):
	'''Code to read velocity data from a 2d LAMMPS output file'''

	if name == 'Nature1997':
		H = 22.82
		rhof = 0.81
		filename = '%s/%s/vel.nemd_H%s_rhof%s_v%s'%(name,thermo,H,rhof,V)
	else:
		H = 22.82
		rhof = 0.5
		filename = 'rhos%s/%s/%s/vel.nemd_H%s_rhof%s_v%s'%(rhos,name,thermo,H,rhof,V)

	f = open(filename,'r')
	data = f.read()
    	data_lines = data.split('\n')
	

    	ycoord = []
    	ycoord_tmp =[]
    	var = []
    	var_tmp = []

	idx = 0
	count  = 0
	count_tmp = 1



	nmbr = 3
	
	for j in range(4,len(data_lines)-1):
       		if len(data_lines[j].split()) == 2:
				# Count number of timesteps collected
			count += 1
       		elif len(data_lines[j].split()) != 2:

            		y_val = float(data_lines[j].split()[1])
            		VAR_val = float(data_lines[j].split()[nmbr])
			if count ==1:
				# Only collect data in channel region
           			var.append(VAR_val)
            			ycoord.append(y_val)

			elif count >1:
				# Average over all timesteps
				if count != count_tmp:
					idx = 0
            			var[idx] = (VAR_val+var[idx])/2
				idx +=1	
				count_tmp = count	


    	ycoord = np.unique(np.array(ycoord))


	# Cut array to get rid of wall region
	
	mid = int(len(ycoord)/2)
	toll = 0.1
	tolr = 0.2
	for k in range(mid,len(ycoord)):
		#print abs(var[::-1][k-1]-var[::-1][k]) 
        	if abs(var[::-1][k-1]-var[::-1][k]) >= toll:
			print 'left'
            
            		LEFT = len(ycoord)-k
            		break
	for k in range(mid,len(ycoord)): 
		#print abs(var[::-1][k-1]-var[::-1][k]), tolr        
        	if abs(var[k-1]-var[k]) >= tolr:
            
			print 'right'            		
			RIGHT = k
			break


	buff = -1

	# redefine LEFT for this setup
	for i in range(len(ycoord)):
		if ycoord[i] > 13:
			LEFT = i
			break
        for i in range(len(ycoord)):
                if ycoord[i] > 36:
                        RIGHT = i
                        break	
	LEFT_minus = LEFT - buff
	RIGHT_plus = RIGHT + buff

	ycoord = ycoord[LEFT_minus:RIGHT_plus]
	height = ycoord[-1]-ycoord[0]

	var = var[LEFT_minus:RIGHT_plus]
	ycoord -=ycoord[0]


	return ycoord, var, height

def derivative(V, y):

	#y = y[-10:]
	#V = V[-10:]
	
	# Interpolate spline to smooth function and find derivative
	f = splrep(y,V,k=3,s=3)
	fitV = splev(y,f)
	dv1 = splev(y,f,der=1)
	b1 = fitV[-1] - dv1[-1]*y[-1]

	'''fig = plt.figure()
	plt.plot(y, fitV, label="fitted")
	#plt.plot(y, dv1, label="1st derivative")
	plt.plot(y, dv1[-1]*y+b1,label='1')
	plt.plot(y,V, label='profile')
	plt.xlim(10,22)
	#plt.ylim(0,6)
	plt.legend(loc='lower left')
	plt.show()'''



	return fitV[-1], dv1[-1]
	
def slipLength(V, y,v_w):
	# velocity at stationary wall
	#U_f = V[0]
	# Velocity at the wall
	U_f = V[-1]

	# wall velocity
	U_w = float(v_w)
	# gradient assuming straight line profile
	grad = (V[-1]-V[0])/(y[-1]-y[0])

	L_s = (U_w - U_f)/grad
	#L_s = U_f/grad
	return L_s

#---------------PLOTTING------------------

matplotlib.rcParams.update({'font.size': 15})
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)


slipL = []
slipL_eval = []
shearRate = []
fig = plt.figure()
for v in vel:
	Y, VelDistr, H = read_data(v)
	vWall, dvWall = derivative(VelDistr,Y)
	#Ls = vWall/dvWall
	Ls = slipLength(VelDistr,Y,v)
	shearRate.append(float(v)/H)
	slipL.append(Ls)	
	if 2 <= float(v) < 8:
		slipL_eval.append(Ls)
	plt.plot(Y, VelDistr, label = '$v = %s$'%(v),lw=2.0)
print plt_name
plt.legend(loc='upper left')
plt.xlabel('Y/$\sigma$')
plt.ylabel('v/ ($\sigma$/$\\tau$)')
plt.xlim(np.min(Y),np.max(Y))
plt.ylim(0,10)
plt.savefig('PLOTS/vel_%s.png'%plt_name)
plt.show()

# Calculate the average slip length with error
slipL_eval = np.array(slipL_eval)
averageLs = np.mean(slipL_eval)
error = np.std(slipL_eval)/np.sqrt(len(slipL_eval))

print 'The slip length is: ', averageLs, ' +/- ', error

# print to file
data_list = [vel, slipL]
data_array = np.transpose(np.array(data_list))
print data_array
np.savetxt('DATA/slip_%s.dat'%(plt_name), data_array,fmt="%s",header = 'Velocity L_s')


fig = plt.figure()
plt.plot(np.log10(shearRate), slipL)
plt.ylabel('$L_s$')
plt.xlabel('$\mathrm{log}_{10}(\dot{\gamma})$')
plt.savefig('PLOTS/slip_%s.png'%plt_name)
plt.show()

fig = plt.figure()
plt.plot(vel, slipL, linestyle = 'None', marker='D')
plt.ylabel('$L_s$/$\sigma$')
plt.xlabel('v/ ($\sigma$/$\\tau$)')
plt.savefig('PLOTS/slip_%s.png'%plt_name)
plt.show()

