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
                    help="Delta x",required=True)


args = parser.parse_args()
vel = args.v



def read_data(V):
	'''Code to read velocity data from a 2d LAMMPS output file'''


	filename = 'temp_rescale_Nature1997/vel.nemd_H22.82_rhof0.96_v%s'%(V)


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
	tol = 0.03

	for k in range(1,mid):
		print abs(var[k-1]-var[k]) 
        	if abs(var[k-1]-var[k]) >= tol:
			print 'left'
            
            		LEFT = k
            		break
	for k in range(1,mid): 
		print abs(var[::-1][k-1]-var[::-1][k]), tol        
        	if abs(var[::-1][k-1]-var[::-1][k]) >= tol:
            
			print 'right'            		
			RIGHT = len(ycoord)-k
			break


	buff = -4
	LEFT_minus = LEFT - buff
	RIGHT_plus = RIGHT + buff

	ycoord = ycoord[LEFT_minus:RIGHT_plus]
	height = ycoord[-1]-ycoord[0]

	var = var[LEFT_minus:RIGHT_plus]


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
	

#---------------PLOTTING------------------

matplotlib.rcParams.update({'font.size': 15})
matplotlib.rc('text', usetex=True)

plt_name = 'Nature1997'

slipL = []
shearRate = []
fig = plt.figure()
for v in vel:
	Y, VelDistr, H = read_data(v)
	vWall, dvWall = derivative(VelDistr,Y)
	Ls = vWall/dvWall

	shearRate.append(float(v)/H)
	slipL.append(Ls)	

	plt.plot(Y, VelDistr, label = '$V = %s$'%(v))
#plt.legend(loc='upper left')
plt.xlim(2,22)
plt.savefig('vel_%s.png'%plt_name)
plt.show()

fig = plt.figure()
plt.plot(np.log10(shearRate), slipL)
plt.ylabel('$L_s$')
plt.xlabel('$\mathrm{log}_{10}(\dot{\gamma})$')
plt.savefig('slip_%s.png'%plt_name)
plt.show()


