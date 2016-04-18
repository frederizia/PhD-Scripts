#!/usr/bin/env python
from __future__ import division
import numpy as np

def read_data(h,f,eps,rhof,md,VAR,extra):
	'''Code to read in density and velocity data from a 2d LAMMPS output file'''

	if extra != 'None':	
		filepath = 'H%s/eps%s/%s'%(h,eps,extra)
	else:
		filepath = 'H%s/eps%s'%(h,eps)


	if md == 'nemd':
		filename = '%s/dens.nemd_xy_H%s_rhof%s_f%s'%(filepath,h,rhof,f)
	if md == 'equ':
		filename = '%s/dens.eq_xy_H%s_rhof%s'%(filepath,h,rhof)

	f = open(filename,'r')
	data = f.read()
    	data_lines = data.split('\n')

    	xcoord = []
    	xcoord_tmp =[]
    	ycoord = []
    	ycoord_tmp =[]
    	var = []
    	var_tmp = []

	idx = 0
	count  = 0
	count_tmp = 1

	lowerY = 4.0

	if VAR == 'rho': nmbr = 4
	if VAR == 'u': nmbr = 5
    
	for j in range(4,len(data_lines)-1):
       		if len(data_lines[j].split()) == 2:
				# Count number of timesteps collected
			count += 1
       		elif len(data_lines[j].split()) != 2:
            		x_val = float(data_lines[j].split()[1])
            		y_val = float(data_lines[j].split()[2])
            		VAR_val = float(data_lines[j].split()[nmbr])

			if count ==1 and x_val > 11 and x_val < 38 and y_val > lowerY and y_val < (lowerY+float(h)+2.5):
				# Only collect data in channel region
            			var.append(VAR_val)
            			xcoord.append(x_val)
            			ycoord.append(y_val)

			elif count >1 and x_val > 11 and x_val < 38 and y_val > lowerY and y_val < (lowerY+float(h)+2.5):
				# Average over all timesteps
				if count != count_tmp:
					idx = 0
            			var[idx] = (VAR_val+var[idx])/2
				idx +=1	
				count_tmp = count	



    	xcoord = np.unique(np.array(xcoord))
    	ycoord = np.unique(np.array(ycoord))


    	VAR_array = np.zeros((len(xcoord),len(ycoord)))

    
    
    	count = 0
    	for i in range(len(xcoord)):
        	for j in range(len(ycoord)):
            		VAR_array[i][j] = var[count]

            		count +=1




	return xcoord, ycoord, VAR_array

def mid_point(X,Y,VAR):
	'''Code to find the mid point in the array'''

	# Find mid point
	VAR_avg_Y = np.average(VAR, axis=1) # use average value to find mid point
	mid = int(len(Y)/2)


	for k in range(1,mid):
        	if abs(VAR_avg_Y[k-1]-VAR_avg_Y[k]) >= 0.03:
            
            		LEFT = k
            		break
	for k in range(1,mid):        
        	if abs(VAR_avg_Y[::-1][k-1]-VAR_avg_Y[::-1][k]) >= 0.03:
            
            		RIGHT = len(Y)-k
			break

	mid_p = int(LEFT+(RIGHT - LEFT)/2)
	return mid_p, LEFT, RIGHT

def symmetry(X,Y,VAR):
	'''Code to take MD data and use the symmetry along the axis to average over more values'''

	mid_p, LEFT, RIGHT = mid_point(X,Y,VAR)

	# Account for averaging used in finding LEFT and RIGHT and criteria
	buff = int(LEFT/2)
	LEFT_minus = LEFT - buff
	RIGHT_plus = RIGHT + buff
	mid_p_new = mid_p-(LEFT_minus)

	VAR_tmp = np.copy(VAR)
	Y_tmp = Y[LEFT_minus:RIGHT_plus]

	VAR_tmp = VAR_tmp[LEFT_minus:RIGHT_plus,:]

	for i in range(len(X)):
		for j in range(mid_p_new):
			VAR_i = VAR_tmp[:,i]
			VAR_tmp[j][i] = (VAR_i[j]+VAR_i[::-1][j])/2
	if len(Y)%2 == 0:
		VARnew = VAR_tmp[:mid_p_new,:]
		Ynew = Y_tmp[:mid_p_new]
	else:
		VARnew = VAR_tmp[:mid_p_new+1,:]
		Ynew = Y_tmp[:mid_p_new+1]
	return Ynew, VARnew, Y[LEFT], Y[RIGHT]


def rho_init(rhovals):
	'''Compute simple profile for computation comparison'''

	# initialise arrays
	rho_init = np.zeros(rhovals.shape) # i=y, j=x
	P_init = np.zeros(rhovals.shape) # i=y, j=x

	D1=0.03233
	D2=6.049

	rho_in = np.mean(rhovals[:,0:10])
	rho_out = np.mean(rhovals[:,-10:-1])

	P_in = P17(rho_in)
	P_out = P17(rho_out)
	print 'P_in: ', P_in
	print 'P_out: ',P_out
	

	for j in range(rhovals.shape[1]):
		Pval = ((P_out-P_in)/rhovals.shape[1])*j + P_in

		rhoval = np.log(Pval/D1)/D2
		rho_init[:,j] = rhoval

	return rho_init

	
def P17(rho):

	D1=0.03233
	D2=6.049
	return D1*np.exp(D2*rho)






