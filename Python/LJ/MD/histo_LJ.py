#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import os


# Parse arguments form command line
parser = argparse.ArgumentParser()

parser.add_argument("-type", type=str, nargs=1, \
                    help="com or xyz",required=True)




args = parser.parse_args()
TYPE = args.type[0]





def hist_data(N):

    cwd = os.getcwd()

    fname = '%s/pmf_pos_N%i.dat'%(cwd, N)
    
    # Define correction value to adjust periodic boundary conditions
    # Find correct position
    N_corr = N/10
	k_spring = 300
    #timesteps = 

    # read in data
    f = open(fname,'r')
    data = f.read()
    data_lines = data.split('\n')

    steps = []

    x_data = []
    x_diff = []

    count=0

                    
    count2 = 0
    for j in range(0,len(data_lines)-1):
   		 #print data_lines[j]
    	pmf_data = data_lines[j].split()
    	#print pmf_data
	#print j, j%10
	
		
	if pmf_data[0] == '1' and len(pmf_data) >1:
		count2 +=1
		if count2 < timesteps:
                    	zlo = 0.0
                    	zhi = 6.37002
                    	dz = zhi - zlo
                    	z = float(pmf_data[3])
                    	delz = z - N_corr

                #print pmf_data[3]
			count += 1
			steps.append(count)                        
                    	if delz > 20: z -= dz         
                    	if delz < -20: z += dz
        		z_data.append(z)
			z_diff.append(float(pmf_data[3])-N_corr)
			x_data.append(float(pmf_data[1]))



    print "Number of timesteps collected: ", count-1
    
    reduced_data_file = np.array([steps, z_data]).transpose()
    reduced_data_file_gwham = np.array([steps, z_diff]).transpose()
    reduced_data_file_2d = np.array([steps, x_data, z_data]).transpose()

    gwham_header = "UMBRELLA  3.0\nComponent selection: 0 0 1\nnSkip 1\nRef. Group 'TestAtom'\nNr. of pull groups 1\nGroup 1 'H2O' Umb. Pos. %i Umb. Cons. %f\n#####" %(N, k_spring)

    np.savetxt('LJ/data_LJ_N%i_%i'%(N,run),reduced_data_file, fmt=['%3i','%.18f'])
    np.savetxt('LJ/gwham_LJ_N%i_%i.pdo'%(N,run),reduced_data_file_gwham, fmt=['%3i','%.18f'], header=gwham_header)
    #np.savetxt('delx%s/data_CH_delx9_N%i_%i_2d'%(dx,N,run),reduced_data_file_2d, fmt=['%3i','%.18f','%.18f'])

    return z_data



# PLOTTING

matplotlib.rcParams.update({'font.size': 19})

data = []

num_list = range(1,23)
num_list.extend(range(24,32))

for i in num_list:#(1,101)
		j = 3

		

		print i, j
		data = hist_data(TYPE,i,j)


		fig = plt.figure()
		ax = fig.add_subplot(111)

		plt.xlabel('$z$ ($\AA$)')
		plt.ylabel('hist')
		#plt.plot(z_bins, z_frames)
		ax.hist(data, 50, facecolor='b', alpha=0.5)
		plt.xlim(0,6.5)
		plt.title('z=%i, run=%i'%(i*0.2,j))

		plt.savefig('LJ/hist_LJ_N%i_%i_%s.pdf'%(i,j,TYPE))
		#plt.show()



def main():
    args = Get_Args()
    filename = args.f
    index = args.i
    label = args.l

    # read file and extract values
    timesteps, values = read_file(filename, index)

    # disentangle values
    pure_values = split_values(values)

    # calculate error
    error = error_calc(pure_values)

    print 'The final averaged value is', values[-1]
    print 'The mean is', pure_values.mean()
    print 'The error is', error

    write_file(timesteps, values, pure_values, label, index)

    return

if __name__ == '__main__':
    main()
