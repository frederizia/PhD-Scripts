#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import sys
'''Script to generate histogram data from Umbrella sampling'''

def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-delx", type=str, nargs=1, \
                        help="Delta x",required=True)
    parser.add_argument("-func", type=str, nargs=1, \
                        help="CH or COO",required=False, default='COO')
    parser.add_argument("-mol", type=str, nargs=1, \
                        help="molecule",required=False, default='200')
    parser.add_argument("-k", type=str, nargs=1, \
                        help="spring constant",required=False, default='1.2')
    parser.add_argument("-dz", type=str, nargs=1, \
                        help="step size",required=True, default='0.2')
    parser.add_argument("-n", type=int, nargs=2, \
                        help="Nmin, Nmax",required=True)
    args = parser.parse_args()
    return args


def hist_data(fn, func, dx, dz, N, run, corr_flag):
    
    N_corr = N*float(dz)
    timesteps=500000

    #if dx == '5':
    #    if dz == '1.0' and k_spring == '0.5':
    #    	timesteps = 3000000
    #    else:
    #    	timesteps = 1000000
    #if dx == '9':
	#   timesteps = 1000000
    
    # read in data
    f = open(fn,'r')
    data = f.read()


    data_lines = data.split('\n')

    steps = []

    z_data = []
    z_diff = []
    x_data = []
    count=0

        

    for j in range(0,len(data_lines)-1):
        pmf_data = data_lines[j].split()
        if pmf_data[0] == 'O':

            zlo = 0.0
            zhi = 100
            d_z = zhi - zlo
            z = float(pmf_data[3])
            delz = z - N_corr

            #print pmf_data[3]
            count += 1
            steps.append(count)                        
            if delz > 20: z -= d_z         
            if delz < -20: z += d_z
            z_data.append(z)
            z_diff.append(float(pmf_data[3])-N_corr)
            x_data.append(float(pmf_data[1]))



    print "Number of timesteps collected: ", count-1

    if corr_flag ==1:
        # Correlation function
        corr = np.correlate(z_data,z_data, mode='same')
        return_val = [z_data, corr]
    else:
        return_val = z_data
    
    reduced_data_file = np.array([steps, z_data]).transpose()

    np.savetxt('HIST/data_%s_delx%s_N%i_%i'%(func,dx,N,run),reduced_data_file, fmt=['%3i','%.18f'])
    #np.savetxt('delx%s/data_CH_delx9_N%i_%i_2d'%(dx,N,run),reduced_data_file_2d, fmt=['%3i','%.18f','%.18f'])

    return return_val


def main():
    args = GetArgs()
    DELX = args.delx[0]
    DX = DELX[0]
    FUNC = args.func[0]
    MOL = args.mol[0]
    K = args.k[0]
    DZ = args.dz[0]
    Nmin = args.n[0]
    Nmax = args.n[1]

    run = 3
    corr_flag = 0

	# PLOTTING


    data = []

    # Loop over all collected data
    for i in range(int(Nmin/float(DZ)),int(Nmax/float(DZ))+1):#(50,151):

        filename = "pmf_positions_%s_delx%s_N%i_%i.xyz" % (FUNC, DELX,i,run)
        try:
            DATA = hist_data(filename, FUNC, DELX,DZ,i,run, corr_flag)

            if corr_flag == 0:

                fig = plt.figure()
                ax = fig.add_subplot(111)

                plt.xlabel('$z$ ($\AA$)')
                plt.ylabel('hist')
                #plt.plot(z_bins, z_frames)
                ax.hist(DATA, 50, facecolor='b', alpha=0.5)
                plt.xlim(0,100)
                plt.title('z=%i, run=%i'%(i,run))

                plt.savefig('hist_delx9_N%i_%i.pdf'%(i, run))
                plt.close(fig)

            else:
                fig = plt.figure()
                ax = fig.add_subplot(111)

                plt.xlabel('time')
                plt.ylabel('Autocorrelation')
                #plt.plot(z_bins, z_frames)
                ax.hist(data, 50, facecolor='b', alpha=0.5)
                #plt.xlim(0,100)
                plt.title('z=%i'%(i))

                plt.savefig('delx%s/mol_%s_dz_%s_k_%s_%s_%s/corr_delx%s_N%i_%i_%s.pdf'%(DELX,MOL,DZ,K,Nmin,Nmax,DELX,i,j,TYPE))
                #plt.show()
        except IOError:
        	print filename
        	print 'File does not exist.'
        	continue
    return

if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS')
    main()


