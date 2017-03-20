#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
from scipy.integrate import simps
import datetime

def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', required=False, default='vel.dat', type=str, action='store',
                       help='Filename')
    parser.add_argument('-t', '--tsteps', required=False, default=5000, type=int, action='store',
                       help='Number of timesteps')
    parser.add_argument('-z', '--dz', required=False, type=int, nargs=2, action='store', default=[0,100],
                       help='zmin and zmax')
    args = parser.parse_args()
    return args


def read_data(fn, zmin,zmax, timesteps):
    print 'Reading in data...'
    filename = fn
    
    # read in data
    f = open(filename,'r')
    data = f.read()
    
    data_lines = data.split('\n')
    
    data_list = []
    data_tmp = []
    flag_list = []
    flag_tmp = []

    A2m = 1e-10
    fs2s = 1e-15

    count = 0
    for i in range(9,len(data_lines)-1): 
        if count < timesteps:     
            try:
                if data_lines[i]=='ITEM: TIMESTEP':
                    if data_tmp != []:
                        data_list.append(np.array(data_tmp))
                        flag_list.append(np.array(flag_tmp))
                    data_tmp = []
                    flag_tmp = []
                    count += 1
                    #print 'Reading in timestep', count
                #if statement for coords
                if len(data_lines[i].split()) == 8:
                    zcoord = float(data_lines[i].split()[4])
                    if zmin < zcoord < zmax:
                        flag = 0
                    else:
                        flag = 1
                    vx = float(data_lines[i].split()[5])*A2m/fs2s
                    vy = float(data_lines[i].split()[6])*A2m/fs2s
                    vz = float(data_lines[i].split()[7])*A2m/fs2s
                    vel = np.array([vx, vy, vz])
                    data_tmp.append(vel)
                    flag_tmp.append(flag)
            except:
                continue

    data_list, flag_list = np.array(data_list), np.array(flag_list)
    print 'Read in', count, 'timesteps.'
    #print flag_list

    return data_list, flag_list

def correlation(data, flag, tsteps):
    vacf = []
    times_range = int(tsteps/10)
    flag_count = np.zeros((1,data[0].shape[0]))
    flag_count_prev = np.zeros((1,data[0].shape[0]))


    print 'Calculating the autocorrelation function for', times_range*2, 'fs.'
    for i in range(data.shape[0]-times_range):
        Cdata = []
        for j in range(times_range+1):
            Cdelta = []
            C_tmp = 0
            count = 0
            for k in range(data[i].shape[0]):
                flag_count[0][k]        += flag[i+j][k]
                #if (flag_count[0][k] < 1:
                if (flag_count_prev[0][k] == 0) and (flag_count[0][k] == 0):
                    C_tmp += np.dot(data[i][k], data[i+j][k])
                    count += 1
                flag_count_prev[0][k]   = flag_count[0][k]
            if count == 0:
                #print 'No eligible points for this origin.'
                C_tmp = 0
            else:
                C_tmp = C_tmp/count
                Cdata.append(C_tmp)
        if Cdata != []:
            vacf.append(Cdata)
    vacf = np.array(vacf)
    print vacf
    mean_vacf = np.mean(vacf, axis=0)
    norm_vacf = mean_vacf/mean_vacf[0]

    return mean_vacf, norm_vacf

def time_data(Cdata):
    print 'Creating corresponding time array.'
    fs2s = 1e-15
    times = []
    times_si = []
    # get corresponding timesteps
    for i in range(len(Cdata)):
        t_tmp = i*2 # data sampled every 1, dt=2fs
        times.append(t_tmp)
        times_si.append(t_tmp*fs2s)

    return times, times_si


def diffusion(vacf, times):
    fs2s = 1e-15
    dt = 2*fs2s

    int_vacf = np.sum(vacf)*dt
    diff = (1/3)*int_vacf

    diff2 = simps(vacf, times)/3

    print 'The self-diffusion coefficient is', diff, 'or', diff2

    return diff

def plotting(tdat, cdat, zmin, zmax, tsteps):
    matplotlib.rcParams.update({'font.size': 19})
    matplotlib.rc('text', usetex=True)

    plt.figure()
    plt.plot(tdat, cdat)
    plt.ylabel('$\Psi(t)$')
    plt.xlabel('$t (fs)$')
    plt.xlim(0,500)
    plt.savefig('vacf_%s_%sA_%its.pdf'%(zmin, zmax,tsteps))
    plt.show()

    return

def main():
    args = GetArgs()
    file = args.filename
    zmin, zmax = args.dz
    tsteps = args.tsteps


    print 'The start time is:', datetime.datetime.now()

    print 'Calculate the VACF between', zmin, 'and', zmax, 'Angstrom.'

    try:
        TDAT        = np.loadtxt('vacf_time_%s_%sA_%its.dat'%(zmin, zmax,tsteps))
        TDAT_SI     = np.loadtxt('vacf_timesi_%s_%sA_%its.dat'%(zmin, zmax,tsteps))
        CDAT        = np.loadtxt('vacf_corr_%s_%sA_%its.dat'%(zmin, zmax,tsteps))
        CDAT_norm   = np.loadtxt('vacf_corrnorm_%s_%sA_%its.dat'%(zmin, zmax,tsteps))
        print 'Data already created and read in from file.'

    except IOError:
        DAT, FLAG       = read_data(file, zmin, zmax, tsteps)
        CDAT, CDAT_norm = correlation(DAT, FLAG, tsteps)
        TDAT, TDAT_SI   = time_data(CDAT)

        # save data
        np.savetxt('vacf_time_%s_%sA_%its.dat'%(zmin, zmax,tsteps), TDAT, fmt="%s")
        np.savetxt('vacf_timesi_%s_%sA_%its.dat'%(zmin, zmax,tsteps), TDAT_SI, fmt="%s")
        np.savetxt('vacf_corr_%s_%sA_%its.dat'%(zmin, zmax,tsteps), CDAT, fmt="%s")
        np.savetxt('vacf_corrnorm_%s_%sA_%its.dat'%(zmin, zmax,tsteps), CDAT_norm, fmt="%s")

    except:
        print 'An error occured.'

    diffusion(CDAT, TDAT_SI)
    plotting(TDAT, CDAT_norm, zmin, zmax, tsteps)

    print 'The end time is:', datetime.datetime.now()


    return


if __name__ == "__main__":
    main()


