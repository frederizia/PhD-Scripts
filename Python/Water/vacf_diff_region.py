#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
from scipy.integrate import simps
import datetime
import sys
#from plotting_params import *

def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', required=False, default='vel.dat', type=str, action='store',
                       help='Filename')
    parser.add_argument('-t', '--tsteps', required=False, default=5000, type=int, action='store',
                       help='Number of timesteps')
    parser.add_argument('-z', '--dz', required=False, type=int, nargs=2, action='store', default=[0,100],
                       help='zmin and zmax')
    parser.add_argument('-u', '--units', required=False, type=str, action='store', default=[0,100],
                       help='metal or real')
    parser.add_argument('-dt', required=False, type=float, action='store', default=2,
                       help='timestep')
    args = parser.parse_args()
    return args


def read_data(fn, zmin,zmax, timesteps, tc):
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

    count = 0
    for i in xrange(9,len(data_lines)-1): 
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
                    vx = float(data_lines[i].split()[5])*A2m/tc
                    vy = float(data_lines[i].split()[6])*A2m/tc
                    vz = float(data_lines[i].split()[7])*A2m/tc
                    vel = np.array([vx, vy, vz])
                    data_tmp.append(vel)
                    flag_tmp.append(flag)
            except:
                continue

    data_list, flag_list = np.array(data_list), np.array(flag_list)
    print 'Read in', count, 'timesteps.'
    data_list = np.ma.array(data_list)
    data_list[flag_list == 1] = np.ma.masked

    return data_list, flag_list

def correlation(data, flag, tsteps, dt, unit):
    vacf = []
    times_range = int(tsteps/2)
    flag_count = np.zeros((1,data[0].shape[0])) # size of number of atoms
    flag_count_prev = np.zeros((1,data[0].shape[0]))


    print 'Calculating the autocorrelation function for', times_range*dt, unit
    # iterate over different time origins
    for i in xrange(data.shape[0]-times_range): 
        Cdata = [] # stores correlation for each time origin
        for j in xrange(times_range+1): # the n in t_0+ndt
            C_tmp = 0
            count = 0
            # cycle through atoms
            for k in xrange(data[i].shape[0]):
                #flag_count[0][k]        += flag[i+j][k]
                # check if both current and previous timestep was in the zone
                #if (flag_count_prev[0][k] == 0) and (flag_count[0][k] == 0):
                C_tmp += np.dot(data[i][k], data[i+j][k])
                count += 1 # counts number of elegible atoms
                #flag_count_prev[0][k]   = flag_count[0][k]
            if count == 0:
                print 'At', (i+j)*dt, unit, 'no atoms were eligible in the desired region.'
                C_tmp = 0
            else:
                C_tmp = C_tmp/count # C_v(t)
            Cdata.append(C_tmp)
        if Cdata != []:
            vacf.append(Cdata)
    vacf = np.array(vacf)
    
    #mean_vacf = np.apply_along_axis(lambda v: np.mean(v[np.nonzero(v)]), 0, vacf)
    #mean_vacf[np.isnan(mean_vacf)]=0.
    mean_vacf = np.mean(vacf, axis=0)
    norm_vacf = mean_vacf/mean_vacf[0]

    return mean_vacf, norm_vacf

def time_data(Cdata,tc, dt):
    print 'Creating corresponding time array.'


    times = []
    times_si = []
    # get corresponding timesteps
    for i in xrange(len(Cdata)):
        t_tmp = i*dt # data sampled every 1, dt=2fs
        times.append(t_tmp)
        times_si.append(t_tmp*tc)

    return times, times_si


def diffusion(vacf, times,tc,dt):

    DT = dt*tc
    

    int_vacf = np.sum(vacf)*DT
    diff = (1/3)*int_vacf
    #print diff


    #diff2 = simps(vacf, times)/3

    print 'The self-diffusion coefficient is', diff#, 'or', diff2

    return diff

def plotting(tdat, cdat, zmin, zmax, tsteps, u):
    matplotlib.rcParams.update({'font.size': 19})
    matplotlib.rc('text', usetex=True)

    plt.figure()
    plt.plot(tdat, cdat)
    plt.ylabel('$\Psi(t)$')
    if u == 'metal':
        plt.xlabel('$t (ps)$')
        plt.xlim(0,0.2)
    if u == 'real':
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
    dt = args.dt
    units = args.units
    if units == 'metal':
        time_conv = 1e-12
        unit = 'ps'
    elif units == 'real':
        time_conv = 1e-15
        unit = 'fs'
    else:
        print 'Invalid unit'
        sys.exit(1)


    print 'The start time is:', datetime.datetime.now()

    print 'Calculate the VACF between', zmin, 'and', zmax, 'Angstrom.'

    try:
        TDAT        = np.loadtxt('vacf_time_%s_%sA_%its.dat'%(zmin, zmax,tsteps))
        TDAT_SI     = np.loadtxt('vacf_timesi_%s_%sA_%its.dat'%(zmin, zmax,tsteps))
        CDAT        = np.loadtxt('vacf_corr_%s_%sA_%its.dat'%(zmin, zmax,tsteps))
        CDAT_norm   = np.loadtxt('vacf_corrnorm_%s_%sA_%its.dat'%(zmin, zmax,tsteps))
        print 'Data already created and read in from file.'

    except IOError:
        DAT, FLAG       = read_data(file, zmin, zmax, tsteps, time_conv)
        CDAT, CDAT_norm = correlation(DAT, FLAG, tsteps, dt, unit)
        TDAT, TDAT_SI   = time_data(CDAT, time_conv, dt)

        # save data
        np.savetxt('vacf_time_%s_%sA_%its.dat'%(zmin, zmax,tsteps), TDAT, fmt="%s")
        np.savetxt('vacf_timesi_%s_%sA_%its.dat'%(zmin, zmax,tsteps), TDAT_SI, fmt="%s")
        np.savetxt('vacf_corr_%s_%sA_%its.dat'%(zmin, zmax,tsteps), CDAT, fmt="%s")
        np.savetxt('vacf_corrnorm_%s_%sA_%its.dat'%(zmin, zmax,tsteps), CDAT_norm, fmt="%s")

    except:
        print 'An error occured.'

    diffusion(CDAT, TDAT_SI, time_conv, dt)
    plotting(TDAT, CDAT_norm, zmin, zmax, tsteps, units)

    print 'The end time is:', datetime.datetime.now()


    return


if __name__ == "__main__":
    sys.path.insert(0,'/home/fred/SCRIPTS/Python')
    main()


