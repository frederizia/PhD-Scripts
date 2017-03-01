#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
from scipy.integrate import simps

def read_data():
    print 'Reading in data...'
    filename = "vel.dat"
    
    # read in data
    f = open(filename,'r')
    data = f.read()
    
    data_lines = data.split('\n')
    
    data_list = []
    data_tmp = []

    A2m = 1e-10
    ps2s = 1e-12

    Natoms = 0
    count = 0
    for i in range(9,len(data_lines)-1): 
        if count < 1000:     
            try:
                if data_lines[i]=='ITEM: TIMESTEP':
                    if data_tmp != 0:
                        data_list.append(np.array(data_tmp))
                    data_tmp = []
                    count += 1
                    #print 'Reading in timestep:', count
                #if statement for coords
                if len(data_lines[i].split()) == 8:
                    Natoms += 1
                    vx = float(data_lines[i].split()[5])*A2m/ps2s
                    vy = float(data_lines[i].split()[6])*A2m/ps2s
                    vz = float(data_lines[i].split()[7])*A2m/ps2s
                    vel = np.array([vx, vy, vz])
                    data_tmp.append(vel)
            except:
                continue

    data_list = np.array(data_list)
    print 'Read in', count, 'timesteps.'

    return data_list

def correlation(data):
    vacf = []
    times_range = 40
    print 'Calculating the autocorrelation function for', times_range*100*0.0005, 'ps.'
    for i in range(data.shape[0]-times_range):
        Cdata = []
        for j in range(times_range+1):
            Cdelta = []
            C_tmp = 0
            count = 0
            for k in range(data[i].shape[0]):
                count += 1
                C_tmp += np.dot(data[i][k], data[i+j][k])
                
            C_tmp = C_tmp/count
            Cdata.append(C_tmp)
        vacf.append(Cdata)
    vacf = np.array(vacf)
    mean_vacf = np.mean(vacf, axis=0)
    norm_vacf = mean_vacf/mean_vacf[0]
    print vacf
    print mean_vacf
    print norm_vacf

    return norm_vacf

def time_data(Cdata):
    print 'Creating corresponding time array.'
    times = []
    # get corresponding timesteps
    for i in range(len(Cdata)):
        t_tmp = i*100*0.0005 # data sampled every 100, dt=0.0005ps
        times.append(t_tmp)

    return times



def find_region(dist, pot):
    dist, pot = np.array(dist), np.array(pot)

    mask = np.logical_and(dist > 0, dist < 4*10**(-10))
    new_dist = dist[mask]
    new_pot = pot[mask]
    new_Eint = sum(new_pot)
    return new_Eint, new_pot, new_dist


def main():
    DAT = read_data()
    CDAT = correlation(DAT)
    TDAT = time_data(CDAT)

    plt.figure()
    plt.plot(TDAT, CDAT)
    plt.show()



    return


if __name__ == "__main__":
    main()


