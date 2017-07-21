#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
import sys
from scipy.integrate import simps


def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', nargs='+', required=False, default='1', action='store',
                       help='File names')
    args = parser.parse_args()
    return args
    
def read_log(filename):
    '''Code to read in log file'''
    print('Reading in data for', filename)

    f = open(filename,'r')

    data = f.read()
    data_lines = data.split('\n')

    DATA = []
    flag=0
    for i in range(len(data_lines)-1):
        if data_lines[i].split() != [] and data_lines[i].split()[0] == 'unfix' and (data_lines[i].split()[1] == 'RESCALE' or data_lines[i].split()[1] == 'all_rescale'):
            flag = 1
        if flag == 1 and data_lines[i].split() != []:
            DATA.append(data_lines[i].split())
    print('Finished.')
    return DATA
    
def props(data):
    print('Reading in thermodynamic properties.')
    temp = []
    press = []
    press_cum = []
    press_cum_new = []
    dP = []
    vol = []
    p_cum = 0
    count = 0
    for i in range(len(data)-1):
        if len(data[i]) == 12 and data[i][0] != 'Step':
            try:
                temp.append(float(data[i][4]))
                press.append(float(data[i][6]))
                press_cum.append(float(data[i][7]))
                if p_cum == 0:
                    p_cum = float(data[i][6])
                    press_cum_new.append(p_cum)
                    dP.append(0)
                    count += 1
                else:
                    count += 1
                    p_cum +=  float(data[i][6])
                    press_cum_new.append(p_cum/count)
                    dP.append(float(data[i][6])-(p_cum/count))
                vol.append(float(data[i][8]))
            except:
                print(data[i])
                print('No data in this line.')
    temp, press, press_cum, press_cum_new, dP, vol = np.array(temp), np.array(press), np.array(press_cum), np.array(press_cum_new), np.array(dP), np.array(vol)
    print('Finished.')
    return temp, press, press_cum, press_cum_new, dP, vol
    
def correlation(data):
    acf = []
    times_range = 20

    print('Calculating the autocorrelation function.')
    # iterate over different time origins
    for i in range(len(data)-times_range): 
        C_tmp = 0 
        Cdata = [] # stores correlation for each time origin
        for j in range(times_range+1): # the n in t_0+ndt
            C_tmp += np.dot(data[i], data[i+j])
            Cdata.append(C_tmp)
        if Cdata != []:
            acf.append(Cdata)
    acf = np.array(acf)
    
    mean_acf = np.mean(acf, axis=0)
    norm_acf = mean_acf/mean_acf[0]

    return mean_acf, norm_acf

def main():
    args    = GetArgs()
    fname  = args.filename[0]


    # Definitions
    fig_size_long = (14,5)
    fig_size_sq   = (9,7)

    EXT = 'PDF'
    ext = 'pdf'


    markers = ['o', 'D', 's', 'v', '^', 'd', '*']
    colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']
    
    DATA = read_log(fname)
    T, P, Pc, Pcnew, dP, V = props(DATA)
    ACF = correlation(P)
    
    Vave = np.mean(V)
    Tave = np.mean(T)
    kB = 1.38e-23
    A2m = 1e-10
    ps2s = 1e-12
    bar2Pa = 1e5
    convert = A2m**3*ps2s*bar2Pa**2
    dt = 0.0005
    sampint = 20*5000
    
    eta_scale = (Vave*sampint*dt*convert)/(Tave*kB)
    ACF_int = simps(ACF[0])
    
    ETA = eta_scale*ACF_int
    
    print('The viscosity is:', ETA)

    # Example figure

    fig1 = plt.figure(figsize=fig_size_sq)
    ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])
    ax1.plot(range(len(P)), P, label = 'P(t)')
    ax1.plot(range(len(P)), Pc, label = 'P cum Lammps')
    ax1.plot(range(len(P)), Pcnew, label = 'P cum here')
    #ax1.set_xlabel('no. of sheets')
    #ax1.set_ylabel('$k$ (10$^{-17}$ m$^2$/sPa)')
    ax1.legend()
    fig1.savefig('pressure.{}'.format( ext))
    
    fig2 = plt.figure(figsize=fig_size_sq)
    ax2  = fig2.add_axes([0.15,0.15,0.75,0.75])
    ax2.plot(range(len(P)), dP, label = 'dP')
    #ax1.set_xlabel('no. of sheets')
    #ax1.set_ylabel('$k$ (10$^{-17}$ m$^2$/sPa)')
    ax2.legend()
    fig2.savefig('dP.{}'.format( ext))
    
    fig3 = plt.figure(figsize=fig_size_sq)
    ax3  = fig3.add_axes([0.15,0.15,0.75,0.75])
    ax3.plot(range(len(ACF[0])), ACF[0], label = 'ACF')
    #ax1.set_xlabel('no. of sheets')
    #ax1.set_ylabel('$k$ (10$^{-17}$ m$^2$/sPa)')
    ax3.legend()
    fig3.savefig('ACF_P.{}'.format(ext))


    return

if __name__ == "__main__":
    #sys.path.append('/home/fred/SCRIPTS/Python')
    #from plotting_params import *
    main()