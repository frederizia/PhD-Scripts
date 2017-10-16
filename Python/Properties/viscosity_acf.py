#!/usr/bin/env python

"""Reads in log files and extracts pressure, shear and bulk viscosity data and compares to data from Meier et al 2004 JCP"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
from matplotlib import rc
import itertools
import pandas as pd
from scipy.optimize import curve_fit

import sys

np.set_printoptions(precision=16)
# Parse arguments form command line

def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-T", type=str, nargs='+', \
                        help="Temperature",required=True)
    parser.add_argument("-rhof", type=str, nargs='+', \
                        help="Fluid density",required=False, default = ['0.1', '0.15','0.2', '0.25','0.3', '0.35', '0.4', '0.45', '0.5', '0.55', '0.6', '0.65', '0.7', '0.75', '0.8', '0.85', '0.9'])#, '0.95', '0.975', '1.0'])

    args = parser.parse_args()
    return args

def fit(x, y, xmin, xmax):
    x = x[:2000]
    y = y[:2000]
    params, cov = curve_fit(f1, np.array(x), np.array(y), p0=(y[0], -1))
    xdat = np.linspace(xmin, xmax, 100)
    fit = []
    print params
    for x in xdat:
        fit.append(params[0]*np.exp(params[1]*x))
    return xdat, fit

def f1(x, A, B):
    return A*np.exp(B*x)

def acf(T,Type, den):
    filename = '/media/fred/My Passport/LJ_relations_results/CORRECT/rc_2.5/T{}/rhof{}/acf{}.rhof{}'.format(T,den,Type,den)
    if Type == 'bv':
        idx = 3
    elif Type == 'sv':
        idx = 3

    count = 0
    with open(filename) as infile:
        for line in infile:
            items = line.split()
            if len(items) == 2:
                ACF = []
                delta = []
                count += 1
                

            elif line[0] != '#' and len(items) != 2 and len(items) != 3:
                Tdelta = float(items[1])
                Acf    = float(items[3])
                delta.append(Tdelta*0.001)
                ACF.append(Acf)

    print count, 'timesteps collected.'
            
    return delta, ACF


def main():
    args = GetArgs()
    Temp = args.T
    dens = args.rhof

    fig1 = plt.figure(figsize=(9,7))
    ax1  = fig1.add_axes([0.1,0.15,0.8,0.75])
    fig2 = plt.figure(figsize=(9,7))
    ax2  = fig2.add_axes([0.1,0.15,0.8,0.75])
    fig3 = plt.figure(figsize=(9,7))
    ax3  = fig3.add_axes([0.1,0.15,0.8,0.75])
    fig4 = plt.figure(figsize=(9,7))
    ax4  = fig4.add_axes([0.1,0.15,0.8,0.75])


    name_plot1 = 'ACF_SV_T'
    name_plot2 = 'ACF_BV_T'
    colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']
    #colours = {'1.35':'#3a42d4', '1.7': '#d43a3a', '2.0':'#60c1dc',\
    #'2.5': '#60dc96','0.8':'#60c1dc',}
    markers = ['D', 's', 'v', '^', 'd', '*']

    nameT = ''
    name_rho = ''
    
    for T in Temp:
        nameT+=T
        idx=0

        for d in dens:
            name_rho +=d

            print 'Reading in data for rho =', d
            Tdelta_sv, ACF_sv = acf(T,'sv',d)
            Td_sv_fit, ACF_sv_fit = fit(Tdelta_sv, ACF_sv,0,1.0)
            Tdelta_bv, ACF_bv = acf(T,'bv',d)
            Td_bv_fit, ACF_bv_fit = fit(Tdelta_bv, ACF_bv,0,1.0)

            # PLOTTING
            #ax1.locator_params(axis='y', nbins=5)
            #ax2.locator_params(axis='y', nbins=5)
            #ax3.locator_params(axis='y', nbins=5)

            ax1.semilogy(Tdelta_sv,ACF_sv,c=colours[idx],label='$\\rho = {}$'.format(d))
            ax1.semilogy(Td_sv_fit,ACF_sv_fit,c=colours[idx], ls='dashed')
            ax1.set_xlabel('$t^*$')
            ax1.set_ylabel('ACF SV')
            ax1.set_xlim(0,1)
            ax1.legend()
            #ax1.set_xlim(0,0.9)
            #ax1.set_ylim(0,4)


            ax2.semilogy(Tdelta_bv,ACF_bv,c=colours[idx],label='$\\rho = {}$'.format(d))
            ax2.semilogy(Td_bv_fit,ACF_bv_fit,c=colours[idx], ls='dashed')
            ax2.set_xlabel('$t^*$')
            ax2.set_ylabel('ACF BV')
            ax2.set_xlim(0,1)
            ax2.legend()
            #ax1.set_xlim(0,0.9)
            #ax1.set_ylim(0,4)

            

            ax3.plot(Tdelta_sv,ACF_sv,c=colours[idx],label='$\\rho = {}$'.format(d))
            ax3.plot(Td_sv_fit,ACF_sv_fit,c=colours[idx], ls='dashed')
            ax3.set_xlabel('$t^*$')
            ax3.set_ylabel('ACF SV')
            ax3.set_xlim(0,1)
            ax3.legend()
            #ax1.set_ylim(0,4)


            ax4.plot(Tdelta_bv,ACF_bv,c=colours[idx],label='$\\rho = {}$'.format(d))
            ax4.plot(Td_bv_fit,ACF_bv_fit,c=colours[idx], ls='dashed')
            ax4.set_xlabel('$t^*$')
            ax4.set_ylabel('ACF BV')
            ax4.set_xlim(0,1)
            ax4.legend()
            #ax1.set_ylim(0,4)

            idx+=1

        plt.show()

        fig1.savefig('PLOTS/{}{}_rhof{}_log.pdf'.format(name_plot1,nameT,name_rho))
        fig2.savefig('PLOTS/{}{}_rhof{}_log.pdf'.format(name_plot2,nameT,name_rho))
        fig3.savefig('PLOTS/{}{}_rhof{}.pdf'.format(name_plot1,nameT,name_rho))
        fig4.savefig('PLOTS/{}{}_rhof{}.pdf'.format(name_plot2,nameT,name_rho))


        


    #plt.show()

    return


if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()
