#!/usr/bin/env python

"""Reads in log files and extracts pressure, shear and bulk viscosity data and compares to data from Meier et al 2004 JCP"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
from matplotlib import rc
import itertools

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



def read_log(press, shear, bulk, tcond,T, den):

    if T=='1.35':
        filename = "/media/fred/My Passport/LJ_relations_results/T%s/rhof%s/log.rel_rhof%s" % (T,den,den)
    elif T=='0.8':
        filename = "/media/fred/My Passport/LJ_relations_results/T%s/log.rel_rhof%s" % (T,den)
    else:
        filename = "/media/fred/My Passport/LJ_relations_results/CORRECT/T%s/rhof%s/log.rel_rhof%s" % (T,den,den)

    # read in data
    f = open(filename,'r')
    data = f.read()


    data_lines = data.split('\n')
    print 'Density: ', den

    for i in range(0,len(data_lines)-1):

        log_data = data_lines[i].split()

        if len(log_data) >1:
            #print log_data	
            if log_data[1] == 'shear' and log_data[2] == 'viscosity' and log_data[3] == ':':
                print 'Shear viscosity: ', log_data[4]
                shear.append(float(log_data[4]))
            if log_data[1] == 'bulk' and log_data[2] == 'pt':
                print 'Bulk viscosity: ', log_data[5]
                bulk.append(float(log_data[5]))
            if log_data[1] == 'thermal' and log_data[2] == 'conductivity':
                print 'Thermal conductivity: ', log_data[4]
                tcond.append(float(log_data[4]))
            if log_data[0] == 'Loop':
                dat = data_lines[i-1].split()
                print 'Pressure: ', dat[5]
                press.append(float(dat[5]))
            continue

    return press, shear, bulk, tcond



def save_data(T, D, P, SH, BU, TC):

    data_list = np.array([D,P,SH,BU,TC])
    data_array = np.transpose(data_list)
    print data_array
    print type(data_array[0][0])
    np.savetxt('relations_T%s.dat'%(T), data_array,fmt="%s",header = 'Density Pressure Shear Bulk TCond')
    return data_array

def read_data(name):
	JCPdata = np.loadtxt('%s'%name)
	JCPrho = JCPdata[:,0]
	JCPeta = JCPdata[:,2]
	return JCPrho, JCPeta



def Woodcock(Rho, Eta_0, T_in):
	#rho_tmp = np.array(map(float,Rho))
	T_tmp = float(T_in)
	term1 = np.sqrt(2)*(1-T_tmp**(-4)-(T_tmp/8))*Rho
	term2 = 3.025*Rho**4/(T_tmp**(1/3))
	#print (2.535+Eta_0)/Eta_0, term1, term2, Rho**4, (1/T_tmp)**(1/3)
	return Eta_0*(1+term1+term2)

	
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

    name_plot1 = 'comp_rho_v_shear_T'
    name_plot2 = 'comp_rho_v_bulk_T'
    name_plot3 = 'comp_rho_v_ratio_T'
    name_plot4 = 'comp_rho_v_tcond_T'
    colours = {'1.35':'#3a42d4', '1.7': '#d43a3a', '2.0':'#60c1dc',\
    '2.5': '#60dc96','0.8':'#60c1dc',}
    markers = ['D', 's', 'v', '^', 'd', '*']
    idx = 0
    for T in Temp:
        name_plot1+=T
        name_plot2+=T
        name_plot3+=T
        name_plot4+=T


        if T=='1.35' or T=='2.5' or T=='0.8':
            # JCP data for comparison
            # ETA
            name = 'JCPdataEtaT{}.dat'.format(T)
            JCP_rho, JCP_eta = read_data(name)
            #sort arrays by density
            sorti = np.argsort(JCP_rho)
            JCP_rho_eta, JCP_eta = JCP_rho[sorti], JCP_eta[sorti]
            eta_0 = JCP_eta[0]

            # KAPPA
            name = 'JCPdataKappaT{}.dat'.format(T)
            JCP_rho, JCP_lda = read_data(name)
            #sort arrays by density
            sorti = np.argsort(JCP_rho)
            JCP_rho_kappa, JCP_lda = JCP_rho[sorti], JCP_lda[sorti]

        #print Woodcock(0.8442, 0.076, 0.722)

        # Our simulation data

        try:
            final_data = np.loadtxt('relations_T{}.dat'.format(T))
            print 'Data already generated.'
            
        except IOError:
            print 'Generating data for T={}....'.format(T)
            dens_new =[]
            press = []
            shear = []
            bulk = []
            tcond = []


            for d in dens:
                try:
                    press, shear, bulk, tcond = read_log(press,shear,bulk,tcond,T,d)
                    dens_new.append(float(d))
                except IOError: 
                    print 'No data for rho=', d
            final_data = save_data(T, dens_new,press,shear,bulk, tcond)

        ratio = []
        for s,b in itertools.izip(final_data[:,2],final_data[:,3]):
            ratio.append(b/s)

        print T

        # PLOTTING

        ax1.plot(final_data[:,0], final_data[:,2], marker=markers[idx], linestyle = 'none', c=colours[T], label='T={}'.format(T))
        #plt.plot(JCP_rho, Woodcock(JCP_rho, eta_0, T), lw = 2.0,linestyle = 'dashed', c='r', label = "Woodcock")
        if T=='1.35' or T=='2.5' or T=='0.8':
            ax1.plot(JCP_rho_eta, JCP_eta, lw = 2.0, marker = 'x', linestyle = 'none', c=colours[T])#, label='T={}$^a$'.format(T))
        ax1.set_ylabel('$\eta^*$')
        ax1.set_xlabel('$\\rho^*$')
        ax1.set_xlim(0.0,0.9)
        ax1.set_ylim(0.0,4.0)
        ax1.legend(loc='upper left')
        #plt.ylim(0,6)

        ax2.plot(final_data[:,0], final_data[:,3], marker=markers[idx], linestyle = 'none', c=colours[T], label='T={}'.format(T))
        if T=='0.8' or T=='1.35' or T=='2.5':
            ax2.plot(JCP_rho_kappa, JCP_lda, lw = 2.0, marker = 'x', linestyle = 'none', c=colours[T])#, label='T={}$^b$'.format(T))
        ax2.set_ylabel('$\kappa^*$')
        ax2.set_xlabel('$\\rho^*$')
        ax2.set_xlim(0.0,0.9)
        ax2.legend(loc='upper left')
        #plt.ylim(0,6)

        ax3.plot(final_data[:,0], ratio, marker=markers[idx], linestyle = 'none', c=colours[T], label='T={}'.format(T))
        #plt.plot(JCP_rho, Woodcock(JCP_rho, eta_0, T), lw = 2.0,linestyle = 'dashed', c='r', label = "Woodcock")
        #if T=='1.35' or T=='2.5':
        #    ax1.plot(JCP_rho, JCP_eta, lw = 2.0, marker = 'x', linestyle = 'none', c=colours[T])#, label='T={}$^a$'.format(T))
        ax3.set_ylabel('$\kappa^*/\eta^*$')
        ax3.set_xlabel('$\\rho^*$')
        ax3.set_xlim(0.0,0.9)
        ax3.set_ylim(0.0,4.0)
        ax3.legend(loc='upper left')
        #plt.ylim(0,6)

        ax4.plot(final_data[:,0], final_data[:,4], marker=markers[idx], linestyle = 'none', c=colours[T], label='T={}'.format(T))
        ax4.set_ylabel('$\lambda^*$')
        ax4.set_xlabel('$\\rho^*$')
        ax4.set_xlim(0.0,0.9)
        #ax4.set_ylim(0.0,4.0)
        ax4.legend(loc='upper left')
        idx+=1


    fig1.savefig('PLOTS/{}.pdf'.format(name_plot1))
    fig2.savefig('PLOTS/{}.pdf'.format(name_plot2))
    fig3.savefig('PLOTS/{}.pdf'.format(name_plot3))
    fig4.savefig('PLOTS/{}.pdf'.format(name_plot4))


    return


if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()

