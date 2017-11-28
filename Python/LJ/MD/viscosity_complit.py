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

import sys

np.set_printoptions(precision=16)
# Parse arguments form command line

def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-T", type=str, nargs='+', \
                        help="Temperature",required=True)
    parser.add_argument("-rhof", type=str, nargs='+', \
                        help="Fluid density",required=False, default = ['0.1', '0.15','0.2', '0.25','0.3', '0.35', '0.4', '0.45', '0.5', '0.55', '0.6', '0.65', '0.7', '0.75', '0.8', '0.85', '0.9'])#, '0.95', '0.975', '1.0'])
    parser.add_argument("-r", type=str, nargs=1, \
                        help="Rerun", default='n')
    parser.add_argument("-rc", type=str, nargs='+', \
                        help="Rerun", default=['5.5'])

    args = parser.parse_args()
    return args



def read_log(T, den,rcut):

    filename = "/media/fred/My Passport/LJ_relations_results/CORRECT/rc_{}/T{}/rhof{}/log.rel_rhof{}".format(rcut,T,den,den)

    # read in data
    f = open(filename,'r')
    data = f.read()


    data_lines = data.split('\n')
    print 'Density: ', den

    for i in range(0,len(data_lines)-1):

        log_data = data_lines[i].split()

        if len(log_data) >1:
            if log_data[0] == 'Loop':
                dat = data_lines[i-1].split()
                print 'Pressure: ', dat[5]
                press = float(dat[5])
            continue

    return press

def blockAverage(data):

    DataLen     = len(data) 
    BlockSize   = 100       # max: 4 blocs (otherwise can't calc variance)
  
    NoBlocks    = int(DataLen/BlockSize)               # total number of such blocks in datastream
    print NoBlocks
    Data        = np.zeros(NoBlocks)                  # container for parcelling block 

    # Loop to chop datastream into blocks
    # and take average
    for i in range(1,NoBlocks+1):
        
        istart = (i-1) * BlockSize
        iend   =  istart + BlockSize
        Data[i-1] = np.mean(data[istart:iend])

    meanVal  = np.mean(Data)
    meanErr  = np.sqrt(np.var(Data)/(NoBlocks - 1))


    return meanVal, meanErr

def prop(PROP,T,den,rcut):

    filename = "/media/fred/My Passport/LJ_relations_results/CORRECT/rc_{}/T{}/rhof{}/viscosity.rhof{}".format(rcut,T,den,den)


    df = pd.read_csv(filename, delimiter=' ', skiprows=2)
    #read first two lines
    with open(filename, 'r') as f:
        _, line2 = f.readline(), f.readline()
    cols = line2.lstrip('#').strip().split(' ')
    df.columns = cols
    Prop = df['{}'.format(PROP)]
    Prop = Prop.tolist()
    time = df['TimeStep']
    time = time.tolist()

    Prop_list = []
    count = 0
    tlim = 100000
    for t,p in itertools.izip(time, Prop):
        if t > tlim:
            Prop_list.append(p)
            count += 1
    Prop_val, Prop_err = blockAverage(Prop_list)
    print PROP, ':', Prop_val, '+/-', Prop_err
    return Prop_val, Prop_err

def save_data(T, D, P, SH, SHe, BU, BUe, TC):

    data_list = np.array([D,P,SH,SHe,BU,BUe,TC])
    data_array = np.transpose(data_list)
    np.savetxt('relations_err_T%s.dat'%(T), data_array,fmt="%s",header = 'Density Pressure Shear ShearErr Bulk BulkErr TCond')
    return data_array

def read_data(prop,T):
    JCPdata = np.loadtxt('JCPdata{}T{}.dat'.format(prop,T))
    JCPrho = JCPdata[:,0]
    JCPprop = JCPdata[:,2]
    return JCPrho, JCPprop

	
def main():
    args = GetArgs()
    Temp = args.T
    dens = args.rhof
    rerun = args.r[0]
    rcut = args.rc
    #rc    = 5.5

    fig1 = plt.figure(figsize=(9,7))
    ax1  = fig1.add_axes([0.1,0.15,0.8,0.75])
    fig2 = plt.figure(figsize=(9,7))
    ax2  = fig2.add_axes([0.1,0.15,0.8,0.75])
    fig3 = plt.figure(figsize=(9,7))
    ax3  = fig3.add_axes([0.1,0.15,0.8,0.75])

    name_plot1 = 'comp_rho_v_shear_T'
    name_plot2 = 'comp_rho_v_bulk_T'
    name_plot3 = 'comp_rho_v_ratio_T'
    colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']
    #colours = {'1.35':'#3a42d4', '1.7': '#d43a3a', '2.0':'#60c1dc',\
    #'2.5': '#60dc96','0.8':'#60c1dc',}
    markers = ['D', 's', 'v', '^', 'd', 'h']

    nameT = ''
    namerc = ''
    idx=0
    for T in Temp:
        nameT+=T
        for rc in rcut:
            namerc+=rc

            try:
                if rerun=='y':
                    raise IOError

                final_data = np.loadtxt('relations_err_T{}.dat'.format(T))
                dens_new, press, shear, shear_err, bulk, bulk_err, tcond = final_data[:,0], final_data[:,1], final_data[:,2], final_data[:,3], final_data[:,4], final_data[:,5], final_data[:,6]
                print 'Data already generated for T={}.'.format(T)
                
            except IOError:
                print 'Generating data for T={}....'.format(T)
                dens_new =[]
                press = []
                shear = []
                shear_err = []
                bulk = []
                bulk_err = []
                tcond = []


                for d in dens:
                    try:
                        press.append(read_log(T,d,rc))
                        dens_new.append(float(d))
                        Shear = prop('v_etas_vi',T,d,rc)
                        Bulk = prop('v_etab_pt',T,d,rc)
                        Tcond = prop('v_thermc',T,d,rc)
                        shear.append(Shear[0])
                        shear_err.append(Shear[1])
                        bulk.append(Bulk[0])
                        bulk_err.append(Bulk[1])
                        tcond.append(Tcond[0])
                    except IOError: 
                        print 'No data for rho=', d
                final_data = save_data(T, dens_new,press,shear,shear_err,bulk, bulk_err,tcond)

            ratio = []
            ratio_err = []
            for s,se,b,be in itertools.izip(final_data[:,2],final_data[:,3],final_data[:,4],final_data[:,5]):
                ratio.append(b/s)
                r_err = (b/2)*np.sqrt((se/s)**2+(be/b)**2)
                ratio_err.append(r_err)

            try:
                JCP_rho, JCP_eta = read_data('Eta',T)
                sorti = np.argsort(JCP_rho)
                JCP_rho_eta, JCP_eta = JCP_rho[sorti], JCP_eta[sorti]
                ax1.plot(JCP_rho_eta, JCP_eta, ls='None', marker='x', c=colours[idx*2])

                JCP_rho, JCP_kappa = read_data('Kappa',T)
                sorti = np.argsort(JCP_rho)
                JCP_rho_kappa, JCP_kappa = JCP_rho[sorti], JCP_kappa[sorti]
                ax2.plot(JCP_rho_kappa, JCP_kappa, ls='None', marker='x', c=colours[idx*2])
            except:
                print 'No literature data available for T=', T

            if len(rcut)>1:
                label = '$r_{\mathrm{cut}}$ = %s'%(rc)
                colidx = idx
            else:
                label = 'T$^*$ = {}'.format(T)
                colidx=idx*2

            # PLOTTING
            ax1.locator_params(axis='y', nbins=5)
            ax2.locator_params(axis='y', nbins=5)
            ax3.locator_params(axis='y', nbins=5)

            ax1.errorbar(dens_new,shear,yerr=shear_err, ls='None', marker=markers[idx], c=colours[colidx],label=label)
            ax1.set_xlabel('$\\rho^*$')
            ax1.set_ylabel('$\eta^*$')
            ax1.set_xlim(0,0.9)
            ax1.set_ylim(0,4)
            ax1.legend()

            ax2.errorbar(dens_new,bulk,yerr=bulk_err, ls='None', marker=markers[idx], c=colours[colidx],label=label)
            ax2.set_xlabel('$\\rho^*$')
            ax2.set_ylabel('$\kappa^*$')
            ax2.set_xlim(0,0.9)
            ax2.set_ylim(0,1.5)
            ax2.legend()


            ax3.errorbar(dens_new,ratio,yerr=ratio_err, ls='None', marker=markers[idx], c=colours[colidx],label=label)
            ax3.set_xlabel('$\\rho^*$')
            ax3.set_ylabel('$\kappa^*/\eta^*$')
            ax3.set_xlim(0,0.9)
            ax3.set_ylim(0,4)
            ax3.legend()

            
            if idx==0 and len(Temp)<4:
                idx+=2
            else:
                idx+=1


    #ax3.legend(loc='upper center', ncol=3, fontsize=18)
    fig1.savefig('PLOTS/{}{}_rc{}.pdf'.format(name_plot1,nameT,namerc))
    fig2.savefig('PLOTS/{}{}_rc{}.pdf'.format(name_plot2,nameT,namerc))
    fig3.savefig('PLOTS/{}{}_rc{}.pdf'.format(name_plot3,nameT,namerc))

    #plt.show()

    return


if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()
