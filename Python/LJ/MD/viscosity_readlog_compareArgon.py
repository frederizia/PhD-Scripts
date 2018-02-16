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
    parser.add_argument("-type", type=str, nargs=1, \
                        help="Fluid density",required=False, default = ['real'])#, '0.95', '0.975', '1.0'])
    parser.add_argument("-r", type=str, nargs=1, \
                        help="Rerun", default='n')

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


	
def main():
    args = GetArgs()
    Temp = args.T
    dens = args.rhof
    units = args.type[0]
    rerun = args.r[0]
    rc    = 5.5


    fig3 = plt.figure(figsize=(9,7))
    ax3  = fig3.add_axes([0.1,0.15,0.8,0.75])
    ax3_inset=fig3.add_axes([0.68,0.36, 0.2,0.3])
    ax3_inset.set_ylim(0,35)
    



    name_plot3 = 'argon_rho_v_ratio_T'
    LIT_colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']
    colours = {'1.35':'#313695', '1.7': '#fdae61', '2.0':'#60c1dc',\
    '2.5': '#d73027','0.8':'#60c1dc',}
    #colours = {'1.35':'#3a42d4', '1.7': '#d43a3a', '2.0':'#60c1dc',\
    #'2.5': '#60dc96','0.8':'#60c1dc',}
    #colours=['#313695', '#4575b4', '#74add1',\
    #'#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']
    markers = ['D', 's', 'v', '^', 'd', 'h']

    # Additional data
    # Madigosky
    Madi_rho = [0.508, 0.694, 0.812, 0.896, 0.958, 1.008]
    Madi_ratio = [0.2465753425, 0.568627451, 0.55, 0.6567164179, 0.6575342466, 0.7]
    Madi_T = 234.55

    # Naugle 1966 2
    Naugle2_138_rho = [1.062, 1.068]
    Naugle2_138_ratio = [2.43, 2.23]
    Naugle2_138_T = 138

    Naugle2_110_rho = [1.28, 1.287, 1.294]
    Naugle2_110_ratio = [1.21, 0.99, 1.07]
    Naugle2_110_T = 110

    Naugle2_89_rho = [1.402, 1.408, 1.41, 1.413]
    Naugle2_89_ratio = [0.86, 0.85, 0.94, 0.89]
    Naugle2_89_T = 89

    Cowan_150_rho = [0.681,0.69,0.726,0.767,0.822,0.852]
    Cowan_150_ratio = [32,24,17,13,9.2,7.4]
    Cowan_150_T = 150

    Cowan_140_rho = [0.942,0.954,0.968,0.981]
    Cowan_140_ratio = [5.5,4.7,4.3,4]
    Cowan_140_T = 140

    LIT_rho = [Madi_rho, Naugle2_89_rho, Naugle2_110_rho, Naugle2_138_rho, Cowan_140_rho, Cowan_150_rho]
    LIT_ratio = [Madi_ratio, Naugle2_89_ratio, Naugle2_110_ratio, Naugle2_138_ratio, Cowan_140_ratio, Cowan_150_ratio]
    LIT_T = [Madi_T, Naugle2_89_T, Naugle2_110_T, Naugle2_138_T, Cowan_140_T, Cowan_150_T]
    LIT_ref = ['a','b','b','b','c','c']

    sigma = 3.4e-8 # cm
    NA = 6.022e23 # mol-1
    M = 39.95 # g/mol
    kB = 1.38e-23
    eps = 1.65e-21 #J
    idx = 0
    T_c = 1.312
    rho_c = 0.316
    T_c_Ar = 150.65
    rho_c_Ar = 0.54
    for T in Temp:
        if units == 'reduced':
            T_new = round(float(T)/T_c,2)
        else:
            T_new = round(float(T)*(eps/kB),1)

        name_plot3+=T


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
        
        dens_argon = []
        ratio = []
        ratio_err = []
        for s,se,b,be in itertools.izip(final_data[:,2],final_data[:,3],final_data[:,4],final_data[:,5]):
            ratio.append(b/s)
            r_err = (b/2)*np.sqrt((se/s)**2+(be/b)**2)
            ratio_err.append(r_err)
        for d in final_data[:,0]:
            if units=='reduced':
                d_new = d/rho_c
                Tname = 'T$^*$'
                Tunit = ''
            else:
                d_new = d*(M/(NA*sigma**3))
                Tname = 'T'
                Tunit = 'K'
            dens_argon.append(d_new)


        print T

        # PLOTTING

        ax3.errorbar(dens_argon, ratio, yerr=ratio_err,marker=markers[idx], linestyle = 'none', c=colours[T], label='{}={}{}'.format(Tname,T_new,Tunit))
        if T=='1.35':
            ax3_inset.errorbar(dens_argon, ratio, yerr=ratio_err,marker=markers[idx], markersize=6,linestyle = 'none', c=colours[T], label='{}={}{}'.format(Tname,T_new,Tunit))
        #plt.plot(JCP_rho, Woodcock(JCP_rho, eta_0, T), lw = 2.0,linestyle = 'dashed', c='r', label = "Woodcock")
        #if T=='1.35' or T=='2.5':
        #    ax1.plot(JCP_rho, JCP_eta, lw = 2.0, marker = 'x', linestyle = 'none', c=colours[T])#, label='T={}$^a$'.format(T))
        ax3.set_ylabel('$\kappa/\eta$')
        #ax3.set_xlim(0.0,0.9)
        ax3.set_ylim(0.0,10)

        
        
        
        #plt.ylim(0,6)

        if idx==0 and len(Temp)<4:
            idx+=2
        else:
            idx+=1

    idx=0
    for i in range(len(LIT_T)):

        if units == 'reduced':
            LIT_rho[i] = np.array(LIT_rho[i])/rho_c_Ar
            LIT_T[i] = np.array(LIT_T[i])/T_c_Ar
            ax3.plot(LIT_rho[i], LIT_ratio[i], marker=markers[i], fillstyle='none',linestyle = 'none', c=LIT_colours[::-1][idx], label='T$^*$={}$^{}$'.format(round(LIT_T[i],2), LIT_ref[i]))
            ax3.set_xlabel('$\\rho^*$')
        else:
            ax3.plot(LIT_rho[i], LIT_ratio[i], marker=markers[i], fillstyle='none',linestyle = 'none', c=LIT_colours[::-1][idx], label='T={}K$^{}$'.format(LIT_T[i], LIT_ref[i]))
            ax3.set_xlabel('$\\rho$ (g/cm$^3$)')
        if LIT_T[i]==150 or LIT_T[i]==140:
            ax3_inset.plot(LIT_rho[i], LIT_ratio[i], marker=markers[i], fillstyle='none',linestyle = 'none', c=LIT_colours[::-1][idx], label='T$^*$={}$^{}$'.format(round(LIT_T[i],2), LIT_ref[i]))
        elif LIT_T[i]==150/T_c_Ar or LIT_T[i]==140/T_c_Ar:
            ax3_inset.plot(LIT_rho[i], LIT_ratio[i], marker=markers[i], fillstyle='none',linestyle = 'none', c=LIT_colours[::-1][idx], label='T={}K$^{}$'.format(LIT_T[i], LIT_ref[i]))
        if idx==0:
            idx+=3
        else:
            idx+=1
    #ax3.legend(loc='upper right')
    ax3.legend(loc='upper center', ncol=3, fontsize=18)
    fig3.savefig('PLOTS/{}_{}.eps'.format(name_plot3,units),format='eps', dpi=1000)


    return


if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()


