#! /usr/bin/env python
# Find phase distribution in two-phase system

import argparse
import matplotlib.pyplot as plt
import subprocess
import matplotlib
import sys
import numpy as np
import re

def GetArgs():
    ## Parse command line arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model', nargs='+', required=False, default='liquid', action='store',
                       help='Model')
    parser.add_argument('-p', '--pressure', nargs='+', required=False, default='liquid', action='store',
                       help='Target pressure of bulk system')
    parser.add_argument('-e', '--ext', nargs=1, type=str, required=False, default=['pdf'], action='store')
    
    args = parser.parse_args()
    return args

def coords(fluid,f):
    filename = '/media/fred/8TB/Bulk_properties/{}/log.{}'.format(fluid,f)
    f = open(filename,'r')

    data = f.read()
    data_lines = data.split('\n')

    coord_line = re.sub('[\(\)]', '', data_lines[25]).split('=')[1].split('to')
    xlo, ylo, zlo = coord_line[0].strip(' ').split(' ')
    xhi, yhi, zhi = coord_line[1].strip(' ').split(' ')

    return float(xlo), float(xhi), float(ylo), float(yhi), float(zlo), float(zhi)

def read_log(fluid,f):
    '''Code to read in log file'''
    filename = '/media/fred/8TB/Bulk_properties/{}/log.{}'.format(fluid,f)
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
    return DATA

def avg_density(fluid, f):
    data = read_log(fluid,f)

    for i in range(len(data)-1):
        if data[i+1][0] == 'Loop':
            dens = float(data[i][2])
    return dens

def main():


    args = GetArgs()

    model = args.model[0]
    press = args.pressure
    ext = args.ext[0]
    d = 6
    T = 300
    fluid = 'CO2'

    Nconf = 1000
    Ndiv = 50

    if fluid == 'CO2':
        lims = (0, 1.4)
    elif fluid == 'Water':
        lims = (0,3.5)



    # Definitions
    fig_size_long = (14,5)
    fig_size_sq   = (9,7)
    fig_size_sq2 = (7,9)

    markers = ['o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P','X', '8', 'p',\
    'o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P','X', '8', 'p']
    colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026',\
    '#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']


    

    for p in press:

        filename = '{}_T{}_P{}_{}'.format(model,T,p,d)
        print '\n\nReading in', filename
        print '\n#-----------------------Finding box dimensions------------------#\n\n'
        xlo, xhi, ylo, yhi, zlo, zhi = coords(fluid,filename)
        offset = -xlo
        boxlength = xhi-xlo
        avg_dens = avg_density(fluid,filename)
        delta_r = 0.04/avg_dens #boxlength/1400
        print 'Box dimensions: ( {}, {}, {}), ({}, {}, {})'.format(xlo, ylo, zlo, xhi, yhi, zhi)
        print 'Box offset:', offset
        print 'Box length:', boxlength
        print 'Average density:', avg_dens



        print '\n#-----------------------Amending files------------------#\n\n'
        # replace box dimensions in files
        inp_cmd = 'sed -e "2s/CONF/{}/g" -e "6s/LEN/{}/g" -e "8s/DIV/{}/g" -e "10s/DR/{}/g" data_tmp.inp > data.inp'.format(Nconf, boxlength,Ndiv,delta_r)
        subprocess.Popen(inp_cmd, shell=True)
        print 'data.inp: Done'

        phase_cmd = 'sed -e "24s/NAME/{}/g" -e "55s/NAME/{}/g" -e "92s/OFF/{}/g" ../source/phase_sph_tmp.f90 > ../source/phase_sph.f90'.format(filename, filename, offset)
        subprocess.Popen(phase_cmd, shell=True)
        print 'phase_sph.f90: Done'

        print '\n#-----------------------Compiling code------------------#\n\n'
        subprocess.call('./make.sh', shell=True)
        print 'Done'

        print '\n#-----------------------Executing code------------------#\n\n'
        subprocess.call('rm freq_{}.dat'.format(filename), shell=True)
        subprocess.call('./composition.exe', shell=True)
        print 'Done'

        print '\n#-----------------------Plotting------------------#\n\n'
        # ------------------- Initialise figures -----------------------

        fig1 = plt.figure(figsize=fig_size_sq)
        ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])


        # Let's read in the file
        HIST = np.loadtxt('freq_{}.dat'.format(filename))
        #print HIST
        Xax  = HIST[:,0]
        Yax  = HIST[:,1]
        RHOax = avg_dens*Xax
        
        width = RHOax[1]-RHOax[0]
        #RHOax = RHOax+width/2


        if fluid == 'CO2' and p == '1':
            lims = (0,0025)

        ax1.bar(RHOax, Yax, width,color=colours[0], alpha=0.9)
        ax1.plot([avg_dens]*len(Yax), Yax, c=colours[6], linestyle='dotted')
        ax1.set_xlabel('$\\rho$ (g/cm$^3$)')
        ax1.set_ylabel('Frequency')
        ax1.set_xlim(lims)
        ax1.set_ylim(0,0.35)
        fig1.savefig('PLOTS/{}/hist_{}_T{}_P{}_{}.{}'.format(ext,model, T, p, d, ext))
        print 'Plot saved'
        #plt.show()

        print '\n#-----------------------Plotting finished------------------\n#'
    print '\n#-----------------------Done------------------\n#'
    return

if __name__ == "__main__":
    sys.path.insert(0,'/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()