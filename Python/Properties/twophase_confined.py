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
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model', nargs='+', required=False, default='liquid', action='store',
                        help='Model')
    parser.add_argument('-s', '--spacing', nargs='+', required=False, default='liquid', action='store',
                        help='Channel height')
    parser.add_argument('-e', '--ext', nargs=1, type=str, required=False, default=['pdf'], action='store')
    parser.add_argument('-d', '--den', nargs=1, type=str, required=False, default=[6], action='store')
    parser.add_argument('-r', '--density', nargs=1, required=False, default=['RHO1'], action='store',
                        help='Avg channel density')

    args = parser.parse_args()
    return args


def coords(fluid, f, r, i):
    filename = ('/media/fred/8TB/Bulk_properties/{}/{}/log.{}'
                .format(fluid, r, f))
    f = open(filename, 'r')

    data = f.read()
    data_lines = data.split('\n')

    coord_line = re.sub('[\(\)]', '', data_lines[i]).split('=')[1].split('to')
    xlo, ylo, zlo = coord_line[0].strip(' ').split(' ')
    xhi, yhi, zhi = coord_line[1].strip(' ').split(' ')

    return float(xlo), float(xhi), float(ylo), float(yhi), float(zlo), float(zhi)


def read_log(fluid, f, r):
    '''Code to read in log file'''
    filename = ('/media/fred/8TB/Bulk_properties/{}/{}/log.{}'
                .format(fluid, r, f))
    f = open(filename, 'r')

    data = f.read()
    data_lines = data.split('\n')

    DATA = []
    flag = 0
    for i in range(len(data_lines)-1):
        if data_lines[i].split() != [] and data_lines[i].split()[0] == 'unfix' and (data_lines[i].split()[1] == 'RESCALE' or data_lines[i].split()[1] == 'all_rescale'):
            flag = 1
        if flag == 1 and data_lines[i].split() != []:
            DATA.append(data_lines[i].split())
    return DATA


def num_water(fluid, f, r):
    filename = ('/media/fred/8TB/Bulk_properties/{}/{}/log.{}'
                .format(fluid, r, f))
    f = open(filename, 'r')
    data = f.read()
    data_lines = data.split('\n')

    for i in range(len(data_lines)-1):
        if data_lines[i].split() != [] and data_lines[i].split()[0] == 'group' and data_lines[i].split()[1] == 'water':
            numwater = int(data_lines[i+1].split()[0])
            break
    return numwater


def avg_density(fluid, f, r):
    filename = ('/media/fred/8TB/Bulk_properties/{}/{}/dens.{}'
                .format(fluid, r, f))
    f = open(filename, 'r')
    _, __, data = (f.readline(), f.readline(), f.readline())

    RHOeff = float(data.split()[1])
    return RHOeff


def main():

    args = GetArgs()

    model = args.model[0]
    DZ = args.spacing
    DIR = args.density[0]
    ext = args.ext[0]
    d = args.den[0]
    T = 298
    fluid = 'Water_Graphene'

    Nconf = 2000
    idx = 45

    # Definitions
    fig_size_long = (14, 5)
    fig_size_sq = (9, 7)
    fig_size_sq2 = (7, 9)

    markers = ['o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P', 'X', '8', 'p',\
               'o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P', 'X', '8', 'p']
    colours = ['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026',\
    '#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']

    for dz in DZ:

        filename = '{}_T{}_z{}_eps1.0_{}'.format(model, T, dz, d)
        print '\n\nReading in', filename
        print '\n#-----------------------Finding box dimensions------------------#\n\n'
        xlo, xhi, ylo, yhi, zlo, zhi = coords(fluid, filename, DIR, idx)
        offset = -xlo
        Lx = xhi-xlo
        Ly = yhi-ylo
        Lz = float(dz)-3.19
        avg_dens = avg_density(fluid, filename, DIR)
        delta_r = 0.04  # 0.02/avg_dens  #boxlength/1400
        Natom = num_water(fluid, filename, DIR)
        Ndiv = 6  # int(Natom/Natom)
        print 'Box dimensions: ( {}, {}, {}), ({}, {}, {})'.format(xlo, ylo, zlo, xhi, yhi, zhi)
        print 'Box offset:', offset
        print 'Box lengths:', (Lx, Ly, Lz)
        print 'Number of atoms:', Natom
        print 'Average density:', avg_dens

        print '\n#-----------------------Amending files------------------#\n\n'
        # replace box dimensions in files
        inp_cmd = ('sed -e "2s/CONF/{}/g" -e "4s/ATOM/{}/g" '
                   '-e "6s/LENX/{}/g" -e "6s/LENY/{}/g" -e "6s/LENZ/{}/g" '
                   '-e "8s/DIV/{}/g" -e "10s/DR/{}/g" data_tmp.inp'
                   ' > data.inp').format(Nconf, Natom, Lx, Ly, Lz, Ndiv, delta_r)
        subprocess.Popen(inp_cmd, shell=True)
        print 'data.inp: Done'

        phase_cmd = ('sed -e "24s/NAME/{}/g" -e "54s/NAME/{}/g" -e "91s/OFF/{}/g" ../source/phase_sph_2d_tmp.f90 > ../source/phase_sph_2d.f90'
                     .format(filename, filename, offset))
        subprocess.Popen(phase_cmd, shell=True)
        print 'phase_sph_2d.f90: Done'

        print '\n#-----------------------Compiling code------------------#\n\n'
        subprocess.call('./make_2d.sh', shell=True)
        print 'Done'

        print '\n#-----------------------Executing code------------------#\n\n'
        subprocess.call('rm freq_{}.dat'.format(filename), shell=True)
        subprocess.call('./composition_2d.exe', shell=True)
        print 'Done'

        print '\n#-----------------------Plotting------------------#\n\n'
        # ------------------- Initialise figures -----------------------

        fig1 = plt.figure(figsize=fig_size_sq)
        ax1 = fig1.add_axes([0.15, 0.15, 0.75, 0.75])

        # Let's read in the file
        HIST = np.loadtxt('freq_{}.dat'.format(filename))
        # print HIST
        Xax = HIST[:, 0]
        Yax = HIST[:, 1]
        RHOax = avg_dens*Xax
        width = RHOax[1]-RHOax[0]
        # RHOax = RHOax+width/2

        ylim = np.max(Yax)+0.1*np.max(Yax)
        Yave = np.linspace(0, ylim, 10)
        lims = (0,1.2)

        ax1.bar(RHOax, Yax, width, color=colours[0], alpha=0.9)
        ax1.plot([avg_dens]*len(Yave), Yave, c=colours[6], linestyle='dotted')
        ax1.set_xlabel('$\\rho$ (g/cm$^3$)')
        ax1.set_ylabel('Frequency')
        ax1.set_xlim(lims)
        #ax1.set_ylim(0, ylim)
        fig1.savefig('PLOTS_C/{}/hist_{}_T{}_z{}_eps1.0_{}.{}'.format(ext, model, T, dz, d, ext))
        print 'Plot saved'
        # plt.show()

        print '\n#-----------------------Plotting finished------------------#\n\n'
    print '\n#-----------------------Done------------------\n#'
    return


if __name__ == "__main__":
    sys.path.insert(0, '/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()
