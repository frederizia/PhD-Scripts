#! /usr/bin/env python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from scipy.spatial import Voronoi, voronoi_plot_2d
import re


def GetArgs():
    ## Parse command line arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', '--width', nargs='+', required=False, default='liquid', action='store',
                       help='Width of channel')
    parser.add_argument('-rho', '--rho', required=False, default='RHO1', action='store',
                       help='Average channel density')
    parser.add_argument('-s', '--start', nargs='+',required=False, default='None', action='store',
                       help='Starting configurations')
    parser.add_argument('-r', '--density', nargs='+',required=False, default=None, action='store',
                       help='Density for constant channel')

    args = parser.parse_args()
    return args


def read_xy(DIR, z_name):
    # Loop over different times

    filename = 'Water_Graphene/{}/xyz.{}.xyz'.format(DIR,z_name)

    steps = 0
    positions = []
    positions_tmp = []
    with open(filename) as infile:
        for line in infile:
            l = line.split()
            if len(l) == 1:
                if steps != 0 and positions_tmp != []:
                    positions_tmp = np.array(positions_tmp)
                    positions.append(positions_tmp)
                steps += 1
                positions_tmp = []
            elif len(l) == 4 and l[0]=='O':
                x_val = float(l[1])
                y_val = float(l[2])

                positions_tmp.append([x_val,y_val])

    return np.array(positions)


def coords(DIR, z_name):
    i = 45
    filename = 'Water_Graphene/{}/log.{}'.format(DIR,z_name)
    f = open(filename, 'r')

    data = f.read()
    data_lines = data.split('\n')

    coord_line = re.sub('[\(\)]', '', data_lines[i]).split('=')[1].split('to')
    xlo, ylo, zlo = coord_line[0].strip(' ').split(' ')
    xhi, yhi, zhi = coord_line[1].strip(' ').split(' ')

    return float(xlo), float(xhi), float(ylo), float(yhi), float(zlo), float(zhi)


def avg_density(DIR, z_name):

    # open density file
    file = open('Water_Graphene/{}/dens.{}'.format(DIR, z_name), 'r')
    _, __, data = (file.readline(), file.readline(),
                   file.readline())

    RHOeff = float(data.split()[1])
    print '\n\n#---------------------Average density = {} --------------------#'.format(RHOeff)
    
    return 

def plotting(DIR, z_name, opos):
    colours = ['#313695', '#4575b4', '#74add1',\
               '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026',\
               '#313695', '#4575b4', '#74add1',\
               '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026',\
               '#313695', '#4575b4', '#74add1',\
               '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']


    xlo, xhi, ylo, yhi, zlo, zhi = coords(DIR, z_name)
    vor = Voronoi(opos)
    fig = voronoi_plot_2d(vor, show_vertices=False, point_size=0.1)

    # add colours
    for region in vor.regions:
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]
            plt.fill(*zip(*polygon), c=colours[2*len(region)+1])
            #plt.fill
    plt.xlim(0, xhi)
    plt.ylim(ylo, yhi)
    plt.xlabel('$x$ (\AA)')
    plt.ylabel('$y$ (\AA)')
    plt.savefig('PLOTS_C/{}/voronoi_2d_{}.pdf'.format(DIR,z_name),
                bbox_inches='tight')
    #plt.show()

    return 


def main():


    args = GetArgs()

    DZ = args.width
    DIR = args.rho
    configs = args.start
    dens = args.density


    compname = ''
    count  = 0



    for dz in DZ:
        for c in configs:
            if dens == None:
                z_name = 'spce_T298_z{}_eps1.0_{}'.format(dz, c)
                print '\n\n#---------------------Reading in dz = {}, DEN = {} --------------------#'.format(dz,c)
                pos = read_xy(DIR, z_name)
                latest_pos = pos[-1]
                plotting(DIR, z_name, latest_pos)
                

            else:
                for r in dens:
                    z_name = 'spce_T298_z{}_r{}_eps1.0_{}'.format(dz, r, c)
                    avg_density(DIR, z_name)
                    pos = read_xy(DIR, z_name)
                    latest_pos = pos[-1]
                    plotting(DIR, z_name, latest_pos)

        count += 1

    return

if __name__ == "__main__":
    sys.path.insert(0,'/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()

