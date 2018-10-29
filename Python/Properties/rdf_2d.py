#! /usr/bin/env python
'''Hexagonal and square order parameter script by Maziar Fayas-Torshizi, amended by Frederike Jaeger'''

from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt
import argparse

def GetArgs():
    ## Parse command line arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-z', '--width', nargs='+', required=False, default='liquid', action='store',
                       help='Width of channel')
    parser.add_argument('-rho', '--rho', required=False, default='RHO1', action='store',
                       help='Average channel density')
    parser.add_argument('-u', '--units', required=False, type=str, action='store', default='metal',
                        help='metal or real')
    parser.add_argument('-s', '--start', nargs='+',required=False, default='None', action='store',
                       help='Starting configurations')
    parser.add_argument('-r', '--density', nargs='+',required=False, default=None, action='store',
                       help='Density for constant channel')

    args = parser.parse_args()
    return args


def apply_pbc(arr):
    Lx = 37.9
    Ly = 32.8


    # prepare arrays
    # plus
    pbcxp = np.copy(arr)
    pbcyp = np.copy(arr)
    pbcxyp = np.copy(arr)
    # minus
    pbcxm = np.copy(arr)
    pbcym = np.copy(arr)
    pbcxym = np.copy(arr)
    # mixed
    pbcxmyp = np.copy(arr)
    pbcxpym = np.copy(arr)


    # apply pbc
    pbcxp[:,0] += Lx
    pbcyp[:,1] += Ly
    pbcxyp[:,0] += Lx
    pbcxyp[:,1] += Ly

    pbcxm[:,0] -= Lx
    pbcym[:,1] -= Ly
    pbcxym[:,0] -= Lx
    pbcxym[:,1] -= Ly

    pbcxpym[:,0] += Lx
    pbcxpym[:,1] -= Ly
    pbcxmyp[:,0] -= Lx
    pbcxmyp[:,1] += Ly

    arr = np.concatenate((arr, pbcxp, pbcyp, pbcxyp, pbcxm, pbcym, pbcxym, pbcxpym, pbcxmyp),axis=0)

    return arr

def cut_and_shift(arr):

    arr = arr[np.all(arr[:,:2]>-32,axis=1) & np.all(arr[:,:2]<64,axis=1)]
    arr[:,0] += 32
    arr[:,1] += 32
    

    return arr


def rdf_calc(DIR, dt, z_name):


    #------------------------Code---------------------------#
    orderArr  = []
    fig_size_sq   = (9,7)
    

    # Loop over different times

    filename = 'Water_Graphene/{}/xyz.{}.xyz'.format(DIR,z_name)

    steps = 0
    positions = []
    rdf_tot = []
    with open(filename) as infile:
        for line in infile:
            l = line.split()
            if len(l) == 1:
                if steps != 0 and positions != []:
                    
                    positions = np.array(positions)
                    positions = apply_pbc(positions)
                    positions = cut_and_shift(positions)
                    xpos = positions[:,0]
                    ypos = positions[:,1]
                    zpos = positions[:,2]

                    extent = 96
                    rMax = 20
                    dr = 0.1
                    rdf = pairCorrelationFunction_2D(xpos, ypos, extent, rMax, dr)
                    print rdf[0]
                    print 'Number of reference atoms:', len(rdf[2])
                    rdf_tot.append(rdf[0])
                    # Plot RDF
                    #fig1 = plt.figure(figsize=fig_size_sq)
                    #ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])
                    
                    #ax1.plot(rdf[1], rdf[0])
                    #ax1.set_xlabel('r (\AA)')
                    #ax1.set_ylabel('g(r)')
                    #ax1.set_ylim( (0, 1) ) 
                    #ax1.legend()

                    #plt.show()
                steps += 1
                no_atoms = int(l[0])
                positions = []
            elif len(l) == 3:
                time = int(l[2])*dt
                print '\nAt step {} the time is {} ps.'.format(steps,time)

            elif len(l) == 4 and l[0]=='O':
                x_val = float(l[1])
                y_val = float(l[2])
                z_val = float(l[3])

                positions.append([x_val,y_val,z_val])


    rdf_mean = np.mean(np.array(rdf_tot), axis=0)
    rdf_xval = rdf[1]

    # save rdf to file
    write_values = np.array([rdf_xval, rdf_mean]).T
    np.savetxt('DATA/{}/rdf_2d_{}.dat'.format(DIR,z_name), write_values, fmt='%.6f, %.6f', header='r, rdf')

    fig1 = plt.figure(figsize=fig_size_sq)
    ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])
    
    ax1.plot(rdf_xval, rdf_mean)
    ax1.set_xlabel('r (\AA)')
    ax1.set_ylabel('g(r)')
    ax1.set_xlim( (0, rMax) ) 
    ax1.legend()
    #plt.show()

    fig1.savefig('PLOTS_C/{}/rdf_2d_{}.pdf'.format(DIR,z_name),bbox_inches='tight')
    return rdf_xval, rdf_mean


def main():


    args = GetArgs()

    DZ = args.width
    DIR = args.rho
    units = args.units
    configs = args.start
    dens = args.density


    if units == 'metal':
        time_conv = 1e-12
        space_conv = 1e-10
        dt = 0.0005
    elif units == 'real':
        time_conv = 1e-15
        space_conv = 1e-10
        dt = 2
    elif units == 'lj':
        time_conv = 1
        space_conv = 1
        dt = 0.001
    else:
        print 'Invalid unit.'
        sys.exit(1)



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


    ls = ['--','-', '--', '-.', ':','-', '--', '-.', ':']

    fig2 = plt.figure(figsize=fig_size_sq)
    ax2  = fig2.add_axes([0.15,0.15,0.75,0.75])

    compname = ''
    count  = 0

    for dz in DZ:
        compname += dz
        for c in configs:
            if dens == None:
                print '\n\n#---------------------Reading in dz = {}, DEN = {} --------------------#'.format(dz,c)

                z_name = 'spce_T298_z{}_eps1.0_{}'.format(dz, c)
                rdf_xval, rdf_mean = rdf_calc(DIR, dt, z_name)
                if count == 0:
                    ccount = 0
                else:
                    ccount = 2*count+4
                if dz=='6':
                    Hlabel = 'liquid'
                elif dz=='7.5':
                    Hlabel='frozen'
                else:
                    Hlabel='H = {} \AA'.format(dz)
                ax2.plot(rdf_xval, rdf_mean, linestyle = ls[count], c=colours[ccount], label=Hlabel)
            else:
                for r in dens:
                    print '\n\n#---------------------Reading in dz = {}, DEN = {}, r= {} --------------------#'.format(dz,c,r)

                    z_name = 'spce_T298_z{}_r{}_eps1.0_{}'.format(dz, r, c)
                    rdf_xval, rdf_mean = rdf_calc(DIR, dt, z_name)


            
            
        count += 1

    rMax = 20
    ax2.set_xlabel('r (\AA)')
    ax2.set_ylabel('g(r)')
    ax2.set_xlim( (0, rMax) ) 
    ax2.legend()
    #plt.show()

    fig2.savefig('PLOTS_C/{}/rdf_2d_comp_{}.pdf'.format(DIR,compname),bbox_inches='tight')

    return

if __name__ == "__main__":
    sys.path.insert(0,'/home/fred/SCRIPTS/Python')
    from plotting_params import *
    sys.path.insert(0,'/home/fred/SCRIPTS/Other/Internet')
    from paircorrelation_ShockSolution import *
    main()