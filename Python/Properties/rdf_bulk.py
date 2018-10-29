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
    parser.add_argument('-u', '--units', required=False, type=str, action='store', default='metal',
                        help='metal or real')
    parser.add_argument('-s', '--start', nargs='+',required=False, default='None', action='store',
                       help='Starting configurations')
    parser.add_argument('-d', '--dim', nargs='+',required=False, default='None', action='store',
                       help='Dimension')

    args = parser.parse_args()
    return args


def apply_pbc(arr):
    Lx = 39.44926
    Ly = 39.44926
    Lz = 39.44926


    # prepare arrays
    # plus
    pbcxp = np.copy(arr)
    pbcyp = np.copy(arr)
    pbczp = np.copy(arr)
    pbcxyp = np.copy(arr)
    pbcxzp = np.copy(arr)
    pbcyzp = np.copy(arr)
    # minus
    pbcxm = np.copy(arr)
    pbcym = np.copy(arr)
    pbczm = np.copy(arr)
    pbcxym = np.copy(arr)
    pbcxzm = np.copy(arr)
    pbcyzm = np.copy(arr)
    # mixed
    pbcxmyp = np.copy(arr)
    pbcxpym = np.copy(arr)
    pbcxmzp = np.copy(arr)
    pbcxpzm = np.copy(arr)
    pbczmyp = np.copy(arr)
    pbczpym = np.copy(arr)


    # apply pbc
    pbcxp[:,0] += Lx
    pbcyp[:,1] += Ly
    pbczp[:,2] += Lz
    pbcxyp[:,0] += Lx
    pbcxyp[:,1] += Ly
    pbcxzp[:,0] += Lx
    pbcxzp[:,2] += Lz
    pbcyzp[:,2] += Lz
    pbcyzp[:,1] += Ly

    pbcxm[:,0] -= Lx
    pbcym[:,1] -= Ly
    pbcym[:,2] -= Lz
    pbcxym[:,0] -= Lx
    pbcxym[:,1] -= Ly
    pbcxzm[:,0] -= Lx
    pbcxzm[:,2] -= Lz
    pbcyzm[:,2] -= Lz
    pbcyzm[:,1] -= Ly

    pbcxpym[:,0] += Lx
    pbcxpym[:,1] -= Ly
    pbcxmyp[:,0] -= Lx
    pbcxmyp[:,1] += Ly
    pbcxpzm[:,0] += Lx
    pbcxpzm[:,2] -= Lz
    pbcxmzp[:,0] -= Lx
    pbcxmzp[:,2] += Lz
    pbczpym[:,2] += Lz
    pbczpym[:,1] -= Ly
    pbczmyp[:,2] -= Lz
    pbczmyp[:,1] += Ly

    arr = np.concatenate((arr, pbcxp, pbcyp, pbczp, 
                          pbcxyp, pbcxzp, pbcyzp,
                          pbcxm, pbcym, pbczm,
                          pbcxym, pbcxzm, pbcyzm,
                          pbcxpym, pbcxmyp, pbcxpzm, pbcxmzp,
                          pbczpym, pbczmyp),axis=0)

    return arr

def cut_and_shift(arr):

    arr = arr[np.all(arr[:,:3]>-38,axis=1) & np.all(arr[:,:3]<76,axis=1)]
    arr[:,0] += 38
    arr[:,1] += 38
    arr[:,2] += 38

    return arr




def main():


    args = GetArgs()

    units = args.units
    configs = args.start
    dims = args.dim


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
    P_name = 1

    for d in dims:
        for c in configs:

            print '\n\n#---------------------Reading in DEN = {} --------------------#'.format(c)


            #------------------------Code---------------------------#
            orderArr  = []
            z_name = 'spce_T298_P1_{}'.format(c)

            # Loop over different times

            filename = 'Water/xyzfull.{}.xyz'.format(z_name)

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

                            extent = 114
                            rMax = 20
                            dr = 0.1
                            if d == '2':
                                rdf = pairCorrelationFunction_2D(xpos, ypos, extent, rMax, dr)
                            elif d == '3':
                                rdf = pairCorrelationFunction_3D(xpos, ypos, zpos, extent, rMax, dr)
                            else:
                                raise ValueError('Incorrect dimension')
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

            fig1 = plt.figure(figsize=fig_size_sq)
            ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])
            
            ax1.plot(rdf_xval, rdf_mean)
            ax1.set_xlabel('r (\AA)')
            ax1.set_ylabel('g(r)')
            ax1.set_xlim( (0, rMax) ) 
            ax1.legend()
            #plt.show()

            fig1.savefig('PLOTS/rdf_{}d_P{}.pdf'.format(d,P_name),bbox_inches='tight')

            if count == 0:
                ccount = 0
            else:
                ccount = 2*count+4
            Hlabel='Bulk'
            ax2.plot(rdf_xval, rdf_mean, linestyle = ls[count], c=colours[ccount], label=Hlabel)
        count += 1

    # ax2.set_xlabel('r (\AA)')
    # ax2.set_ylabel('g(r)')
    # ax2.set_xlim( (0, rMax) ) 
    # ax2.legend()
    # plt.show()

    # fig2.savefig('PLOTS/rdf_3d_comp_{}.pdf'.format(DIR,compname),bbox_inches='tight')

    return

if __name__ == "__main__":
    sys.path.insert(0,'/home/fred/SCRIPTS/Python')
    from plotting_params import *
    sys.path.insert(0,'/home/fred/SCRIPTS/Other/Internet')
    from paircorrelation_ShockSolution import *
    main()