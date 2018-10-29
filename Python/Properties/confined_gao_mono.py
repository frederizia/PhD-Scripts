#! /usr/bin/env python
# Plots for confined water data

import argparse
import matplotlib.pyplot as plt
import matplotlib
import csv
import sys
import numpy as np


def GetArgs():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-rho', '--rho', nargs='+', required=False, default=['RHO1'],
                        action='store',
                        help='Average channel density')
    args = parser.parse_args()
    return args


def main():

    args = GetArgs()

    rho = args.rho

    # Definitions
    fig_size_long = (14, 5)
    fig_size_sq = (9, 7)
    fig_size_sq2 = (7, 9)

    markers = ['o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P', 'X', '8', 'p',
               'o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P', 'X', '8', 'p',
               'o', 'D', 's', 'v', '^', 'd', '*', '>', '<', 'P', 'X', '8', 'p']
    colours = ['#313695', '#4575b4', '#74add1',
               '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026',
               '#313695', '#4575b4', '#74add1',
               '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026',
               '#313695', '#4575b4', '#74add1',
               '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']
    LABEL = {'RHO1': '$\\rho$ = 1 g/cm$^3$', 'RHO1.1': '$\\rho$ = 1.1 g/cm$^3$', 'RHOtest': 'H=6.8 \AA'}

    ls = ['--', '-', '--', '-.', ':', '-', '--', '-.', ':']

    NAVO = 6.02214086*1e23  # mol^-1
    MW = 18.01528  # g/mol

    # ------------------- Initialise figures -----------------------

    # rho (z)
    fig1 = plt.figure(figsize=fig_size_sq)
    ax1 = fig1.add_axes([0.15, 0.15, 0.75, 0.75])

    ax1.set_xlabel('H (\AA)')
    ax1.set_ylabel('$\\rho$ (g/cm$^3$)')


    # read GAO data
    print '\n#-----------------------Reading in Gao------------------\n#'

    for i in [1, 2]:
        H_GAO = []
        RHO_GAO = []
        with open('DATA/gao2018mono{}.csv'.format(i), 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:

                H = float(row[0])*10
                H_eff_tmp = (H-3.19)
                H_GAO.append(H)

                RHO_tmp = float(row[1])*MW/(NAVO*H_eff_tmp*1e-8*1e-14)
                RHO_GAO.append(RHO_tmp)
                print H_eff_tmp, RHO_tmp

        ax1.plot(H_GAO, RHO_GAO, marker='o', c=colours[0])
        



    # Read sim data
    DIR_name = ''

    count = 0
    for DIR in rho:
        DIR_name += DIR
        H_list, RHO_list = [], []
        print '\n#-----------------------Reading in {}------------------\n#'.format(DIR)

        with open('DATA/{}/rho_H_mono.csv'.format(DIR), 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                if row[0] != 'rho':
                    H = float(row[1])
                    H_list.append(H)

                    RHO_tmp = float(row[0])
                    RHO_list.append(RHO_tmp)
        ax1.plot(H_list, RHO_list, linestyle='None', marker='x', label=LABEL[DIR], c=colours[5+count])
        count += 1

    ax1.legend()
    fig1.savefig('PLOTS_C/gao_comp_{}.pdf'.format(DIR_name),
                  box_inches='tight')

    print '\n#-----------------------Done------------------\n#'

    return


if __name__ == "__main__":
    sys.path.insert(0, '/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()
