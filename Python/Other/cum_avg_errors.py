#!/usr/bin/python
'''Script to extract values and their errors from cumultative averages from \
LAMMPS correlative functions'''
from __future__ import division
import argparse
import os
import numpy as np
  

def Get_Args():
    ## Parse arguments form command line
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", type=str, \
                        help="Data filename",required=True)
    parser.add_argument("-i", type=int, \
                        help="Index of quantity to determine",required=True) 
    parser.add_argument("-l", type=str, \
                        help="Label of file",required=False, default='values') 
    args = parser.parse_args()
    return args

def read_file(fname, idx):
    cwd = os.getcwd()

    f = open('%s/%s'%(cwd,fname),'r')
    data = f.read()
    data_lines = data.split('\n')
    
    vals = []
    t_steps = []

    for i in range(len(data_lines)):
        if '#' in data_lines[i]:
            continue
        else:
            line_items = data_lines[i].split()
            if line_items == []:
                continue
            else:
                t_steps.append(float(line_items[0]))
                vals.append(float(line_items[idx-1]))
    return np.array(t_steps), np.array(vals)

def split_values(vals):

    N = 1
    single_vals = []
    single_vals.append(vals[0])
    for i in range(1,len(vals)):
        N += 1
        true_v = N*vals[i] - (N-1)*vals[i-1]
        single_vals.append(true_v)

    return np.array(single_vals)

def error_calc(p_vals):
    N = len(p_vals)
    stdev = np.std(p_vals)
    err = stdev/N
    return err

def write_file(t_steps, vals, p_vals, lbl, idx):
    cwd = os.getcwd()

    f = open('%s/%s_%i.dat'%(cwd,lbl, idx),'w+')
    
    write_values = np.array([t_steps, vals, p_vals]).T
    np.savetxt(f, write_values, header='cum. avg. values | single values')
 
    return 


def main():
    args = Get_Args()
    filename = args.f
    index = args.i
    label = args.l

    # read file and extract values
    timesteps, values = read_file(filename, index)

    # disentangle values
    pure_values = split_values(values)

    # calculate error
    error = error_calc(pure_values)

    print 'The final averaged value is', values[-1]
    print 'The mean is', pure_values.mean()
    print 'The error is', error

    write_file(timesteps, values, pure_values, label, index)

    return

if __name__ == '__main__':
    main()

