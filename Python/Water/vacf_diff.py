#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
from scipy.integrate import simps

def data(A_dict):
    print 'Reading in data...'
    filename = "pair.local"
    
    # read in data
    f = open(filename,'r')
    data = f.read()
    
    data_lines = data.split('\n')
    
    type1 = []
    type2 = []
    Int = []
    Dist = []
    Dist2 = []
    eV2J = 1.602176565e-19
    A2m = 1e-10

    Natoms = 0
    count = 0
    for i in range(9,len(data_lines)-1): 
        if count < 100:     
            try:
                if data_lines[i]=='ITEM: TIMESTEP':
                    count += 1
                    print 'Reading in timestep:', count
                type1 = data_lines[i].split()[1]
                type2 = data_lines[i].split()[2]
                Dist.append(float(data_lines[i].split()[3]))
                if (A_dict[type1]=='C' and A_dict[type2] == 'O') or (A_dict[type1]=='O' and A_dict[type2] == 'C'):
                    Dist2.append(float(data_lines[i].split()[3])*A2m)
                    Int.append(float(data_lines[i].split()[4])*eV2J)
                    Natoms += 1
            except:
                continue

    
    Eint = sum(Int)/count
    print 'Finished reading in data. %i data points collected.' % Natoms


    return Int, Eint, Dist2, count

def Atom_lookup2():
    filename = "data.graphene"
    # read in data
    f = open(filename,'r')
    data = f.read()

    atom_dict = {}
    
    data_lines = data.split('\n')
    blank_count = 0

    for i in range(0,len(data_lines)-1):
        if blank_count == 5 and len(data_lines[i].split()) > 1:
            atom_dict['%s'%data_lines[i].split()[0]] = data_lines[i].split()[2]
        if data_lines[i].split() == []:
            blank_count +=1

    return atom_dict

def Atom_lookup():
    filename = "water.xyz"
    # read in data
    f = open(filename,'r')
    data = f.read()

    atom_dict = {}
    
    data_lines = data.split('\n')
    count = 0
    atom_count = 0

    for i in range(0,len(data_lines)-1):
        if count == 2 and len(data_lines[i].split()) > 1:
            atom_count += 1
            atom_dict['%i'%atom_count] = data_lines[i].split()[0]
        if data_lines[i].split()[0]!= ('C' or 'H' or 'O'):
            count +=1

    return atom_dict

def Bins(dist, pot, no_bins):
    dist = np.array(dist)
    pot = np.array(pot)
    hist, bins = np.histogram(dist, bins=no_bins, density=False)

    # Find number density
    density = []
    for i in range(len(hist)):
        density.append(hist[i]/(bins[i+1]-bins[i]))

    

    regular_dist = bins # Regular heights every 100m
    regular_pot = []

    # Find average potential values
    for i in range(len(regular_dist)-1):
        mask = np.logical_and(dist > regular_dist[i], dist < regular_dist[i+1])
        mean = np.mean(pot[mask])
        regular_pot.append(mean)

    regular_pot = np.hstack((regular_pot))

    # find average bin values
    avg_bins = []
    for i in range(len(bins)-1):
        avg_bins.append((bins[i+1]+bins[i])/2)

    return density, avg_bins, regular_pot

def Area():

    filename = "data.graphene"#"graphene_slit_9_9_40_1.20.dat"
    # read in data
    f = open(filename,'r')
    data = f.read()
    
    data_lines = data.split('\n')
    xlo, xhi = float(data_lines[11].split()[0]), float(data_lines[11].split()[1])
    ylo, yhi = float(data_lines[12].split()[0]), float(data_lines[12].split()[1])

    A = (xhi-xlo)*(yhi-ylo)*10**(-20) # in m^2

    return A

def find_region(dist, pot):
    dist, pot = np.array(dist), np.array(pot)

    mask = np.logical_and(dist > 0, dist < 4*10**(-10))
    new_dist = dist[mask]
    new_pot = pot[mask]
    new_Eint = sum(new_pot)
    return new_Eint, new_pot, new_dist


def main():
    print 'Creating atom type dictionary'
    Atom_dict = Atom_lookup()
    Int, Eint, Dist, ts = data(Atom_dict)
    A = Area()
    dens, bins, pot = Bins(Dist, Int, 100)

    Eint_layer, int_layer, dist_layer = find_region(Dist, Int)
    dens_layer, bins_layer, pot_layer = Bins(dist_layer, int_layer, 20)

    WA = Eint/(2*A)
    WA_layer = Eint_layer/(2*A)
    print 'The work of adhesion is: ', WA
    print 'The work of adhesion in the layer is: ', WA_layer

    # integrate V(z)*n(z)
    Vn = [] 
    Vn_layer = []
    for i in range(len(pot)):
        Vn.append(pot[i]*dens[i])
    WA_int = simps(Vn, bins)/(2*ts*A)

    for i in range(len(pot_layer)):
        Vn_layer.append(pot_layer[i]*dens_layer[i])
    WA_int_layer = simps(Vn_layer, bins_layer)/(2*ts*A)

    print 'The integrated value of WA is:', WA_int
    print 'The integrated value of WA in the layer is:', WA_int_layer


    plt.plot(bins, Vn)
    #plt.show()

    plt.plot(Dist, Int, linestyle='None', marker='x')
    plt.plot(bins, pot, linestyle='None', marker='D')
    #plt.show()

    plt.plot(bins, dens, linestyle='None', marker='x')
    #plt.show()



if __name__ == "__main__":
    main()


