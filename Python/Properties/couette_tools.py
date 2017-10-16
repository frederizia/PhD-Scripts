#! /usr/bin/env python
# Quick diffusion coefficient calculation given VACF

import argparse
import matplotlib.pyplot as plt
import matplotlib
import sys
import re
import numpy as np
from scipy.integrate import simps

def read_log(f):
    '''Code to read in log file'''
    #print 'Reading in data for', f

    file = open('Water_Graphene/log.{}_1'.format(f),'r')

    data = file.read()
    data_lines = data.split('\n')

    DATA = []
    flag=0
    for i in range(len(data_lines)-1):
        if data_lines[i].split() != [] and data_lines[i].split()[0] == 'unfix' and (data_lines[i].split()[1] == 'RESCALE' or data_lines[i].split()[1] == 'all_rescale'):
            flag = 1
        if flag == 1 and data_lines[i].split() != []:
            DATA.append(data_lines[i].split())
    #print 'Finished.' 
    return DATA



def mid_point(Y,VAR, tol=0.03):
    '''Code to find the mid point in the array'''

    # Find mid point
    mid = int(len(Y)/2)



    for k in range(1,mid):
        if abs(VAR[k-1]-VAR[k]) >= tol:           
            LEFT = k
            break
    for k in range(1,mid):        
        if abs(VAR[::-1][k-1]-VAR[::-1][k]) >= tol:           
            RIGHT = len(Y)-k-1
            break

    mid_p = int(LEFT+(RIGHT - LEFT)/2)
    return mid_p, LEFT, RIGHT

def rho_wall_max(RHO, LEFT):
    '''Code to find the mid point in the array'''
    rho_wall = RHO[LEFT:LEFT+30]
    rho_max = np.max(rho_wall)
    return rho_max

def read_velprof(folder, pre, f):
    '''Code to read in density  data from a 1d LAMMPS output file'''

    filename = '{}/{}.{}'.format(folder,pre,f)
    f = open(filename,'r')

    data = f.read()
    data_lines = data.split('\n')

    zcoord = []
    zcoord_tmp =[]
    dens = []
    dens_tmp = []
    vx = []
    vx_tmp = []

    idx = 0
    count  = 0
    count_tmp = 1

    for j in range(4,len(data_lines)-1):
        if len(data_lines[j].split()) == 3:
            # Count number of timesteps collected
            count += 1
            zcoord.append(zcoord_tmp)
            dens.append(dens_tmp)
            vx.append(vx_tmp)
            zcoord_tmp = []
            dens_tmp = []
            vx_tmp = []

        elif len(data_lines[j].split()) != 3:
            z_val = float(data_lines[j].split()[1])
            dens_val = float(data_lines[j].split()[3])
            vx_val = float(data_lines[j].split()[4])
            zcoord_tmp.append(z_val)
            dens_tmp.append(dens_val)
            vx_tmp.append(vx_val)


    zcoord = np.array(zcoord[0])
    dens   = np.mean(np.array(dens), axis=0)
    vx   = np.mean(np.array(vx), axis=0)

    return zcoord, dens, vx

def derivative(X, Y):
    diff = np.diff(Y)
    dx = X[1]-X[0]
    der = diff/dx
    return der


def coords(folder,pre,f):
    filename = '{}/{}.{}'.format(folder,pre,f)
    f = open(filename,'r')

    data = f.read()
    data_lines = data.split('\n')


    coord_line = re.sub('[\(\)]', '', data_lines[41]).split('=')[1].split('to')
    xlo, ylo, zlo = coord_line[0].strip(' ').split(' ')
    xhi, yhi, zhi = coord_line[1].strip(' ').split(' ')

    return float(xlo), float(xhi), float(ylo), float(yhi), float(zlo), float(zhi)

def geometry(folder):
    xlo, xhi, ylo, yhi, zlo, zhi = coords(folder,'log', 'lammps')
    area_xy = (xhi-xlo)*(yhi-ylo)
    area_flow = (yhi-ylo)*(zhi-zlo)
    volume = (xhi-xlo)*(yhi-ylo)*(zhi-zlo)
    return area_xy, area_flow, volume


def stress_prof(folder,pre,f):

    filename = '{}/{}.{}'.format(folder,pre,f)
    file = open(filename,'r')
    data_lines = file.readlines()
    file.close()

    area = geometry(folder)[0]

    count_tmp = 0
    dz = 0
    for line in data_lines:
        if line[0] != '#' and len(line.split()) != 2 and len(line.split()) != 3: 
            dz = float(line.split()[1]) - dz
            count_tmp +=1
            if count_tmp ==2:
                break
    vol = area*dz

    coords, pressure, pxy, pxz, pyz, pxx, pyy, pzz = [], [], [], [], [], [], [], []
    coords_tot, pressure_tot, pxy_tot, pxz_tot, pyz_tot, pxx_tot, pyy_tot, pzz_tot = [], [], [], [], [], [], [], []
    count_tmp = -1
    for line in data_lines:
        items = line.split()
        if len(items) == 3:
            count_tmp += 1
            if count_tmp == 0:
                continue
            else:
                coords_tot.append(coords)
                pressure_tot.append(pressure)
                pxy_tot.append(pxy)
                pxz_tot.append(pxz)
                pyz_tot.append(pyz)
                pxx_tot.append(pxx)
                pyy_tot.append(pyy)
                pzz_tot.append(pzz)
                coords, pressure, pxy, pxz, pyz, pxx, pyy, pzz = [], [], [], [], [], [], [], []
            

        elif line[0] != '#' and len(items) != 2 and len(items) != 3:
            Coord = float(items[1])
            Ncount = float(items[2])
            xx = float(items[3])*Ncount/vol
            yy = float(items[4])*Ncount/vol
            zz = float(items[5])*Ncount/vol
            xy = float(items[6])*Ncount/vol
            xz = float(items[7])*Ncount/vol
            yz = float(items[8])*Ncount/vol  
            Pxx = -xx
            Pyy = -yy
            Pzz = -zz
            Pxy = -xy
            Pxz = -xz
            Pyz = -yz
            pres = (Pxx+Pyy+Pzz)*0.101325/3

            coords.append(Coord)
            pressure.append(pres)
            pxy.append(Pxy)
            pxz.append(Pxz)
            pyz.append(Pyz)
            pxx.append(Pxx)
            pyy.append(Pyy)
            pzz.append(Pzz)

    coords_tot.append(coords)
    pressure_tot.append(pressure)
    pxy_tot.append(pxy)
    pxz_tot.append(pxz)
    pyz_tot.append(pyz)
    pxx_tot.append(pxx)
    pyy_tot.append(pyy)
    pzz_tot.append(pzz)
    coords = np.mean(np.array(coords_tot),axis=0)
    pressure = np.mean(np.array(pressure_tot),axis=0)
    pxx_tot = np.mean(np.array(pxx_tot),axis=0)
    pyy_tot = np.mean(np.array(pyy_tot),axis=0)
    pzz_tot = np.mean(np.array(pzz_tot),axis=0)
    pxy_tot = np.mean(np.array(pxy_tot),axis=0)
    pxz_tot = np.mean(np.array(pxz_tot),axis=0)
    pyz_tot = np.mean(np.array(pyz_tot),axis=0)
    zmax = np.max(coords)
  
    return coords, pressure, np.array(pressure_tot), np.array(pxy_tot),np.array(pxz_tot),np.array(pyz_tot), np.array(pxx_tot), np.array(pyy_tot), np.array(pzz_tot), dz  

def props(data):
    #print 'Reading in thermodynamic properties.' 
    temp = []
    press = []
    press_cum = []
    press_cum_new = []
    dP = []
    vol = []
    p_cum = 0
    count = 0
    for i in range(len(data)-1):
        if len(data[i]) == 12 and data[i][0] != 'Step':
            try:
                temp.append(float(data[i][4]))
                press.append(float(data[i][6]))
                press_cum.append(float(data[i][7]))
                if p_cum == 0:
                    p_cum = float(data[i][6])
                    press_cum_new.append(p_cum)
                    dP.append(0)
                    count += 1
                else:
                    count += 1
                    p_cum +=  float(data[i][6])
                    press_cum_new.append(p_cum/count)
                    dP.append(float(data[i][6])-(p_cum/count))
                vol.append(float(data[i][8]))
            except:
                pass
                #print('No data in this line.')
    temp, press, press_cum, press_cum_new, dP, vol = np.array(temp), np.array(press), np.array(press_cum), np.array(press_cum_new), np.array(dP), np.array(vol)
    #print 'Finished.' 
    return temp, press, press_cum, press_cum_new, dP, vol


def correlation(data, clen, sint):
    acf = []
    dt = 0.0005
    sampint = sint
    times_range = clen # want this to be between 100-1000 clen

    print 'Calculating the autocorrelation function for', dt*sampint*clen, 'ps.'
    # iterate over different time origins
    Cdata = np.zeros(times_range+1) # stores correlation for each time origin
    for i in range(0, len(data)-times_range,100): 
        C_tmp = 0 
        count = 0  
        for j in range(times_range+1): # the n in t_0+ndt
            C_tmp = np.dot(data[i], data[i+j])
            if count == 0:
                Cdata[j] = C_tmp
            else:
                Cdata[j] = (Cdata[j]+C_tmp)/2
        count += 1
            #Cdata.append(C_tmp)
        #if Cdata != []:
            #acf.append(Cdata)
    #acf = np.array(acf)
    print Cdata
    
    #mean_acf = np.mean(acf, axis=0)
    mean_acf = Cdata
    norm_acf = mean_acf/mean_acf[0]

    return mean_acf, norm_acf

def eta(f, data, delz):
    corrlen = 2000
    sampint = 10 #20
    dt = 0.0005
    acf, acf_norm  = correlation(data,corrlen, sampint)

    Log_data = read_log(f)
    T, P, Pc, Pcnew, dP, V = props(Log_data)

    Vave = np.mean(V)
    Tave = np.mean(T)
    kB = 1.38e-23
    A2m = 1e-10
    ps2s = 1e-12
    bar2Pa = 1e5
    convert = A2m**3*ps2s*bar2Pa**2
    Vnew = area = geometry(f)[0]*delz

    eta_scale = (Vnew*sampint*dt*convert)/(Tave*kB)
    ACF_int = simps(acf)
    ETA = eta_scale*ACF_int

    return ETA

def pmf_prof(f):

    filename = ilename = 'Water_Graphene/pmf.{}_1'.format(f)
    
    # read in data
    f = open(filename,'r')
    data = f.read()
    
    data_lines = data.split('\n')
    
    
    zbin = []
    pmf = []


    for j in range(1,len(data_lines)-1):
        if data_lines[j].split()[0][0] != '#':
            bin_val = float(data_lines[j].split()[0])
            pmf_val = float(data_lines[j].split()[1])

            zbin.append(bin_val)
            pmf.append(pmf_val)


    zbin = np.array(zbin)
    pmf = np.array(pmf)

    return zbin, pmf