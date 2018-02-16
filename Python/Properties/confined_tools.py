#! /usr/bin/env python
# Quick diffusion coefficient calculation given VACF

import argparse
import matplotlib.pyplot as plt
import matplotlib
import sys
import re
import numpy as np
from scipy.integrate import simps
from scipy import stats
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
import pandas as pd

def read_log(f):
    '''Code to read in log file'''
    #print 'Reading in data for', f

    file = open('Water_Graphene/log.{}'.format(f),'r')

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

def diffusion(vacf, DT, t_conv, s_conv, pre):

    int_vacf = simps(vacf, dx=DT)

    diff = int_vacf/3
    diff_2d = int_vacf/2

    conv = s_conv**2/t_conv

    diff_si = diff*conv
    diff_2d_si = diff_2d*conv

    if pre == 'C_vv':
        print 'The self-diffusion coefficient in 3D is', diff_si
        diff_ret = diff_si
    else:
        print 'The self-diffusion coefficient in 2D is', diff_2d_si
        diff_ret = diff_2d_si

    return diff_ret

def mid_point(Y,VAR):
    '''Code to find the mid point in the array'''

    # Find mid point
    mid = int(len(Y)/2)
    tol = 0.01

    for k in range(1,mid):
        if abs(VAR[k-1]-VAR[k]) >= tol:           
            LEFT = k
            break
    for k in range(1,mid):        
        if abs(VAR[::-1][k-1]-VAR[::-1][k]) >= tol:           
            RIGHT = len(Y)-k
            break

    mid_p = int(LEFT+(RIGHT - LEFT)/2)
    return mid_p, LEFT, RIGHT

def rho_wall_max(RHO, LEFT):
    '''Code to find the mid point in the array'''
    rho_wall = RHO[LEFT:LEFT+30]
    rho_max = np.max(rho_wall)
    return rho_max

def read_densprof(f):
    '''Code to read in density  data from a 1d LAMMPS output file'''

    filename = 'Water_Graphene/densprof.{}'.format(f)
    f = open(filename,'r')

    data = f.read()
    data_lines = data.split('\n')

    zcoord = []
    zcoord_tmp =[]
    dens = []
    dens_tmp = []

    idx = 0
    count  = 0
    count_tmp = 1

    for j in range(4,len(data_lines)-1):
        if len(data_lines[j].split()) == 3:
            # Count number of timesteps collected
            count += 1
            zcoord.append(zcoord_tmp)
            dens.append(dens_tmp)
            zcoord_tmp = []
            dens_tmp = []

        elif len(data_lines[j].split()) != 3:
            z_val = float(data_lines[j].split()[1])
            dens_val = float(data_lines[j].split()[4])
            zcoord_tmp.append(z_val)
            dens_tmp.append(dens_val)


    zcoord = np.array(zcoord[0])
    dens   = np.mean(np.array(dens), axis=0)

    return zcoord, dens

def read_densprof_2d(f,atom):

    if atom==None:
        filename = 'Water_Graphene/dens.{}'.format(f)
    else:
        filename = 'Water_Graphene/dens{}.{}'.format(atom,f)
    
    
    ybin = []
    zbin = []
    dx = 0.2
    ylo = -0.5
    yhi = 37.3
    zlo = -0.5
    zhi = 32.3

    ylist = list(np.arange(ylo,yhi+dx,dx))
    zlist = list(np.arange(zlo, zhi+dx, dx))

    # round values in list so that they agree with data
    ylist = [round(n, 1) for n in ylist]
    zlist = [round(n, 1) for n in zlist]


    ydim = len(ylist)
    zdim = len(zlist)

    density  = np.zeros((ydim,zdim))
    coords_y = np.zeros((ydim,zdim))
    coords_z = np.zeros((ydim,zdim))


    count = 1
    count_tmp = 0
    count_tmp2 = 0
    with open(filename) as infile:
        for line in infile:
            if len(line.split()) != 2 and len(line.split()) != 3 and line.split()[0] != '#':
                y_val = float(line.split()[1])
                z_val = float(line.split()[2])
                den_val = float(line.split()[5])


                # only select values within slit
                if y_val in ylist and z_val in zlist:
                    y_ind = ylist.index(y_val)
                    z_ind = zlist.index(z_val)
                    #print z_ind, z_val
                    ybin.append(y_val)
                    zbin.append(z_val)

                    # average over all time steps
                    if count_tmp == 0:
                        density[y_ind][z_ind] = den_val
                    else:
                        density[y_ind][z_ind] = (density[y_ind][z_ind]+den_val)/2
                    coords_y[y_ind][z_ind] = ylist[y_ind]
                    coords_z[y_ind][z_ind] = zlist[z_ind]
                    count_tmp += 1
            count =2
      

    #density_z = np.average(density, axis = 0)
    #density_z_err = stats.sem(density, axis = 0)
    coords_z = coords_z[0,:] #np.average(coords_z, axis = 0)
    #density_y = np.average(density, axis = 1)
    coords_y = coords_y[:,0] #np.average(coords_x, axis = 1)


    return density, coords_y, coords_z

def coords(f):
    filename = 'Water_Graphene/log.{}'.format(f)
    f = open(filename,'r')

    data = f.read()
    data_lines = data.split('\n')

    coord_line = re.sub('[\(\)]', '', data_lines[45]).split('=')[1].split('to')
    xlo, ylo, zlo = coord_line[0].strip(' ').split(' ')
    xhi, yhi, zhi = coord_line[1].strip(' ').split(' ')

    return float(xlo), float(xhi), float(ylo), float(yhi), float(zlo), float(zhi)

def geometry(f):
    xlo, xhi, ylo, yhi, zlo, zhi = coords(f)
    area_xy = (xhi-xlo)*(yhi-ylo)
    area_flow = (yhi-ylo)*(zhi-zlo)
    volume = (xhi-xlo)*(yhi-ylo)*(zhi-zlo)
    return area_xy, area_flow, volume


def stress_prof(f,xtra='None'):
    if xtra == 'None':
        filename = 'Water_Graphene/stress.{}'.format(f)
    else:
        filename = 'Water_Graphene/stress{}.{}'.format(xtra,f)
    file = open(filename,'r')
    data_lines = file.readlines()
    file.close()

    area = geometry(f)[0]

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
            pres = (Pxx+Pyy+Pzz)*0.1/3 # in MPa

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
    print data.shape
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
    ETA = eta_scale*ACF_int # in Pas

    return ETA

def eta_diff(diff):
    kB = 1.38*1e-23
    a = 1.7*1e-10
    T = 298
    return (kB*T)/(3*np.pi*a*diff)

def visc_gk(file, H, visc_type, eta_type):
    filename = 'Water_Graphene/visc.{}'.format(file)
    df = pd.read_csv(filename, delimiter=' ', skiprows=2)
    #read first two lines
    with open(filename, 'r') as f:
        _, line2 = f.readline(), f.readline()
    cols = line2.lstrip('#').strip().split(' ')
    df.columns = cols
    Shear = float(df['v_{}_{}'.format(visc_type,eta_type)].tolist()[-1])
    #area, area_f, vol = geometry(file)
    #print 'vol = ', vol
    #print 'volf = ', area*H
    #Shear_channel = Shear * (area*H)/vol 
    #print Shear, Shear_channel
    return Shear #_channel

def eta_gk_acf(file, H):
    corrlen = 2000
    sampint = 10 #20
    dt = 0.0005
    tlim = 16000000
    kB = 1.38e-23

    # find volume and T
    area, area_f, vol = geometry(file)
    volf = area*H

    Log_data = read_log(file)
    T, P, Pc, Pcnew, dP, V = props(Log_data)
    Tave = np.mean(T)

    # conversion
    A2m = 1e-10
    ps2s = 1e-12
    bar2Pa = 1e5
    convert = A2m**3*ps2s*bar2Pa**2

    filename = 'Water_Graphene/acfsv.{}'.format(file)
    #read first two lines
    with open(filename, 'r') as f:
        _, line2, line3 = f.readline(), f.readline(), f.readline()
    cols = line3.lstrip('#').strip().split(' ')
    df = pd.read_csv(filename, delimiter=' ', names=cols,skiprows=4)
    #df.columns = cols
    try:
        final_slice = df[df.Index == tlim].index[0]
    except:
        final_slice = df[df.Index == 14000000].index[0]
    df = df[final_slice+1:] # only use final average
    acfsv = df['v_pxy*v_pxy'].tolist()
    acfsv_normal = df['v_pxz*v_pxz'].tolist()
    int_acfsv = simps(acfsv) #, dx=dt*sampint)
    

    # Calculate eta
    eta_scale = (volf*sampint*dt*convert)/(Tave*kB)
    Shear = int_acfsv*eta_scale

    # plot ACF parrallel and normal to compare behaviour
    fig1 = plt.figure(figsize=(9,7)) 
    ax1  = fig1.add_axes([0.1,0.15,0.8,0.75])
    ax1.plot(np.arange(len(acfsv))*20, acfsv/acfsv[0], label = 'parallel')
    ax1.plot(np.arange(len(acfsv))*20, acfsv_normal/acfsv_normal[0], label = 'normal')
    ax1.set_xlabel('$\Delta$ t')
    ax1.set_ylabel('Normalised ACF')
    ax1.set_xlim(0,5000)
    ax1.legend()
    fig1.savefig('PLOTS_C/ACF_{}.pdf'.format(file),bbox_inches='tight')
    fig1.clear()

    return Shear 

def pmf_prof(f):

    filename = ilename = 'Water_Graphene/pmf.{}'.format(f)
    
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

def exp_fit(x, y, xmin, xmax, guess=[1,1,1]):
    params, cov = curve_fit(expf, np.array(x), np.array(y), p0=guess)
    xdat = np.linspace(xmin, xmax, 100)
    fit = expf(xdat,*params)
    return xdat, fit

def expf(x, A,B,C):
    return A*np.exp(-B*x)+C

def averaging(k,v):
    averages = {}
    errors = {}
    counts = {}
    for name, value in zip(k, v):
        if name in averages:
            averages[name].append(value)
            counts[name] += 1
        else:
            averages[name] = [value]
            counts[name] = 1
    for name in averages:
        ave_array = np.array(averages[name])
        ave_val = np.mean(ave_array)
        ave_stdev = np.std(ave_array)
        ave_error = ave_stdev/np.sqrt(float(counts[name]))
        averages[name] = ave_val#averages[name]/float(counts[name])
        errors[name] = ave_error 

    k = map(float, averages.keys())
    v = averages.values()
    e = errors.values()
    sorti = np.argsort(k)
    k, v, e = np.array(k)[sorti], np.array(v)[sorti], np.array(e)[sorti]
    return k,v,e

def diff_msd(file):
    filename = 'Water_Graphene/msd.{}'.format(file)
    tlim = 8*1e6
    dt = 0.0005
    space_convert = 1e-20
    time_convert = 1e-12
    df = pd.read_csv(filename, delimiter=' ', skiprows=2)
    #read first two lines
    with open(filename, 'r') as f:
        _, line2 = f.readline(), f.readline()
    cols = line2.lstrip('#').strip().split(' ')
    df.columns = cols
    df   = df[df.TimeStep>=tlim]
    MSDx = df['c_rmsd[1]'].tolist()
    MSDy = df['c_rmsd[2]'].tolist()
    MSDz = df['c_rmsd[3]'].tolist()
    T = df['TimeStep']*dt
    dT = np.max(T)-np.min(T)

    # fit
    T_fit, msdfit_x, slope_x, slope_x_err = straight_fit(T, MSDx, np.min(T), np.max(T))
    T_fit, msdfit_y, slope_y, slope_y_err = straight_fit(T, MSDy, np.min(T), np.max(T))
    T_fit, msdfit_z, slope_z, slope_z_err = straight_fit(T, MSDz, np.min(T), np.max(T))

    # diffusion coefficients
    diff_x = (slope_x/2)*(space_convert/time_convert)
    diff_y = (slope_y/2)*(space_convert/time_convert)
    diff_z = (slope_z/2)*(space_convert/time_convert)
    # errors
    diff_x_err = (slope_x_err/2)*(space_convert/time_convert)
    diff_y_err = (slope_y_err/2)*(space_convert/time_convert)
    diff_z_err = (slope_z_err/2)*(space_convert/time_convert)

    diff_ave = (diff_x+diff_y)/2
    diff_ave_err = error_avg(diff_x_err, diff_y_err)

    print 'Average 2D diffusion coefficient xy:', diff_ave

    fig1 = plt.figure(figsize=(9,7)) 
    ax1  = fig1.add_axes([0.1,0.15,0.8,0.75])
    ax1.plot(T, MSDx, label = 'x')
    ax1.plot(T_fit, msdfit_x, linestyle='dashed', label = 'x fit')
    ax1.plot(T, MSDy, label = 'y')
    ax1.plot(T_fit, msdfit_y, linestyle='dashed', label = 'y fit')
    #ax1.plot(T, MSDz, label = 'z')
    #ax1.plot(T_fit, msdfit_z, linestyle='dashed', label = 'z fit')
    ax1.set_xlabel('t (ps)')
    ax1.set_ylabel('MSD (\AA$^2$)')
    #ax1.set_xlim(0,5000)
    ax1.legend()
    fig1.savefig('PLOTS_C/MSD_{}.pdf'.format(file),bbox_inches='tight')
    fig1.clear()

    return diff_ave, diff_ave_err

def straight_fit(x, y, xmin, xmax):
    params , cov = curve_fit(f1, np.array(x), np.array(y), bounds=(0, np.inf))
    slope, inter = params[0], params[1]
    slope_err = np.sqrt(np.diag(cov))[0]
    xdat = np.linspace(xmin, xmax, 100)
    fit = []
    for x in xdat:
        fit.append(slope*x+inter)
    return xdat, fit, slope, slope_err

def f1(x, A, B):
    return A*x+B

def error_avg(err1, err2):
    return 0.5*np.sqrt(err1**2+err2**2)

def interpolate_derivative(xdat, ydat):
    # make data less exact/more crude
    spl = UnivariateSpline(xdat, ydat, k=4,s=0)
    roots = spl.derivative().roots()
    yfit = spl(xdat)
    return yfit, roots