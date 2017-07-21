'''Tools for GO calculations'''
from __future__ import division
import numpy as np
import sys, string
import re
import itertools
from scipy import stats
from scipy.optimize import curve_fit

class FLOW:
    def __init__(self, n,o,dx,f):
        self.n = n
        self.o = o
        self.dx = dx
        self.f = f

    def read_log(self):
        '''Code to read in log file'''

        filename = 'DATA/log.n{}_o{}_delx{}_F{}'.format(self.n,self.o,self.dx,self.f)
        f = open(filename,'r')

        data = f.read()
        data_lines = data.split('\n')

        DATA = []
        flag=0
        for i in range(len(data_lines)-1):
            if data_lines[i].split() != [] and data_lines[i].split()[0] == 'run' and data_lines[i].split()[1] == '${t_run}':
                flag = 1
            if flag == 1 and data_lines[i].split() != []:
                DATA.append(data_lines[i].split())

        return DATA

    def coords(self):
        filename = 'DATA/log.n{}_o{}_delx{}_F{}'.format(self.n,self.o,self.dx,self.f)
        f = open(filename,'r')

        data = f.read()
        data_lines = data.split('\n')

        coord_line = re.sub('[\(\)]', '', data_lines[3]).split('=')[1].split('to')
        xlo, ylo, zlo = coord_line[0].strip(' ').split(' ')
        xhi, yhi, zhi = coord_line[1].strip(' ').split(' ')

        return float(xlo), float(xhi), float(ylo), float(yhi), float(zlo), float(zhi)

    def geometry(self):
        xlo, xhi, ylo, yhi, zlo, zhi = self.coords()
        area_xy = (xhi-xlo)*(yhi-ylo)
        area_flow = (yhi-ylo)*float(self.dx)
        volume = (xhi-xlo)*(yhi-ylo)*(zhi-zlo)
        return area_xy, area_flow, volume

    def stress_prof(self):
        filename = 'DATA/stress.n{}_o{}_delx{}_F{}'.format(self.n,self.o,self.dx,self.f)
        f = open(filename,'r')
        data_lines = f.readlines()
        f.close()

        area = self.geometry()[0]

        count_tmp = 0
        dz = 0
        for line in data_lines:
            if line[0] != '#' and len(line.split()) != 2 and len(line.split()) != 3: 
                dz = float(line.split()[1]) - dz
                count_tmp +=1
                if count_tmp ==2:
                    break
        vol = area*dz

        coords, pressure = [], []
        coords_tot, pressure_tot = [], []
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
                    coords, pressure = [], []
                

            elif line[0] != '#' and len(items) != 2 and len(items) != 3:
                Coord = float(items[1])
                Ncount = float(items[2])
                xx = float(items[3])*Ncount/vol
                yy = float(items[4])*Ncount/vol
                zz = float(items[5])*Ncount/vol 
                Pxx = -xx
                Pyy = -yy
                Pzz = -zz
                pres = (Pxx+Pyy+Pzz)*0.101325/3

                coords.append(Coord)
                pressure.append(pres)

        coords_tot.append(coords)
        #print coords_tot
        pressure_tot.append(pressure)
        coords = np.mean(np.array(coords_tot),axis=0)
        pressure = np.mean(np.array(pressure_tot),axis=0)
        zmax = np.max(coords)

        # Pleft and Pright
        pleft, pright = [], []
        for c,p in itertools.izip(coords,pressure):
            if 15 < c < 35:
                pleft.append(p)
            if zmax-20 < c < zmax:
                pright.append(p)
        Pleft = np.mean(np.array(pleft))
        Pright = np.mean(np.array(pright))
        deltaP = abs(Pright-Pleft)

        print 'The pressure drop for F =', self.f, 'is', deltaP, 'MPa'

                
        return coords, pressure, deltaP

    def density_prof(self, xlo):

        filename = filename = 'DATA/dens.n{}_o{}_delx{}_F{}'.format(self.n,self.o,self.dx,self.f)
        
        # read in data
        f = open(filename,'r')
        data = f.read()        
        data_lines = data.split('\n')
        
        
        xbin = []
        zbin = []

        if self.n == '1':
            zhi = 99.75
            if self.dx == '6':
                xhi = 17.1
            else:
                xhi = 17.1 #21.3
            if xlo == 'None':
                xlist = list(np.arange(xhi,xhi+float(self.dx)+0.2,0.2))
            else:
                xlist = list(np.arange(xlo,xhi+float(self.dx)+0.2,0.2))
        else:
            zhi = 109.9
            if xlo == 'None':
                xlist = list(np.arange(17.1,17.1+float(self.dx)+0.2,0.2))
                xlist.extend(list(np.arange(5.1,5.1+float(self.dx)+0.2,0.2)))
            else:
                xlist = list(np.arange(xlo,17.1+float(self.dx),0.2))
        zlist = list(np.arange(0.1, zhi, 0.2))
        #zlist_mem = list(np.arange(48.75, 51.25, 0.5))

        # round values in list so that they agree with data
        xlist = [round(n, 1) for n in xlist]
        zlist = [round(n, 1) for n in zlist]
        #zlist_mem = [round(n, 2) for n in zlist_mem]


        xdim = len(xlist)
        zdim = len(zlist)
        #zdim_mem = len(zlist_mem)
        density  = np.zeros((xdim,zdim))
        velocity = np.zeros((xdim,zdim))
        #density_mem = np.zeros((xdim,zdim_mem))
        coords_x = np.zeros((xdim,zdim))
        coords_z = np.zeros((xdim,zdim))


        count = 1
        count_tmp = 0
        count_tmp2 = 0
        for j in range(4,len(data_lines)-1):

            if len(data_lines[j].split()) != 2 and len(data_lines[j].split()) != 3:
                x_val = float(data_lines[j].split()[1])
                z_val = float(data_lines[j].split()[2])
                den_val = float(data_lines[j].split()[4])
                vel_val = float(data_lines[j].split()[7])

                # only select values within slit
                if x_val in xlist and z_val in zlist:
                    x_ind = xlist.index(x_val)
                    z_ind = zlist.index(z_val)
                    #print z_ind, z_val
                    xbin.append(x_val)
                    zbin.append(z_val)

                    # average over all time steps
                    if count_tmp == 0:
                        density[x_ind][z_ind] = den_val
                        velocity[x_ind][z_ind] = vel_val
                    else:
                        density[x_ind][z_ind] = (density[x_ind][z_ind]+den_val)/2
                        velocity[x_ind][z_ind] = (velocity[x_ind][z_ind]+vel_val)/2
                    coords_x[x_ind][z_ind] = xlist[x_ind]
                    coords_z[x_ind][z_ind] = zlist[z_ind]
                    count_tmp += 1
            count =2
      

        density_z = np.average(density, axis = 0)
        density_z_err = stats.sem(density, axis = 0)
        coords_z = coords_z[0,:] #np.average(coords_z, axis = 0)
        density_x = np.average(density, axis = 1)
        velocity_z = np.average(velocity, axis = 0)
        velocity_z_err = stats.sem(velocity, axis = 0)
        coords_x = coords_x[:,0] #np.average(coords_x, axis = 1)

        #print coords_z

        rho_mean = np.mean(density)
        zmax = np.max(coords_z)

        # Rholeft and Rhoright
        rholeft, rhoright = [], []
        for c,r in itertools.izip(coords_z,density_z):
            if 15 < c < 35:
                rholeft.append(r)
            if zmax-20 < c < zmax:
                rhoright.append(r)
        Rholeft = np.mean(np.array(rholeft))
        Rhoright = np.mean(np.array(rhoright))
        deltaRho = abs(Rhoright-Rholeft)


        print 'The density drop for F =', self.f, 'is', deltaRho, 'g/cm^3'

        tol = 0.03
        tmp = 0
        for k in range(25, 250, 1):
            if density_z[k] > 1.10:    
                L_left = coords_z[k]-4
                break
        for k in range(25, 250, 1):
            if density_z[::-1][k] > 1.08:
                L_right = coords_z[::-1][k]+4
                break

    
        L_diff = L_right - L_left
        print 'Membrane width:', L_diff


        return coords_z, density_z, coords_x, density_x, deltaRho, velocity_z, L_diff, rho_mean, density_z_err, velocity_z_err

    def flux(self):
        mol_mass_water = 18.01
        convert = 10**8 # in 10^3 mol/m^2s
        if self.n == '1':
            data = self.density_prof(0.25)
        else:
            data = self.density_prof(0.1)
        rho = data[1]
        vel = data[5]
        rho_err = data[8]
        vel_err = data[9]


        Jz = (1/rho.shape[0])*sum(rho*vel)

        Jz_mol = convert*Jz/mol_mass_water 

        Jz_err = err_Jz(rho, rho_err, vel, vel_err)*convert/mol_mass_water

         
        return Jz_mol, Jz_err

def straight_fit(x, y, xmin, xmax):
    slope, cov = curve_fit(f, np.array(x), np.array(y))
    slope_err = np.sqrt(np.diag(cov))
    xdat = np.linspace(xmin, xmax, 100)
    fit = []
    for x in xdat:
        fit.append(slope*x)
    return xdat, fit, slope, slope_err

def poly_fit(x, y, xmin, xmax):
    params, cov = curve_fit(f2, np.array(x), np.array(y))
    #slope_err = np.sqrt(np.diag(cov))
    xdat = np.linspace(xmin, xmax, 100)
    fit = []
    for x in xdat:
        fit.append(params[0]*x**2+params[1]*x)
    return xdat, fit, params

def exp_fit(x, y, xmin, xmax):
    params, cov = curve_fit(f3, np.array(x), np.array(y))
    #slope_err = np.sqrt(np.diag(cov))
    xdat = np.linspace(xmin, xmax+1, 100)
    fit = []
    for x in xdat:
        fit.append(f3(x,params[0],params[1]))
    return xdat, fit, params

def exp_fit2(x, y, xmin, xmax):
    params, cov = curve_fit(f5, np.array(x), np.array(y))
    #slope_err = np.sqrt(np.diag(cov))
    xdat = np.linspace(xmin, xmax+1, 100)
    fit = []
    for x in xdat:
        fit.append(f5(x,params[0],params[1],params[2]))
    return xdat, fit, params

def f(x, A):
    return A*x

def f2(x, A, B):
    return A*x**2+B*x

def f3(x,A,B):
    return A*np.exp(B*x)-A

def f4(x, a, b, c):
    return a * np.exp(b * x) + c

def f5(x, a, b, c):
    return a*b**x + c

def perm_units(permeance, width, rho_mean, w_err, p_err):
    width_si = 1e-10*width
    mol2vol = (18.01/rho_mean)*1e-6
    permeance_si = permeance*1e-3*mol2vol
    permeability_si = permeance_si*width_si
    perm_err = permeability_si*np.sqrt((p_err/permeance)**2+(w_err/width)**2)
    return permeability_si, perm_err, permeance_si

def err_Jz (rho_list, rho_err, v_list, v_err):
    # pv_err
    pv_err = []
    for i in range(len(rho_list)):
        pv = rho_list[i]*v_list[i]
        rho_bracket = (rho_err[i]/rho_list[i])**2
        v_bracket = (v_err[i]/v_list[i])**2
        sq_pv = np.sqrt(rho_bracket+v_bracket)
        pv_err.append((pv*sq_pv)**2)

    # Jz_err
    Jz_err = (1/len(rho_list))*np.sqrt(sum(pv_err))
    return Jz_err
