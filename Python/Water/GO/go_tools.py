'''Tools for GO calculations'''
from __future__ import division
import numpy as np
import sys, string
import re
import itertools

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

        for line in data_lines:
            if line[0] != '#' and len(line.split()) != 2: 
                dz = float(line.split()[1])
                break
        vol = area*dz

        coords, pressure = [], []
        for line in data_lines:
            items = line.split()
            if line[0] != '#' and len(items) != 2: 
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

        # Pleft and Pright
        pleft, pright = [], []
        for c,p in itertools.izip(coords,pressure):
            if 15 < c < 35:
                pleft.append(p)
            if 75 < c < 95:
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

        if xlo == 'None':
            xlist = list(np.arange(21.25,21.25+float(self.dx),0.5))
        else:
            xlist = list(np.arange(xlo,21.25+float(self.dx),0.5))
        zlist = list(np.arange(0.25, 99.75, 0.5))
        #zlist_mem = list(np.arange(48.75, 51.25, 0.5))

        # round values in list so that they agree with data
        xlist = [round(n, 2) for n in xlist]
        zlist = [round(n, 2) for n in zlist]
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

            if len(data_lines[j].split()) != 2:
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
                # only select values within slit
                #if x_val in xlist and z_val in zlist_mem:
                #    x_ind = xlist.index(x_val)
                #    z_ind = zlist_mem.index(z_val)
                    # average over all time steps
                #    if count_tmp2 == 0:
                #        density_mem[x_ind][z_ind] = den_val
                #        velocity_mem[x_ind][z_ind] = vel_val
                #    else:
                #        density_mem[x_ind][z_ind] = (density_mem[x_ind][z_ind]+den_val)/2
                #        velocity_mem[x_ind][z_ind] = (velocity_mem[x_ind][z_ind]+vel_val)/2
                #    count_tmp2 += 1
            count =2
      

        density_z = np.average(density, axis = 0)
        coords_z = np.average(coords_z, axis = 0)
        density_x = np.average(density, axis = 1)
        velocity_mem_z = np.average(velocity, axis = 0)
        coords_x = np.average(coords_x, axis = 1)

        # Rholeft and Rhoright
        rholeft, rhoright = [], []
        for c,r in itertools.izip(coords_z,density_z):
            if 15 < c < 35:
                rholeft.append(r)
            if 75 < c < 95:
                rhoright.append(r)
        Rholeft = np.mean(np.array(rholeft))
        Rhoright = np.mean(np.array(rhoright))
        deltaRho = abs(Rhoright-Rholeft)

        print 'Rholeft:', Rholeft
        print 'Rhoright:', Rhoright

        print 'The density drop for F =', self.f, 'is', deltaRho, 'g/cm^3'

        tol = 0.03
        tmp = 0
        for k in range(25, 100, 1):
            if density_z[k] > 1.10:    
                L_left = coords_z[k]-4
                break
        for k in range(25, 100, 1):
            if density_z[::-1][k] > 1.08:
                L_right = coords_z[::-1][k]+4
                break
                
    
    
        L_diff = L_right - L_left
        print 'Membrane width:', L_diff

        vel_max = velocity_mem_z[int(len(velocity_mem_z)/2)]

        return coords_z, density_z, coords_x, density_x, deltaRho, velocity_mem_z, vel_max

    def flux(self):
        mol_mass_water = 18.01
        convert = 10**8 
        data = self.density_prof(0.25)
        rho = data[1]
        vel = data[5]


        Jz = (1/rho.shape[0])*sum(rho*vel)

        Jz_mol = convert*Jz/mol_mass_water
        return Jz_mol
