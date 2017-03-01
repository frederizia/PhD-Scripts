#!/usr/bin/env python
from __future__ import division
import numpy as np
import csv



def read_log(f,m,T,P,idx):
    '''Code to read in log file'''

    try:
        filename = '%s/log.%s_T%s_P%s_%s'%(f,m,T,P,idx)
        f = open(filename,'r')
    except:
        filename = '%s/log.%s_T%s_z%s_eps%s'%(f,m,T,P,idx)
        f = open(filename,'r')
    data = f.read()
    data_lines = data.split('\n')

    DATA = []
    flag=0
    for i in range(len(data_lines)-1):
        if data_lines[i].split() != [] and data_lines[i].split()[0] == 'unfix' and data_lines[i].split()[1] == 'RESCALE':
            flag = 1
        if flag == 1 and data_lines[i].split() != []:
            DATA.append(data_lines[i].split())

    return DATA

def read_msd(f,m,T,P, idx):
    '''Code to read in msd file'''
    try:
        filename = '%s/msd.%s_T%s_P%s_%s'%(f,m,T,P, idx)
        f = open(filename,'r')
    except:
        filename = '%s/msd.%s_T%s_z%s_eps%s'%(f,m,T,P,idx)
        f = open(filename,'r')
    data = f.read()
    data_lines = data.split('\n')

    ts = []
    msd_x = []
    msd_y = []
    msd_z = []
    msd_tot = []
    flag=0
    for i in range(2,len(data_lines)-1):
        ts.append(float(data_lines[i].split()[0]))
        msd_x.append(float(data_lines[i].split()[1]))
        msd_y.append(float(data_lines[i].split()[2]))
        msd_z.append(float(data_lines[i].split()[3]))
        msd_tot.append(float(data_lines[i].split()[4]))

    ts, msd_x, msd_y, msd_z, msd_tot = np.array(ts), np.array(msd_x), np.array(msd_y), np.array(msd_z), np.array(msd_tot)

    return ts, msd_x, msd_y, msd_z, msd_tot

def read_dens(f,m,T,P, eps):
    '''Code to read in msd file'''
    try:
        filename = '%s/dens.%s_T%s_P%s'%(f,m,T,P)
        f = open(filename,'r')
    except:
        filename = '%s/dens.%s_T%s_z%s_eps%s'%(f,m,T,P,eps)
        f = open(filename,'r')
    data = f.read()
    data_lines = data.split('\n')

    rho = data_lines[-2].split()[1]
    return rho

def read_nist(T,f):
    '''Code to read in nist data file'''

    filename = '%s/nist_data_T%s.csv'%(f,T)

    f = open(filename,'r')
    reader = csv.DictReader(f)

    rho = []
    press = []
    shear = []
    diff = []

    for row in reader:
        rho.append(float(row['Rho']))
        press.append(float(row['P']))
        shear.append(float(row['Shear']))
        diff.append(float(row['Diff']))

    return rho, press, shear, diff

def fluid_vol(f,m,T,P,eps):
    '''Code to read in msd file'''
    try:
        filename = '%s/data.%s_T%s_z%s_eps%s'%(f,m,T,P,eps)
        file = open(filename,'r')
    except:
        filename = '%s/data.%s_T%s_P%s'%(f,m,T,P)
        file = open(filename,'r')

    data = file.read()
    data_lines = data.split('\n')
    blank_count = 0

    zlo = 15
    zhi = 5

    if f == 'Water_Graphene':
        indx, indy = 11, 12
        blank_idx = 6
    else:
        indx, indy = 5, 6
        blank_idx = 5

    xlo, xhi = float(data_lines[indx].split()[0]), float(data_lines[indx].split()[1])
    ylo, yhi = float(data_lines[indy].split()[0]), float(data_lines[indy].split()[1])

    count = 0
    for i in range(0,len(data_lines)-1):
        if blank_count == blank_idx and len(data_lines[i].split()) > 1:
            if data_lines[i].split()[2]== '1':
                z_coord = float(data_lines[i].split()[6])
                count += 1
                if z_coord > zhi:
                    zhi = z_coord
                if z_coord < zlo:
                    zlo = z_coord
        if data_lines[i].split() == []:
            blank_count +=1

    # account for sizes
    zhi -= 1.5
    zlo += 1.5
    vol = (xhi-xlo)*(yhi-ylo)*(zhi-zlo)
    fluid_count = count/3
    sep = zhi-zlo

    return vol, fluid_count, sep

# create class

class bulk_properties:

    def __init__(self,f,m,T,P,idx):
        self.f = f
        self.m = m
        self.T = T
        self.P = P
        self.idx = idx
        self.data = read_log(self.f,self.m,self.T,self.P, self.idx)
        self.msd = read_msd(self.f,self.m,self.T,self.P, self.idx)

    def press(self):
        for i in range(len(self.data)-1):
            if self.data[i+1][0] == 'Loop':
                P = float(self.data[i][6])

        return P

    def temp(self):
        for i in range(len(self.data)-1):
            if self.data[i+1][0] == 'Loop':
                T = float(self.data[i][4])

        return T

    def rho(self):
        for i in range(len(self.data)-1):
            if self.data[i+1][0] == 'Loop':
                r = float(self.data[i][2])

        return r

    def shear(self):
        for i in range(len(self.data)-1):
            try:
                if self.data[i][1] == 'shear' and self.data[i][2] == 'viscosity':
                    sv = 1000*float(self.data[i][4])
            except:
                continue
        return sv

    def bulk(self):
        for i in range(len(self.data)-1):
            try:
                if self.data[i][1] == 'bulk' and self.data[i][2] == 'nvt' and self.data[i][3] == 'viscosity':
                    bv = 1000*float(self.data[i][5])
            except:
                continue
        return bv

    def visc_ratio(self):
        ratio = self.bulk()/self.shear()
        return ratio

    def diff(self,tstamp, dt):
        data = self.msd
        idx = list(data[0]).index(tstamp)
        convert = 10**(-5)

        delta = []
        diff  = []
        for i in [1,2,3,4]:
            delta.append(data[i][-1]-data[i][idx])
            diff.append(convert*(data[i][-1]-data[i][idx])/dt)
        diff[3] = 1e9*diff[3]/3 # As this is averaging over all directions

        return diff

class confined_properties:

    def __init__(self,f,m,T,z,eps):
        self.f = f
        self.m = m
        self.T = T
        self.z = z
        self.eps = eps
        self.data = read_log(self.f,self.m,self.T,self.z, self.eps)
        self.msd = read_msd(self.f,self.m,self.T,self.z,self.eps)
        self.dens = read_dens(self.f,self.m,self.T,self.z,self.eps)

    def press(self):
        for i in range(len(self.data)-1):
            if self.data[i+1][0] == 'Loop':
                P = self.data[i][7]

        return P

    def temp(self):
        for i in range(len(self.data)-1):
            if self.data[i+1][0] == 'Loop':
                T = self.data[i][5]

        return T

    def rho(self):
        #V, cnt, sep = fluid_vol(self.f,self.m,self.T,self.z)
        #mass = 18.0153/6.0221415e23
        #r = mass*cnt/(V*10**(-24)) # ?/cm^3
        r = self.dens
        return r

    def shear(self):
        for i in range(len(self.data)-1):
            try:
                if self.data[i][1] == 'shear' and self.data[i][2] == 'viscosity':
                    sv = 1000*float(self.data[i][4])
            except:
                continue
        return sv

    def shear2(self, T, tstamp, dt):
        Ds = self.diff(tstamp, dt)[2]/1e9
        alpha = 1.7*10**(-10)
        kB = 1.38*10**(-23)
        sv = (kB*T)/(3*np.pi*alpha*Ds)

        return sv

    def bulk(self):
        for i in range(len(self.data)-1):
            try:
                if self.data[i][1] == 'bulk' and self.data[i][2] == 'nvt' and self.data[i][3] == 'viscosity':
                    bv = 1000*float(self.data[i][5])
            except:
                continue
        return bv

    def visc_ratio(self):
        ratio = self.bulk()/self.shear()
        return ratio

    def diff(self,tstamp, dt):
        data = self.msd
        idx = list(data[0]).index(tstamp)
        convert = 10**(-5)

        delta = []
        diff  = []
        for i in [1,2,3,4]:
            delta.append(data[i][-1]-data[i][idx])
            diff.append(1e9*convert*(data[i][-1]-data[i][idx])/dt)
        diff[3] = diff[3]/3 # As this is averaging over all directions

        return diff

    def separation(self):
        sep = fluid_vol(self.f,self.m,self.T,self.z, self.eps)[2]
        return sep

    def wa(self):
        for i in range(len(self.data)-1):
            try:
                if self.data[i][0] == 'WA' and self.data[i][1] == 'layer':
                    wa_val = 0.5*float(self.data[i][3])
            except:
                continue
        return wa_val
