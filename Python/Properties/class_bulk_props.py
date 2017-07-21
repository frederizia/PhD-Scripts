#!/usr/bin/env python
from __future__ import division
import numpy as np
import csv
import pandas as pd
import itertools



def read_log(f,m,T,P,idx,rhos):
    '''Code to read in log file'''

    try:
        if f == 'LJ_channel':
            filename = '%s/log.%s_T%s_z%s_eps%s_rhos%s'%(f,m,T,P, idx, rhos)
        else:    
            filename = '%s/log.%s_T%s_P%s_%s'%(f,m,T,P,idx)
        f = open(filename,'r')
    except:
        filename = '%s/log.%s_T%s_z%s_eps%s_%s'%(f,m,T,P,idx,rhos)
        f = open(filename,'r')

    data = f.read()
    data_lines = data.split('\n')

    DATA = []
    flag=0
    for i in range(len(data_lines)-1):
        if data_lines[i].split() != [] and data_lines[i].split()[0] == 'unfix' and (data_lines[i].split()[1] == 'RESCALE' or data_lines[i].split()[1] == 'all_rescale'):
            flag = 1
        if flag == 1 and data_lines[i].split() != []:
            DATA.append(data_lines[i].split())

    return DATA

def read_msd(f,m,T,P, idx, rhos):
    '''Code to read in msd file'''
    try:
        if rhos != 'None':
            filename = '%s/msd.%s_T%s_z%s_eps%s_rhos%s'%(f,m,T,P, idx, rhos)
        else:    
            filename = '%s/msd.%s_T%s_P%s_%s'%(f,m,T,P,idx)
        f = open(filename,'r')
    except:
        filename = '%s/msd.%s_T%s_z%s_eps%s_%s'%(f,m,T,P,idx,rhos)
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

def read_dens(f,m,T,P, eps, rhos):
    '''Code to read in density file'''

    try:
        if rhos != 'None':
            filename = '%s/dens.%s_T%s_z%s_eps%s_rhos%s'%(f,m,T,P, eps, rhos)
        else:    
            filename = '%s/dens.%s_T%s_P%s_%s'%(f,m,T,P,idx)
        f = open(filename,'r')
    except:
        filename = '%s/dens.%s_T%s_z%s_eps%s_%s'%(f,m,T,P,eps,rhos)
        f = open(filename,'r')
    data = f.read()
    data_lines = data.split('\n')

    rho = data_lines[-2].split()[1]
    return rho

def read_densprof(f,m,T,P, eps, rhos):
    '''Code to read in density  data from a 1d LAMMPS output file'''

    try:
        if rhos != 'None':
            filename = '%s/densprof.%s_T%s_z%s_eps%s_rhos%s'%(f,m,T,P, eps, rhos)
        else:    
            filename = '%s/densprof.%s_T%s_P%s_%s'%(f,m,T,P,idx)
        f = open(filename,'r')
    except:
        filename = '%s/densprof.%s_T%s_z%s_eps%s_%s'%(f,m,T,P,eps,rhos)
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

    zcoord = zcoord[0]
    dens   = np.mean(np.array(dens), axis=0)

    return zcoord, dens

def read_nist(T,f):
    '''Code to read in nist data file'''
    diff = []

    filename = '%s/nist_data_T%s.txt'%(f,T)

    df = pd.read_csv(filename, delimiter='\t')

    rho = df['Density (kg/m3)']/1000
    rho = rho.tolist()
    press = df['Pressure (bar)'].tolist()
    shear = df['Viscosity (Pa*s)']*1000
    shear = shear.tolist()


    return rho, press, shear, diff

def press_fluct(p):
    filename = 'thermo_p{}.dat'.format(p)

    df = pd.read_csv(filename, delimiter=' ', skiprows=2)
    #read first two lines
    with open(filename, 'r') as f:
        _, line2 = f.readline(), f.readline()
    cols = line2.lstrip('#').strip().split(' ')
    df.columns = cols
    time = df['TimeStep']
    time = time.tolist()
    vol  = df['v_vol_one']*1e-30
    vol  = vol.tolist()
    temp = df['v_tem_one']
    temp = temp.tolist()

    temp_list = []
    vol_list = []
    for t,v,temp in itertools.izip(time, vol, temp):
        if t > 70000:
            vol_list.append(v)
            temp_list.append(temp)


    vol_mean   = np.mean(np.array(vol_list))
    vol_var    = np.var(np.array(vol_list))
    temp_mean  = np.mean(np.array(temp_list))


    return temp_mean, vol_mean, vol_var

def mean_vals(p):
    filename = 'thermo_p{}.dat'.format(p)
    #try:
        #filename = '%s/visc.%s_T%s_P%s_%s'%(f,m,T,P,idx)
        #f = open(filename,'r')
    #except:
        #filename = '%s/visc.%s_T%s_z%s_eps%s'%(f,m,T,P,idx)
        #f = open(filename,'r')

    df = pd.read_csv(filename, delimiter=' ', skiprows=2)
    #read first two lines
    with open(filename, 'r') as f:
        _, line2 = f.readline(), f.readline()
    cols = line2.lstrip('#').strip().split(' ')
    df.columns = cols
    vol  = df['v_vol_one']*1e-30
    vol  = vol.tolist()
    press = df['v_pre_one']*1e5
    press = press.tolist()

    vol_mean   = np.mean(np.array(vol))
    press_mean  = np.mean(np.array(press))

    return vol_mean, press_mean

def fluid_vol(f,m,T,P,eps,rhos):
    '''Code to read in msd file'''
    try:
        filename = '%s/data.%s_T%s_z%s_eps%s_%s'%(f,m,T,P,eps,rhos)
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
        if m=='spce':
            blank_idx = 6
        else:
            blank_idx = 4
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
        self.data = read_log(self.f,self.m,self.T,self.P, self.idx,'None')
        self.msd = read_msd(self.f,self.m,self.T,self.P, self.idx,'None')

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
        sv = 0
        for i in range(len(self.data)-1):
            try:
                if self.data[i][1] == 'shear' and self.data[i][2] == 'viscosity':
                    sv = 1000*float(self.data[i][4])
            except:
                continue
        return sv

    def shear2(self):
        try:
            filename = '%s/visc.%s_T%s_P%s_%s'%(self.f,self.m,self.T,self.P, self.idx)
            #f = open(filename,'r')
        except:
            filename = '%s/visc.%s_T%s_z%s_eps%s'%(self.f,self.m,self.T,self.z,self.idx)
            #f = open(filename,'r')

        df = pd.read_csv(filename, delimiter=' ', skiprows=2)
        #read first two lines
        with open(filename, 'r') as f:
            _, line2 = f.readline(), f.readline()
        cols = line2.lstrip('#').strip().split(' ')
        df.columns = cols
        shear = df['v_etas']/1000
        shear = shear.tolist()
        time = df['TimeStep']
        time = time.tolist()

        shear_list = []
        count = 0
        if self.m == 'tip4p' or self.m == 'TraPPErigid' or self.m == 'EPM2rigid' or self.m == 'EPM2flex':
            tlim = 4000000
        else:
            tlim = 8000000
        for t,s in itertools.izip(time, shear):
            if t > tlim:
                shear_list.append(1e6*s)
                count += 1
        shear_val = np.mean(np.array(shear_list))
        shear_err = np.std(np.array(shear_list))/np.sqrt(len(shear_list))

        return shear_val, shear_err

    def bulk(self):
        bv = 0
        for i in range(len(self.data)-1):
            try:
                if self.data[i][1] == 'bulk' and self.data[i][2] == 'nvt' and self.data[i][3] == 'viscosity':
                    bv = 1000*float(self.data[i][5])
            except:
                continue
        return bv

    def bulk2(self):
        try:
            filename = '%s/visc.%s_T%s_P%s_%s'%(self.f,self.m,self.T,self.P, self.idx)
            #f = open(filename,'r')
        except:
            filename = '%s/visc.%s_T%s_z%s_eps%s_%s'%(self.f,self.m,self.T,self.P,self.idx)
            #f = open(filename,'r')

        df = pd.read_csv(filename, delimiter=' ', skiprows=2)
        #read first two lines
        with open(filename, 'r') as f:
            _, line2 = f.readline(), f.readline()
        cols = line2.lstrip('#').strip().split(' ')
        df.columns = cols
        bulk = df['v_etab_nvt']/1000
        bulk = bulk.tolist()
        time = df['TimeStep']
        time = time.tolist()

        bulk_list = []
        count = 0
        if self.m == 'tip4p' or self.m == 'TraPPErigid' or self.m == 'EPM2rigid' or self.m == 'EPM2flex':
            tlim = 4000000
        else:
            tlim = 8000000
        for t,s in itertools.izip(time, bulk):
            if t > tlim:
                bulk_list.append(1e6*s)
                count += 1
        bulk_val = np.mean(np.array(bulk_list))
        bulk_err = np.std(np.array(bulk_list))/np.sqrt(len(bulk_list))

        return bulk_val, bulk_err


    def visc_ratio(self):
        ratio = self.bulk()/self.shear()
        return ratio

    def diff(self,tstamp, dt):
        data = self.msd
        if self.m == 'tip4p' or self.m == 'TraPPErigid' or self.m == 'EPM2rigid' or self.m == 'EPM2flex':
            tstamp = 2000000
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

    def __init__(self,f,m,T,z,eps,rhos):
        self.f = f
        self.m = m
        self.T = T
        self.z = z
        self.eps = eps
        self.rhos = rhos
        self.data = read_log(self.f,self.m,self.T,self.z, self.eps, self.rhos)
        self.msd = read_msd(self.f,self.m,self.T,self.z,self.eps, self.rhos)
        #self.dens = read_dens(self.f,self.m,self.T,self.z,self.eps, self.rhos)
        self.densprof = read_densprof(self.f,self.m,self.T,self.z,self.eps, self.rhos)
    def press(self):
        for i in range(len(self.data)-1):
            if self.data[i+1][0] == 'Loop':
                P = float(self.data[i][7])

        return P

    def temp(self):
        for i in range(len(self.data)-1):
            if self.data[i+1][0] == 'Loop':
                T = float(self.data[i][5])

        return T

    def rho(self):
        #V, cnt, sep = fluid_vol(self.f,self.m,self.T,self.z)
        #mass = 18.0153/6.0221415e23
        #r = mass*cnt/(V*10**(-24)) # ?/cm^3
        try:
            r = read_dens(self.f,self.m,self.T,self.z,self.eps, self.rhos)
        except:
            r = np.mean(np.array(self.densprof[1]))

        return r

    def shear(self):
        for i in range(len(self.data)-1):
            try:
                if self.data[i][1] == 'shear' and self.data[i][2] == 'viscosity':
                    sv = 1000*float(self.data[i][4])
            except:
                continue
        return sv

    def shear3(self, T, tstamp, dt):
        Ds = self.diff(tstamp, dt)[2]/1e9
        alpha = 1.7*10**(-10)
        kB = 1.38*10**(-23)
        sv = (kB*T)/(3*np.pi*alpha*Ds)

        return sv

    def shear2(self):
        try:
            filename = '%s/visc.%s_T%s_P%s_%s'%(self.f,self.m,self.T,self.P, self.idx)
            #f = open(filename,'r')
        except:
            filename = '%s/visc.%s_T%s_z%s_eps%s_%s'%(self.f,self.m,self.T,self.z,self.eps, self.rhos)
            #f = open(filename,'r')

        df = pd.read_csv(filename, delimiter=' ', skiprows=2)
        #read first two lines
        with open(filename, 'r') as f:
            _, line2 = f.readline(), f.readline()
        cols = line2.lstrip('#').strip().split(' ')
        df.columns = cols
        shear = df['v_etas']/1000
        shear = shear.tolist()
        time = df['TimeStep']
        time = time.tolist()

        shear_list = []
        count = 0
        tlim = 10000000
        for t,s in itertools.izip(time, shear):
            if t > tlim:
                shear_list.append(1e6*s)
                count += 1
        shear_val = np.mean(np.array(shear_list))
        shear_val *= (float(self.z)-1.0)/(float(self.z)+1.5)
        shear_err = np.std(np.array(shear_list))/np.sqrt(len(shear_list))

        return shear_val, shear_err
    

    def bulk(self):
        bv = 0
        for i in range(len(self.data)-1):
            try:
                if self.data[i][1] == 'bulk' and self.data[i][2] == 'nvt' and self.data[i][3] == 'viscosity':
                    bv = 1000*float(self.data[i][5])
            except:
                continue
        return bv

    def bulk2(self):
        try:
            filename = '%s/visc.%s_T%s_P%s_%s'%(self.f,self.m,self.T,self.P, self.idx)
            #f = open(filename,'r')
        except:
            filename = '%s/visc.%s_T%s_z%s_eps%s_%s'%(self.f,self.m,self.T,self.z,self.eps, self.rhos)
            #f = open(filename,'r')

        df = pd.read_csv(filename, delimiter=' ', skiprows=2)
        #read first two lines
        with open(filename, 'r') as f:
            _, line2 = f.readline(), f.readline()
        cols = line2.lstrip('#').strip().split(' ')
        df.columns = cols
        bulk = df['v_etab_nvt']/1000
        bulk = bulk.tolist()
        time = df['TimeStep']
        time = time.tolist()

        bulk_list = []
        count = 0
        tlim = 10000000
        for t,s in itertools.izip(time, bulk):
            if t > tlim:
                bulk_list.append(1e6*s)
                count += 1
        bulk_val = np.mean(np.array(bulk_list))
        bulk_val *= (float(self.z)+1.5)/(float(self.z)-1.0)
        bulk_err = np.std(np.array(bulk_list))/np.sqrt(len(bulk_list))

        return bulk_val, bulk_err


    def fric(self, coord):
        try:
            filename = '%s/fric.%s_T%s_z%s_eps%s_%s'%(self.f,self.m,self.T,self.z, self.eps, self.rhos)
            if self.f == 'LJ_channel':
                filename = '%s/fric.%s_T%s_z%s_eps%s_rhos%s'%(self.f,self.m,self.T,self.z, self.eps, self.rhos)
            df = pd.read_csv(filename, delimiter=' ', skiprows=2)
        except:
            filename = '%s/visc.%s_T%s_z%s_eps%s'%(self.f,self.m,self.T,self.z, self.eps)
            df = pd.read_csv(filename, delimiter=' ', skiprows=2)

        
        #read first two lines
        with open(filename, 'r') as f:
            _, line2 = f.readline(), f.readline()
        cols = line2.lstrip('#').strip().split(' ')
        df.columns = cols
        fric = df['v_ffin_%s'%coord]
        fric = fric.tolist()
        time = df['TimeStep']
        time = time.tolist()

        fric_list = []
        count = 0
        for t,f in itertools.izip(time, fric):
            if t > 800000:
                fric_list.append(f/10000)
                count += 1
        fric_val = np.mean(np.array(fric_list))
        fric_err = np.std(np.array(fric_list))/np.sqrt(len(fric_list))

        return fric_val

    def visc_ratio(self):
        ratio = self.bulk()/self.shear()
        return ratio

    def diff(self,tstamp, dt):
        data = self.msd
        idx = list(data[0]).index(tstamp)
        convert = 10**(-5) # A2m^2/(2*timestep*ps2s)

        delta = []
        diff  = []
        for i in [1,2,3,4]:
            delta.append(data[i][-1]-data[i][idx])
            diff.append(1e9*convert*(data[i][-1]-data[i][idx])/dt)
        diff[3] = diff[3]/3 # As this is averaging over all directions

        return diff

    def separation(self):
        #sep2 = fluid_vol(self.f,self.m,self.T,self.z, self.eps)[2]
        sep = self.wall_pos()
        
        return sep

    def wa(self):

        filename = '%s/wa.%s_T%s_z%s_eps%s_%s'%(self.f,self.m,self.T,self.z, self.eps, self.rhos)
        if self.f == 'LJ_channel':
            filename = '%s/wa.%s_T%s_z%s_eps%s_rhos%s'%(self.f,self.m,self.T,self.z, self.eps, self.rhos)
            #f = open(filename,'r')

        df = pd.read_csv(filename, delimiter=' ', skiprows=2)
        #read first two lines
        with open(filename, 'r') as f:
            _, line2 = f.readline(), f.readline()
        cols = line2.lstrip('#').strip().split(' ')
        df.columns = cols
        wa = df['v_WA_layer']
        wa = wa.tolist()
        time = df['TimeStep']
        time = time.tolist()

        for i in range(len(self.data)-1):
            try:
                
                if self.data[i][0] == 'WA' and self.data[i][1] == 'layer':
                    wa_tmp = float(self.data[i][3])
            except:
                continue
        factor = wa[-1]/wa_tmp

        wa_list = []
        count = 0
        for t,w in itertools.izip(time, wa):
            if t > 800000:
                wa_list.append(w/factor)
                count += 1
        wa_val = np.mean(np.array(wa_list))
        wa_err = np.std(np.array(wa_list))/np.sqrt(len(wa_list))

        return wa_val

    def profile(self):
        z, rho = self.densprof
        return z, rho

    def wall_pos(self):
        Y, RHO = self.densprof
        lower = 0
        lower_idx = 0
        upper = len(RHO)
        upper_idx = len(RHO)

        ref = round((np.average(np.array(RHO[:10]))+np.average(np.array(RHO[-10:])))/2,2)

        for i in range(len(RHO)):
            if round(RHO[i],2) == ref and round(RHO[i+1],2) > ref:
                lower_idx = i
            if round(RHO[i],2) > ref and round(RHO[i],2) > round(RHO[i+1],2):
                lowerpeak_idx = i
                break
        print 'Lower peak:', Y[lowerpeak_idx]
        LOWER = lower_idx+int((lowerpeak_idx-lower_idx)/4)
        ylower = Y[LOWER]
        
        count = 0
        for j in list(reversed(range(len(RHO)))):
            if round(RHO[j],2) == ref and round(RHO[j-1],2) > ref:
                upper_idx = j
                count +=1
            if round(RHO[j],2) > ref and round(RHO[j],2) > round(RHO[j-1],2) and count > 0:
                upperpeak_idx = j
                break

        print 'Upper peak:', Y[upperpeak_idx]
        UPPER = upper_idx-int((upper_idx-upperpeak_idx)/4)
        yupper = Y[UPPER]
        

        sep = yupper-ylower
        print 'Separation:', sep

        return sep