'''Tools for Poiseuille flow calculations'''
from __future__ import division
import numpy as np
import sys, string
import matplotlib.pyplot as plt
from matplotlib import cm
import re
import itertools
from scipy import stats
from scipy.optimize import curve_fit
import pandas as pd
from scipy.integrate import simps


def blockAverage(data):

    DataLen     = len(data) 
    BlockSize   = 100       # max: 4 blocs (otherwise can't calc variance)
  
    NoBlocks    = int(DataLen/BlockSize)               # total number of such blocks in datastream
    Data        = np.zeros(NoBlocks)                  # container for parcelling block 

    # Loop to chop datastream into blocks
    # and take average
    for i in range(1,NoBlocks+1):
        
        istart = (i-1) * BlockSize
        iend   =  istart + BlockSize
        Data[i-1] = np.mean(data[istart:iend])

    meanVal  = np.mean(Data)
    meanErr  = np.sqrt(np.var(Data)/(NoBlocks - 1))


    return meanVal, meanErr

def read_log(f,w,m,e,s):
    '''Code to read in log file'''

    filename = '{}_{}/log.eq_{}_eps{}_s{}_1'.format(f,w,m,e,s)
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

def read_log_nemd(f,w,m,e,s,force):
    '''Code to read in log file'''

    filename = '{}_{}/log.f{}_{}_eps{}_s{}_1'.format(f,w,force,m,e,s)
    f = open(filename,'r')
    data = f.read()
    data_lines = data.split('\n')

    DATA = []
    flag=0
    for i in range(len(data_lines)-1):

        if data_lines[i].split() != [] and data_lines[i].split()[0] == 'run' and data_lines[i].split()[1] == "${t_run}":
            flag = 1
        if flag == 1 and data_lines[i].split()[0] == 'SHAKE':
            flag = 2
        if flag == 2 and len(data_lines[i].split()) == 10 and data_lines[i].split()[0] != 'Step':
            DATA.append(map(float,data_lines[i].split()))
        if flag == 2 and data_lines[i].split()[0]=='Loop':
            break
    DATA = np.array(DATA)
    return DATA

def prop(PROP,file,factor,f,w,m,e,s):

    filename = '{}_{}/{}.eq_{}_eps{}_s{}_1'.format(f,w,file,m,e,s)

    df = pd.read_csv(filename, delimiter=' ', skiprows=2)
    #read first two lines
    with open(filename, 'r') as f:
        _, line2 = f.readline(), f.readline()
    cols = line2.lstrip('#').strip().split(' ')
    df.columns = cols
    Prop = df['{}'.format(PROP)]
    Prop = Prop.tolist()
    time = df['TimeStep']
    time = time.tolist()

    Prop_list = []
    count = 0
    tlim = 800000
    for t,p in itertools.izip(time, Prop):
        if t > tlim:
            Prop_list.append(p/factor) # correct dimensions etc
            count += 1
    Prop_val, Prop_err = blockAverage(Prop_list)
    print PROP, ':', Prop_val, '+/-', Prop_err
    return Prop_val, Prop_err, Prop_list

def fric(f,w,m,e,s):
    Factor = 1 # 10000 if we want 10^4
    Fric_x = prop('v_ffin_x','fric',Factor,f,w,m,e,s)
    Fric_y = prop('v_ffin_y','fric',Factor,f,w,m,e,s)

    # Might need to consider direction of flow and just use that one
    Fric_avg = (Fric_x[0]+Fric_y[0])/2
    #Fric_avg_err = 0.5*np.sqrt(Fric_x[1]**2+Fric_y[1]**2)
    Fric_avg_err = np.std(np.array([Fric_x[0], Fric_y[0]]))/np.sqrt(2)
    
    return Fric_avg, Fric_avg_err


def wa(f,w,m,e,s):

    WA_list = prop('v_WA_layer','wa',1,f,w,m,e,s)[2]
    Log_Data = read_log(f,w,m,e,s)
    for i in range(len(Log_Data)-1):
        try:
            if Log_Data[i][0] == 'WA' and Log_Data[i][1] == 'layer':
                wa_tmp = float(Log_Data[i][3])
        except:
            continue
    Factor = WA_list[-1]/wa_tmp

    WA_val, WA_err, WA_list = prop('v_WA_layer','wa',Factor,f,w,m,e,s)

    return WA_val, WA_err

def diff(f,w,m,e,s):
    try:
        C_vv_array = np.loadtxt("{}_{}/C_vv_y_eq_{}_eps{}_s{}_1_1_5000_z0_5.dat".format(f,w,m,e,s))
        dt = 0.001
    except:
        C_vv_array = np.loadtxt("../Bulk_properties/{}_{}/C_vv_y_{}_T298_z{}_eps{}_Ds.dat".format(f,w,m,s,e))
        dt = 0.0005
    times =  C_vv_array[:,0]
    C_vv_ave = C_vv_array[:,1]

    time_conv = 1e-12
    space_conv = 1e-10

    int_vacf = simps(C_vv_ave, dx=dt)

    diff_2d = int_vacf/2
    conv = space_conv**2/time_conv
    diff_2d_si = diff_2d*conv

    print 'D_s in 2D:', diff_2d_si

    return diff_2d_si


def wall_pos(Y, RHO):
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
    return sep, ylower, yupper

def dimensions(f,w,m,e,s, rerun):
    name = '{}_{}_{}_eps{}_s{}_1'.format(f,w,m,e,s)
    
    try:
        if rerun=='y':
            raise IOError  
        DEN = np.loadtxt('DATA/den_eq_{}.dat'.format(name))
        VEL = np.loadtxt('DATA/vel_eq_{}.dat'.format(name))
        CY  = np.loadtxt('DATA/cy_eq_{}.dat'.format(name))
        CZ  = np.loadtxt('DATA/cz_eq_{}.dat'.format(name))
    except IOError:
        print 'Generating density data...'
        DEN, VEL, CY, CZ = density_prof(f,w,m,e,s,'eq',-1.9, 66.1)
        np.savetxt('DATA/den_eq_{}.dat'.format(name),DEN)
        np.savetxt('DATA/vel_eq_{}.dat'.format(name),VEL)
        np.savetxt('DATA/cy_eq_{}.dat'.format(name),CY)
        np.savetxt('DATA/cz_eq_{}.dat'.format(name),CZ)

    # Length
    for i in range(DEN.shape[0]-1):
        if DEN[i,0] > 0 and DEN[i+1,0] == 0:
            idxlower = i+1
        if DEN[i,0] ==0 and DEN[i+1,0] > 0:
            idxupper = i
    L = CY[idxupper]-CY[idxlower]
    # new density array
    DEN = DEN[idxlower:idxupper,:]
    CY  = CY[idxlower:idxupper]
    DENz = np.average(DEN, axis = 0)
    DENave = np.average(DEN[DEN!=0])


    H, Zlo, Zhi = wall_pos(CZ, DENz)

    xlo,xhi,ylo,yhi,zlo,zhi = coords(f,w,m,e,s)
    xWidth = xhi-xlo

    print 'The channel dimensions are: L=', L, ', H=', H, 'A.'


    return L, H, CY[0], CY[-1], Zlo, Zhi, xWidth, DENave

def coords(f,w,m,e,s):
    filename = '{}_{}/log.eq_{}_eps{}_s{}_1'.format(f,w,m,e,s)
    f = open(filename,'r')

    data = f.read()
    data_lines = data.split('\n')


    coord_line = re.sub('[\(\)]', '', data_lines[45]).split('=')[1].split('to')
    xlo, ylo, zlo = coord_line[0].strip(' ').split(' ')
    xhi, yhi, zhi = coord_line[1].strip(' ').split(' ')

    return float(xlo), float(xhi), float(ylo), float(yhi), float(zlo), float(zhi)

def density_prof(f,w,m,e,s,mdtype, ylo, yhi):

    if mdtype == 'eq':
        filename = '{}_{}/dens.eq_{}_eps{}_s{}_1'.format(f,w,m,e,s)
    else:
        filename = '{}_{}/dens.f{}_{}_eps{}_s{}_1'.format(f,w,mdtype,m,e,s)
    
    
    ybin = []
    zbin = []
    dx = 0.2

    ylist = list(np.arange(ylo,yhi+0.2,0.2))
    zlist = list(np.arange(-0.5, 63.9, 0.2))

    # round values in list so that they agree with data
    ylist = [round(n, 1) for n in ylist]
    zlist = [round(n, 1) for n in zlist]


    ydim = len(ylist)
    zdim = len(zlist)

    density  = np.zeros((ydim,zdim))
    velocity = np.zeros((ydim,zdim))
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
                den_val = float(line.split()[4])
                vel_val = float(line.split()[6])

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
                        velocity[y_ind][z_ind] = vel_val
                    else:
                        density[y_ind][z_ind] = (density[y_ind][z_ind]+den_val)/2
                        velocity[y_ind][z_ind] = (velocity[y_ind][z_ind]+vel_val)/2
                    coords_y[y_ind][z_ind] = ylist[y_ind]
                    coords_z[y_ind][z_ind] = zlist[z_ind]
                    count_tmp += 1
            count =2
      

    #density_z = np.average(density, axis = 0)
    #density_z_err = stats.sem(density, axis = 0)
    coords_z = coords_z[0,:] #np.average(coords_z, axis = 0)
    #density_y = np.average(density, axis = 1)
    #velocity_z = np.average(velocity, axis = 0)
    #velocity_z_err = stats.sem(velocity, axis = 0)
    coords_y = coords_y[:,0] #np.average(coords_x, axis = 1)


    return density, velocity, coords_y, coords_z

def Qmolcount(f,w,m,e,s,force,rhoave,xwidth):

    filename = '{}_{}/xyzvel.f{}_{}_eps{}_s{}_1'.format(f,w,force,m,e,s)
    
    
    Ybarrier = 40
    mH2O = 2.9e-26 # kg
    dt = 0.001e-12 # s
    rhoave *= 1000 # kg/m^3
    xwidth *= 1e-10 # m
    fac = 100 #1500
    steps_thresh = 0 #4000


    steps = 0
    molcount = 0
    with open(filename) as infile:
        for line in infile:
            if len(line.split()) == 1:
                steps += 1
            if len(line.split()) == 4 and steps==1 and line.split()[0]=='O':
                molcount += 1
            if steps > 1:
                break

    flag_prev1 = np.zeros(molcount)
    flag_prev2 = np.zeros(molcount)
    flag_prev3 = np.zeros(molcount)
    flag_prev4 = np.zeros(molcount)
    flag_prev5 = np.zeros(molcount)
    flag_prev6 = np.zeros(molcount)
    flag_new = np.zeros(molcount)
    steps = 0
    flowcount = 0

    with open(filename) as infile:
        for line in infile:
            if len(line.split()) == 1:
                steps += 1
                flag_prev1 = flag_prev2
                flag_prev2 = flag_prev3
                flag_prev3 = flag_prev4
                flag_prev4 = flag_prev5
                flag_prev5 = flag_prev6
                flag_prev6 = flag_new
                flag_new = np.zeros(molcount)
                idx = 0
            if steps > steps_thresh and len(line.split()) == 4 and line.split()[0]=='O':
                y_val = float(line.split()[2])

                if y_val > Ybarrier and steps == 1:
                    flag_prev6[idx] = 1
                elif y_val > Ybarrier and steps > 1:
                    flag_new[idx] = 1
                if flag_prev6[idx]==0 and flag_prev5[idx]==0 and flag_prev4[idx]==0 and flag_prev3[idx]==0 and flag_prev2[idx]==0 and flag_prev1[idx]==0 and flag_new[idx]==1:
                    flowcount +=1
                #print y_val > Ybarrier, idx, y_val, flag_prev[idx], flag_new[idx]

                idx += 1    
#            if steps > 10:
#                break
    print steps, 'timesteps collected.'
    print flowcount, 'molecules transferred the barrier.'
    Q = (flowcount*mH2O)/(xwidth*rhoave*steps*fac*dt)
    print 'Qmc=',Q


    return Q

def delP(f,w,m,e,s,force):

    NEMD = read_log_nemd(f,w,m,e,s,force)

    # for now only pick the last 2 nanoseconds
    # should rerun with longer equilibration
    Pleft = NEMD[:,8]
    Pright = NEMD[:,9]

    # plot pressure to check convergence
    fig1 = plt.figure(figsize=(9,7)) 
    ax1  = fig1.add_axes([0.1,0.15,0.8,0.75])
    ax1.plot(NEMD[:,0][::10], (NEMD[:,8][::10]-NEMD[:,9][::10])*1e-1)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('$\Delta$P (MPa)')
    fig1.savefig('PLOTS/PDF/dP_{}_{}_{}_eps{}_s{}_f{}.pdf'.format(f,w,m,e,s,force))
    fig1.clear()


    Pleft_val = np.mean(Pleft)
    Pleft_err = stats.sem(Pleft)

    Pright_val = np.mean(Pright)
    Pright_err = stats.sem(Pright)

    deltaP = (Pright_val - Pleft_val)*1e5
    deltaP_err = 1e5*np.sqrt(Pleft_err**2+Pright_err**2)
    print 'Pressure drop: ', deltaP*1e-6, '+/-', deltaP_err*1e-6, 'MPa.'
    return deltaP, deltaP_err

def flux(f,w,m,e,s,force,ylo,yhi,zlo,zhi,rerun):
    name = '{}_{}_{}_eps{}_s{}_1'.format(f,w,m,e,s)
    try:
        if rerun=='y':
            raise IOError  
        DEN = np.loadtxt('DATA/den_f{}_{}.dat'.format(force,name))
        VEL = np.loadtxt('DATA/vel_f{}_{}.dat'.format(force,name))
        CY  = np.loadtxt('DATA/cy_f{}_{}.dat'.format(force,name))
        CZ  = np.loadtxt('DATA/cz_f{}_{}.dat'.format(force,name))
    except IOError:
        print 'Generating density and velocity data...'
        DEN, VEL, CY, CZ = density_prof(f,w,m,e,s,force,ylo,yhi)
        np.savetxt('DATA/den_f{}_{}.dat'.format(force,name),DEN)
        np.savetxt('DATA/vel_f{}_{}.dat'.format(force,name),VEL)
        np.savetxt('DATA/cy_f{}_{}.dat'.format(force,name),CY)
        np.savetxt('DATA/cz_f{}_{}.dat'.format(force,name),CZ)


    DZ = 0.2
    convert = 1e-8 # in m^2/s

    VEL_avg = np.average(VEL,axis=0)
    intVEL = convert*simps(VEL_avg)*DZ
    print 'Qmd = ', intVEL, 'm^2/s.'

    den_max = np.max(DEN)
    den_min = np.min(DEN)
    fig1 = plt.figure(figsize=(14,5)) 
    ax1  = fig1.add_axes([0.1,0.15,0.8,0.75])
    #fig1.text(0.44, 0.025, '$y$ (\AA)', ha='center', va='center', fontsize=26)
    ax1.set_xlim(ylo,yhi)
    ax1.set_ylim(zlo,zhi)
    ax1.set_xlabel('$y$ (\AA)')
    ax1.set_ylabel('$z$ (\AA)')
    #plt.tick_params(pad=7)
    ctest=ax1.contourf(CY, CZ, np.transpose(DEN), cmap=cm.RdBu, levels=np.linspace(den_min,den_max,500))
    fig1.colorbar(ctest)
    fig1.savefig('PLOTS/PNG/den_2d_f{}_{}.png'.format(force,name))
    fig1.clear()

    vel_max = np.max(VEL)
    vel_min = np.min(VEL)
    fig2 = plt.figure(figsize=(14,5)) 
    ax2  = fig2.add_axes([0.1,0.15,0.8,0.75])
    #fig1.text(0.44, 0.025, '$y$ (\AA)', ha='center', va='center', fontsize=26)
    ax2.set_xlim(ylo,yhi)
    ax2.set_ylim(zlo,zhi)
    ax2.set_xlabel('$y$ (\AA)')
    ax2.set_ylabel('$z$ (\AA)')
    #plt.tick_params(pad=7)
    ctest=ax2.contourf(CY, CZ, np.transpose(VEL), cmap=cm.RdBu, levels=np.linspace(vel_min,vel_max,500))
    fig2.colorbar(ctest)
    fig2.savefig('PLOTS/PNG/vel_2d_f{}_{}.png'.format(force,name))
    fig2.clear()

    fig3 = plt.figure(figsize=(9,7)) 
    ax3  = fig3.add_axes([0.1,0.15,0.8,0.75])
    #ax3.plot(CZ, np.average(DEN,axis=0))
    ax3.plot(CZ, VEL_avg)
    ax3.set_xlim(zlo,zhi)
    #ax3.set_ylim(-1,1)
    ax3.set_xlabel('$z$ (\AA)')
    ax3.set_ylabel('$v_y$ (\AA/ps)')

    fig3.savefig('PLOTS/PDF/vel_z_f{}_{}.pdf'.format(force,name))
    fig3.clear()

    #Jz_mol = convert*Jz/mol_mass_water 

    #Jz_err = err_Jz(rho, rho_err, vel, vel_err)*convert/mol_mass_water

     
    return intVEL

def vol_flow_P(dP, eta, L, D):
    Q = (-2/3)*(dP/(eta*L))*D**3
    print 'QP = ', Q, 'm^2/s.'
    return Q



def straight_fit(x, y, xmin, xmax):
    slope, cov = curve_fit(f, np.array(x), np.array(y))
    slope_err = np.sqrt(np.diag(cov))
    xdat = np.linspace(xmin, xmax, 100)
    fit = []
    for x in xdat:
        fit.append(slope*x)
    return xdat, fit, slope, slope_err

def straight_fit2(X, Y, xmin, xmax):
    Xnew, Ynew = [],[]
    for x,y in itertools.izip(X,Y):
        if y > 0.2:
            Xnew.append(x)
            Ynew.append(y)
    print Xnew, Ynew
    params, cov = curve_fit(f1, np.array(Xnew), np.array(Ynew))
    #slope_err = np.sqrt(np.diag(cov))
    xdat = np.linspace(xmin, xmax, 100)
    fit = []
    for x in xdat:
        fit.append(params[0]*x+params[1])
    return xdat, fit

def f(x, A):
    return A*x

def f1(x, A, B):
    return A*x+B


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

