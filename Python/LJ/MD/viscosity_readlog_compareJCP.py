#!/usr/bin/env python

"""Reads in log files and extracts pressure, shear and bulk viscosity data and compares to data from Meier et al 2004 JCP"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
from matplotlib import rc


np.set_printoptions(precision=16)
# Parse arguments form command line


parser = argparse.ArgumentParser()


parser.add_argument("-T", type=str, nargs='+', \
                    help="Temperature",required=True)
parser.add_argument("-rhof", type=str, nargs='+', \
                    help="Fluid density",required=False, default = ['0.1', '0.15','0.2', '0.25','0.3', '0.35', '0.4', '0.45', '0.5', '0.55', '0.6', '0.65', '0.7', '0.75', '0.8', '0.85', '0.9'])#, '0.95', '0.975', '1.0'])



args = parser.parse_args()
T = args.T
dens = args.rhof

if len(T) ==1:
	T = T[0]
'''
press = []
shear = []
bulk = []
density = []


def read_log(T, den):

    if T=='1.7':
    	filename = "/media/fred/My Passport/LJ_relations_results/T%s/SMALLandEQN/rhof%s/log.lammps_rel" % (T, den)

    else:
    	filename = "/media/fred/My Passport/LJ_relations_results/T%s/log.rel_rhof%s" % (T,den)


    # read in data
    f = open(filename,'r')
    data = f.read()


    data_lines = data.split('\n')
    density.append(float(den))
    print 'Density: ', den

    for i in range(0,len(data_lines)-1):

        log_data = data_lines[i].split()

	if len(log_data) >1:
			
		if log_data[1] == 'shear':
			print 'Shear viscosity: ', log_data[4]
			shear.append(float(log_data[4]))
		if log_data[1] == 'bulk' and log_data[2] == 'pt':
			print 'Bulk viscosity: ', log_data[5]
			bulk.append(float(log_data[5]))
		if log_data[0] == 'Loop':
			dat = data_lines[i-1].split()
			print 'Pressure: ', dat[5]
			press.append(float(dat[5]))
		continue
		
    return 



def save_data(D, P, SH, BU):

	data_list = [D,P,SH,BU]
	data_array = np.transpose(np.array(data_list))
	print data_array
	print type(data_array[0][0])
	np.savetxt('relations_T%s.dat'%(T), data_array,fmt="%s",header = 'Density Pressure Shear Bulk')
	return data_array'''

def read_data(name):
	JCPdata = np.loadtxt('%s'%name)
	JCPrho = JCPdata[:,0]
	JCPeta = JCPdata[:,2]
	return JCPrho, JCPeta


def P_17(rho):
	rho = np.array(map(float,rho))
	P = 0.03233*np.exp(6.049*rho)
	return P

def P_20(rho):
	rho = np.array(map(float,rho))
	P = 0.07503*np.exp(5.318*rho)
	return P

def P_25(rho):
	rho = np.array(map(float,rho))
	return 0.1568*np.exp(4.773*rho)

def shear_17(rho):
	rho = np.array(map(float,rho))
	return 5.426*rho**(5.344) + 0.3384

def shear_20(rho):
	rho = np.array(map(float,rho))
	return 4.196*rho**(4.336) + 0.4596

def shear_25(rho):
	rho = np.array(map(float,rho))
	return 2.675*rho**(2.302) + 0.1549

def bulk_17(rho):
	rho = np.array(map(float,rho))
	return 2.926*rho**(2.114)

def bulk_20(rho):
	rho = np.array(map(float,rho))
	return 0.01689*np.exp(5.463*rho)

def bulk_25(rho):
	rho = np.array(map(float,rho))
	return 2.882*rho**(1.74)
	#return 12.18*rho**(104.6)+1.167

def Woodcock(rho, eta_0, T):
	rho = np.array(map(float,rho))
	T = float(T)
	term1 = np.sqrt(2)*(1-(1/T)**4-(T/8))*rho
	term2 = 3.025*(1/T)**(1/3)*rho**4
	return eta_0*(1+term1+term2)

	
# PLOTTING

name = 'JCPdataEtaT1.35.dat'
JCP_rho, JCP_eta = read_data(name)
#sort arrays by density
sorti = np.argsort(JCP_rho)
JCP_rho, JCP_eta = JCP_rho[sorti], JCP_eta[sorti]


'''dens_new =[]

for d in dens:
	try:
		read_log(T,d)
		dens_new.append(d)
	except IOError: pass'''
eta_0 = JCP_eta[0]

#final_data = save_data(dens_new,press,shear,bulk)


matplotlib.rcParams.update({'font.size': 19})
rc('text', usetex=True)



fig = plt.figure()
#plt.plot(final_data[:,0], final_data[:,2], marker ='D', linestyle = 'none')
plt.plot(JCP_rho, Woodcock(JCP_rho, eta_0, T), lw = 2.0,linestyle = 'dashed', c='r', label = "Woodcock")
plt.plot(JCP_rho, JCP_eta, lw = 2.0,linestyle = 'dashed', marker = 'D',c='b', label = "JCP 2004")
plt.ylabel('$\mathrm{Shear}$ $\mathrm{viscosity}$')
plt.xlabel('$\\rho$')
plt.xlim(0.1,0.9)
plt.legend()
#plt.ylim(0,6)
plt.savefig('PLOTS/comp_rho_v_shear_T%s.pdf'%(T))
plt.show()




