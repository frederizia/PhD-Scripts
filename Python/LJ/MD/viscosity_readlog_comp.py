#!/usr/bin/env python

"""Reads in log files and extracts pressure, shear and bulk viscosity data"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
from matplotlib import rc
import re

np.set_printoptions(precision=16)
# Parse arguments form command line
parser = argparse.ArgumentParser()


parser.add_argument("-var", type=str, nargs=1, \
                    help="Variable",required=True)
parser.add_argument("-T", type=str, nargs='+', \
                    help="Temperature",required=True)
parser.add_argument("-es", type=str, nargs='+', \
                    help="epsilon or sigma",required=True)
parser.add_argument("-rhof", type=str, nargs='+', \
                    help="Fluid density",required=False, default = ['0.1', '0.15','0.2', '0.25','0.3', '0.35', '0.4', '0.45', '0.5', '0.55', '0.6', '0.65', '0.7', '0.75', '0.8', '0.85', '0.9'])#, '0.95', '0.975', '1.0'])



args = parser.parse_args()
var = args.var[0]
dens = args.rhof

if var == 'T':
	T = args.T

if var == 'eps' or var == 'sigma':
	T_tmp = args.T[0]
	ES = args.es
	T = []
	lgd_items = []
	for es in ES:
		T.append(T_tmp+'_'+es)
		test_str = '%s(.*)'%var
		lgd_items.append(re.search(test_str,es).group(1))
	if var == 'eps':
		var_lgd = '$\epsilon$'
        if var == 'sigma':
                var_lgd = '$\sigma$'


def read_log(Tval, den, density, press, shear, bulk):

    filename = "/media/fred/My Passport/LJ_relations_results/CORRECT/T%s/rhof%s/log.rel_rhof%s" % (Tval,den,den)

    # read in data
    f = open(filename,'r')
    data = f.read()


    data_lines = data.split('\n')
  
    density.append(float(den))
    print 'Density: ', den

    for i in range(0,len(data_lines)-1):

        log_data = data_lines[i].split()

	if len(log_data) >1:
			
                if log_data[1] == 'shear' and log_data[2] == 'viscosity' and log_data[3] == ':':
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
		
    return density, shear, bulk, press 



def save_data(Tval, D, P, SH, BU):

	data_list = [D,P,SH,BU]
	data_array = np.transpose(np.array(data_list))
	print data_array
	np.savetxt('relations_T%s.dat'%(Tval), data_array,fmt="%s",header = 'Density Pressure Shear Bulk')
	return data_array

def P_17(rho):
	rho = np.array(map(float,rho))
	P = 0.163*np.exp(-1.078*rho)+0.02561*np.exp(6.305*rho)
	return P

def shear_17(rho):
	rho = np.array(map(float,rho))
	return 0.167*np.exp(1.715*rho)+0.008242*np.exp(6.292*rho)

def bulk_17(rho):
	rho = np.array(map(float,rho))
	return 14.41*rho**2.805*np.exp(-2.936*rho)


	
# PLOTTING


final_list = [] 
for t in T:
	print t
	Dens =[]
	Press = []
	Shear = []
	Bulk = []
	for d in dens:
		try:
			read_log(t,d, Dens, Press, Shear, Bulk)
		except IOError: 
			print t,d
			pass

	final_data = save_data(t,Dens,Press,Shear,Bulk)
	final_list.append(final_data)

rhovals = np.arange(0.05,0.91,0.01)

matplotlib.rcParams.update({'font.size': 19})
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

fig = plt.figure()
for i in range(len(T)):
	plt.plot(final_list[i][:,0], final_list[i][:,1], marker ='D' , linestyle = 'none', label = '%s = %s'%(var_lgd, lgd_items[i]))
plt.ylabel('$\mathrm{P}$')
plt.xlabel('$\\rho$')
plt.xlim(0.1,0.9)
plt.legend(loc='upper left')
plt.savefig('PLOTS/rho_v_P_%s.pdf'%(var))
plt.show()

fig = plt.figure()
for i in range(len(T)):
        plt.plot(final_list[i][:,0], final_list[i][:,2], marker ='D' , linestyle = 'none', label = '%s = %s'%(var_lgd, lgd_items[i]))
plt.ylabel('$\eta$')
plt.xlabel('$\\rho$')
plt.xlim(0.1,0.9)
plt.legend(loc='upper left')
#plt.ylim(0,6)
plt.savefig('PLOTS/rho_v_shear_%s.pdf'%(var))
plt.show()

fig = plt.figure()
for i in range(len(T)):
        plt.plot(final_list[i][:,0], final_list[i][:,3], marker ='D' , linestyle = 'none', label = '%s = %s'%(var_lgd, lgd_items[i]))
plt.ylabel('$\kappa$')
plt.xlabel('$\\rho$')
plt.xlim(0.1,0.9)
#plt.ylim(-2,14)
plt.legend(loc='upper left')
plt.savefig('PLOTS/rho_v_bulk_%s.pdf'%(var))
plt.show()



# COMPARE EXPRESSIONS FOR T
'''

fig = plt.figure()
plt.plot(rhovals, P_135(rhovals), label = '$T=1.35$', lw = 2.0)
plt.plot(rhovals, P_17(rhovals), label = '$T=1.7$', lw = 2.0)
plt.plot(rhovals, P_20(rhovals), label = '$T=2.0$', lw = 2.0)
plt.plot(rhovals, P_25(rhovals), label = '$T=2.5$', lw = 2.0)
plt.ylabel('$\mathrm{P}$')
plt.xlabel('$\\rho$')
plt.legend()
plt.xlim(0.1,0.9)
#plt.ylim(0,20)
plt.savefig('PLOTS/rho_v_P.pdf')
#plt.show()

fig = plt.figure()
plt.plot(rhovals, shear_135(rhovals), label = '$T=1.35$', lw = 2.0)
plt.plot(rhovals, shear_17(rhovals), label = '$T=1.7$', lw = 2.0)
plt.plot(rhovals, shear_20(rhovals), label = '$T=2.0$', lw = 2.0)
plt.plot(rhovals, shear_25(rhovals), label = '$T=2.5$', lw = 2.0)
plt.ylabel('$\eta$')
plt.xlabel('$\\rho$')
plt.legend()
plt.xlim(0.1,0.9)
#plt.ylim(0,20)
plt.savefig('PLOTS/rho_v_shear.pdf')
#plt.show()

fig = plt.figure()
plt.plot(rhovals, bulk_135(rhovals), label = '$T=1.35$', lw = 2.0)
plt.plot(rhovals, bulk_17(rhovals), label = '$T=1.7$', lw = 2.0)
plt.plot(rhovals, bulk_20(rhovals), label = '$T=2.0$', lw = 2.0)
plt.plot(rhovals, bulk_25(rhovals), label = '$T=2.5$', lw = 2.0)
plt.ylabel('$\kappa$')
plt.xlabel('$\\rho$')
plt.legend()
plt.xlim(0.1,0.9)
#plt.ylim(0,20)
plt.savefig('PLOTS/rho_v_bulk.pdf')
#plt.show()


'''
