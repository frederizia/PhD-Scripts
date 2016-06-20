from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import ast
from matplotlib import cm
import pylab as p
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
import argparse
import re

from scipy.optimize import curve_fit

# Script to 1) Fit data to a curve 2) Convert inbetween units metal (lammps) and our scaled units


# define mu = a1 + a2rho^(80.01)
# define lambda = b1 exp((-rho+b2)^2/b3)+b4 exp(-(rho+b5)/b6)
# define rho = c1 p^2 + c2 p + c3


# FUNCTIONS

def mu_fn(rho,A1,A2):
	mu = A1 + A2*rho**(80.01)
	return mu

def lmda_fn(rho,B1,B2,B3,B4,B5,B6):
	lmda = B1*np.exp(-((rho-B2)/B3)**2)+B4*np.exp(-((rho-B5)/B6)**2)
	return lmda

def rho_fn(p,C1,C2,C3):
	rho = C1*p**2 + C2*p + C3
	return rho

def Lfunc(x, a, b, c):
    return a * np.exp(-b * x) + c

def P(rho,D1,D2,D3):
	return D1 + np.sqrt(D2+D3*rho)

# FIT DATA FOR LAMBDA

den_test = np.array([0.999842,0.999164,0.995594,0.992099,0.988268,0.982537,0.976562,
0.970808,0.963678,0.956278])

lmda_test = np.array([3.292,2.553,1.981,1.318,1.252,1.198,0.836,0.724,0.764,0.739])


# fit curve for lambda
popt, pcov = curve_fit(Lfunc, den_test, lmda_test,maxfev=100000)


# Conversion

viscosityToSI = 10**(-3) # mPas -> Pas
densityToSI = 10**(3) # g/cm^3 -> kg/m^3
pressureToSI = 10**(3)# kPa -> Pa

viscosityToA = 10**(-12) # mPas -> kg/As
densityToA = 10**(-21) # g/cm^3 -> g/A^3
pressureToA = 10**(6)# kPa -> kg/As^2


# in units metal
# mu: mPas
# lambda: mPas
# rho: g/cm^3
# P: kPa

a1 = 0.2859 # mPas
a2 = 0.7016 # (g/cm^3)^(-79.01)*mPas

'''
b1 = 1.529e18#1.54e15#5.559*10**(-9) # mPas
b2 = 1.109#1.243#-1.508 # g/cm^3
b3 = 0.01691#0.04141#0.1085 # g/cm^3
b4 = 6864#7.055#8.574e23 # mPas
b5 = 1.632#1.199#72.55 # g/cm^3
b6 = 0.2204#0.1573#9.876 # g/cm^3'''

b1 = popt[0] # mPas
b2 = popt[1] # 1/(g/cm^3)
b3 = popt[2] # mPas


c1 = 4.45*10**(-6)/densityToSI # (kPa)^-2*g/cm^3
c2 = -0.000134/densityToSI # (kPa)^-1*g/cm^3
c3 = 990.2/densityToSI # g/cm^3

d1 = -c2/(2*c1) # kPa
d2 = d1**2-(c3/c1) # kPa^2
d3 = (1/c1) # kPa^2/(g/cm^3)

# convert these to SI




a1_SI = a1*viscosityToSI
a2_SI = a2*pressureToSI**(-79.01)*viscosityToSI

'''b1_SI = b1*viscosityToSI
b2_SI = b2*densityToSI
b3_SI = b3*densityToSI
b4_SI = b4*viscosityToSI
b5_SI = b5*densityToSI
b6_SI = b6*densityToSI'''

b1_SI = b1*viscosityToSI
b2_SI = b2/densityToSI
b3_SI = b3*viscosityToSI

c1_SI = c1*pressureToSI**(-2)*densityToSI
c2_SI = c2*pressureToSI**(-1)*densityToSI
c3_SI = c3*pressureToSI

d1_SI = d1*pressureToSI
d2_SI = d2 * pressureToSI**2
d3_SI = pressureToSI**2/densityToSI




a1_A = a1*viscosityToA
a2_A = a2*pressureToA**(-79.01)*viscosityToA


b1_A = b1*viscosityToA
b2_A = b2/densityToA
b3_A = b3*viscosityToA

c1_A = c1*pressureToA**(-2)*densityToA
c2_A = c2*pressureToA**(-1)*densityToA
c3_A = c3*pressureToA

d1_A = d1*pressureToA
d2_A = d2 * pressureToA**2
d3_A = pressureToA**2/densityToA



'''
# formula on drive not complete
# go by graph

lmda_exp = np.log(-np.log((0.75-b1*np.exp(-((0.95+b2)/b3)**2))/b4))/np.log((0.95+b5)/b6)
print lmda_exp

print lmda_fn(0.95,b1,b2,b3,b4,b5,b6)
print lmda_fn(950,b1_SI,b2_SI,b3_SI,b4_SI,b5_SI,b6_SI)'''


print 'P(1)', P(1,d1,d2,d3), 'P(1.05)', P(1.05,d1,d2,d3)

# scale to initial values

P_0 = P(1,d1,d2,d3)#100000#1e-7
rho_0 = rho_fn(P_0,c1,c2,c3)
mu_0 = mu_fn(rho_0,a1,a2)
#lmda_0 = lmda_fn(b1,b2,b3,b4,b5,b6,rho_0)
lmda_0 = Lfunc(rho_0,b1,b2,b3)

a1_sc = a1/mu_0
a2_sc = (a2/mu_0)*rho_0**(80.01)

'''
b1_sc = (b1/lmda_0)
b2_sc = (b2/rho_0)
b3_sc = (b3/rho_0)
b4_sc = (b4/lmda_0)
b5_sc = (b5/rho_0)
b6_sc = (b6/rho_0)'''

b1_sc = b1/lmda_0
b2_sc = b2*rho_0
b3_sc = b3/lmda_0


c1_sc = (c1/rho_0)*P_0**2 
c2_sc = (c2/rho_0)*P_0
c3_sc = (c3/rho_0)

d1_sc = d1/P_0
d2_sc = d2/P_0**2
d3_sc = (d3*rho_0)/P_0**2


print P_0, rho_0, mu_0

print "mu"
print a1, a1_SI, a1_sc, a1_A
print a2, a2_SI, a2_sc, a2_A

print "lambda"
print b1, b1_SI, b1_sc, b1_A
print b2, b2_SI, b2_sc, b2_A
print b3, b3_SI, b3_sc, b3_A
#print b4, b4_SI, b4_sc
#print b5, b5_SI, b5_sc
#print b6, b6_SI

print "rho"
print c1, c1_SI, c1_sc, c1_A
print c2, c2_SI, c2_sc, c2_A
print c3, c3_SI, c3_sc, c3_A

print "P"
print d1, d1_SI, d1_sc, d1_A
print d2, d2_SI, d2_sc, d2_A
print d3, d3_SI, d3_sc, d3_A

print 'scaled'
print rho_0, 'mu_sc(rho_0)', mu_fn(rho_0,a1_sc,a2_sc)
print 'P(1)', P(1,d1_sc,d2_sc,d3_sc), 'P(1.05)', P(1.05,d1_sc,d2_sc,d3_sc)

den_array = np.linspace(0.950,1.050,100)
den_array_SI = np.linspace(950,1050,1000)
den_array_sc = den_array/rho_0
P_array = np.arange(0,600,100)
P_array_SI = np.arange(0,600000,100)
P_array_sc = P_array/P_0


#print lmda_fn(b1,b2,b3,b4,b5,b6,den_array)

#----------------------------PLOTTING-----------------
'''

plt.figure()
#plt.plot(den_array, lmda_fn(b1,b2,b3,b4,b5,b6,den_array),label='metal')
#plt.plot(den_test,lmda_test,linestyle='None', marker = 'D')
plt.plot(den_array,P(den_array,d1,d2,d3))
#plt.plot(den_array, mu(a1_SI,a2_SI,den_array_SI),label='SI')
#plt.plot(P_array, rho_fn(P_array,c1,c2,c3))
plt.xlim(0.95,1.01)
#plt.ylim(0,5)
plt.legend()
plt.show()

plt.figure()
#plt.plot(den_array, mu(a1,a2,den_array),label='metal')
#plt.plot(den_array_SI, lmda_fn(b1_SI,b2_SI,b3_SI,b4_SI,b5_SI,b6_SI,den_array_sc),label='SI')
#plt.plot(den_array, P(den_array_SI,d1_SI,d2_SI,d3_SI))
plt.legend()
#plt.xlim(960,1010)
plt.show()

'''

