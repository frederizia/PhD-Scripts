#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import argparse
from class_bulk_props import *
from plotting_params import *
import sys

def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-T', '--temp', nargs='+', required=False, default='298', action='store',
                       help='Temp')
    parser.add_argument('-P', '--press', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-s', '--sep', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-f', '--fluid', nargs='+', required=False, default='Water', action='store')
    parser.add_argument('-m', '--model', nargs='+', required=False, default='spce', action='store')
    parser.add_argument('-e', '--eps', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-r', '--rhos', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-d', '--den', nargs='+', required=False, default=[5,6], action='store')
    parser.add_argument('-png', '--png', nargs=1, type=str, required=False, default='n', action='store')
    args = parser.parse_args()
    return args

def poly1(xlist,m,b):
    result = []
    for x in xlist:
        result.append(m*x+b)
    return result

def averaging2(k,v):
    averages = {}
    counts = {}
    for name, value in zip(k, v):
        if name in averages:
            averages[name] += value
            counts[name] += 1
        else:
            averages[name] = value
            counts[name] = 1
    for name in averages:
        averages[name] = averages[name]/float(counts[name]) 
        print averages[name]
    k = averages.keys()
    v = averages.values()
    return k,v

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

    k = averages.keys()
    v = averages.values()
    e = errors.values()
    return k,v,e

def averaging3(k,v,e):
    averages = {}
    errors = {}
    counts = {}
    for name, value, error in zip(k, v, e):
        if name in averages:
            averages[name].append(value)
            errors[name].append(error)
            counts[name] += 1
        else:
            averages[name] = [value]
            errors[name] = [error]
            counts[name] = 1
    for name in averages:
        ave_array = np.array(averages[name])
        ave_val = np.mean(ave_array)
        err_array = np.array(errors[name])
        ave_error = np.sqrt(np.sum(np.square(err_array)))/counts[name]
        averages[name] = ave_val
        errors[name] = ave_error 

    k = averages.keys()
    v = averages.values()
    e = errors.values()
    return k,v,e

def main():
    args = GetArgs()
    temp = args.temp
    press = args.press
    fluid = args.fluid
    model = args.model
    den = args.den
    png = args.png[0]
    fig_size = (9,7)



    print 'Evaluating properties for', fluid, 'using the', model, 'model at', temp, 'K for pressures', press
    
    # Define plots
    fig1 = plt.figure(figsize=fig_size)
    ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])
    fig2 = plt.figure(figsize=fig_size)
    ax2  = fig2.add_axes([0.1,0.15,0.8,0.75])
    fig3 = plt.figure(figsize=fig_size)
    ax3  = fig3.add_axes([0.1,0.15,0.8,0.75])
    fig4 = plt.figure(figsize=fig_size)
    ax4  = fig4.add_axes([0.1,0.15,0.8,0.75])
    fig5 = plt.figure(figsize=fig_size)
    ax5  = fig5.add_axes([0.1,0.15,0.8,0.75])
    fig6 = plt.figure(figsize=fig_size)
    ax6  = fig6.add_axes([0.1,0.15,0.8,0.75])
    fig7 = plt.figure(figsize=fig_size)
    ax7  = fig7.add_axes([0.1,0.15,0.8,0.75])
    fig8 = plt.figure(figsize=fig_size)
    ax8  = fig8.add_axes([0.1,0.15,0.8,0.75])


    markers = ['D', 's', 'v', '^', 'd', '*']
    legend_names = {'spce':'SPC/E', 'tip4p':'TIP4P/2005', 'TraPPE': 'TraPPE', 'SAFT': 'SAFT Dimer', \
    'SAFT1': 'SAFT Monomer', 'TraPPEnc': 'TraPPE (no charges)', 'OPLS': 'OPLS',\
    'SAFT1_rc2.5': 'SAFT (rc=2.5)', '12_6':'12-6 model', '12_6_rcsaft': '12-6 model (rc from SAFT)',\
    '12_6_23_666': '12-6 model (SAFT exponents)',
    'SAFTflex': 'SAFT Dimer (flexible II)', 'EPM2':'EPM2', 'SAFTrigid': 'SAFT (rigid)',\
    'SAFT1vle': 'SAFT-vle', 'SAFT1ift': 'SAFT-ift', 'SAFTflex2': 'SAFT Dimer (flexible)',\
    'OPLSmod': 'L-OPLS', 'TraPPEflex': 'TraPPE (flexible)'}

    cmap = matplotlib.cm.get_cmap('viridis')
    # Using contourf to provide my colorbar info, then clearing the figure
    Z = [[0,0],[0,0]]
    Min, Max = (0.9, 1.1)
    step = 50
    levels = np.linspace(Min,Max,step)
    CS3 = plt.contourf(Z, levels, cmap=cmap)
    plt.clf()

    try:
        NIST = read_nist(temp, fluid)

        ax1.plot(NIST[1], NIST[0], linestyle = 'dashed', marker='+', label='NIST', c='k')

        ax2.plot(NIST[0], NIST[2], linestyle = 'dashed', marker='+', label='NIST', c='k')

        ax7.plot(NIST[1], NIST[2], linestyle = 'dashed', marker='+', label='NIST', c='k')

        #ax4.plot(NIST[0], NIST[3], linestyle = 'dashed', label='NIST', c='k')
    except IOError:
        print 'No NIST data available'
    except:
        print 'Some error. Check' 

    count = 0
    name_plot = ''
    timestamp = 8000000
    deltat = 6000000
    for f in fluid:
        for m in model:
            for t in temp:
                tmp=0
                print '\n\n++++++++++++++++++++++++ Fluid:%s, model: %s, Temp: %s +++++++++++++++++++++\n' %(f,m,t)
                press_dat = []
                temp_dat =[]
                rho_dat_init = []
                shear_dat = []
                shear_err = []
                ratio_dat = []
                bulk_dat = []
                bulk_err = []
                diff_dat = []

                if f=='LJ' and m=='lj':
                    try:
                        print 'Reading in data..'
                        final_data = np.loadtxt('LJ/relations_T{}.dat'.format(t))
                        rho_dat = list(final_data[:,0])
                        press_dat = list(final_data[:,1])
                        shear_dat = list(final_data[:,2])
                        bulk_dat = list(final_data[:,3])
                        ratio_dat=[]
                        for s,b in itertools.izip(final_data[:,2],final_data[:,3]):
                            ratio_dat.append(b/s)
                        tmp+=1
                    except:
                        print 'ERROR: No data for %s, %s and T=%s exists.\n'%(f,m,t)

                else:
                    for p in press:
                        for i in den:
                            print 'Evaluating %s, DEN=%s and P=%s.'%(m,i,p)                           
                            try:
                                results = bulk_properties(f,m,t,p,i)

                                press_dat.append(results.press())
                                temp_dat.append(results.temp())
                                rho_dat_init.append(results.rho())
                                shear_dat.append(results.shear2()[0])
                                shear_err.append(results.shear2()[1])
                                bulk_dat.append(results.bulk2()[0])
                                bulk_err.append(results.bulk2()[1])
                                ratio_dat.append(results.visc_ratio())
                                diff_dat.append(results.diff(timestamp,deltat)[-1])
                                print 'Data collection successful.\n'
                                tmp+=1
                            except:
                                print 'ERROR: Maybe no data for %s, DEN=%s and P=%s exists.\n'%(m,i,p)

                if tmp != 0:
                    name_plot += f
                    if f !='LJ':
                        rho_dat, press_dat, press_err   = averaging(rho_dat_init, press_dat)
                        rho_dat, shear_dat, shear_err   = averaging(rho_dat_init, shear_dat)#, shear_err)
                        rho_dat, bulk_dat, bulk_err     = averaging(rho_dat_init, bulk_dat)#, bulk_err)
                        rho_dat, ratio_dat, ratio_err   = averaging(rho_dat_init, ratio_dat)
                    else:
                        shear_err = np.zeros(len(shear_dat))
                        bulk_err = np.zeros(len(shear_dat))
                        press_err = np.zeros(len(shear_dat))
                        ratio_err = np.zeros(len(shear_dat))

                    if f=='Water':
                        Tcrit = 667
                        rhocrit = 0.283
                    elif f=='CO2':
                        Tcrit = 311.13
                        rhocrit = 0.48
                        Pcrit = 80
                    elif f=='LJ':
                        Tcrit = 1.312
                        rhocrit = 0.316
                    else:
                        print 'Unknown fluid'
                        sys.exit(1)

                    rho_dat[:] = [x / rhocrit for x in rho_dat]
                    Tnew = float(t)/Tcrit

                    Tarray = np.linspace(0.9,1.1,50)
                    Carray = np.linspace(0.1,1.0,50)
                    c_tmp  = cmap(Carray[np.digitize(Tnew,Tarray)])

                    # Plotting
                    ax1.errorbar(press_dat, rho_dat, xerr=press_err, linestyle = 'None', c=c_tmp, marker=markers[count], label='%s, $T^*$=%.2f'%(f,Tnew))
                    ax1.set_xlabel('P (bar)')
                    ax1.set_ylabel('$\\rho$/$\\rho_{crit}$')
                    ax1.legend()

                    ax2.errorbar(rho_dat, shear_dat, yerr=shear_err, linestyle = 'None', c=c_tmp, marker=markers[count], label='%s, $T^*$=%.2f'%(f,Tnew))
                    ax2.set_xlabel('$\\rho$/$\\rho_{crit}$')
                    ax2.set_ylabel('$\eta$ (mPa.s)')
                    ax2.legend()

                    ax3.errorbar(rho_dat, bulk_dat, yerr=bulk_err,linestyle = 'None', c=c_tmp, marker=markers[count], label='%s, $T^*$=%.2f'%(f,Tnew))
                    ax3.set_xlabel('$\\rho$/$\\rho_{crit}$')
                    ax3.set_ylabel('$\kappa$ (mPa.s)')
                    ax3.legend()

                    

                    ax5.errorbar(rho_dat, ratio_dat, yerr=ratio_err, linestyle = 'None', c=c_tmp, marker=markers[count], label='%s, $T^*$=%.2f'%(f,Tnew))
                    ax5.set_xlabel('$\\rho$/$\\rho_{crit}$')
                    ax5.set_ylabel('$\kappa$/$\eta$')
                    #cbar = fig5.colorbar(CS3, ticks=[-1, 0, 1])
                    #cbar.set_yticklabels(['0.9', '1.0', '1.1'])
                    ax5.legend()

                    ax7.errorbar(press_dat, shear_dat, yerr=shear_err, linestyle = 'None', marker=markers[count], label='%s, $T^*$=%.2f'%(f,Tnew))
                    ax7.set_xlabel('P (bar)')
                    ax7.set_ylabel('$\eta$ (mPa.s)')
                    ax7.legend()

                    ax8.errorbar(press_dat, bulk_dat, yerr=bulk_err,linestyle = 'None', marker=markers[count], label='%s, $T^*$=%.2f'%(f,Tnew))
                    ax8.set_xlabel('P (bar)')
                    ax8.set_ylabel('$\kappa$ (mPa.s)')
                    ax8.legend()

                    count +=1
                else:
                    print 'No data collected for:', f,m,t
    cbar = fig1.colorbar(CS3, ticks=[0.9,1.0,1.1])
    cbar = fig2.colorbar(CS3, ticks=[0.9,1.0,1.1])
    cbar = fig3.colorbar(CS3, ticks=[0.9,1.0,1.1])
    cbar = fig5.colorbar(CS3, ticks=[0.9,1.0,1.1])

    if png == 'y':
        fig1.savefig('PLOTS/PNG/crit_P_%s.png'%(name_plot))
        fig2.savefig('PLOTS/PNG/crit_SV_%s.png'%(name_plot))
        fig3.savefig('PLOTS/PNG/crit_BV_%s.png'%(name_plot))
        fig5.savefig('PLOTS/PNG/crit_BS_%s.png'%(name_plot))
        fig7.savefig('PLOTS/PNG/crit_SV_vP_%s.png'%(name_plot))
        fig8.savefig('PLOTS/PNG/crit_BV_vP_%s.png'%(name_plot))
    else:
        fig1.savefig('PLOTS/crit_P_%s.pdf'%(name_plot))
        fig2.savefig('PLOTS/crit_SV_%s.pdf'%(name_plot))
        fig3.savefig('PLOTS/crit_BV_%s.pdf'%(name_plot))
        fig5.savefig('PLOTS/crit_BS_%s.pdf'%(name_plot))
        fig7.savefig('PLOTS/crit_SV_vP_%s.pdf'%(name_plot))
        fig8.savefig('PLOTS/crit_BV_vP_%s.pdf'%(name_plot))

    return




if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS')
    main()
