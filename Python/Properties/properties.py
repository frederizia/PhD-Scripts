#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
from class_bulk_props import *
from plotting_params import *
import sys

def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--properties', required=False, default ='all', help='Properiese')
    parser.add_argument('-T', '--temp', required=False, default='298', action='store',
                       help='Temp')
    parser.add_argument('-P', '--press', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-s', '--sep', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-f', '--fluid', required=False, default='Water', action='store')
    parser.add_argument('-m', '--model', nargs='+', required=False, default='298', action='store')
    parser.add_argument('-e', '--eps', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-d', '--den', nargs='+', required=False, default=[1,2,3,4,5,6], action='store')
    args = parser.parse_args()
    return args

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


def main():
    args = GetArgs()
    props = args.properties
    temp = args.temp
    press = args.press
    fluid = args.fluid
    model = args.model
    sep = args.sep
    eps = args.eps
    den = args.den

    #matplotlib.rcParams.update({'font.size': 18})
    #matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    #matplotlib.rc('text', usetex=True)

    if sep == 'None':

        print 'Evaluating', props, 'for', fluid, 'using the', model, 'model at', temp, 'K for pressures', press
        
        # Define plots
        fig1 = plt.figure(figsize=(9,6))
        ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])
        fig2 = plt.figure(figsize=(9,6))
        ax2  = fig2.add_axes([0.1,0.15,0.8,0.75])
        fig3 = plt.figure(figsize=(9,6))
        ax3  = fig3.add_axes([0.1,0.15,0.8,0.75])
        fig4 = plt.figure(figsize=(9,6))
        ax4  = fig4.add_axes([0.1,0.15,0.8,0.75])
        fig5 = plt.figure(figsize=(9,6))
        ax5  = fig5.add_axes([0.1,0.15,0.8,0.75])
        fig6 = plt.figure(figsize=(9,6))
        ax6  = fig6.add_axes([0.1,0.15,0.8,0.75])


        markers = ['D', 's', 'v', '^', 'd', '*']
        legend_names = {'spce':'spce', 'tip4p':'TIP4P/2005', 'TraPPE': 'TraPPE', 'SAFT': 'SAFT Dimer', 'SAFT1': 'SAFT Monomer'}

        try:
            NIST = read_nist(temp, fluid)

            ax1.plot(NIST[1], NIST[0], linestyle = 'dashed', marker='+', label='NIST', c='k')

            ax2.plot(NIST[0], NIST[2], linestyle = 'dashed', marker='+', label='NIST', c='k')

            ax4.plot(NIST[0], NIST[3], linestyle = 'dashed', label='NIST', c='k')
        except IOError:
            print 'No NIST data available'
        except:
            print 'Some error. Check' 

        count = 0
        name_plot = ''
        timestamp = 8000000
        deltat = 6000000
        for m in model:
            name_plot += m
            press_dat = []
            temp_dat =[]
            rho_dat_init = []
            shear_dat = []
            ratio_dat = []
            bulk_dat = []
            diff_dat = []
            for p in press:
                for i in den:
                    try:
                        results = bulk_properties(fluid,m,temp,p,i)

                        press_dat.append(results.press())
                        temp_dat.append(results.temp())
                        rho_dat_init.append(results.rho())
                        shear_dat.append(results.shear())
                        bulk_dat.append(results.bulk())
                        ratio_dat.append(results.visc_ratio())
                        diff_dat.append(results.diff(timestamp,deltat)[-1])
                    except:
                        print 'There was an error. Maybe no data for %s and P=%s exists.'%(m,p)

            
    
            rho_dat, press_dat, press_err   = averaging(rho_dat_init, press_dat)
            rho_dat, temp_dat, temp_err     = averaging(rho_dat_init, temp_dat)
            rho_dat, shear_dat, shear_err   = averaging(rho_dat_init, shear_dat)
            rho_dat, bulk_dat, bulk_err     = averaging(rho_dat_init, bulk_dat)
            rho_dat, ratio_dat, ratio_err   = averaging(rho_dat_init, ratio_dat)
            rho_dat, diff_dat, diff_err     = averaging(rho_dat_init, diff_dat)
            

            print 'Diffusion for %s' % m, np.mean(diff_dat)


            # Plotting
            ax1.errorbar(press_dat, rho_dat, xerr=press_err, linestyle = 'None', marker=markers[count], label='%s'%(legend_names[m]))
            ax1.set_xlabel('P (bar)')
            ax1.set_ylabel('$\\rho$ (g/cm$^3$)')
            ax1.legend(loc='upper left')

            ax2.errorbar(rho_dat, shear_dat, yerr=shear_err, linestyle = 'None', marker=markers[count], label='%s'%(legend_names[m]))
            ax2.set_xlabel('$\\rho$ (g/cm$^3$)')
            ax2.set_ylabel('$\eta$ (mPa.s)')
            ax2.legend(loc='upper left')

            ax3.errorbar(rho_dat, bulk_dat, yerr=bulk_err,linestyle = 'None', marker=markers[count], label='%s'%(legend_names[m]))
            ax3.set_xlabel('$\\rho$ (g/cm$^3$)')
            ax3.set_ylabel('$\kappa$ (mPa.s)')
            ax3.legend(loc='upper left')

            ax4.errorbar(rho_dat, diff_dat, yerr=diff_err,linestyle = 'None', marker=markers[count], label='%s'%(legend_names[m]))
            #ax4.set_ylim(0,100)
            ax4.set_xlabel('$\\rho$ (g/cm$^3$)')
            ax4.set_ylabel('$D_s$ ($10^{-9}m^2/s$)')
            ax4.legend()
            #ax4_inset=fig4.add_axes([0.22,0.22, 0.34,0.34])
            #ax4_inset.plot(rho_dat, diff_dat, linestyle = 'None', marker=markers[count], label='%s'%(legend_names[m]))
            #ax4_inset.set_xlim(0,0.12)
            #ax4_inset.set_ylim(100,)

            ax5.errorbar(rho_dat, ratio_dat, yerr=ratio_err, linestyle = 'None', marker=markers[count], label='%s'%(legend_names[m]))
            ax5.set_xlabel('$\\rho$ (g/cm$^3$)')
            ax5.set_ylabel('$\kappa$/$\eta$')
            ax5.legend(loc='upper left')

            ax6.errorbar(rho_dat, temp_dat, yerr=temp_err,linestyle = 'None', marker=markers[count], label='%s'%(legend_names[m]))
            ax6.set_xlabel('$\\rho$ (g/cm$^3$)')
            ax6.set_ylabel('T (K)')
            ax6.legend(loc='upper left')

            count +=1


        fig1.savefig('PLOTS/P_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig2.savefig('PLOTS/SV_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig3.savefig('PLOTS/BV_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig4.savefig('PLOTS/Ds_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig5.savefig('PLOTS/BS_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig6.savefig('PLOTS/T_%s_T%s_%s.pdf'%(fluid, temp,name_plot))

    elif press == 'None':

        print 'Evaluating', props, 'for', fluid, 'using the', model, 'model at', temp, 'K for separations', sep
        
        # Define plots
        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()
        fig3, ax3 = plt.subplots()
        fig4, ax4 = plt.subplots()
        fig5, ax5 = plt.subplots()
        fig6, ax6 = plt.subplots()
        fig7, ax7 = plt.subplots()
        fig8, ax8 = plt.subplots()
        fig9, ax9 = plt.subplots()
        fig10, ax10 = plt.subplots()
        fig11, ax11 = plt.subplots()

        markers = ['D', 's', 'v', '^']

        try:
            NIST = read_nist(temp, fluid)

            ax1.plot(NIST[1], NIST[0], linestyle = 'dashed', marker='+', label='NIST', c='k')

            ax2.plot(NIST[0], NIST[2], linestyle = 'dashed', marker='+', label='NIST', c='k')

            ax4.plot(NIST[0], NIST[3], linestyle = 'dashed', label='NIST', c='k')
        except IOError:
            print 'No NIST data available'
        except:
            print 'Some error. Check' 

        count = 0
        name_plot = ''
        for m in model:
            #name_plot += m
            sep_dat = []
            temp_dat =[]
            rho_dat = []
            shear_dat = []
            ratio_dat = []
            #shear2_dat = []
            bulk_dat = []
            diff_dat = []
            wa_dat = []
            for e in eps:
                name_plot += e
                sep_dat = []
                temp_dat =[]
                rho_dat = []
                shear_dat = []
                ratio_dat = []
                #shear2_dat = []
                bulk_dat = []
                diff_dat = []
                wa_dat = []
                for s in sep:
                    try:
                        results = confined_properties(fluid,m,temp,s, e)
                        timestamp = 800000
                        deltat = 600000
                        
                        sep_dat.append(results.separation())
                        temp_dat.append(results.temp())
                        rho_dat.append(results.rho())
                        shear_dat.append(results.shear())
                        bulk_dat.append(results.bulk())
                        ratio_dat.append(results.visc_ratio())
                        diff_dat.append(results.diff(timestamp,deltat)[-1])
                        wa_dat.append(results.wa())
                        #shear2_dat.append(results.shear2(float(temp),timestamp,deltat))
                    except:
                        print 'There was an error. Maybe no data for %s and s=%s exists.'%(m,s)

                #print sep_dat, rho_dat, wa_dat, temp_dat
                # Plotting
                ax1.plot(sep_dat, rho_dat, linestyle = 'None', marker=markers[count], label='$\epsilon$ = %s'%(e))
                ax1.legend()

                ax2.plot(sep_dat, shear_dat, linestyle = 'None', marker=markers[count], label='$\epsilon$ = %s'%(e))
                #ax2.plot(sep_dat, shear2_dat, linestyle = 'None', marker=markers[count], c='g', label='%s'%(m))
                ax2.set_xlabel('$\Delta z$ $(\AA)$')
                #ax2.set_xlim(0,18)
                ax2.set_ylabel('$\eta$ (mPa.s)')
                ax2.legend(loc='upper right')

                ax3.plot(sep_dat, bulk_dat, linestyle = 'None', marker=markers[count], label='$\epsilon$ = %s'%(e))
                ax3.set_xlabel('$\Delta z$ $(\AA)$')
                ax3.set_ylabel('$\kappa$ (mPa.s)')
                ax3.legend(loc='upper left')

                ax4.plot(sep_dat, diff_dat, linestyle = 'None', marker=markers[count], label='$\epsilon$ = %s'%(e))
                #ax4.set_xlim(0,18)
                ax4.set_xlabel('$\Delta z$ $(\AA)$')
                ax4.set_ylabel('$D_s$ ($10^{-9}m^2/s$)')
                ax4.legend(loc='upper left')

                ax5.plot(rho_dat, shear_dat, linestyle = 'None', marker=markers[count], label='$\epsilon$ = %s'%(e))
                ax5.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax5.set_ylabel('$\eta$ (mPa.s)')
                ax5.legend(loc='upper left')

                ax6.plot(rho_dat, bulk_dat, linestyle = 'None', marker=markers[count], label='$\epsilon$ = %s'%(e))
                ax6.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax6.set_ylabel('$\kappa$ (mPa.s)')
                ax6.legend(loc='upper left')

                ax7.plot(rho_dat, diff_dat, linestyle = 'None', marker=markers[count], label='$\epsilon$ = %s'%(e))
                #ax7.set_ylim(0,100)
                ax7.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax7.set_ylabel('$D_s$ ($10^{-9}m^2/s$)')
                ax7.legend()

                ax8.plot(sep_dat, ratio_dat, linestyle = 'None', marker=markers[count], label='$\epsilon$ = %s'%(e))
                ax8.set_xlabel('$\Delta z$ $(\AA)$')
                ax8.set_ylabel('$\kappa$/$\eta$')
                ax8.legend(loc='upper left')

                ax9.plot(rho_dat, ratio_dat, linestyle = 'None', marker=markers[count], label='$\epsilon$ = %s'%(e))
                ax9.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax9.set_ylabel('$\kappa$/$\eta$')
                ax9.legend(loc='upper left')

                ax10.plot(sep_dat, wa_dat, linestyle = 'None', marker=markers[count], label='$\epsilon$ = %s'%(e))
                ax10.set_xlabel('$\Delta z$ $(\AA)$')
                ax10.set_ylabel('$W_A$')
                ax10.legend(loc='upper left')

                ax11.plot(rho_dat, wa_dat, linestyle = 'None', marker=markers[count], label='$\epsilon$ = %s'%(e))
                ax11.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax11.set_ylabel('$W_A$')
                ax11.legend(loc='upper left')

                count +=1


        fig1.savefig('PLOTS/Graph_S_rho_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig2.savefig('PLOTS/Graph_SV_s_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig3.savefig('PLOTS/Graph_BV_s_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig4.savefig('PLOTS/Graph_Ds_s_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig5.savefig('PLOTS/Graph_SV_rho_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig6.savefig('PLOTS/Graph_BV_rho_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig7.savefig('PLOTS/Graph_Ds_rho_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig8.savefig('PLOTS/Graph_BS_s_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig9.savefig('PLOTS/Graph_BS_rho_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig10.savefig('PLOTS/Graph_WA_s_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
        fig11.savefig('PLOTS/Graph_WA_rho_%s_T%s_%s.pdf'%(fluid, temp,name_plot))


    return




if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS')
    main()
