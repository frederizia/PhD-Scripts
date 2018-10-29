#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import argparse
from class_bulk_props import *
import sys
import csv
from itertools import izip,imap

def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-T', '--temp', required=False, default='298', action='store',
                       help='Temp')
    parser.add_argument('-P', '--press', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-s', '--sep', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-f', '--fluid', required=False, default='Water', action='store')
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

    k = map(float, averages.keys())
    v = averages.values()
    ke = map(float, errors.keys())
    e = errors.values()
    sorti = np.argsort(k)
    sortj = np.argsort(ke)
    k, v, e = np.array(k)[sorti], np.array(v)[sorti], np.array(e)[sortj]
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
    sep = args.sep
    eps = args.eps
    den = args.den
    rhos = args.rhos[0]
    png = args.png[0]
    fig_size = (9,7)

    colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027','#a50026',\
    '#313695','#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']


    if png == 'y':
        EXT = 'png'
    else:
        EXT = 'pdf'


    if sep == 'None':

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
        fig8a = plt.figure(figsize=fig_size)
        ax8a  = fig8a.add_axes([0.1,0.15,0.8,0.75])
        fig9 = plt.figure(figsize=fig_size)
        ax9  = fig9.add_axes([0.1,0.15,0.8,0.75])
        fig10 = plt.figure(figsize=fig_size)
        ax10  = fig10.add_axes([0.1,0.15,0.8,0.75])

        if fluid == 'CO2':
            saft_term = 'dimer'
            #Kvib = kappa_vib(temp)
            #print Kvib
            
        elif fluid== 'Decane':
            saft_term = 'three-bead'
        else:
            saft_term = ''
        markers = ['D', 'd', 'v', '^', 's', '*','d', 'v']
        legend_names = {'spce':'SPC/E', 'tip4p':'TIP4P/2005', 'TraPPE': 'TraPPE', 'SAFT': 'SAFT %s'%saft_term, \
        'SAFT1': 'SAFT monomer', 'TraPPEnc': 'TraPPE (no charges)', 'OPLS': 'OPLS',\
        'SAFT1_rc2.5': 'SAFT (rc=2.5)', '12_6':'12-6 model', '12_6_rcsaft': '12-6 model (rc from SAFT)',\
        '12_6_23_666': '12-6 model (SAFT exponents)',
        'SAFTflex': 'SAFT Dimer (flexible II)', 'EPM2':'EPM2', 'SAFTrigid': 'SAFT (rigid)',\
        'SAFT1vle': 'SAFT-vle', 'SAFT1ift': 'SAFT-ift', 'SAFTflex2': 'SAFT dimer (flexible)',\
        'OPLSmod': 'L-OPLS', 'TraPPEflex': 'TraPPE (flexible)',\
        'EPM2rigid': 'EPM2', 'EPM2flex': 'EPM2 (flexible)', 'TraPPErigid': 'TraPPE',\
        'EPM2nc': 'EPM2 (no charges)', 'spcenc': 'SPC/E (no charges)', 'EPM2vib': 'EPM2 (fully flexible)',\
        'EPM2angle1': 'EPM2 ($k_{\\theta}=0.5$)','EPM2angle2': 'EPM2 ($k_{\\theta}=1.0$)','EPM2angle3': 'EPM2 ($k_{\\theta}=3.0$)',\
        'EPM2angle4': 'EPM2 ($k_{\\theta}=6.405$)','EPM2bond1': 'EPM2 ($k_{b}=0.5$)','EPM2bond2': 'EPM2 ($k_{b}=1.0$)',\
        'EPM2bond3': 'EPM2 ($k_{b}=3.0$)','EPM2bond4': 'EPM2 ($k_{b}=6.405$)',\
        'OPLSrigid': 'OPLS (rigid)','OPLSmodrigid': 'L-OPLS (rigid)', 'spcerattle': 'SPC/E (RATTLE)'}

        try:
            NIST = read_nist(temp, fluid)
            nist = 0

            ax1.plot(NIST[1], NIST[0], linestyle = 'dashed', marker='None', label='NIST', c='k')
            #shear, density
            #ax2.plot(NIST[0], NIST[2], linestyle = 'dashed', marker='None', label='NIST', c='k')
            # shear, press
            ax7.plot(NIST[1], NIST[2], linestyle = 'dashed', marker='None', label='NIST', c='k')

            #ax4.plot(NIST[0], NIST[3], linestyle = 'dashed', label='NIST', c='k')
        except IOError:
            print 'No NIST data available'
            nist = 1
        except:
            nist = 1
            print 'Some error. Check' 

        count = 0
        name_plot = ''
        timestamp = 8000000
        deltat = 6000000
        steps = 3000

        kVib, kRot = kappa_int_est(temp,fluid)
        print kVib, kRot
        for m in model:
            print '\n\n++++++++++++++++++++++++ %s +++++++++++++++++++++\n' %m
            name_plot += m
            press_dat = []
            temp_dat =[]
            rho_dat_init = []
            shear_dat = []
            shear_msd_dat = []
            shear_vacf_dat = []
            shear_err = []
            ratio_dat = []
            ratio_plus_dat = []
            bulk_dat = []
            bulk_err = []
            bulk_plus_dat = []
            diff_dat = []
            diff_vacf_dat = []
            for p in press:
                for i in den:
                    print 'Evaluating %s, DEN=%s and P=%s.'%(m,i,p)
                    try:
                        results = bulk_properties(fluid,m,temp,p,i)
                        Pval = results.press()
                        press_dat.append(Pval)
                        temp_dat.append(results.temp())
                        rho_dat_init.append(results.rho())
                        shear_dat.append(results.shear())
                        shear_err.append(results.shear2()[1])
                        shear_msd_dat.append(results.shear_diff('msd', 'None'))
                        bulk_dat.append(results.bulk())
                        bulk_err.append(results.bulk2()[1])
                        ratio_dat.append(results.visc_ratio())
                        diff_dat.append(results.diff()[0]*1e9)#diff(timestamp,deltat)[-1])
                        if m=='SAFT1':
                            bulk_plus_dat.append(results.bulk2()[0]+(kVib+kRot))
                            ratio_plus_dat.append((results.bulk2()[0]+(kVib+kRot))/results.shear2()[0])
                        else:
                            bulk_plus_dat.append(results.bulk2()[0]+kVib)
                            ratio_plus_dat.append((results.bulk2()[0]+kVib)/results.shear2()[0])
                        if (m=='tip4p' or m=='spce') and temp=='300' and i=='5':
                            print 'Doing extra evaluation...'
                            shear_vacf_dat.append(results.shear_diff('vacf', 1000))
                            diff_vacf_dat.append(results.diff2(steps=1000)*1e9)
                        elif m=='spce' and temp=='298':
                            print 'Doing extra evaluation...'
                            shear_vacf_dat.append(results.shear_diff('vacf', steps))
                            diff_vacf_dat.append(results.diff2(steps=steps)*1e9)
                        #print results.rho(), results.bulk2()[0], results.bulk2()[1]
                        print 'Data collection successful.\n'
                    except:
                        print 'ERROR: Maybe no data for %s, DEN=%s and P=%s exists.\n'%(m,i,p)

            rho_dat, press_dat, press_err   = averaging(rho_dat_init, press_dat)
            rho_dat, temp_dat, temp_err     = averaging(rho_dat_init, temp_dat)
            rho_dat, shear_msd_dat, shear_msd_err   = averaging(rho_dat_init, shear_msd_dat)
            if len(den)==1:
                rho_dat, shear_dat, shear_err   = averaging3(rho_dat_init, shear_dat, shear_err)
                rho_dat, bulk_dat, bulk_err     = averaging3(rho_dat_init, bulk_dat, bulk_err)
            else:
                rho_dat, shear_dat, shear_err   = averaging(rho_dat_init, shear_dat)#, shear_err)
                rho_dat, bulk_dat, bulk_err     = averaging(rho_dat_init, bulk_dat)#, bulk_err)
            rho_dat, ratio_dat, ratio_err   = averaging(rho_dat_init, ratio_dat)
            rho_dat, diff_dat, diff_err     = averaging(rho_dat_init, diff_dat)
            rho_dat, bulk_plus_dat, bulk_plus_err   = averaging(rho_dat_init, bulk_plus_dat)
            rho_dat, ratio_plus_dat, ratio_plus_err   = averaging(rho_dat_init, ratio_plus_dat)
            
            if m=='spce' and temp=='298':
                rho_dat, diff_vacf_dat, diff_vacf_err     = averaging(rho_dat_init, diff_vacf_dat)
                rho_dat, shear_vacf_dat, shear_vacf_err     = averaging(rho_dat_init, shear_vacf_dat)

            if fluid=='CO2':
                rho_gas = 0.275
                rho_liquid = 0.679
                test = np.ma.masked_where((rho_dat > rho_gas) & (rho_dat < rho_liquid), rho_dat)
                co2_tp_mask = test.mask
                #if nist == 0:
                #    test2 = np.ma.masked_where((NIST[0] > rho_gas) & (NIST[0] < rho_liquid), NIST[0])
                #    co2_tp_mask2 = test2.mask
            #if (m=='tip4p' or m=='spce') and temp=='300' and shear_vacf_dat != []:
            #    rho_dat, shear_vacf_dat, shear_vacf_err   = averaging(rho_dat_init, shear_vacf_dat)#, shear_err)
            #    rho_dat, diff_vacf_dat, diff_vacf_err   = averaging(rho_dat_init, diff_vacf_dat)
            
            print 'Diffusion for %s' % m,'at P=1bar:', diff_dat[0], '+/-', diff_err[0]
            print 'Shear viscosity for %s' % m,'at P=1bar:', shear_dat[0], '+/-', shear_err[0]
            print 'Bulk viscosity for %s' % m,'at P=1bar:', bulk_dat[0], '+/-', bulk_err[0]
            print 'Viscosity ratio for %s' % m,'at P=1bar:', ratio_dat[0], '+/-', ratio_err[0]
            print 'Final viscosity ratio for %s' % m, ratio_dat[-1], '+/-', ratio_err[-1]
            
            print 'rho'
            print rho_dat
            print 'Ds (MSD)'
            print diff_dat
            print 'Ds (VACF)'
            print diff_vacf_dat
            print 'kappa/eta (plus):'
            print ratio_plus_dat, ratio_plus_err
            print 'kappa/eta:'
            print ratio_dat, ratio_err
            print 'Pressure:'
            print press_dat
            print 'eta:'
            print shear_dat, shear_err
            print 'eta from Ds (MSD):'
            print shear_msd_dat
            print 'eta from Ds (VACF):'
            print shear_vacf_dat
            print 'kappa:'
            print bulk_plus_dat, bulk_plus_err
            print "-------------------------------"
            print 'Average diffusion for %s' % m, np.mean(diff_dat), '+/-', np.std(np.array(diff_dat))/np.sqrt(len(diff_dat))
            print 'Average shear viscosity for %s' % m, np.mean(shear_dat), '+/-', np.std(np.array(shear_dat))/np.sqrt(len(shear_dat))
            print 'Average bulk viscosity for %s' % m, np.mean(bulk_plus_dat), '+/-', np.std(np.array(bulk_plus_dat))/np.sqrt(len(bulk_plus_dat))
            print 'Average viscosity ratio for %s' % m, np.mean(ratio_plus_dat), '+/-', np.std(np.array(ratio_plus_dat))/np.sqrt(len(ratio_plus_dat))


            # SAVE DATA
            print #--------------Saving data ....-----------------#
            with open('DATA/viscosity_{}_{}_T{}.csv'.format(fluid,m,temp), 'wb') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(("density", "pressure", "shear","bulk", "ratio"))
                if temp=='323' and fluid=='CO2':
                    writer.writerows(imap(lambda a,b,c,d,e,f,g,h: (round(a,3), round(b,1), "{} \pm {}".format(round(c,3), round(d,3)), \
                        "{} \pm {}".format(round(e,3), round(f,3)),"{} \pm {}".format(round(g,3), round(h,3))), \
                        rho_dat, press_dat, shear_dat, shear_err, bulk_dat, bulk_err, ratio_dat, ratio_err))#izip(dZ,round(Hlist,2)))
                elif temp=='613' and fluid=='Water':
                    print 'Saving file.'
                    writer.writerows(imap(lambda a,b,c,d,e,f,g,h: (round(a,4), round(b,1), "{} \pm {}".format(round(c,3), round(d,3)), \
                        "{} \pm {}".format(round(e,3), round(f,3)),"{} \pm {}".format(round(g,3), round(h,3))), \
                        rho_dat, press_dat, shear_dat, shear_err, bulk_dat, bulk_err, ratio_dat, ratio_err))
                elif temp=='393' and fluid=='Water':
                    print 'Saving file.'
                    writer.writerows(imap(lambda a,b,c,d,e,f,g,h: (round(a,3), round(b,1), "{} \pm {}".format(round(c,3), round(d,3)), \
                        "{} \pm {}".format(round(e,3), round(f,3)),"{} \pm {}".format(round(g,3), round(h,3))), \
                        rho_dat, press_dat, shear_dat, shear_err, bulk_dat, bulk_err, ratio_dat, ratio_err))
                else:    
                    writer.writerows(imap(lambda a,b,c,d,e,f,g,h: (round(a,3), round(b,1), "{} \pm {}".format(round(c,3), round(d,3)), \
                        "{} \pm {}".format(round(e,3), round(f,3)),"{} \pm {}".format(round(g,3), round(h,3))), \
                        rho_dat, press_dat, shear_dat, shear_err, bulk_plus_dat, bulk_plus_err, ratio_plus_dat, ratio_plus_err))#izip(dZ,round(Hlist,2)))



            # Plotting
            if m=='SAFT1' and fluid=='Water':
                if temp=='300':
                    label_tmp = 'SAFT CGW1-ift'
                    modnum = 6
                else:
                    label_tmp = 'SAFT CGW1-vle'
                    modnum = 7
            elif m=='SAFT1vle':
                label_tmp = 'SAFT CGW1-vle'
                modnum = 7
            elif m=='SAFT1ift':
                label_tmp = 'SAFT CGW1-ift'
                modnum = 6
            elif m=='spce':
                label_tmp = legend_names[m]
                modnum = 0
            elif m=='tip4p':
                label_tmp = legend_names[m]
                modnum = 4
            else:
                label_tmp = legend_names[m]
                modnum = 1

            if fluid=='Water':
                ax1.errorbar(press_dat, rho_dat, xerr=press_err, linestyle = 'None', marker=markers[modnum], c=colours[modnum], label='%s'%(label_tmp))
            else:
                ax1.errorbar(press_dat, rho_dat, xerr=press_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(label_tmp))
            ax1.set_xlabel('P (bar)')
            ax1.set_ylabel('$\\rho$ (g/cm$^3$)')
            ax1.legend()

            if fluid=='CO2' and temp=='300':
                if count==0:
                    nist_rho = np.array(NIST[0])
                    nist_shear = np.array(NIST[2])
                    lower_mask = np.where(nist_rho<rho_gas)
                    upper_mask = np.where(nist_rho>rho_liquid)

                    ax2.plot(nist_rho[lower_mask], nist_shear[lower_mask], linestyle = 'dashed', marker='None', label='NIST', c='k')
                    ax2.plot(nist_rho[upper_mask], nist_shear[upper_mask], linestyle = 'dashed', marker='None', c='k')
                ax2.errorbar(rho_dat[co2_tp_mask], shear_dat[co2_tp_mask], yerr=shear_err[co2_tp_mask], linestyle = 'None', marker=markers[count], markerfacecolor='None', c=colours[2*count])
                ax2.errorbar(rho_dat[~co2_tp_mask], shear_dat[~co2_tp_mask], yerr=shear_err[~co2_tp_mask], linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            else:
                if nist == 0 and count==0:
                    ax2.plot(NIST[0], NIST[2], linestyle = 'dashed', marker='None', label='NIST', c='k')
                ax2.errorbar(rho_dat, shear_dat, yerr=shear_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            ax2.set_xlabel('$\\rho$ (g/cm$^3$)')
            ax2.set_ylabel('$\eta$ (mPa.s)')
            if m=='EPM2rigid' and temp=='300':
                #rho_gas = 0.275
                #rho_liquid = 0.679
                ax2.plot([rho_gas,rho_gas], [0, 0.26], linestyle='dotted',c='k', lw=1.5)
                ax2.plot([rho_liquid,rho_liquid], [0, 0.26], linestyle='dotted',c='k', lw=1.5)
                ax2.plot([rho_gas,rho_liquid],[0.04,0.04], linestyle='dashed',c='k', lw=1.5)
                ax2.arrow(rho_gas, 0.04, -0.1, 0, head_width=0.007, head_length=0.02, fc='k', ec='k')
                ax2.arrow(rho_liquid, 0.04, 0.1, 0, head_width=0.007, head_length=0.02, fc='k', ec='k')
                ax2.text(rho_gas-0.1,0.05,'gas')
                ax2.text(rho_gas+0.05,0.05,'two-phase')
                ax2.text(rho_liquid+0.01,0.02,'liquid')
                ax2.set_ylim(0,0.26)
            ax2.legend()

            

            ax3.errorbar(rho_dat, bulk_dat, yerr=bulk_err,linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            ax3.set_xlabel('$\\rho$ (g/cm$^3$)')
            ax3.set_ylabel('$\kappa_{\mathrm{conf}}$ (mPa.s)')
            if m=='EPM2rigid':
                rho_min = np.min(rho_dat)
                rho_max = np.max(rho_dat)
                rho_fit, bulk_fit, params = poly_fit(rho_dat, bulk_dat, rho_min, rho_max)
                ax3.plot(rho_fit, bulk_fit, linestyle='dotted', lw=3, c='k')#, label='$\\rho^2$ fit')
            #if fluid == 'CO2' and 'TraPPE' in model:
                #ax3.set_yscale("log", nonposy='clip')
                #ax3.set_ylim(ymin=0.1)
                #ax3.set_xlim(0.78,1.1)
                #ax3_inset=fig3.add_axes([0.67,0.3, 0.2,0.2])
                #ax3_inset.plot(rho_dat, bulk_dat, linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
                #ax3_inset.set_xlim(0,0.35)
                #ax3_inset.set_ylim(0,1)
            ax3.legend()

            ax4.errorbar(rho_dat, diff_dat, yerr=diff_err,linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            #ax4.set_ylim(0,100)
            ax4.set_xlabel('$\\rho$ (g/cm$^3$)')
            ax4.set_ylabel('$D_s$ ($10^{-9}m^2/s$)')
            ax4.legend()
            #ax4_inset=fig4.add_axes([0.22,0.22, 0.34,0.34])
            #ax4_inset.plot(rho_dat, diff_dat, linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            #ax4_inset.set_xlim(0,0.12)
            #ax4_inset.set_ylim(100,)

            if len(model) == 1:
                fig2a = plt.figure(figsize=fig_size)
                ax2a  = fig2a.add_axes([0.1,0.15,0.8,0.75])
                ax2a.errorbar(rho_dat, shear_dat, yerr=shear_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label='GK')
                ax2a.plot(rho_dat, shear_msd_dat, linestyle = 'None', marker=markers[count+1], c=colours[2*count+1], label='$D_s$ (MSD)')
                ax2a.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax2a.set_ylabel('$\eta$ (mPa.s)')

                
                if ((m=='tip4p' or m=='spce') and temp=='300' and shear_vacf_dat != []) or (m=='spce' and temp=='298' and len(model)==1):
                    ax2a.plot(rho_dat, shear_vacf_dat, linestyle = 'None', marker=markers[count+2], c=colours[2*count+2], label='$D_s$ (VACF)')
                    
                    # straight fit for both VACF and GK eta
                    xfit_gk, shear_fit_gk, slope_gk, slope_gk_err = straight_fit(1e9/np.array(diff_vacf_dat), 1e-3*np.array(shear_dat), 1e9/np.max(diff_vacf_dat), 1e9/np.min(diff_vacf_dat))
                    xfit_vacf, shear_fit_vacf, slope_vacf, slope_vacf_err = straight_fit(1e9/np.array(diff_vacf_dat), 1e-3*np.array(shear_vacf_dat), 1e9/np.max(diff_vacf_dat), 1e9/np.min(diff_vacf_dat))

                    # effect diameter
                    kB = 1.38*1e-23
                    T=int(temp)
                    convert = 1e10
                    factor = (kB*T)/(3*np.pi)

                    alpha_vacf = convert*factor/slope_vacf # should be 1.7
                    alpha_vacf_err = convert*factor*(-1)*(slope_vacf_err/slope_vacf**2)
                    alpha_gk = convert*factor/slope_gk
                    alpha_gk_err = convert*factor*(-1)*(slope_gk_err/slope_gk**2)
                    xfit_gk, shear_fit_gk = 1e-9*np.array(xfit_gk), 1e3*np.array(shear_fit_gk)
                    xfit_vacf, shear_fit_vacf = 1e-9*np.array(xfit_vacf), 1e3*np.array(shear_fit_vacf)
                    print 'The effective diameters are:', alpha_vacf, '+/-',alpha_vacf_err, 'A (VACF) and', alpha_gk, '+/-',alpha_gk_err, 'A (GK)'

                    fig4a = plt.figure(figsize=fig_size)
                    ax4a  = fig4a.add_axes([0.1,0.15,0.8,0.75])
                    
                    ax4a.plot(1/np.array(diff_vacf_dat), shear_vacf_dat,linestyle = 'None', marker=markers[count+1], c=colours[2*count+1], label='Stokes-Einstein')
                    ax4a.plot(1/np.array(diff_vacf_dat), shear_dat,linestyle = 'None', marker=markers[count], c=colours[2*count], label='Green-Kubo')
                    ax4a.plot(xfit_vacf, shear_fit_vacf, linestyle='dashed',c=colours[2*count+1], label='fit (SE)')
                    ax4a.plot(xfit_gk, shear_fit_gk, linestyle='dashed',c=colours[2*count],label='fit (GK)')
                    if m=='spce' and temp=='298':
                        xfit_msd, shear_fit_msd, slope_msd, slope_msd_err = straight_fit(1e9/np.array(diff_dat), 1e-3*np.array(shear_msd_dat), 1e9/np.max(diff_dat), 1e9/np.min(diff_dat))
                        xfit_msd, shear_fit_msd = 1e-9*np.array(xfit_msd), 1e3*np.array(shear_fit_msd)
                        print 'Average shear viscosity (MSD) for %s' % m, np.mean(shear_msd_dat), '+/-', np.std(np.array(shear_msd_dat))/np.sqrt(len(shear_msd_dat)), '\n'
                        ax4a.plot(1/np.array(diff_dat), shear_msd_dat,linestyle = 'None', marker=markers[count+1], c=colours[2*count+2], label='Stokes-Einstein (MSD)')
                        ax4a.plot(xfit_msd, shear_fit_msd, linestyle='dashed',c=colours[2*count+2], label='fit (SE) (MSD)')
                    ax4a.set_xlabel('$1/D_s$ ($10^{9}$s/m$^2$)')
                    ax4a.set_ylabel('$\eta$ (mPa.s)')
                    ax4a.legend()
                    fig4a.savefig('PLOTS/Stokes_Einstein_%s_T%s_%s.pdf'%(fluid, temp,name_plot),bbox_inches='tight')
                    

                    fig4b = plt.figure(figsize=fig_size)
                    ax4b  = fig4b.add_axes([0.1,0.15,0.8,0.75])
                    ax4b.errorbar(rho_dat, diff_dat, yerr=diff_err,linestyle = 'None', marker=markers[count], c=colours[2*count], label='MSD')
                    if temp=='298':
                        ax4b.errorbar(rho_dat, diff_vacf_dat, yerr=diff_vacf_err,linestyle = 'None', marker=markers[count+1], c=colours[2*count+1], label='VACF')
                    else:
                        ax4b.plot(rho_dat, diff_vacf_dat,linestyle = 'None', marker=markers[count+1], c=colours[2*count+1], label='VACF')
                    #ax4.set_ylim(0,100)
                    ax4b.set_xlabel('$\\rho$ (g/cm$^3$)')
                    ax4b.set_ylabel('$D_s$ ($10^{-9}$m$^2$/s)')
                    ax4b.legend()
                    fig4b.savefig('PLOTS/Ds_comp_%s_T%s_%s.pdf'%(fluid, temp,name_plot),bbox_inches='tight')

                    
                    print 'Average shear viscosity (VACF) for %s' % m, np.mean(shear_vacf_dat), '+/-', np.std(np.array(shear_vacf_dat))/np.sqrt(len(shear_vacf_dat)), '\n'
                ax2a.legend()
                fig2a.savefig('PLOTS/SV_comp_%s_T%s_%s.pdf'%(fluid, temp,name_plot),bbox_inches='tight')

            if fluid=='CO2' and temp=='300':
                ax5.errorbar(rho_dat[co2_tp_mask], ratio_dat[co2_tp_mask], yerr=ratio_err[co2_tp_mask], linestyle = 'None', markerfacecolor='None', marker=markers[count], c=colours[2*count])
                ax5.errorbar(rho_dat[~co2_tp_mask], ratio_dat[~co2_tp_mask], yerr=ratio_err[~co2_tp_mask], linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            else:
                ax5.errorbar(rho_dat, ratio_dat, yerr=ratio_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            ax5.set_xlabel('$\\rho$ (g/cm$^3$)')
            ax5.set_ylabel('$\kappa_{\mathrm{conf}}$/$\eta$')
            if m=='EPM2rigid':
                #rho_gas = 0.25
                #rho_liquid = 0.679
                ylimupper = 8.8
                ax5.plot([rho_gas,rho_gas], [-0.4,ylimupper], linestyle='dotted',c='k', lw=1.5)
                ax5.plot([rho_gas,rho_liquid],[1.0, 1.0], linestyle='dashed',c='k', lw=1.5)
                ax5.plot([rho_liquid,rho_liquid], [-0.4,ylimupper], linestyle='dotted',c='k', lw=1.5)
                ax5.arrow(rho_gas, 1.0, -0.1, 0, head_width=0.16, head_length=0.02, fc='k', ec='k')
                ax5.arrow(rho_liquid, 1.0, 0.1, 0, head_width=0.16, head_length=0.02, fc='k', ec='k')
                ax5.text(rho_gas-0.09,1.25,'gas')
                ax5.text(rho_gas+0.11,1.25,'two-phase')
                ax5.text(rho_liquid+0.01,0.55,'liquid')
                ax5.set_ylim(-0.4,ylimupper)
                ax5.set_xlim(0,1.2)
            ax5.legend(loc='upper right')

            ax6.errorbar(rho_dat, temp_dat, yerr=temp_err,linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            ax6.set_xlabel('$\\rho$ (g/cm$^3$)')
            ax6.set_ylabel('T (K)')
            ax6.legend(loc='upper left')

            ax7.errorbar(press_dat, shear_dat, yerr=shear_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            ax7.set_xlabel('P (bar)')
            ax7.set_ylabel('$\eta$ (mPa.s)')
            ax7.legend()

            ax8.errorbar(press_dat, bulk_dat, yerr=bulk_err,linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            ax8.set_xlabel('P (bar)')
            ax8.set_ylabel('$\kappa_{\mathrm{conf}}$ (mPa.s)')
            ax8.legend()

            ax8a.errorbar(press_dat, bulk_plus_dat, yerr=bulk_plus_err,linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            ax8a.set_xlabel('P (bar)')
            ax8a.set_ylabel('$\kappa$ (mPa.s)')
            ax8a.legend()

            if fluid=='CO2' and temp=='300':
                ax9.errorbar(rho_dat[co2_tp_mask], bulk_plus_dat[co2_tp_mask], yerr=bulk_plus_err[co2_tp_mask], linestyle = 'None', marker=markers[count], c=colours[2*count], markerfacecolor='None')
                ax9.errorbar(rho_dat[~co2_tp_mask], bulk_plus_dat[~co2_tp_mask], yerr=bulk_plus_err[~co2_tp_mask], linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            else:
                ax9.errorbar(rho_dat, bulk_plus_dat, yerr=bulk_plus_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            if m=='EPM2rigid':
                bulk_plus_fit = np.array(bulk_fit)+np.min(bulk_plus_dat)
                ax9.plot(rho_fit, bulk_plus_fit, linestyle='dashed', lw=3, c='g')#, label='$\\rho^2$ fit')
                #rho_gas = 0.25
                #rho_liquid = 0.679
                ax9.plot([rho_gas,rho_gas], [32, 32.6], linestyle='dotted',c='k', lw=1.5)
                ax9.plot([rho_gas,rho_liquid],[32.05, 32.05], linestyle='dashed',c='k', lw=1.5)
                ax9.plot([rho_liquid,rho_liquid], [32, 32.6], linestyle='dotted',c='k', lw=1.5)
                ax9.arrow(rho_gas, 32.05, -0.1, 0, head_width=0.014, head_length=0.02, fc='k', ec='k')
                ax9.arrow(rho_liquid, 32.05, 0.1, 0, head_width=0.014, head_length=0.02, fc='k', ec='k')
                ax9.text(rho_gas-0.09,32.07,'gas')
                ax9.text(rho_gas+0.11,32.07,'two-phase')
                ax9.text(rho_liquid+0.01,32.07,'liquid')
                ax9.set_ylim(32, 32.6)

            ax9.set_xlabel('$\\rho$ (g/cm$^3$)')
            ax9.set_ylabel('$\kappa$ (mPa.s)')
            ax9.legend(loc='upper left')

            if fluid=='CO2' and temp=='300':
                ax10.errorbar(rho_dat[co2_tp_mask], ratio_plus_dat[co2_tp_mask], yerr=ratio_plus_err[co2_tp_mask], linestyle = 'None', marker=markers[count], c=colours[2*count], markerfacecolor='None')
                ax10.errorbar(rho_dat[~co2_tp_mask], ratio_plus_dat[~co2_tp_mask], yerr=ratio_plus_err[~co2_tp_mask], linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            else:
                ax10.errorbar(rho_dat, ratio_plus_dat, yerr=ratio_plus_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label='%s'%(legend_names[m]))
            ax10.set_xlabel('$\\rho$ (g/cm$^3$)')
            ax10.set_ylabel('$\kappa$/$\eta$')
            if m=='EPM2rigid':
                #ax9.plot(rho_fit, bulk_plus_fit, linestyle='dashed', lw=3, c='g')#, label='$\\rho^2$ fit')
                #rho_gas = 0.25
                #rho_liquid = 0.679
                ax10.plot([rho_gas,rho_gas], [0,4000], linestyle='dotted',c='k', lw=1.5)
                ax10.plot([rho_gas,rho_liquid],[630, 630], linestyle='dashed',c='k', lw=1.5)
                ax10.plot([rho_liquid,rho_liquid], [0,4000], linestyle='dotted',c='k', lw=1.5)
                ax10.arrow(rho_gas, 630, -0.1, 0, head_width=100, head_length=0.02, fc='k', ec='k')
                ax10.arrow(rho_liquid, 630, 0.1, 0, head_width=100, head_length=0.02, fc='k', ec='k')
                ax10.text(rho_gas-0.09,710,'gas')
                ax10.text(rho_gas+0.01,710,'two-phase')
                ax10.text(rho_liquid+0.01,710,'liquid')
                ax10.set_ylim(0,4000)
            ax10.legend()

            if count == 0 and fluid!='CO2':
                count+= 2
            else:
                count +=1

        if fluid=='Water' and (temp=='298' or temp=='300'):
            lower_bound = [2.7,2.7]
            upper_bound = [2.8,2.8]
            rho_min = np.min(rho_dat)
            rho_max = np.max(rho_dat)

            ax5.set_xlim(rho_min-0.001*rho_min, rho_max+0.001*rho_max)
            ax5.fill_between([rho_min-0.2*rho_min,rho_max+0.2*rho_max], lower_bound, upper_bound, facecolor='#ffd6cc', zorder=1)
            ax10.set_xlim(rho_min-0.001*rho_min, rho_max+0.001*rho_max)
            ax10.fill_between([rho_min-0.2*rho_min,rho_max+0.2*rho_max], lower_bound, upper_bound, facecolor='#ffd6cc', zorder=1)
        elif fluid=='Decane' and temp=='300':
            lower_bound = [3.23,3.23] #[3.1,3.1]
            upper_bound = [3.29,3.29] #[3.35,3.35]
            rho_min = np.min(rho_dat)
            rho_max = np.max(rho_dat)
            ax5.set_xlim(rho_min-0.001*rho_min, rho_max+0.001*rho_max)
            ax5.fill_between([rho_min-0.2*rho_min,rho_max+0.2*rho_max], lower_bound, upper_bound, facecolor='#ffd6cc', zorder=1)
            ax10.set_xlim(rho_min-0.001*rho_min, rho_max+0.001*rho_max)
            ax10.fill_between([rho_min-0.2*rho_min,rho_max+0.2*rho_max], lower_bound, upper_bound, facecolor='#ffd6cc', zorder=1)
        elif fluid=='CO2' and temp=='300':
            lower_bound = [843,843]
            upper_bound = [3849,3849]
            rho_min = np.min(rho_dat)
            rho_lim = np.max(rho_dat)
            rho_max = 0.15
            # add constant kint to plot
            rho_tmp = np.linspace(rho_min-0.1*rho_min, rho_max+0.01*rho_max,10)
            kInt_list_tmp = [(kVib+kRot)]*len(rho_tmp)
            kInt_list = [(kVib+kRot)]*len(NIST[0])
            kInt_ratio = np.divide(np.array(kInt_list),NIST[2])
            
            ax9.plot(rho_tmp, kInt_list_tmp, linestyle='dotted', lw=3, c='k')
            ax9.set_xlim(rho_min-0.1*rho_min, rho_lim+0.01*rho_lim)
            #ax10.plot(NIST[0], kInt_ratio, linestyle='dotted', lw=3, c='k')
            ax10.set_xlim(rho_min-0.001*rho_min, rho_lim+0.001*rho_lim)
            ax10.fill_between([rho_min-0.2*rho_min,rho_max+0.2*rho_max], lower_bound, upper_bound, facecolor='#ffd6cc', zorder=1)




        if png == 'y':
            fig1.savefig('PLOTS/PNG/P_%s_T%s_%s.png'%(fluid, temp,name_plot))
            fig2.savefig('PLOTS/PNG/SV_%s_T%s_%s.png'%(fluid, temp,name_plot))
            fig3.savefig('PLOTS/PNG/BV_%s_T%s_%s.png'%(fluid, temp,name_plot))
            fig4.savefig('PLOTS/PNG/Ds_%s_T%s_%s.png'%(fluid, temp,name_plot))
            fig5.savefig('PLOTS/PNG/BS_%s_T%s_%s.png'%(fluid, temp,name_plot))
            fig6.savefig('PLOTS/PNG/T_%s_T%s_%s.png'%(fluid, temp,name_plot))
            fig7.savefig('PLOTS/PNG/SV_vP_%s_T%s_%s.png'%(fluid, temp,name_plot))
            fig8.savefig('PLOTS/PNG/BV_vP_%s_T%s_%s.png'%(fluid, temp,name_plot))
            fig8a.savefig('PLOTS/PNG/BVplus_vP_%s_T%s_%s.png'%(fluid, temp,name_plot))
            fig9.savefig('PLOTS/PNG/BVplus_%s_T%s_%s.png'%(fluid, temp,name_plot))
            fig10.savefig('PLOTS/PNG/BSplus_%s_T%s_%s.png'%(fluid, temp,name_plot))
        elif png == 'eps':
            fig1.savefig('PLOTS/EPS/P_%s_T%s_%s.eps'%(fluid, temp,name_plot),format='eps', dpi=1000)
            fig2.savefig('PLOTS/EPS/SV_%s_T%s_%s.eps'%(fluid, temp,name_plot), format='eps', dpi=1000, bbox_inches='tight')
            fig3.savefig('PLOTS/EPS/BV_%s_T%s_%s.eps'%(fluid, temp,name_plot),format='eps', dpi=1000,bbox_inches='tight')
            fig4.savefig('PLOTS/EPS/Ds_%s_T%s_%s.eps'%(fluid, temp,name_plot),format='eps', dpi=1000)
            fig5.savefig('PLOTS/EPS/BS_%s_T%s_%s.eps'%(fluid, temp,name_plot),format='eps', dpi=1000,bbox_inches='tight')
            fig6.savefig('PLOTS/EPS/T_%s_T%s_%s.eps'%(fluid, temp,name_plot),format='eps', dpi=1000)
            fig7.savefig('PLOTS/EPS/SV_vP_%s_T%s_%s.eps'%(fluid, temp,name_plot),format='eps', dpi=1000)
            fig8.savefig('PLOTS/EPS/BV_vP_%s_T%s_%s.eps'%(fluid, temp,name_plot),format='eps', dpi=1000)
            fig8a.savefig('PLOTS/EPS/BVplus_vP_%s_T%s_%s.eps'%(fluid, temp,name_plot),format='eps', dpi=1000,bbox_inches='tight')
            fig9.savefig('PLOTS/EPS/BVplus_%s_T%s_%s.eps'%(fluid, temp,name_plot),format='eps', dpi=1000,bbox_inches='tight')
            fig10.savefig('PLOTS/EPS/BSplus_%s_T%s_%s.eps'%(fluid, temp,name_plot),format='eps', dpi=1000,bbox_inches='tight')
        else:
            fig1.savefig('PLOTS/P_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
            fig2.savefig('PLOTS/SV_%s_T%s_%s.pdf'%(fluid, temp,name_plot),bbox_inches='tight')
            fig3.savefig('PLOTS/BV_%s_T%s_%s.pdf'%(fluid, temp,name_plot),bbox_inches='tight')
            fig4.savefig('PLOTS/Ds_%s_T%s_%s.pdf'%(fluid, temp,name_plot),bbox_inches='tight')
            fig5.savefig('PLOTS/BS_%s_T%s_%s.pdf'%(fluid, temp,name_plot),bbox_inches='tight')
            fig6.savefig('PLOTS/T_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
            fig7.savefig('PLOTS/SV_vP_%s_T%s_%s.eps'%(fluid, temp,name_plot))
            fig8.savefig('PLOTS/BV_vP_%s_T%s_%s.pdf'%(fluid, temp,name_plot))
            fig8a.savefig('PLOTS/BVplus_vP_%s_T%s_%s.pdf'%(fluid, temp,name_plot),bbox_inches='tight')
            fig9.savefig('PLOTS/BVplus_%s_T%s_%s.pdf'%(fluid, temp,name_plot),bbox_inches='tight')
            fig10.savefig('PLOTS/BSplus_%s_T%s_%s.pdf'%(fluid, temp,name_plot),bbox_inches='tight')
        print 'The vibrational bulk viscosity contribution is:', kVib
        print 'The rotational bulk viscosity contribution is:', kRot

    elif press == 'None':
        if fluid=='Water_Graphene':
            extra_plot_name = 'Graph'
            fric_coords = ['x','y']
        elif fluid=='LJ_channel':
            extra_plot_name = 'LJ'
            fric_coords = ['x','z']
        else:
            print 'Unknown channel.'
            sys.exit(1)

        print 'Evaluating properties for', fluid, 'using the', model, 'model at', temp, 'K for separations', sep
        
        # Define plots
        fig1  = plt.figure(figsize=fig_size)
        ax1   = fig1.add_axes([0.15,0.15,0.75,0.75])
        fig2  = plt.figure(figsize=fig_size)
        ax2   = fig2.add_axes([0.15,0.15,0.75,0.75])
        fig3  = plt.figure(figsize=fig_size)
        ax3   = fig3.add_axes([0.15,0.15,0.75,0.75])
        fig4  = plt.figure(figsize=fig_size)
        ax4   = fig4.add_axes([0.15,0.15,0.75,0.75])
        fig5  = plt.figure(figsize=fig_size)
        ax5   = fig5.add_axes([0.15,0.15,0.75,0.75])
        fig6  = plt.figure(figsize=fig_size)
        ax6   = fig6.add_axes([0.15,0.15,0.75,0.75])
        fig7  = plt.figure(figsize=fig_size)
        ax7   = fig7.add_axes([0.15,0.15,0.75,0.75])
        fig8  = plt.figure(figsize=fig_size)
        ax8   = fig8.add_axes([0.15,0.15,0.75,0.75])
        fig9  = plt.figure(figsize=fig_size)
        ax9   = fig9.add_axes([0.15,0.15,0.75,0.75])
        fig10 = plt.figure(figsize=fig_size)
        ax10  = fig10.add_axes([0.15,0.15,0.75,0.75])
        fig11 = plt.figure(figsize=fig_size)
        ax11  = fig11.add_axes([0.15,0.15,0.75,0.75])
        fig12 = plt.figure(figsize=fig_size)
        ax12  = fig12.add_axes([0.15,0.15,0.75,0.75])
        fig13 = plt.figure(figsize=fig_size)
        ax13  = fig13.add_axes([0.15,0.15,0.75,0.75])
        fig14 = plt.figure(figsize=fig_size)
        ax14  = fig14.add_axes([0.15,0.15,0.75,0.75])
        fig15 = plt.figure(figsize=fig_size)
        ax15  = fig15.add_axes([0.15,0.15,0.75,0.75])

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
            name_plot += m
            sep_dat = []
            temp_dat =[]
            rho_dat = []
            shear_dat = []
            ratio_dat = []
            #shear2_dat = []
            bulk_dat = []
            diff_dat = []
            wa_dat = []
            fric1_dat = []
            fric2_dat = []
            for e in eps:
                print '\n\n++++++++++++++++++++++++ eps = %s +++++++++++++++++++++\n' %e
                name_plot += e
                z_dat = []
                sep_dat = []
                temp_dat =[]
                rho_dat = []
                rho_dat_init = []
                shear_dat = []
                ratio_dat = []
                #shear2_dat = []
                bulk_dat = []
                diff_dat = []
                wa_dat = []
                fric_dat = []
                fric1_dat = []
                fric2_dat = []
                fricinv_dat = []
                fw_dat = []
                name_plot2 = name_plot
                for s in sep:
                    print '\nSeparation = ', s
                    for i in den:
                        try:
                            name_plot2 += '_%s'% s
                            if fluid == 'LJ_channel':
                                results = confined_properties(fluid,m,temp,s, e, rhos)
                            else:
                                results = confined_properties(fluid,m,temp,s, e, i)
                            timestamp = 800000
                            deltat = 600000
                            
                            
                            sep_dat.append(results.separation())
                            temp_dat.append(results.temp())
                            rho_dat.append(results.rho())
                            shear_dat.append(results.shear2()[0])
                            bulk_dat.append(results.bulk2()[0])
                            ratio_dat.append(results.visc_ratio())
                            diff_dat.append(1e9*results.diff2())#(timestamp,deltat)[-1])
                            wa_dat.append(results.wa())
                            f1 = results.fric(fric_coords[0])
                            f2 = results.fric(fric_coords[1])
                            fric1_dat.append(f1)
                            fric2_dat.append(f2)
                            fric_dat.append((f1+f2)/2)
                            fricinv_dat.append(2/(f1+f2))
                            zcoord, densprof = results.profile()
                            fw_dat.append((results.wa()*results.separation())/(2*results.diff(timestamp,deltat)[-1]))
                            z_dat.append(int(s))

                            ax14.plot(zcoord, densprof, label='$\epsilon$ = %s, $\Delta z$ =%s\AA'%(e,s))
                            ax14.set_xlabel('z (\AA)')
                            ax14.set_ylabel('$\\rho$ (g/cm$^3$)')
                            ax14.legend()
                                #shear2_dat.append(results.shear2(float(temp),timestamp,deltat))
                        except:
                            print 'There was an error. Maybe no data for %s and s=%s exists for eps=%s.'%(m,s,e)


                znew_dat, sep_dat, sep_err       = averaging(z_dat, sep_dat)
                znew_dat, temp_dat, temp_err     = averaging(z_dat, temp_dat)
                znew_dat, shear_dat, shear_err   = averaging(z_dat, shear_dat)#, shear_err)
                znew_dat, bulk_dat, bulk_err     = averaging(z_dat, bulk_dat)#, bulk_err)
                znew_dat, ratio_dat, ratio_err   = averaging(z_dat, ratio_dat)
                znew_dat, diff_dat, diff_err     = averaging(z_dat, diff_dat)
                znew_dat, wa_dat, wa_err         = averaging(z_dat, wa_dat)
                znew_dat, fric_dat, fric_err     = averaging(z_dat, fric_dat)
                znew_dat, fw_dat, fw_err         = averaging(z_dat, fw_dat)
                znew_dat, rho_dat, rho_err       = averaging(z_dat, rho_dat)

            
                print znew_dat, sep_dat, rho_dat, fric_dat

                #print sep_dat, rho_dat, wa_dat, temp_dat
                legend_names = {'spce':'SPC/E', 'tip4p':'TIP4P/2005', 'TraPPE': 'TraPPE', 'SAFT': 'SAFT Dimer', \
                    'SAFT1': 'SAFT Monomer', 'TraPPEnc': 'TraPPE (no charges)', 'OPLS': 'OPLS',\
                    'SAFT1_rc2.5': 'SAFT (rc=2.5)', '12_6':'12-6 model', '12_6_rcsaft': '12-6 model (rc from SAFT)',\
                    '12_6_23_666': '12-6 model (SAFT exponents)', 'lj':'LJ'}
                legend_name = '$\epsilon$ = %s (%s)'%(e,legend_names[m])
                # Plotting
                ax1.plot(sep_dat, rho_dat, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                ax1.legend()

                ax2.errorbar(sep_dat, shear_dat, yerr=shear_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                #ax2.plot(sep_dat, shear2_dat, linestyle = 'None', marker=markers[count], c=colours[2*count], c='g', label='%s'%(m))
                ax2.set_xlabel('$\Delta z$ (\AA)')
                #ax2.set_xlim(0,18)
                ax2.set_ylabel('$\eta$ (mPa.s)')
                ax2.legend(loc='upper right')

                ax3.errorbar(sep_dat, bulk_dat, yerr=bulk_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                ax3.set_xlabel('$\Delta z$ (\AA)')
                ax3.set_ylabel('$\kappa$ (mPa.s)')
                ax3.legend(loc='upper left')

                ax4.errorbar(sep_dat, diff_dat, yerr=diff_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                #ax4.set_xlim(0,18)
                ax4.set_xlabel('$\Delta z$ (\AA)')
                ax4.set_ylabel('$D_s$ ($10^{-9}m^2/s$)')
                ax4.legend(loc='upper left')

                ax5.plot(rho_dat, shear_dat, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                ax5.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax5.set_ylabel('$\eta$ (mPa.s)')
                ax5.legend(loc='upper left')

                ax6.plot(rho_dat, bulk_dat, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                ax6.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax6.set_ylabel('$\kappa$ (mPa.s)')
                ax6.legend(loc='upper left')

                ax7.plot(rho_dat, diff_dat, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                #ax7.set_ylim(0,100)
                ax7.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax7.set_ylabel('$D_s$ ($10^{-9}m^2/s$)')
                ax7.legend()

                ax8.errorbar(sep_dat, ratio_dat, yerr=ratio_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                ax8.set_xlabel('$\Delta z$ (\AA)')
                ax8.set_ylabel('$\kappa$/$\eta$')
                ax8.legend(loc='upper left')

                ax9.plot(rho_dat, ratio_dat, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                ax9.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax9.set_ylabel('$\kappa$/$\eta$')
                ax9.legend(loc='upper left')

                ax10.errorbar(sep_dat, wa_dat, yerr=wa_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                ax10.set_xlabel('$\Delta z$ (\AA)')
                ax10.set_ylabel('$W_A$ (N/m)')
                ax10.legend()

                ax11.plot(rho_dat, wa_dat, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                ax11.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax11.set_ylabel('$W_A$ ($10^{-3}$ N/m)')
                ax11.legend()

                ax12.errorbar(sep_dat, fric_dat, yerr=fric_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                ax12.set_xlabel('$\Delta z$ (\AA)')
                ax12.set_ylabel('$\lambda$ ($10^4$ Ns/m$^3$)')
                ax12.legend()

                ax13.plot(rho_dat, fric_dat, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                ax13.set_xlabel('$\\rho$ (g/cm$^3$)')
                ax13.set_ylabel('$\lambda$ ($10^4$ Ns/m$^3$)')
                ax13.legend()

                ax15.errorbar(fw_dat, fric_dat, yerr=fric_err, xerr=fw_err, linestyle = 'None', marker=markers[count], c=colours[2*count], label=legend_name)
                ax15.set_xlabel('$W_A$D/$D_s$')
                ax15.set_ylabel('$\lambda$ ($10^4$ Ns/m$^3$)')
                ax15.legend()

                count +=1

                # Calculate predicted channel length
                slope, b = np.polyfit(np.array(fw_dat), np.array(fric_dat),1)
                xlength = 1/(slope*1e-5)
                print 'The channel has a predicted length of ', xlength,'A.'
                max_val = np.max(np.array(fw_dat))
                min_val = np.min(np.array(fw_dat))
                xdat = np.linspace(min_val, max_val, 100)
                fit_data = poly1(xdat, slope, b)
                ax15.plot(xdat, fit_data, linestyle = 'dashed' , c='k')



        fig1.savefig('PLOTS/%s_S_rho_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot, EXT))
        fig2.savefig('PLOTS/%s_SV_s_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot, EXT))
        fig3.savefig('PLOTS/%s_BV_s_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot, EXT))
        fig4.savefig('PLOTS/%s_Ds_s_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot, EXT))
        fig5.savefig('PLOTS/%s_SV_rho_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot, EXT))
        fig6.savefig('PLOTS/%s_BV_rho_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot, EXT))
        fig7.savefig('PLOTS/%s_Ds_rho_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot, EXT))
        fig8.savefig('PLOTS/%s_BS_s_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot, EXT))
        fig9.savefig('PLOTS/%s_BS_rho_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot, EXT))
        fig10.savefig('PLOTS/%s_WA_s_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot, EXT))
        fig11.savefig('PLOTS/%s_WA_rho_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot, EXT))
        fig12.savefig('PLOTS/%s_fric_s_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot, EXT))
        fig13.savefig('PLOTS/%s_fric_rho_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot, EXT))
        fig14.savefig('PLOTS/%s_prof_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot2, EXT))
        fig15.savefig('PLOTS/%s_fw_%s_T%s_%s.%s'%(extra_plot_name,fluid, temp,name_plot2, EXT))


    return




if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()
