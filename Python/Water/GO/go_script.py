#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
from go_tools import *
import sys


def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--sheets', nargs='+', required=False, default='1', action='store',
                       help='Number of sheets')
    parser.add_argument('-dx', '--delx', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-o', '--off', nargs='+', required=False, default='None', action='store')
    parser.add_argument('-f', '--force', required=False, nargs='+', default='Water', action='store')
    parser.add_argument('-func', '--func', required=False, nargs='+', default=['CH'], action='store')
    parser.add_argument('-png', '--png', nargs=1, type=str, required=False, default='n', action='store')
    args = parser.parse_args()
    return args

def main():
    args    = GetArgs()
    sheets  = args.sheets
    delx    = args.delx
    offset  = args.off
    force   = args.force
    func    = args.func
    png     = args.png[0]

    # Definitions
    fig_size_long = (14,5)
    fig_size_sq   = (9,7)
    if png == 'n':
        EXT = 'PDF'
        ext = 'pdf'
    else:
        EXT = 'PNG'
        ext = 'png'



    

    markers = ['o', 'D', 's', 'v', '^', 'd', '*']
    #colours = ['#3a42d4', '#d43a3a', '#60c1dc',\
    #'#60dc96', '#dc6060', '#addc60', '#dc6079', '#60dcd4']
    colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']

    offset_label = {('1','0.0'):'Single sheet', ('2','0.0'):'$c=0.0$', ('2','0.25'):'$c=0.25$',('2','0.5'):'$c=0.5$',
    ('3','0.0'):'n=3, c=0.0',('4','0.0'):'n=4, c=0.0',('5','0.0'):'n=5, c=0.0',('6','0.0'):'n=6, c=0.0',\
    ('3','0.5'):'n=3, c=0.5',('4','0.5'):'n=4, c=0.5',('5','0.5'):'n=5, c=0.5',('6','0.5'):'n=6, c=0.5'}

    WIDTH = {('CH','5'):3.29, ('CH','6'): 4, ('CH','7'): 4.83, ('CH','9'):6.64, ('CH','11'):8.43,\
    ('CH','13'):10.35, ('CH', '15'):12.26,\
    ('COH','5'):3.04, ('COH','7'):4.49, ('COH','9'):6.25, ('COH','11'):8.08, \
    ('CCOO','5'):1.74, ('CCOO','7'):3.14, ('CCOO','9'):4.63, ('CCOO','11'):6.40, \
    ('CCOOH','5'):1.55, ('CCOOH','7'):2.82, ('CCOOH','9'):4.30, ('CCOOH','11'):6.01}
    LABEL = {'CH':'H', 'CCOO':'COO$^-$', 'CCOOH': 'COOH', 'COH': 'OH'}

    fig10 = plt.figure(figsize=fig_size_sq)
    ax10  = fig10.add_axes([0.15,0.15,0.75,0.75])

    name_func = ''
    count_fn = 0
    for fn in func:
        print '\n#----------------------- func = {} ---------------------------#\n'.format(fn)
        name_func += fn
        d_eff = []
        

        fig4 = plt.figure(figsize=fig_size_sq)
        ax4  = fig4.add_axes([0.15,0.15,0.75,0.75])
        fig6 = plt.figure(figsize=fig_size_sq)
        ax6  = fig6.add_axes([0.15,0.15,0.75,0.75])
        fig7 = plt.figure(figsize=fig_size_sq)
        ax7  = fig7.add_axes([0.15,0.15,0.75,0.75])
        fig8 = plt.figure(figsize=fig_size_sq)
        ax8  = fig8.add_axes([0.15,0.15,0.75,0.75])
        fig9 = plt.figure(figsize=fig_size_sq)
        ax9  = fig9.add_axes([0.15,0.15,0.75,0.75])

        name_sheets = '{}_delx{}_o'.format(fn,delx[0])
        for o in offset:
            name_sheets += o
        name_sheets += '_n'
        perm_sheets, perm_sheets_err = [], []
        count_8 = 0
        for n in sheets:
            print '\n#----------------------- n = {} ---------------------------#\n'.format(n)
            for o in offset:
                name_plot = '{}_n{}_o{}_delx{}_F'.format(fn,n,o,delx)
                name_comb = '{}_n{}_o{}_delx'.format(fn,n,o)
                name_sheets += n
                count = 0
                perm = []
                perm_err = []
                for dx in delx:
                    if int(n)==1 and len(sheets)==1:
                        d_eff.append(WIDTH[(fn,dx)])
                    vel_max = 0

                    print '\n#----------------------- delx = {} ------------------------#\n'.format(dx)
                    name_plot = '{}_n{}_o{}_delx{}_F'.format(fn,n,o,dx)
                    name_comb += dx
                    delP = []
                    delRho = []
                    Jz = []
                    Jz_err = []
                    width = []
                    rho_mean = []

                    fig1 = plt.figure(figsize=fig_size_long)
                    ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])
                    fig2 = plt.figure(figsize=fig_size_long)
                    ax2  = fig2.add_axes([0.15,0.15,0.75,0.75])
                    fig3 = plt.figure(figsize=fig_size_sq)
                    ax3  = fig3.add_axes([0.15,0.15,0.75,0.75])
                    fig5 = plt.figure(figsize=fig_size_sq)
                    ax5  = fig5.add_axes([0.15,0.15,0.75,0.75])

                    try:
                        for f in force:
                            name_plot += f

                            try:
                                flow = FLOW(n,o,dx,f,fn)
                                C, P, dP = flow.stress_prof()
                                CZ, RHOZ, CX, RHOX, dRho, VELZ, LD, RHOmean, RHOerr, VELerr = flow.density_prof('None')
                                JZ_MOL, JZ_ERR = flow.flux()

                                delP.append(dP)
                                delRho.append(dRho)
                                Jz.append(JZ_MOL)
                                Jz_err.append(JZ_ERR)
                                width.append(LD)
                                rho_mean.append(RHOmean)


                                ax1.plot(C,P, label='$A={}$'.format(f))
                                ax2.plot(CZ,RHOZ, label='A={}'.format(f))
                                ax3.plot(CX,RHOX, label='A={}'.format(f))
                                #ax5.plot(CX,VELX, label='F={}'.format(f))
                            except:
                                print '\nThere was an error. Maybe no file for n={}, delx={}, o={} and F={} exists.\n'.format(n,dx,o,f)

                        delP_fit, Jz_fit, permeance, permeance_err = straight_fit(delP, Jz, 0, 300)
                        width_avg, width_err = np.mean(np.array(width)), stats.sem(np.array(width))
                        rho_avg = np.mean(rho_mean)
                        permeability_si, permeability_err, permeance_si = perm_units(permeance[0], width_avg, rho_avg, width_err, permeance_err[0])
                        print '\nThe permeability for delx =', dx, 'is', permeability_si, '+/-', permeability_err, 'm^2/sPa'
                        print 'The permeance for delx =', dx, 'is', permeance_si, 'm/sPa'
                        print 'The average width is', width_avg, 'A \n'
                        perm.append(permeability_si/1e-17)
                        perm_err.append(permeability_err/1e-17)
                        perm_sheets.append(permeability_si/1e-17)
                        perm_sheets_err.append(permeability_err/1e-17)
                        


                        ax4.plot(delRho, delP, marker=markers[count], label='$d_{\mathrm{slit}}$ = %s'%(dx))
                        ax5.errorbar(delP, Jz, yerr=Jz_err, marker=markers[count], c=colours[count], linestyle='None', label='$d_{\mathrm{slit}}$ = %s'%(dx))
                        ax5.plot(delP_fit, Jz_fit, c=colours[count],linestyle='dashed')
                        ax6.errorbar(delP, Jz, yerr=Jz_err, marker=markers[count], c=colours[count],linestyle='None', label='$d_{\mathrm{slit}}$ = %s'%(dx))
                        ax6.plot(delP_fit, Jz_fit, c=colours[count],linestyle='dashed')
                        ax8.errorbar(delP, Jz, yerr=Jz_err, marker=markers[count_8], c=colours[count_8*2],linestyle='None', label=offset_label[(n, o)])
                        ax8.plot(delP_fit, Jz_fit, c=colours[count_8*2],linestyle='dashed')
                        count_8 += 1


                        ax1.set_xlabel('z (\AA)')
                        ax1.set_ylabel('P (MPa)')
                        ax1.set_ylim(-100,400)
                        ax1.set_xlim(0,100)
                        ax1.legend()

                        ax2.set_xlabel('z (\AA)')
                        ax2.set_ylabel('$\\rho$ (g/cm$^3$)')
                        #ax2.set_ylim(-50,300)
                        ax2.set_xlim(0,100)
                        ax2.legend(loc='upper left')

                        ax3.set_xlabel('x (\AA)')
                        ax3.set_ylabel('$\\rho$ (g/cm$^3$)')
                        #ax3.set_ylim(-50,300)
                        #ax3.set_xlim(0,100)
                        ax3.legend()

                        ax4.set_xlabel('$\Delta\\rho$ (g/cm$^3$)')
                        ax4.set_ylabel('$\Delta$P (MPa)')
                        #ax3.set_ylim(-50,300)
                        #ax3.set_xlim(0,100)
                        ax4.legend()

                        ax5.set_xlabel('$\Delta$P (MPa)')
                        ax5.set_ylabel('$J_z$ (10$^3$ mol/m$^2$s)')
                        #ax5.set_ylim(0,vel_max)
                        ax5.set_xlim(0, 150)
                        ax5.legend()

                        ax6.set_xlabel('$\Delta$P (MPa)')
                        ax6.set_ylabel('$J_z$ (10$^3$ mol/m$^2$s)')
                        ax6.set_ylim(0, 150)
                        ax6.set_xlim(0, 230)
                        ax6.legend()

                        ax8.set_xlabel('$\Delta$P (MPa)')
                        ax8.set_ylabel('$J_z$ (10$^3$ mol/m$^2$s)')
                        ax8.set_ylim(0,30)
                        ax8.set_xlim(0, 230)
                        ax8.legend()

                        fig1.savefig('PLOTS/{}/press_{}.{}'.format(EXT, name_plot, ext))
                        fig2.savefig('PLOTS/{}/rhoz_{}.{}'.format(EXT, name_plot, ext))
                        fig3.savefig('PLOTS/{}/rhox_{}.{}'.format(EXT, name_plot, ext))
                        fig5.savefig('PLOTS/{}/Jz_{}.{}'.format(EXT, name_plot, ext))

                        plt.close(fig1)
                        plt.close(fig2)
                        plt.close(fig3)
                        plt.close(fig5)
                        count += 1
                    except:
                        print '\nMaybe no data exists for this set of parameters.\n'
        
        if len(sheets) == 1 and len(delx) >1:
            delx_int = map(int,delx)
            perm_dx_fit = exp_fit(delx_int, perm, 0, 15.5)
            A_lower = 0.0610023944183
            B_lower = 0.42515763262
            A_higher = 1.66724563937
            B_higher =  0.13260576073
            #print 'k fit parameters for delx dependence: A=', perm_dx_fit[2][0],', B=', perm_dx_fit[2][1] 
            ax7.errorbar(delx_int, perm, yerr=perm_err, marker='o', markersize=11, linestyle='None', c=colours[0])
            ax7.plot(perm_dx_fit[0], perm_dx_fit[1], linestyle='dashed', c=colours[0], label='Exponential fit')
            ax7.plot(perm_dx_fit[0], f3(perm_dx_fit[0], A_lower, B_lower), linestyle='dotted', c=colours[1], label='Fit (5-9 \AA)')
            #ax7.plot(perm_dx_fit[0], f3(perm_dx_fit[0], A_higher, B_higher), linestyle='dotted', c=colours[2], label='Fit (11-15 \AA)')
            ax7.set_xlabel('$d_{\mathrm{slit}}$ (\AA)')
            ax7.set_ylabel('$k$ (10$^{-17}$ m$^2$s$^{-1}$Pa$^{-1}$)')
            ax7.set_ylim(-2,12)
            ax7.set_xlim(0,16)
            ax7.legend(loc='upper left')

            if len(sheets)==1 and int(sheets[0])==1:
                perm_deff_fit = poly_fit(d_eff, perm, 0, 15.5) #straight_fit2(d_eff, perm, 0, 15.5)
                ax10.errorbar(d_eff, perm, yerr=perm_err, marker='o', markersize=11, linestyle='None', c=colours[count_fn*2], label=LABEL[fn])
                ax10.plot(perm_deff_fit[0], perm_deff_fit[1], linestyle='dashed',c=colours[count_fn*2])
                ax10.set_xlabel('$d_{\mathrm{eff}}$ (\AA)')
                ax10.set_ylabel('$k$ (10$^{-17}$ m$^2$s$^{-1}$Pa$^{-1}$)')
                ax10.set_ylim(0,7)
                ax10.set_xlim(0,9)
                ax10.legend(loc='upper left')

            ## inset
            #a7 = fig7.add_axes([0.25, 0.35, 0.3, 0.3])
            #a7.errorbar(delx_int, perm, yerr=perm_err, marker='o', markersize=11, linestyle='None', c=colours[0])
            #a7.plot(perm_dx_fit[0], perm_dx_fit[1], linestyle='dashed', c=colours[0])
            #a7.plot(perm_dx_fit[0], f3(perm_dx_fit[0], A_lower, B_lower), linestyle='dotted', c=colours[1])
            #a7.set_ylim(0,3)
            #a7.set_xlim(4.5,9.5)

            fig7.savefig('PLOTS/{}/perm_{}.{}'.format(EXT, name_comb, ext))
        if len(offset)==1 and len(sheets) > 1:
            if offset[0]=='0.0':
                ax9.set_xlim(0,7)
                ax9.set_ylim(0,1.4)
            sheets_int = map(int,sheets)
            print sheets_int, perm_sheets
            perm_sheets_fit = exp_fit2(map(float,sheets_int), perm_sheets, 0, 7)
            ax9.errorbar(sheets_int, perm_sheets, yerr=perm_sheets_err, marker='o', linestyle='None')
            ax9.plot(perm_sheets_fit[0], perm_sheets_fit[1], linestyle='dashed', c=colours[0])
            ax9.set_xlabel('no. of sheets')
            ax9.set_ylabel('$k$ (10$^{-17}$ m$^2$s$^{-1}$Pa$^{-1}$)')
            fig9.savefig('PLOTS/{}/perm_{}.{}'.format(EXT, name_sheets, ext))


        fig4.savefig('PLOTS/{}/drop_{}.{}'.format(EXT, name_comb, ext))
        fig6.savefig('PLOTS/{}/Jz_{}.{}'.format(EXT, name_comb, ext))
        fig8.savefig('PLOTS/{}/Jz_{}.{}'.format(EXT, name_sheets, ext))
        count_fn +=1
    # func combined plots
    fig10.savefig('PLOTS/{}/perm_{}.{}'.format(EXT, name_func, ext))


    #plt.show()


    

    return

if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()