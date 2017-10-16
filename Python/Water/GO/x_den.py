#!/usr/bin/python

from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
import sys

def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-delx", type=str, nargs='+', \
                    help="Delta x",required=True)
    parser.add_argument("-func", type=str, nargs='+', \
                    help="Functionalisation",required=True)
    parser.add_argument('-png', '--png', nargs=1, type=str, required=False, default='n', action='store')
    args = parser.parse_args()
    return args


def den_data_x(func,dx):

    filename = "DATA/dens.{}_n1_o0.0_delx{}_F1.0".format(func,dx)

    # read in data
    f = open(filename,'r')
    data = f.read()

    data_lines = data.split('\n')


    xbin = []
    zbin = []


    xlist = list(np.arange(17.1,17.1+float(dx)+0.2,0.2))
    zlist = list(np.arange(48.7, 51.5+0.2, 0.2))

    # round values in list so that they agree with data
    xlist = [round(n, 1) for n in xlist]
    zlist = [round(n, 1) for n in zlist]

    #print xlist, zlist

    xdim = len(xlist)
    zdim = len(zlist)
    density = np.zeros((xdim,zdim))
    coords_x = np.zeros((xdim,zdim))
    coords_z = np.zeros((xdim,zdim))


    count = 1
    for j in range(4,len(data_lines)-1):
        if len(data_lines[j].split()) != 3:
            x_val = float(data_lines[j].split()[1])
            z_val = float(data_lines[j].split()[2])
            den_val = float(data_lines[j].split()[4])

            # only select values within slit
            if x_val in xlist and z_val in zlist:
                x_ind = xlist.index(x_val)
                z_ind = zlist.index(z_val)
                xbin.append(x_val)
                zbin.append(z_val)


                # average over all time steps

                density[x_ind][z_ind] = (density[x_ind][z_ind]+den_val)/2
                coords_x[x_ind][z_ind] = xlist[x_ind]
                coords_z[x_ind][z_ind] = zlist[z_ind]

                count += 1
    print count-1, 'timesteps collected'


    density = np.average(density, axis = 1)
    coords_x = np.average(coords_x, axis = 1)


    return coords_x, density


#--------------------------------------------------------------------#

def main():
    args    = GetArgs()
    DELX    = args.delx
    FUNC    = args.func
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


    # PLOTTING

    den_vals = []
    x_vals = []
    delx_vals = []

    colours=['#313695', '#4575b4', '#74add1',\
    '#abd9e9', '#fdae61', '#f46d43', '#d73027', '#a50026']
    LABEL = {'CH':'H', 'CCOO':'COO$^-$', 'CCOOH': 'COOH', 'COH': 'OH'}



    fig1 = plt.figure(figsize=fig_size_sq)
    ax1  = fig1.add_axes([0.15,0.15,0.75,0.75])

    
    func_name = ''
    delx_name = ''
    c_count=0
    for f in FUNC:
        print '############################### {} ########################'.format(f)
        func_name += f
        for dx in DELX:
            print '############################### delx = {} ########################'.format(dx)
            try:
                if f=='CH' and dx=='9':
                    continue
                else:

                    X, DEN = den_data_x(f,dx)
                    if float(dx)%2 == 0:
                    	centre = 17.1+(float(dx)/2.)-0.6
                    else:
                    	centre = 17.1+(float(dx)/2.)-0.7


                    for i in range(len(X)):
                    	if np.allclose(X[i], centre) == True:
                    		cen_pos = i

                    den_vals.append(DEN)
                    x_vals.append(X-X[cen_pos])
                    ax1.plot(X-X[cen_pos],DEN, label="{}".format(LABEL[f]),c=colours[c_count*2],lw =2.0)
                    c_count +=1
                    delx_name += dx
            except:
                print 'File not found.'

    plotname = 'xden_'+func_name+'_'+delx_name

    ax1.plot([-6,6], [1,1],linestyle='dashed', c='k')
    #plt.title('$(a)$')
    #plt.legend(prop={'size':13})
    plt.xlim(-6,6)
    plt.ylim(0,4)
    #plt.setp(a, plt.ylim(-4,4))#, plt.legend(), 
    plt.ylabel('$\\rho$ $(g/\mathrm{cm}^3)$')
    plt.xlabel('Distance from center (\AA)')
    ax1.legend()


    fig1.savefig('PLOTS/{}/{}.{}'.format(EXT, plotname, ext))
    plt.show()

    return


if __name__ == "__main__":
    sys.path.append('/home/fred/SCRIPTS/Python')
    from plotting_params import *
    main()
