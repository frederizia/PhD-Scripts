from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#lda,sx,hf,gj,fj4,fj5,exp
x_vals = np.array([1,2,3,4,5,6])

gamma = -1*np.array([-2.577405261,-3.792932647,-8.847171971,-2.670535066,
                  -2.695278393,-3.040965219, -3.37])
L_point = -1*np.array([-2.792312303,-4.124349772,-9.831080881,-2.693446481,
                    -2.699703336,-2.93720423, -3.3])

X_point = -1*np.array([-3.5269563,-4.899066432,-10.41396872,-3.569778736,
                    -3.594804796,-3.947183324,-4.03])
BG = np.array([0.5,1.5,6.3,0.7,0.7,1.2,1.17])

#print len(gamma)
Names = ('LDA','SX','HF','GJ','$f_2$', '$f_5$', 'Exp.')
N = 7
ind = np.arange(1,N+1)
width = 0.2

matplotlib.rcParams.update({'font.size': 18})
fig = plt.figure()
ax = fig.add_subplot(111)
gam_bar = ax.bar(ind, gamma, width, color = 'b')

L_bar = ax.bar(ind+width, L_point, width, color='r')

X_bar = ax.bar(ind+2*width, X_point, width, color='g')

ax.set_ylabel('Energy Gap (eV)')#, fontsize = 20)
plt.tick_params(pad=9)
ax.set_xticks(ind+width)
ax.set_xticklabels(Names)

ax.legend((gam_bar[0],L_bar[0],X_bar[0]),('$\Gamma$','L','X'))
plt.savefig('bg_points.pdf')
plt.show()



