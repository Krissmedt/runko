from math import sqrt, fsum, pi, exp, cos, sin, floor
from decimal import Decimal
import io
import pickle as pk
import matplotlib.pyplot as plt
import numpy as np
import cmath as cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import h5py as h5
import sys
import traceback
from collections import OrderedDict

from pyhack.py_runko_aux import *

data_root = "./"
fig_name = "test"

filenames = {}

filenames["Boris-SDC M5K5"] = "sdc55_wp_gyro10.h5"
filenames["Leapfrog-Py"] = "lf_wp_gyro10.h5"
filenames["Velocity-Verlet-Py"] = "vv_wp_gyro10.h5"
filenames["Boris-SDC M3K3"] = "sdc_wp_gyro10.h5"


plot_params = {}
plot_params['legend.fontsize'] = 16
plot_params['figure.figsize'] = (12,8)
plot_params['axes.labelsize'] = 20
plot_params['axes.titlesize'] = 20
plot_params['xtick.labelsize'] = 16
plot_params['ytick.labelsize'] = 16
plot_params['lines.linewidth'] = 4
plot_params['axes.titlepad'] = 5
plot_params['legend.loc'] = 'upper right'
plt.rcParams.update(plot_params)

for key,value in filenames.items():
    file = h5.File(data_root+value,'r')
    Nt = file["fields/Nt"]
    dumptime = file["fields/t"]
    x = file["fields/x"]
    y = file["fields/y"]
    vx = file["fields/vx"]
    vy = file["fields/vy"]

    if key == "Boris-SDC M5K5":
        solx = np.copy(x[-1,0])

    errors = np.abs(x[:,0]-solx)/np.abs(solx)

    print(solx)
    print(x)
    print(errors)

    label = key
    rhs = Nt

    ##Order Plot w/ rhs
    fig_rhs = plt.figure(1)
    ax_rhs = fig_rhs.add_subplot(1, 1, 1)
    ax_rhs.plot(rhs,errors,marker="o",label=label)

    ##Order Plot w/ Nt
    fig_nt = plt.figure(2)
    ax_nt = fig_nt.add_subplot(1, 1, 1)
    ax_nt.plot(Nt,errors,marker="o",label=label)

handles, labels = fig_rhs.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax_rhs.legend(by_label.values(), by_label.keys(),loc='best')

handles, labels = fig_nt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax_nt.legend(by_label.values(), by_label.keys(),loc='best')


axnl_list = []
axnl_list.append(ax_rhs)
axnl_list.append(ax_nt)

i = 0
for ax in axnl_list:
    i +=1
    if i == 1:
        orderSlope = -1
        ax.set_xlabel('RHS evaluations')
    else:
        ax.set_xlabel(r'$N t$')
        orderSlope = -1

    ax.set_xscale('log')
    #ax_rhs.set_xlim(10**3,10**5)
    ax.set_yscale('log')
    ax.set_ylim(10**(-15),10**(-1))
    ax.set_ylabel(r'$\Delta x^{\mathrm{rel}}$')

    xRange = ax.get_xlim()
    yRange = ax.get_ylim()

    ax.plot(xRange,orderLines(1*orderSlope,xRange,yRange),
                ls='dashdot',c='0.2')
    ax.plot(xRange,orderLines(2*orderSlope,xRange,yRange),
                ls='dotted',c='0.4')
    ax.plot(xRange,orderLines(4*orderSlope,xRange,yRange),
                ls='dashed',c='0.6')

    fig_rhs.savefig(data_root + 'gyro_'+ fig_name + '_rhs.pdf', dpi=150, facecolor='w', edgecolor='w',orientation='portrait',pad_inches=0.0,bbox_inches = 'tight')
    fig_nt.savefig(data_root + 'gyro_' + fig_name + '_nt.pdf', dpi=150, facecolor='w', edgecolor='w',orientation='portrait',pad_inches=0.0,bbox_inches = 'tight')
