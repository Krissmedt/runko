from math import sqrt, fsum, pi, exp, cos, sin, floor, isclose
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
import warnings
warnings.filterwarnings("ignore")

from pyhack.py_runko_aux import *

data_root = "./"
fig_name = "wp"

filenames = {}

filenames["Ref"] = "coll_M5_pen_ref.h5"
# filenames["Collocation M5"] = "coll_M5_pen.h5"
filenames["Velocity-Verlet"] = "vv_pen.h5"

filenames["Boris-SDC M2K1"] = "sdc_M2K1_pen.h5"
filenames["Boris-SDC M2K2"] = "sdc_M2K2_pen.h5"
filenames["Boris-SDC M3K2"] = "sdc_M3K2_pen.h5"
filenames["Boris-SDC M3K3"] = "sdc_M3K3_pen.h5"
filenames["Boris-SDC M3K4"] = "sdc_M3K4_pen.h5"
filenames["Boris-SDC M5K4"] = "sdc_M5K4_pen.h5"
filenames["Boris-SDC M5K5"] = "sdc_M5K5_pen.h5"
filenames["Boris-SDC M5K6"] = "sdc_M5K6_pen.h5"

# filenames["Boris-SDC M2K1"] = "sdc_M2K1_pen.h5"
# filenames["Leapfrog-Py"] = "lf_pen.h5"


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
r = 1
g = 1
b = 1

for key,value in filenames.items():
    file = h5.File(data_root+value,'r')
    Nt = file["fields/Nt"]
    dumptime = file["fields/t"]
    x = file["fields/x"]
    y = file["fields/y"]
    z = file["fields/z"]
    vx = file["fields/vx"]
    vy = file["fields/vy"]
    vz = file["fields/vz"]

    if key == "Ref":
        solx = np.copy(x[-1,0])
        soly = np.copy(y[-1,0])
        solz = np.copy(z[-1,0])
        solvx = np.copy(vx[-1,0])
        solvy = np.copy(vy[-1,0])
        solvz = np.copy(vz[-1,0])
        ref_time = np.copy(dumptime)
        continue

    assert np.all(np.isclose(dumptime[:],dumptime[0]))
    assert np.isclose(dumptime[0],ref_time)

    x_errors = np.abs(x[:,0]-solx)/np.abs(solx)
    y_errors = np.abs(y[:,0]-soly)/np.abs(soly)
    z_errors = np.abs(z[:,0]-solz)/np.abs(solz)
    vx_errors = np.abs(vx[:,0]-solvx)/np.abs(solvx)
    vy_errors = np.abs(vy[:,0]-solvy)/np.abs(solvy)
    vz_errors = np.abs(vz[:,0]-solvz)/np.abs(solvz)

    errors = (x_errors+y_errors+z_errors)/3
    v_errors = (vx_errors+vy_errors+vz_errors)/3
    # errors = np.abs(x[:-1,0]-x[-1,0])/np.abs(x[-1,0])
    # v_errors = np.abs(vx[:-1,0]-vx[-1,0])/np.abs(vx[-1,0])

    xfactors = np.log2(errors[:-1]/errors[1:])
    vfactors = np.log2(v_errors[:-1]/v_errors[1:])
    print(key+" x order: {0}".format(xfactors))
    print(key+" v order: {0}".format(vfactors))

    label = key
    # Nt = Nt[1:]
    rhs = Nt


    if key == "Velocity-Verlet":
        c = "black"
        rhs = Nt[:]
    if "Boris-SDC M2" in key:
        sims = 2
        c = (0,g,0)
        g -= 1/sims
        M = int(key[key.find("M")+1])
        K = int(key[key.find("K")+1])
        rhs = (M-1)*K*Nt[:]
    if "Boris-SDC M3" in key:
        sims = 3
        c = (0,0,b)
        b -= 1/sims
        M = int(key[key.find("M")+1])
        K = int(key[key.find("K")+1])
        rhs = (M-1)*K*Nt[:]
    if "Boris-SDC M5" in key:
        sims = 3
        c = (r,0,0)
        r -= 1/sims
        M = int(key[key.find("M")+1])
        K = int(key[key.find("K")+1])
        rhs = (M-1)*K*Nt[:]


    ## x order Plot w/ rhs
    fig_rhs = plt.figure(1)
    ax_rhs = fig_rhs.add_subplot(1, 1, 1)
    ax_rhs.plot(rhs,errors,marker="o",color=c,label=label)

    ## x order Plot w/ Nt
    fig_nt = plt.figure(2)
    ax_nt = fig_nt.add_subplot(1, 1, 1)
    ax_nt.plot(Nt,errors,marker="o",color=c,label=label)

    ## v order Plot w/ Nt
    fig_v_nt = plt.figure(3)
    ax_v_nt = fig_v_nt.add_subplot(1, 1, 1)
    ax_v_nt.plot(Nt,v_errors,marker="o",color=c,label=label)

handles, labels = fig_rhs.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax_rhs.legend(by_label.values(), by_label.keys(),loc='lower left')

handles, labels = fig_nt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax_nt.legend(by_label.values(), by_label.keys(),loc='upper right')

handles, labels = fig_v_nt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax_v_nt.legend(by_label.values(), by_label.keys(),loc='upper right')

axnl_list = []
axnl_list.append(ax_rhs)
axnl_list.append(ax_nt)
axnl_list.append(ax_v_nt)

i = 0
for ax in axnl_list:
    i +=1
    if i == 1:
        ax.set_xlabel('RHS evaluations')
    elif i == 3:
        ax.set_ylabel(r'$\Delta v^{\mathrm{rel}}$')
    else:
        ax.set_ylabel(r'$\Delta x^{\mathrm{rel}}$')
        ax.set_xlabel(r'$N t$')

    orderSlope = -1
    ax.set_xscale('log')
    #ax_rhs.set_xlim(10**3,10**5)
    ax.set_yscale('log')
    ax.set_ylim(10**(-10),10**(0))

    xRange = ax.get_xlim()
    yRange = ax.get_ylim()

    ax.plot(xRange,orderLines(1*orderSlope,xRange,yRange),
                ls='dashdot',c='0.2')
    ax.plot(xRange,orderLines(2*orderSlope,xRange,yRange),
                ls='dotted',c='0.4')
    ax.plot(xRange,orderLines(4*orderSlope,xRange,yRange),
                ls='dashed',c='0.6')
    ax.plot(xRange,orderLines(8*orderSlope,xRange,yRange),
                ls=(0,(5,1)),c='0.8')

fig_rhs.savefig(data_root + 'pen_'+ fig_name + '_rhs.pdf', dpi=150, facecolor='w', edgecolor='w',orientation='portrait',pad_inches=0.05,bbox_inches = 'tight')
fig_nt.savefig(data_root + 'pen_' + fig_name + '_nt.pdf', dpi=150, facecolor='w', edgecolor='w',orientation='portrait',pad_inches=0.05,bbox_inches = 'tight')
fig_v_nt.savefig(data_root + 'pen_' + fig_name + '_v_nt.pdf', dpi=150, facecolor='w', edgecolor='w',orientation='portrait',pad_inches=0.05,bbox_inches = 'tight')
