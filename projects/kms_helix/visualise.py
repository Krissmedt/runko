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

from pyhack.py_runko_aux_3d import *

data_root = "./"
fig_name = ""
plim = 1

filenames = {}

filenames["Collocation M3"] = "coll_M3_helix_data.h5"
filenames["Velocity-Verlet"] = "vv_helix_data.h5"
filenames["Boris-SDC M3K4"] = "sdc_M3K4_helix_data.h5"
filenames["Runko Leapfrog"] = "lf_helix_data.h5"


plot_params = {}
plot_params['legend.fontsize'] = 22
plot_params['figure.figsize'] = (12,8)
plot_params['axes.labelsize'] = 24
plot_params['axes.titlesize'] = 24
plot_params['xtick.labelsize'] = 24
plot_params['ytick.labelsize'] = 24
plot_params['lines.linewidth'] = 4
plot_params['axes.titlepad'] = 5
plot_params['legend.loc'] = 'upper right'
plt.rcParams.update(plot_params)
r = 1
g = 1
b = 1

fig_isotraj = plt.figure(1)
ax_iso = fig_isotraj.gca(projection='3d')

fig_traj = plt.figure(2)
ax_traj = fig_traj.add_subplot(111)

linewidths = [5,5,5,5]
linestyles = ["-","--","-.",":"]

i = 0
for key,value in filenames.items():
    file = h5.File(data_root+value,'r')
    Nt = file["fields/Nt"]
    t = file["fields/t"]
    x = file["fields/x"]
    y = file["fields/y"]
    z = file["fields/z"]
    vx = file["fields/vx"]
    vy = file["fields/vy"]
    vz = file["fields/vz"]
    if "SDC" in key:
        xres = file["fields/xres"]
        vres = file["fields/vres"]

    ################ Plot isometric trajectory ################################
    for pii in range(0,plim):
        ax_iso.plot3D(x[:,pii],
                  y[:,pii],
                  zs=z[:,pii],
                  label=key,linewidth=linewidths[i],linestyle=linestyles[i])

    ax_iso.set_xlim([2,8])
    ax_iso.set_ylim([2,8])
    ax_iso.set_zlim([1,6])
    ax_iso.set_xlabel('\n$x$')
    ax_iso.set_ylabel('\n$y$')
    ax_iso.set_zlabel('\n$z$')
    ax_iso.legend()


    ################ Plot xy-plane trajectory ################################
    ax_traj.plot(x,y,label=key,linewidth=linewidths[i],linestyle=linestyles[i])
    ax_traj.set_xlim([2,8])
    ax_traj.set_ylim([2,8])
    ax_traj.set_aspect('equal')
    ax_traj.set_xlabel('$x$')
    ax_traj.set_ylabel('$y$')
    ax_traj.legend()

    ################ Plot SDC residual  ######################################
    if "SDC" in key:
        pass
        # fig_xres = plt.figure(3)
        # ax_xres = fig_traj.add_subplot(111)
        # ax_xres.plot(t,xres,label=key,linewidth=linewidths[i],linestyle=linestyles[i])
        # ax_xres.set_xlabel('$x$')
        # ax_xres.set_ylabel('$y$')
        # ax_xres.legend()
        #
        # fig_xres.savefig(data_root + 'helix_xres'+ fig_name + '.pdf', dpi=150, facecolor='w', edgecolor='w',orientation='portrait',pad_inches=0.2,bbox_inches = 'tight')
        # fig_vres.savefig(data_root + 'helix_vres'+ fig_name + '.pdf', dpi=150, facecolor='w', edgecolor='w',orientation='portrait',pad_inches=0.2,bbox_inches = 'tight')
        #

    i+=1

fig_isotraj.savefig(data_root + 'helix_isotraj'+ fig_name + '.pdf', dpi=150, facecolor='w', edgecolor='w',orientation='portrait',pad_inches=0.2,bbox_inches = 'tight')
fig_traj.savefig(data_root + 'helix_traj'+ fig_name + '.pdf', dpi=150, facecolor='w', edgecolor='w',orientation='portrait',pad_inches=0.2,bbox_inches = 'tight')
