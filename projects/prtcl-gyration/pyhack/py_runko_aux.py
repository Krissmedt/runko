import numpy as np
import matplotlib.pyplot as plt
from pyhack.analytical_gyro import *

def py_pos(cont):
    # Retrieve position vectors from particle container and rearrange to numpy Nx3
    locx = cont.loc(0)
    locy = cont.loc(1)
    locz = cont.loc(2)

    # loc is meant to be passed by reference, but container doesnt change for some reason
    # locx += dt * vel[:,0]*g

    pos = np.array([locx,locy,locz]).transpose()

    return pos

def py_vel(cont):
    # Retrieve velocity vectors from particle container and rearrange to numpy Nx3
    velx = cont.vel(0)
    vely = cont.vel(1)
    velz = cont.vel(2)
    vel = np.array([velx,vely,velz]).transpose()

    return vel

def py_em(cont):
    # Retrieve EM vectors from particle container and rearrange to numpy Nx3
    locx = cont.loc(0)
    nq = len(locx)

    E = np.array([[cont.ex(i),cont.ey(i),cont.ez(i)] for i in range(0,nq)])
    B = np.array([[cont.bx(i),cont.by(i),cont.bz(i)] for i in range(0,nq)])

    return E,B


def ginv(c,v):
    gamma = c/np.sqrt(c**2 + np.sum(v**2,axis=1))
    return gamma


def output_stag(t,x,y,vx,vy,conf,name):
    t = np.array(t)
    x = np.array(x)
    y = np.array(y)
    vx = np.array(vx)
    vy = np.array(vy)

    plot_traj(t,x,y,conf,name)
    plot_vel(t-conf.cfl/2,vx,vy,conf,name)


def output_sync(t,x,y,vx,vy,conf,name):
    t = np.array(t)
    x = np.array(x)
    y = np.array(y)
    vx = np.array(vx)
    vy = np.array(vy)

    plot_traj(t,x,y,conf,name)
    plot_vel(t,vx,vy,conf,name)


def plot_traj(t,x,y,conf,name):
    xref,yref,dummy1,dummy2 = analytical_gyro_full(t,conf)

    fig_traj = plt.figure(1)
    ax_traj = fig_traj.add_subplot(111)
    ax_traj.plot(x,y,label="sim")
    ax_traj.plot(xref,yref,label="sol")
    # ax_traj.set_xlim([0,90])
    # ax_traj.set_ylim([0,90])
    ax_traj.set_aspect('equal')
    ax_traj.legend()

    fig_traj.savefig(name+'_trajectory_{0}.png'.format(conf.Nt))

    x_error = abs(xref[-1]-x[-1])/abs(xref[-1])
    rl1 = (np.max(x[:,0])-np.min(x[:,0]))/2
    l_error = abs(conf.larmor - rl1)/abs(rl1)

    print(l_error)

def plot_vel(t,vx,vy,conf,name):
    dummy1,dummy2,vxref,vyref = analytical_gyro_full(t,conf)

    fig_xvel = plt.figure(2)
    ax_xvel = fig_xvel.add_subplot(111)
    ax_xvel.plot(t,vx,label="sim")
    ax_xvel.plot(t,vxref,label="sol")
    ax_xvel.set_xlim([0,t[-1]])
    ax_xvel.legend()

    fig_yvel = plt.figure(3)
    ax_yvel = fig_yvel.add_subplot(111)
    ax_yvel.plot(t,vy,label="sim")
    ax_yvel.plot(t,vyref,label="sol")
    ax_yvel.set_xlim([0,t[-1]])
    ax_yvel.legend()

    fig_xvel.savefig(name+'_xvelocity_{0}.png'.format(conf.Nt))
    fig_yvel.savefig(name+'_yvelocity_{0}.png'.format(conf.Nt))

    vx_error = abs(vxref[-1]-vx[-1])/abs(vxref[-1])
