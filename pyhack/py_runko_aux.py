import numpy as np
import matplotlib.pyplot as plt
from pyhack.analytical_gyro import *
import h5py as h5

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

def gvi(c,v):
    ginv = np.sqrt(1 - np.sum(v**2/c**2,axis=1))
    return ginv

def gui(c,u):
    ginv = 1/np.sqrt(1 + np.sum(u**2/c**2,axis=1))
    return ginv

def F(vel,E,B,c,q):
    F = q/1 * (E + np.cross(vel*(gui(c,vel)[:,np.newaxis]/c),B))

    return F


def output_lf(t,x,y,vx,vy,conf,name):
    if conf.plot == 1:
        t = np.array(t)
        x = np.array(x)
        y = np.array(y)
        vx = np.array(vx)
        vy = np.array(vy)

        plot_traj(t,x,y,conf,name)
        t[1:] -= conf.cfl/2
        plot_vel(t,vx,vy,conf,name)


def output_vv(t,x,y,vx,vy,conf,name):
    if conf.plot == 1:
        t = np.array(t)
        x = np.array(x)
        y = np.array(y)
        vx = np.array(vx)
        vy = np.array(vy)

        plot_traj(t,x,y,conf,name)
        plot_vel(t,vx,vy,conf,name)

def output_sdc(t,x,y,vx,vy,xres,vres,conf,name):
    if conf.plot == 1:
        t = np.array(t)
        x = np.array(x)
        y = np.array(y)
        vx = np.array(vx)
        vy = np.array(vy)
        xres = np.array(xres,ndmin=2)
        vres = np.array(vres,ndmin=2)

        plot_traj(t,x,y,conf,name)
        plot_xres(t,xres,name)
        plot_vres(t,vres,name)



def plot_traj(t,x,y,conf,name):
    fig_traj = plt.figure(1)
    ax_traj = fig_traj.add_subplot(111)
    ax_traj.plot(x,y,label="sim")
    ax_traj.set_xlim([0,conf.NxMesh])
    ax_traj.set_ylim([0,conf.NyMesh])
    ax_traj.set_aspect('equal')
    ax_traj.legend()

    fig_traj.savefig(name+'_trajectory_{0}.png'.format(conf.Nt))


def plot_vel(t,vx,vy,conf,name):
    fig_xvel = plt.figure(2)
    ax_xvel = fig_xvel.add_subplot(111)
    ax_xvel.plot(t,vx,label="sim")
    ax_xvel.set_xlim([0,t[-1]])
    ax_xvel.legend()

    fig_yvel = plt.figure(3)
    ax_yvel = fig_yvel.add_subplot(111)
    ax_yvel.plot(t,vy,label="sim")
    ax_yvel.set_xlim([0,t[-1]])
    ax_yvel.legend()

    fig_xvel.savefig(name+'_xvelocity_{0}.png'.format(conf.Nt))
    fig_yvel.savefig(name+'_yvelocity_{0}.png'.format(conf.Nt))


def plot_xres(t,xres,name,label=""):
    fig_xres = plt.figure(5)
    ax_xres = fig_xres.add_subplot(111)
    for k in range(0,xres.shape[1]):
        ax_xres.plot(t[1:],xres[1:,k],label=label+" K={0}".format(k))
    ax_xres.set_xlim([0,t[-1]])
    ax_xres.set_yscale('log')
    ax_xres.legend()

    fig_xres.savefig(name+'_xres.png')


def plot_vres(t,vres,name,label=""):
    fig_vres = plt.figure(6)
    ax_vres = fig_vres.add_subplot(111)
    for k in range(0,vres.shape[1]):
        ax_vres.plot(t[1:],vres[1:,k],label=label+" K={0}".format(k))
    ax_vres.set_xlim([0,t[-1]])
    ax_vres.set_yscale('log')
    ax_vres.legend()

    fig_vres.savefig(name+'_vres.png')

def wp_dump(t,x,y,vx,vy,conf,filename):
    Nt = conf.Nt
    x = np.array(x[-1])
    y = np.array(y[-1])
    vx = np.array(vx[-1])
    vy = np.array(vy[-1])
    t = t[-1]

    if conf.ref_sim == 1:
        try:
            file = h5.File(conf.data_root+filename,'w')
        except OSError:
            file.close()
            file = h5.File(conf.data_root+filename,'w')

        grp = file.create_group('fields')
        grp.create_dataset('Nt',data=np.array(Nt)[np.newaxis],maxshape=(None,))
        grp.create_dataset('x',data=x[np.newaxis,:],maxshape=(None,None))
        grp.create_dataset('y',data=y[np.newaxis,:],maxshape=(None,None))
        grp.create_dataset('vx',data=vx[np.newaxis,:],maxshape=(None,None))
        grp.create_dataset('vy',data=vy[np.newaxis,:],maxshape=(None,None))
        grp.create_dataset('t',data=np.array(t)[np.newaxis],maxshape=(None,))
        file.close()

    elif conf.workprec == 1:
        try:
            file = h5.File(conf.data_root+filename,'r+')
        except OSError:
            file.close()
            file = h5.File(conf.data_root+filename,'r+')

        file["fields/Nt"].resize((file["fields/Nt"].shape[0]+1),axis=0)
        file["fields/Nt"][-1] = conf.Nt

        file["fields/x"].resize((file["fields/x"].shape[0]+1),axis=0)
        file["fields/x"][-1,:] = x

        file["fields/y"].resize((file["fields/y"].shape[0]+1),axis=0)
        file["fields/y"][-1,:] = y

        file["fields/vx"].resize((file["fields/vx"].shape[0]+1),axis=0)
        file["fields/vx"][-1,:] = vx

        file["fields/vy"].resize((file["fields/vy"].shape[0]+1),axis=0)
        file["fields/vy"][-1,:] = vy

        file["fields/t"].resize((file["fields/t"].shape[0]+1),axis=0)
        file["fields/t"][-1] = t

        file.close()


def orderLines(order,xRange,yRange):
    if order < 0:
        a = yRange[1]/xRange[0]**order
    else:
        a = yRange[0]/xRange[0]**order

    oLine = [a*xRange[0]**order,a*xRange[1]**order]
    return oLine
