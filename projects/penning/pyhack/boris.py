import numpy as np
import time
from pyhack.py_runko_aux import *


def boris_rp(vel, E, B, c, qm, ck=0,dtf=1):
    dims = vel.shape[1]

    cinv = 1./c

    E0 = dtf*0.5*(qm*E + ck)
    B0 = dtf*0.5*qm*cinv*B

    u1 = vel + E0

    g1 = gui(c,u1)

    for i in range(0,dims):
        B0[:,i] *= g1

    f = 2./(1. + np.sum(B0*B0,axis=1))

    u2 = u1 + np.cross(u1,B0)

    for i in range(0,dims):
        u2[:,i] *= f

    u3 = u1 + np.cross(u2,B0) + E0
    vel = u3

    return vel


def boris_iter(u,E,B,dt,c,ck,ginv,q):
    t = 0.5*dt*q*B*ginv/c
    s = 2.0*t/(1.0 + np.linalg.norm(t**2,axis=1)[:,np.newaxis])
    u_min = u + 0.5*dt*q*E + 0.5*ck
    u_star  = u_min + np.cross(u_min, t)
    u_plus  = u_min + np.cross(u_star, s)
    return u_plus + 0.5*dt*q*E + 0.5*ck


def boris_rp2(vel, E, B, c, qm, ck=0):
    for pii in range(0,vel.shape[0]):
        cinv = 1./c

        E0 = 0.5*(qm*E[pii,:] + ck)
        B0 = 0.5*qm*cinv*B[pii,:]

        u0 = c*vel[pii,:] + E0

        g1 = c/np.sqrt(c*c + u0[0]*u0[0] + u0[1]*u0[1] + u0[2]*u0[2])

        B0 *= g1

        f = 2./(1. + B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2])

        u1 = (u0 + np.cross(u0,B0))*f

        u0 = u0 + np.cross(u1,B0) + E0
        vel[pii,:] = u0*cinv

    return vel
