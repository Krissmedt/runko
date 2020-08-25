import numpy as np
from pyhack.py_runko_aux import *
from pyhack.boris import *

def boris_staggered(tile,dtf=1):
    c = tile.cfl
    cont = tile.get_container(0)

    pos = py_pos(cont)
    vel = py_vel(cont)
    E,B = py_em(cont)

    nq = pos.shape[0]
    dims = pos.shape[1]

    vel = boris_rp(vel,E,B,c,cont.q,dtf=dtf)

    g = ginv(c,vel*c)

    for i in range(0,dims):
        pos[:,i] += dtf*c*vel[:,i]*g


    tile.delete_all_particles()

    for i in range(0,nq):
        cont.add_particle(pos[i,:],vel[i,:],1.0)


def boris_staggered_first(tile,dtf=1):
    c = tile.cfl
    cont = tile.get_container(0)

    pos = py_pos(cont)
    vel = py_vel(cont)
    E,B = py_em(cont)

    nq = pos.shape[0]
    dims = pos.shape[1]

    vel = boris_rp(vel,E,B,c,cont.q,dtf=0.5*dtf)

    g = ginv(c,vel*c)

    for i in range(0,dims):
        pos[:,i] += dtf*c*vel[:,i]*g


    tile.delete_all_particles()

    for i in range(0,nq):
        cont.add_particle(pos[i,:],vel[i,:],1.0)
