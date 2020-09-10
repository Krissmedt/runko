import numpy as np
from pyhack.py_runko_aux import *
from pyhack.boris import boris_rp

def boris_synced_pos(tile,dtf=1):
    c = tile.cfl
    cont = tile.get_container(0)

    pos = py_pos(cont)
    vel = py_vel(cont)
    E,B = py_em(cont)

    nq = pos.shape[0]
    dims = pos.shape[1]

    g = ginv(c,vel*c)

    vg = vel*g[:,np.newaxis]
    pos = pos + dtf*c*vg + dtf**2*g[:,np.newaxis]/2*cont.q*(E + np.cross(vg,B))

    tile.delete_all_particles()

    for i in range(0,nq):
        cont.add_particle(pos[i,:],vel[i,:],1.0)


def boris_synced_vel(tile,dtf=1):
    c = tile.cfl
    cont = tile.get_container(0)

    pos = py_pos(cont)
    vel = py_vel(cont)
    E,B = py_em(cont)

    nq = pos.shape[0]
    dims = pos.shape[1]

    vel = boris_rp(vel,E,B,c,cont.q,dtf=dtf)

    tile.delete_all_particles()

    for i in range(0,nq):
        cont.add_particle(pos[i,:],vel[i,:],1.0)
