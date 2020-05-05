import numpy as np
from py_runko_aux import *
from boris import boris_rp

def boris_synced_pos(tile):
    dt = tile.cfl
    cont = tile.get_container(0)

    pos = py_pos(cont)
    vel = py_vel(cont)
    E,B = py_em(cont)

    nq = pos.shape[0]
    dims = pos.shape[1]

    g = ginv(dt,vel*dt)

    vg = np.zeros(pos.shape,dtype=np.float)
    for i in range(0,dims):
        vg[:,i] = vel[:,i]*g
        pos[:,i] += dt*vg[:,i]

    pos += dt*dt/2 * (1*E + np.cross(vg,B))

    tile.delete_all_particles()

    for i in range(0,nq):
        cont.add_particle(pos[i,:],vel[i,:],1.0)


def boris_synced_vel(tile):
    dt = tile.cfl
    cont = tile.get_container(0)

    pos = py_pos(cont)
    vel = py_vel(cont)
    E,B = py_em(cont)

    nq = pos.shape[0]
    dims = pos.shape[1]

    vel = boris_rp(vel,E,B,dt,cont.q)

    tile.delete_all_particles()

    for i in range(0,nq):
        cont.add_particle(pos[i,:],vel[i,:],1.0)
