import numpy as np
from pyhack.py_runko_aux import *
from pyhack.boris import boris_rp, boris_iter

def vv_pos(tile,dtf=1):
    c = tile.cfl
    cont = tile.get_container(0)

    pos = py_pos(cont)
    vel = py_vel(cont)*c
    E,B = py_em(cont)

    nq = pos.shape[0]
    dims = pos.shape[1]

    g = gui(c,vel)[:,np.newaxis]
    v_half = vel + dtf/2 * cont.q * (E+np.cross(vel*g/c,B))

    pos = pos + dtf*v_half*gui(c,v_half)[:,np.newaxis]

    tile.delete_all_particles()

    for i in range(0,nq):
        cont.add_particle(pos[i,:],vel[i,:]/c,1.0)

    E_old = np.copy(E)

    return E_old


def vv_vel(tile,dtf=1,E_old=0):
    c = tile.cfl
    cont = tile.get_container(0)

    pos = py_pos(cont)
    vel = py_vel(cont)*c
    E,B = py_em(cont)

    Eh = (E+E_old)/2
    nq = pos.shape[0]
    dims = pos.shape[1]
    # vel = boris_rp(vel,Eh,B,c,cont.q,dtf=dtf)
    ginv = gui(c,vel+dtf*0.5*(cont.q*Eh))
    vel = boris_iter(vel,Eh,B,dtf,c,0,ginv,cont.q)

    tile.delete_all_particles()

    for i in range(0,nq):
        cont.add_particle(pos[i,:],vel[i,:]/c,1.0)
