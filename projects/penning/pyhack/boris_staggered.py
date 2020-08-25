import numpy as np

def boris_staggered(tile):
    c = tile.cfl
    cont = tile.get_container(0)

    pos = py_pos(cont)
    vel = py_vel(cont)
    E,B = py_em(cont)

    nq = pos.shape[0]
    dims = pos.shape[1]

    vel = boris_rp(vel,E,B,dt,1)

    g = ginv(c,vel*c)

    pos[:,i] += vel[:,i]*g[:,np.newaxis]*c


    tile.delete_all_particles()

    for i in range(0,nq):
        cont.add_particle(pos[i,:],vel[i,:],1.0)
