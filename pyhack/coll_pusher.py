import numpy as np
from pyhack.py_runko_aux_3d import *
import scipy.optimize as scop

def implicit_coll(tile,coll,fintp,timer):
    ## Retrieve runko stored data #############################################
    c = tile.cfl
    cont = tile.get_container(0)
    q = cont.q

    pos = py_pos(cont)
    vel = py_vel(cont)*c
    ###########################################################################

    M = coll.M
    nq = coll.nq
    #Remap collocation weights from [0,1] to [tn,tn+1]
    weights =  coll.weights

    for m in range(0,M+1):
        coll.x[m,:,:] = pos
        coll.u[m,:,:] = vel

    Id = coll.Id
    Ix = coll.Ix
    Iv = coll.Iv

    U0 = np.append(coll.x[1:,:,:].ravel(),coll.u[1:,:,:].ravel())
    FU = FXV(U0,coll,tile,fintp,timer)
    sol = scop.root(rootF,U0,args=(coll,U0,tile,fintp,timer),tol=10**-14,jac=False)
    U = sol.x
    md = M*coll.nq*3
    coll.x[1:,:,:] = U[0:md].reshape((M,coll.nq,3))
    coll.u[1:,:,:] = U[md:].reshape((M,coll.nq,3))

    pos = coll.x[-1,:,:]
    vel = coll.u[-1,:,:]

    ## Write to runko ##############################
    tile.delete_all_particles()
    for i in range(0,nq):
        cont.add_particle(pos[i,:],vel[i,:]/c,1.0)
    ################################################

    return pos, vel, coll


def rootF(U,*args):
    coll = args[0]
    U0 = args[1]
    tile = args[2]
    fintp = args[3]
    timer = args[4]

    f = U - coll.Q @ FXV(U,coll,tile,fintp,timer) - U0

    return f


def FXV(U,coll,tile,fintp,timer):
    M = coll.M
    Id = coll.Id
    Ix = coll.Ix
    Iv = coll.Iv
    md = M*coll.nq*3
    x = U[0:md].reshape((M,coll.nq,3))
    v = U[md:].reshape((M,coll.nq,3))
    Fx = np.zeros((M,coll.nq,3),dtype=np.float)
    Fv = np.zeros((M,coll.nq,3),dtype=np.float)
    for m in range(0,M):
        Fx[m,:,:] = v[m,:,:]*gui(coll.c,v[m,:,:])[:,np.newaxis]

        coll.E[m,:,:],coll.B[m,:,:] = interpolate(x[m,:,:],v[m,:,:],coll,tile,fintp,timer)
        Fv[m,:,:] = F(v[m,:,:],coll.E[m,:,:],coll.B[m,:,:] ,c=coll.c,q=coll.q)

    coll.F[1:,:,:] = Fv
    FXV = np.append(Fx.ravel(),Fv.ravel())
    return FXV


def interpolate(x,v,coll,tile,fintp,timer):
    timer.start_comp("interp_em")
    ## Write to runko ##############################
    pos = np.copy(x)
    vel = np.copy(v)
    tile.delete_all_particles()
    cont = tile.get_container(0)
    for i in range(0,coll.nq):
        cont.add_particle(pos[i,:],vel[i,:]/coll.c,1.0)
    ################################################
    fintp.solve(tile)
    cont = tile.get_container(0)
    E,B = py_em(cont)
    timer.stop_comp("interp_em")

    return E,B
