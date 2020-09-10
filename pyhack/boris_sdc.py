import numpy as np
from pyhack.py_runko_aux_3d import *
from pyhack.boris import boris_iter

def boris_sdc_init(tile,coll):
    ## Retrieve runko stored data #############################################
    c = tile.cfl
    cont = tile.get_container(0)
    q = cont.q

    pos = py_pos(cont)
    vel = py_vel(cont)*c
    E,B = py_em(cont)
    ###########################################################################

    ## Retrieve Python stored data
    M = coll.M
    K = coll.K
    nq = coll.nq

    #Remap collocation weights from [0,1] to [tn,tn+1]
    weights =  coll.weights
    Q =  coll.Qmat
    dm =  coll.delta_m

    ## Populate node solutions with x0, v0, F0 ##
    coll.x[0,:,:] = pos
    coll.u[0,:,:] = vel
    coll.E[0,:,:] = E
    coll.B[0,:,:] = B

    coll.xn[0,:,:] = coll.x[0,:,:]
    coll.un[0,:,:] = coll.u[0,:,:]
    coll.En[0,:,:] = coll.E[0,:,:]
    coll.Bn[0,:,:] = coll.B[0,:,:]

    coll.IV = 0
    coll.IF = 0

    if coll.predictor == False:
        for m in range(0,M):
            coll.x[m+1,:,:] = coll.x[m,:,:]
            coll.u[m+1,:,:] = coll.u[m,:,:]
            coll.E[m+1,:,:] = coll.E[m,:,:]
            coll.B[m+1,:,:] = coll.B[m,:,:]

    # if coll.predictor == True:
    #     ###### Initial Step #########################
    #     v_half = vel + 0.5*dm[0]*F(vel,E(pos,q=qe),B(pos,q=qe),c=c)
    #     coll.x[1,:,:] = pos + dm[0]*G(v_half,c=c)
    #
    #     En         = 0.5*(E(pos) + E(coll.x[1,:,:]))*qe
    #     Bn         = B(coll.x[1,:,:])*qe
    #     gamma      = gu(coll.u[0,:,:],c=c)
    #     c_1        = 0.5*dm[0]*np.cross(G(coll.u[0,:,:],c=c), B(coll.x[0,:,:]))*qe
    #     c_2        = -(0.5*dm[0]/gamma)*np.cross(coll.u[0,:,:], Bn) + c_1
    #     coll.u[1,:,:] = boris_daniel(coll.u[0,:,:],En,Bn,dm[0],c_2,gamma,q=1)
    #     coll.F[1,:,:] = F(coll.u[0,:,:],E(coll.x[0,:,:]),B(coll.x[0,:,:]),c=c)
    #
    #
    #         coll.x[1,:,:] = coll.x[0,:,:]
    #         coll.u[1,:,:] = coll.u[0,:,:]
    #         coll.F[1,:,:] = coll.F[0,:,:]
    #
    #     ############################################
    #     ######## Predictor Step ####################
    #     for m in range(1,M):
    #         v_half = coll.u[m,:,:] + 0.5*dm[m]*coll.F[m,:,:]
    #         coll.x[m+1,:,:] = coll.x[m,:,:] + dm[m]*G(v_half,c=c)
    #
    #         En         = 0.5*(E(coll.x[m,:,:]) + E(coll.x[m+1,:,:]))*qe
    #         Bn         = B(coll.x[m,:,:])*qe
    #         gamma      = gu(coll.u[m,:,:],c=c)
    #         c_1        = 0.5*dm[m]*np.cross(G(coll.u[m,:,:],c=c), B(coll.x[m,:,:]))*qe
    #         c_2        = -(0.5*dm[m]/gamma)*np.cross(coll.u[m,:,:], Bn) + c_1
    #         coll.u[m+1,:,:] = boris_daniel(coll.u[m,:,:],En,Bn,dm[m],c_2,gamma,q=1)
    #         coll.F[m+1,:,:] = F(coll.u[m+1,:,:],E(coll.x[m+1,:,:]),B(coll.x[m+1,:,:]),c=c)

    coll.calc_residual(0)


def boris_sdc_pos(tile,coll,m):
    ## Retrieve runko stored data #############################################
    c = tile.cfl
    cont = tile.get_container(0)
    q = cont.q

    pos = py_pos(cont)
    vel = py_vel(cont)*c
    coll.En[m,:,:], coll.Bn[m,:,:] = py_em(cont)

    ###########################################################################
    M = coll.M
    K = coll.K

    nq = coll.nq
    c = coll.c
    q = coll.q

    #Remap collocation weights from [0,1] to [tn,tn+1]
    weights =  coll.weights
    Q =  coll.Qmat
    dm =  coll.delta_m

    # Calculate collocation terms required for pos update
    coll.IV = 0
    for j in range(1,M+1):
        coll.IV += (Q[m+1,j]-Q[m,j])*coll.u[j,:,:]*gui(c,coll.u[j,:,:])[:,np.newaxis]

    v_half = coll.u[m,:,:] + 0.5*dm[m]*F(coll.u[m,:,:],coll.E[m,:,:],coll.B[m,:,:],c,q)
    vn_half = coll.un[m,:,:] + 0.5*dm[m]*F(coll.un[m,:,:],coll.En[m,:,:],coll.Bn[m,:,:],c,q)

    ### POSITION UPDATE FOR NODE m/SWEEP k ###
    coll.xn[m+1,:,:] = coll.xn[m,:,:]
    coll.xn[m+1,:,:] += dm[m]* (vn_half*gui(c,vn_half)[:,np.newaxis]-v_half*gui(c,v_half)[:,np.newaxis])
    coll.xn[m+1,:,:] += coll.IV

    ## Write to runko ##############################
    pos = np.copy(coll.xn[m+1,:,:])
    tile.delete_all_particles()
    for i in range(0,nq):
        cont.add_particle(pos[i,:],vel[i,:]/c,1.0)
    ################################################

    return pos, vel, coll



def boris_sdc_vel(tile,coll,m):
    ## Retrieve runko stored data #############################################
    c = tile.cfl
    cont = tile.get_container(0)
    q = cont.q

    pos = py_pos(cont)
    vel = py_vel(cont)*c
    coll.En[m+1,:,:], coll.Bn[m+1,:,:] = py_em(cont)
    ###########################################################################

    M = coll.M
    K = coll.K

    nq = coll.nq
    c = coll.c
    q = coll.q

    #Remap collocation weights from [0,1] to [tn,tn+1]
    weights =  coll.weights
    Q =  coll.Qmat
    dm =  coll.delta_m

    # Calculate collocation terms required for pos update
    coll.IF = 0
    for j in range(1,M+1):
        coll.IF += (Q[m+1,j]-Q[m,j])*F(coll.u[j,:,:],coll.E[j,:,:],coll.B[j,:,:],c,q)

    Eh         = 0.5*(coll.En[m+1,:,:] + coll.En[m,:,:])
    ginv      = gui(c,coll.u[m+1,:,:])[:,np.newaxis]

    c_1        = 0.5*q*dm[m]/c * np.cross(coll.un[m,:,:]*gui(c,coll.un[m,:,:])[:,np.newaxis], coll.Bn[m,:,:])
    c_1       += -0.5*dm[m]* (F(coll.u[m,:,:],coll.E[m,:,:],coll.B[m,:,:],c,q)
                            + F(coll.u[m+1,:,:],coll.E[m+1,:,:],coll.B[m+1,:,:],c,q))
    c_1       += coll.IF
    ck        = c_1 -(0.5*dm[m]*q*ginv/c)*np.cross(coll.un[m,:,:], coll.Bn[m+1,:,:])
    coll.un[m+1,:,:] = boris_iter(coll.un[m,:,:],Eh,coll.Bn[m+1,:,:],dm[m],c,ck,ginv,q)

    ## Write to runko ##############################
    vel = np.copy(coll.un[m+1,:,:])
    tile.delete_all_particles()
    for i in range(0,nq):
        cont.add_particle(pos[i,:],vel[i,:]/c,1.0)
    ################################################

    return pos, vel, coll
