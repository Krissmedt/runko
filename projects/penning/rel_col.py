import numpy as np
import scipy.optimize as scop
from rvv_functions import *
from rvv_fields import *
from rvv_pushers import *
from gauss_legendre import CollGaussLegendre
from gauss_lobatto import CollGaussLobatto

class coll:
    def __init__(self,collclass,dt,nq,M=3,K=3,q=-1,**kwargs):
        self.collclass = collclass
        coll = self.collclass(M,0,1)

        self.K = K
        self.M = M

        self.nodes = coll._getNodes
        self.weights = coll._getWeights(coll.tleft,coll.tright) #Get M  nodes and weights

        self.Qmat = coll._gen_Qmatrix           #Generate q_(m,j), i.e. the large weights matrix
        self.Smat = coll._gen_Smatrix           #Generate s_(m,j), i.e. the large node-to-node weights matrix

        self.delta_m = coll._gen_deltas         #Generate vector of node spacings

        self.Qmat *= dt
        self.Smat *= dt
        self.delta_m *= dt

        self.ssi = 1

        self.nq = nq
        self.qe = q

        #Define required calculation matrices
        QE = np.zeros((M+1,M+1),dtype=np.float)
        QI = np.zeros((M+1,M+1),dtype=np.float)
        QT = np.zeros((M+1,M+1),dtype=np.float)

        SX = np.zeros((M+1,M+1),dtype=np.float)

        for i in range(0,M):
            QE[(i+1):,i] = self.delta_m[i]
            QI[(i+1):,i+1] = self.delta_m[i]

        QT = 1/2 * (QE + QI)
        QX = QE @ QT + (QE*QE)/2
        SX[:,:] = QX[:,:]
        SX[1:,:] = QX[1:,:] - QX[0:-1,:]

        self.SX = SX
        self.SQ = self.Smat @ self.Qmat

        self.x0 = np.zeros((M+1,nq,3),dtype=np.float)
        self.u0 = np.zeros((M+1,nq,3),dtype=np.float)

        self.xn = np.zeros((M+1,nq,3),dtype=np.float)
        self.un = np.zeros((M+1,nq,3),dtype=np.float)

        self.F = np.zeros((M+1,nq,3),dtype=np.float)
        self.Fn = np.zeros((M+1,nq,3),dtype=np.float)

        self.x_con = np.zeros((K,M))
        self.x_res = np.zeros((K,M))
        self.u_con = np.zeros((K,M))
        self.u_res = np.zeros((K,M))

        self.Rx = np.zeros((K,M),dtype=np.float)
        self.Rv = np.zeros((K,M),dtype=np.float)


    def calc_residual_2018(self,k):
        s = self
        q =  self.Qmat
        M = s.M

        for m in range(1,M+1):
            qvsum = 0
            qfsum = 0
            for j in range(1,M+1):
                qvsum += q[m,j] * s.u[j,:,:]
                qfsum += q[m,j] * s.F[j,:,:]

            s.Rx[k-1,m-1] = np.linalg.norm(s.x[0,:,:] + qvsum - s.x[m,:,:])
            s.Rv[k-1,m-1] = np.linalg.norm(s.u[0,:,:] + qfsum - s.u[m,:,:])


def implicit_coll(pos,vel,coll):
    M = coll.M
    K = coll.K
    nq = coll.nq
    #Remap collocation weights from [0,1] to [tn,tn+1]
    weights =  coll.weights
    q =  coll.Qmat

    x0 = np.ravel(pos)
    v0 = np.ravel(vel)

    Ix = np.array([1,0])
    Iv = np.array([0,1])
    Id = np.identity(nq*3)

    u0 = np.kron(Id,Ix).transpose() @ x0 + np.kron(Id,Iv).transpose() @ v0
    print(x0)
    print()
    print(v0)
    print()
    print(u0)

    # sol = scop.root_scalar(rootF,args=(coll,u0))
    return pos, vel, coll


def G(um):
    gamma = gu(um)
    vm = um/gamma[:,np.newaxis]
    return vm


def rootF(U):
    coll = args[0]
    U = args[1]

    f = U[um:] - coll.Ccoll @ U0 + coll.Qcoll @ FL(U)

    return f


def FL(U):
    FLU = np.zeros(U.shape,dtype=np.float)

    return FLU
