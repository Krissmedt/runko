from __future__ import print_function
from configSetup import Configuration

from numpy import sqrt, pi
import numpy as np


class Configuration_Gyro(Configuration):

    def __init__(self, *file_names, do_print=False):
        Configuration.__init__(self, *file_names)

        #--------------------------------------------------
        # problem specific initializations
        if do_print:
            print("Initializing gyration setup...")

        self.threeD = False
        self.twoD = True

        self.M = 5
        self.K = 5

        #---------cold plasma-----------
        self.bphi=90.0  #Bfield z angle (bphi=0  bz, bphi=90 -> x-y plane)
        self.btheta=0.0   #Bfield x-y angle: btheta=0 -> parallel

        #--------------------------------------------------
        # particle initialization

        #local variables just for easier/cleaner syntax
        c   = self.cfl

	#plasma reaction & subsequent normalization
        self.omp=c/self.c_omp
        self.qe = 1

        #------------------------------------------------------------------
        # problem setup

        #particle 1, initialised to make perfect larmor orbit corresponding
        #to speed gamma.
        self.beta = sqrt(1-1/self.gamma**2.)
        self.vy = self.beta*c
        self.cfreq = (self.qe*self.binit)/(self.gamma*abs(self.qe)*c**2)
        # self.larmor = self.vy*self.gamma*c**2/self.binit
        self.larmor = self.vy/self.cfreq

        #particle 1, initialised to make half Larmor radius orbit of particle 1
        self.vy2 = 0.5*self.vy
        self.beta2 = self.vy2/c
        self.gamma2 = 1/np.sqrt(1-self.beta2**2)
        self.cfreq2 = (self.qe*self.binit)/(self.gamma2*abs(self.qe)*c**2)
        self.larmor2 = self.vy2/self.cfreq2

        self.uy = self.vy * self.gamma
        self.uy2 = self.vy2 * self.gamma2

        #--------------------------------------------------
        # field initialization

        self.Nx = 1
        self.Ny = 1
        self.Nz = 1
        self.NxMesh = 1
        self.NyMesh = 1
        self.NzMesh =1

        self.NxMesh = np.maximum(int(self.larmor*2.5),5)
        self.NyMesh = np.maximum(int(self.larmor*2.5),5)
        # self.NzMesh =1

        self.dx=1.0
        self.dy=1.0
        self.dz=1.0

        self.xmin = 0.0
        self.xmax = self.Nx*self.NxMesh*self.dx #XXX scaled length
        self.ymin = 0.0
        self.ymax = self.Ny*self.NyMesh*self.dy
        self.zmin = 0.0
        self.zmax = self.Nz*self.NzMesh*self.dz

        self.x_start = self.NxMesh/2. + self.larmor
        self.x_start2 = self.NxMesh/2. + self.larmor2
