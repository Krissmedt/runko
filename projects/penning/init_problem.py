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

        self.twoD = False
        self.threeD = True

        self.M = 3
        self.K = 4

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
        self.vx = self.beta_x*c
        self.vy = self.beta_y*c
        self.vz = self.beta_z*c

        self.gamma = np.sqrt(1/(1-(self.vx**2+self.vy**2+self.vz**2)/c**2))

        self.ux = self.vx*self.gamma
        self.uy = self.vy*self.gamma
        self.uz = self.vz*self.gamma

        self.cfreq = (self.qe*self.binit)/(self.gamma*abs(self.qe)*c**2)
        self.larmor = self.vy/self.cfreq

        #--------------------------------------------------
        # field initialization

        self.Nx = 1
        self.Ny = 1
        self.Nz = 1
        # self.NxMesh = 10
        # self.NyMesh = 10
        # self.NzMesh = 10

        # self.NxMesh = np.maximum(int(self.larmor*2.25),5)
        # self.NyMesh = np.maximum(int(self.larmor*2.25),5)
        # self.NzMesh = np.maximum(int(self.larmor*2.25),5)

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
