# -*- coding: utf-8 -*-

# system libraries
from __future__ import print_function
from mpi4py import MPI
import numpy as np
import sys, os
import matplotlib.pyplot as plt

# runko + auxiliary modules
import pytools  # runko python tools

# Runko-Python functionality by Krissmedt
from pyhack.py_runko_aux import *
from pyhack.analytical_gyro import *


# problem specific modules
np.random.seed(1)
from init_problem import Configuration_Gyro as Configuration

debug = False

def py_init(conf):
    xr,yr,vxr,vyr = analytical_gyro_single(-conf.cfl/2,conf)

    t = [0]
    x = [np.array([conf.x_start,conf.x_start2])]
    y = [np.array([conf.NyMesh/2,conf.NyMesh/2])]
    vx = [vxr]
    vy = [vyr]

    return t,x,y,vx,vy

def debug_print(n, msg):
    if debug:
        print("{}: {}".format(n.rank(), msg))
        sys.stdout.flush()


def filler(xloc, ispcs, conf):
    # perturb position between x0 + RUnif[0,1)

    #electrons
    if ispcs == 0:
        delgam  = conf.delgam #* np.abs(conf.mi / conf.me) * conf.temp_ratio

        xx = xloc[0] + np.random.rand(1)
        yy = xloc[1] + np.random.rand(1)
        #zz = xloc[2] + np.random.rand(1)
        zz = 0.5

    #positrons/ions/second species
    if ispcs == 1:
        delgam  = conf.delgam

        #on top of electrons
        xx = xloc[0]
        yy = xloc[1]
        zz = 0.5

    gamma = conf.gamma
    direction = -1
    ux, uy, uz, uu = boosted_maxwellian(delgam, gamma, direction=direction, dims=3)

    x0 = [xx, yy, zz]
    u0 = [ux, uy, uz]
    return x0, u0


def direct_inject(grid, conf):
    cid = grid.id(0,0)
    c = grid.get_tile(cid)
    container = c.get_container(0)

    xr,yr,vxr,vyr = analytical_gyro_single(-conf.cfl/2,conf)

    x = conf.x_start
    x2 = conf.x_start2
    y = conf.NyMesh/2.
    z = 0.5
    x01 = [x,y,z]
    x02 = [x2,y,z]

    vx = vxr[0]
    vx2 = vxr[1]
    vy = vyr[0]
    vy2 = vyr[1]
    vz = 0
    u01 = [vx,vy,vz]
    u02 = [vx2,vy2,vz]

    container.add_particle(x01,u01,1.0)
    container.add_particle(x02,u02,1.0)

    x0 = [x01,x02]
    u0 = [u01,u02]

    return x0,u0


# Field initialization (guide field)
def insert_em(grid, conf):

    #into radians
    btheta = conf.btheta/180.*np.pi
    bphi   = conf.bphi/180.*np.pi
    beta   = conf.beta

    kk = 0
    for cid in grid.get_tile_ids():
        tile = grid.get_tile(cid)
        yee = tile.get_yee(0)

        ii,jj,kk = tile.index if conf.threeD else (*tile.index, 0)

        for n in range(conf.NzMesh):
            for m in range(-3, conf.NyMesh+3):
                for l in range(-3, conf.NxMesh+3):
                    # get global coordinates
                    iglob, jglob, kglob = pytools.ind2loc((ii, jj, kk), (l, m, n), conf)

                    yee.bx[l,m,n] = 0. #conf.binit*np.cos(btheta)
                    yee.by[l,m,n] = 0. #conf.binit*np.sin(btheta)*np.sin(bphi)
                    yee.bz[l,m,n] = conf.binit #conf.binit*np.sin(btheta)*np.cos(bphi)

                    yee.ex[l,m,n] = 0.0
                    yee.ey[l,m,n] = 0. #-beta*yee.bz[l,m,n]
                    yee.ez[l,m,n] = 0. #beta*yee.by[l,m,n]


if __name__ == "__main__":

    x = []
    y = []
    vx = []
    vy = []
    t = []

    do_plots = True
    do_print = False

    if MPI.COMM_WORLD.Get_rank() == 0:
        do_print = False

    if do_print:
        print("Running with {} MPI processes.".format(MPI.COMM_WORLD.Get_size()))

    ##################################################
    # set up plotting and figure
    try:
        if do_plots:
            pass
    except:
        #print()
        pass


    # Timer for profiling
    timer = pytools.Timer()
    timer.start("total")
    timer.start("init")

    timer.do_print = do_print


    # parse command line arguments
    # parser = argparse.ArgumentParser(description='Simple PIC-Maxwell simulations')
    # parser.add_argument('--conf', dest='conf_filename', default=None,
    #                    help='Name of the configuration file (default: None)')
    args = pytools.parse_args()
    if args.conf_filename == None:
        conf = Configuration('gyration.ini', do_print=do_print)
    else:
        if do_print:
            print("Reading configuration setup from ", args.conf_filename)
        conf = Configuration(args.conf_filename, do_print=do_print)

    if conf.threeD:
        # 3D modules
        import pycorgi.threeD as pycorgi  # corgi ++ bindings
        import pyrunko.pic.threeD as pypic  # runko pic c++ bindings
        import pyrunko.fields.threeD as pyfld  # runko fld c++ bindings

    elif conf.twoD:
        # 2D modules
        import pycorgi.twoD as pycorgi  # corgi ++ bindings
        import pyrunko.pic.twoD as pypic  # runko pic c++ bindings
        import pyrunko.fields.twoD as pyfld  # runko fld c++ bindings

    grid = pycorgi.Grid(conf.Nx, conf.Ny, conf.Nz)

    xmin = 0.0
    xmax = conf.Nx*conf.NxMesh #XXX scaled length
    ymin = 0.0
    ymax = conf.Ny*conf.NyMesh
    grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax)

    # compute initial mpi ranks using Hilbert's curve partitioning
    pytools.balance_mpi(grid, conf)

    # load pic tiles into grid
    pytools.pic.load_tiles(grid, conf)

    ##################################################
    # create output folders
    if grid.master():
        pytools.create_output_folders(conf)

    # get current restart file status
    io_stat = pytools.check_for_restart(conf)

    if io_stat["do_initialization"]:
        if do_print:
            print("initializing simulation...")
        lap = 1

        np.random.seed(1)  # sync rnd generator seed for different mpi ranks

        # initialising solution arrays
        t,x,y,vx,vy = py_init(conf)
        # injecting plasma particles
        prtcl_stat = direct_inject(grid,conf) #inject plasma particles individually by loc,vel
        if do_print:
            print("injected:")
            print("    e- prtcls: {}".format(prtcl_stat[0]))
            print("    e+ prtcls: {}".format(prtcl_stat[1]))

        # inserting em grid
        insert_em(grid, conf)

    else:
        if do_print:
            print("restarting simulation from lap {}...".format(io_stat["lap"]))

        # read restart files
        pyfld.read_yee(grid, io_stat["read_lap"], io_stat["read_dir"])
        pypic.read_particles(grid, io_stat["read_lap"], io_stat["read_dir"])

        # step one step ahead
        lap = io_stat["lap"] + 1



    #static load balancing setup; communicate neighbor info once
    debug_print(grid, "analyze bcs")
    grid.analyze_boundaries()
    debug_print(grid, "send tiles")
    grid.send_tiles()
    debug_print(grid, "recv tiles")
    grid.recv_tiles()
    MPI.COMM_WORLD.barrier()

    #sys.exit()

    debug_print(grid, "init virs")
    pytools.pic.load_virtual_tiles(grid, conf)


    timer.stop("init")
    timer.stats("init")


    # end of initialization
    ##################################################
    debug_print(grid, "solvers")


    # visualize initial condition
    if do_plots:
        try:
            plotNode( axs[0], grid, conf)
            #plotXmesh(axs[1], grid, conf, 0, "x")
            saveVisz(-1, grid, conf)
        except:
            pass


    Nsamples = conf.Nt
    #pusher   = pypic.BorisPusher()
    pusher   = pypic.VayPusher()


    fldprop  = pyfld.FDTD2()
    # fldprop  = pyfld.FDTD4()
    fintp    = pypic.LinearInterpolator()
    currint  = pypic.ZigZag()
    # analyzer = pypic.Analyzator()
    flt      = pyfld.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh)

    #enhance numerical speed of light slightly to suppress numerical Cherenkov instability
    fldprop.corr = 1.02

    debug_print(grid, "mpi_e")
    grid.send_data(1)
    grid.recv_data(1)
    grid.wait_data(1)

    debug_print(grid, "mpi_b")
    grid.send_data(2)
    grid.recv_data(2)
    grid.wait_data(2)

    for tile in pytools.tiles_all(grid):
        tile.update_boundaries(grid)

    ##################################################
    sys.stdout.flush()



#-----------------------------------------------------------------------------#
################################ Simulation Loop ##############################
#-----------------------------------------------------------------------------#

    time = lap*(conf.cfl/conf.c_omp)
    for lap in range(lap, conf.Nt+1):
        debug_print(grid, "lap_start")


        #--------------------------------------------------
        # comm B
        timer.start_comp("mpi_b1")
        debug_print(grid, "mpi_b1")

        grid.send_data(2)
        grid.recv_data(2)
        grid.wait_data(2)

        timer.stop_comp("mpi_b1")

        #--------------------------------------------------
        #update boundaries
        timer.start_comp("upd_bc0")
        debug_print(grid, "upd_bc0")

        for tile in pytools.tiles_all(grid):
            tile.update_boundaries(grid)

        timer.stop_comp("upd_bc0")

        #--------------------------------------------------
        #push B half
        timer.start_comp("push_half_b1")
        debug_print(grid, "push_half_b1")

        for tile in pytools.tiles_all(grid):
            fldprop.push_half_b(tile)

        timer.stop_comp("push_half_b1")

        #--------------------------------------------------
        # comm B
        timer.start_comp("mpi_b2")
        debug_print(grid, "mpi_b2")

        grid.send_data(2)
        grid.recv_data(2)
        grid.wait_data(2)

        timer.stop_comp("mpi_b2")

        #--------------------------------------------------
        #update boundaries
        timer.start_comp("upd_bc1")
        debug_print(grid, "upd_bc1")

        for tile in pytools.tiles_all(grid):
            tile.update_boundaries(grid)

        timer.stop_comp("upd_bc1")


        ##################################################
        # move particles (only locals tiles)

        #--------------------------------------------------
        #interpolate fields (can move to next asap)
        timer.start_comp("interp_em")
        debug_print(grid, "interp_em")

        for tile in pytools.tiles_local(grid):
            fintp.solve(tile)

        timer.stop_comp("interp_em")
        #--------------------------------------------------

        #--------------------------------------------------
        #push particles in x and u
        timer.start_comp("push")
        debug_print(grid, "push")

        for tile in pytools.tiles_local(grid):
            pusher.solve(tile)

        timer.stop_comp("push")


        # advance B half

        #--------------------------------------------------
        #push B half
        timer.start_comp("push_half_b2")
        debug_print(grid, "push_half_b2")

        for tile in pytools.tiles_all(grid):
            fldprop.push_half_b(tile)

        timer.stop_comp("push_half_b2")


        #--------------------------------------------------
        # comm B
        timer.start_comp("mpi_e1")
        debug_print(grid, "mpi_e1")

        grid.send_data(1)
        grid.recv_data(1)
        grid.wait_data(1)

        timer.stop_comp("mpi_e1")

        #--------------------------------------------------
        #update boundaries
        timer.start_comp("upd_bc2")
        debug_print(grid, "upd_bc2")

        for tile in pytools.tiles_all(grid):
            tile.update_boundaries(grid)

        timer.stop_comp("upd_bc2")


        ##################################################
        # advance E

        #--------------------------------------------------
        #push E
        timer.start_comp("push_e")
        debug_print(grid, "push_e")

        # for cid in grid.get_tile_ids():
        #     tile = grid.get_tile(cid)
        for tile in pytools.tiles_all(grid):
            fldprop.push_e(tile)

        timer.stop_comp("push_e")


        #Current calculations
        # # --------------------------------------------------
        # # current calculation; charge conserving current deposition
        # t1 = timer.start_comp("comp_curr")
        # for tile in pytools.tiles_local(grid):
        #     currint.solve(tile)
        # timer.stop_comp(t1)
        #
        # # --------------------------------------------------
        # # clear virtual current arrays for boundary addition after mpi
        # t1 = timer.start_comp("clear_vir_cur")
        # for tile in pytools.tiles_virtual(grid):
        #     tile.clear_current()
        # timer.stop_comp(t1)
        #
        # # --------------------------------------------------
        # # mpi send currents
        # t1 = timer.start_comp("mpi_cur")
        # grid.send_data(0)
        # grid.recv_data(0)
        # grid.wait_data(0)
        # timer.stop_comp(t1)
        #
        # # --------------------------------------------------
        # # exchange currents
        # t1 = timer.start_comp("cur_exchange")
        # for tile in pytools.tiles_all(grid):
        #     tile.exchange_currents(grid)
        # timer.stop_comp(t1)



        ##################################################
        # particle communication (only local/boundary tiles)

        #--------------------------------------------------
        #local particle exchange (independent)
        timer.start_comp("check_outg_prtcls")
        debug_print(grid, "check_outg_prtcls")

        for tile in pytools.tiles_local(grid):
            tile.check_outgoing_particles()

        timer.stop_comp("check_outg_prtcls")

        #--------------------------------------------------
        # global mpi exchange (independent)
        timer.start_comp("pack_outg_prtcls")
        debug_print(grid, "pack_outg_prtcls")

        for tile in pytools.tiles_boundary(grid):
            tile.pack_outgoing_particles()

        timer.stop_comp("pack_outg_prtcls")

        # --------------------------------------------------
        # MPI global particle exchange
        # transfer primary and extra data
        t1 = timer.start_comp("mpi_prtcls")
        grid.send_data(3)
        grid.recv_data(3)
        grid.wait_data(3)

        # orig just after send3
        grid.send_data(4)
        grid.recv_data(4)
        grid.wait_data(4)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # global unpacking (independent)
        t1 = timer.start_comp("unpack_vir_prtcls")
        for tile in pytools.tiles_virtual(grid):
            tile.unpack_incoming_particles()
            tile.check_outgoing_particles()
        timer.stop_comp(t1)

        # --------------------------------------------------
        # transfer local + global
        t1 = timer.start_comp("get_inc_prtcls")
        for tile in pytools.tiles_local(grid):
            tile.get_incoming_particles(grid)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # delete local transferred particles
        t1 = timer.start_comp("del_trnsfrd_prtcls")
        for tile in pytools.tiles_local(grid):
            tile.delete_transferred_particles()
        timer.stop_comp(t1)

        # --------------------------------------------------
        # delete all virtual particles (because new prtcls will come)
        t1 = timer.start_comp("del_vir_prtcls")
        for tile in pytools.tiles_virtual(grid):
            tile.delete_all_particles()
        timer.stop_comp(t1)

        ##################################################
        # filter
        timer.start_comp("filter")

        #sweep over npasses times
        for fj in range(conf.npasses):

            #update global neighbors (mpi)
            grid.send_data(0)
            grid.recv_data(0)
            grid.wait_data(0)

            #get halo boundaries
            for cid in grid.get_local_tiles():
                tile = grid.get_tile(cid)
                tile.update_boundaries(grid)

            #filter each tile
            for cid in grid.get_local_tiles():
                tile = grid.get_tile(cid)
                flt.solve(tile)

            MPI.COMM_WORLD.barrier() # sync everybody


        # --------------------------------------------------
        timer.stop_comp("filter")

        # --------------------------------------------------
        # add current to E
        # t1 = timer.start_comp("add_cur")
        # for tile in pytools.tiles_all(grid):
        #     tile.deposit_current()
        # timer.stop_comp(t1)

        # comm E
        t1 = timer.start_comp("mpi_e2")
        grid.send_data(1)
        grid.recv_data(1)
        grid.wait_data(1)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # comm B
        t1 = timer.start_comp("mpi_b1")
        grid.send_data(2)
        grid.recv_data(2)
        grid.wait_data(2)
        timer.stop_comp(t1)

        # --------------------------------------------------
        # update boundaries
        t1 = timer.start_comp("upd_bc0")
        for tile in pytools.tiles_all(grid):
            tile.update_boundaries(grid)
        timer.stop_comp(t1)

        ##################################################
        # data reduction and I/O

        cid = grid.id(0,0)
        c = grid.get_tile(cid)
        container = c.get_container(0)

        x.append(container.loc(0))
        y.append(container.loc(1))
        vx.append(container.vel(0))
        vy.append(container.vel(1))
        t.append(time)


        timer.lap("step")
        if (lap % conf.interval == 0):
            debug_print(grid, "io")
            if do_print:
                print("--------------------------------------------------")
                print("------ lap: {} / t: {}".format(lap, time))

            print("------------------------------------------------------")
            print("x-position:" + str(x[lap]))
            print("y-position:" + str(y[lap]))
            print("x-vel:" + str(vx[lap]))
            print("y-vel:" + str(vy[lap]))
            print("------------------------------------------------------")

            #for cid in grid.get_tile_ids():
            #    tile = grid.get_tile(cid)
            #    tile.erase_temporary_arrays()

            timer.stats("step")
            timer.comp_stats()
            timer.purge_comps()

            #analyze (independent)
            timer.start("io")


            #--------------------------------------------------
            #2D plots
            if do_plots:
                try:
                    pass

                except:
                    #print()
                    pass
            timer.stop("io")


            timer.stats("io")
            timer.start("step") #refresh lap counter (avoids IO profiling)

            sys.stdout.flush()

        #next step
        time += conf.cfl/conf.c_omp
    #end of loop

    timer.stop("total")
    timer.stats("total")

    output_stag(t,x,y,vx,vy,conf,'rk_' + conf.name + '_')

    print("")
    print("------------------------------------- END ------------------------------------")
