# -*- coding: utf-8 -*-

# system libraries
from __future__ import print_function
from mpi4py import MPI
import numpy as np
import sys, os
import matplotlib.pyplot as plt
import time

# runko + auxiliary modules
import pytools  # runko python tools

# Runko-Python functionality by Krissmedt
from pyhack.coll_setup import coll
from pyhack.coll_pusher import *
from pyhack.py_runko_aux_3d import *

# problem specific modules
np.random.seed(1)
from init_problem import Configuration_Gyro as Configuration

debug = False

def py_init(conf):
    t = [0]
    x = [np.array([conf.x_start])]
    y = [np.array([conf.NyMesh/2.])]
    z = [np.array([conf.NzMesh/10.])]
    vx = [np.array([conf.ux])]
    vy = [np.array([conf.uy])]
    vz = [np.array([conf.uz])]
    xres = [np.zeros(1,dtype=np.float)]
    vres = [np.zeros(1,dtype=np.float)]

    return t,x,y,z,vx,vy,vz,xres,vres


def debug_print(n, msg):
    if debug:
        print("{}: {}".format(n.rank(), msg))
        sys.stdout.flush()

def direct_inject(grid, conf):
    cid = grid.id(0,0,0)
    c = grid.get_tile(cid)
    container = c.get_container(0)

    x = conf.x_start
    y = conf.NyMesh/2.
    z = conf.NzMesh/10.
    x01 = [x,y,z]

    vx = conf.ux
    vy = conf.uy
    vz = conf.uz
    u01 = [vx,vy,vz]

    container.add_particle(x01,u01,1.0)
    # container.add_particle(x02,u02,1.0)

    x0 = [x01]
    u0 = [u01]

    return x0,u0


# Field initialization (guide field)
def insert_em(grid, conf):

    #into radians
    btheta = conf.btheta/180.*np.pi
    bphi   = conf.bphi/180.*np.pi

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

                    yee.ex[l,m,n] = 0.
                    yee.ey[l,m,n] = 0. #-beta*yee.bz[l,m,n]
                    yee.ez[l,m,n] = conf.einit #beta*yee.by[l,m,n]


if __name__ == "__main__":

    do_plots = True
    do_print = True

    if MPI.COMM_WORLD.Get_rank() == 0:
        do_print = True

    if do_print:
        print("")
        print("")
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
    grid.set_grid_lims(conf.xmin, conf.xmax, conf.ymin, conf.ymax,conf.zmin,conf.zmax)

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
        t,x,y,z,vx,vy,vz,xres,vres = py_init(conf)
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
    pushloc   = pypic.VerletLocPusher()
    pushvel   = pypic.VerletVelPusher()


    fldprop  = pyfld.FDTD2()
    # fldprop  = pyfld.FDTD4()
    fintp    = pypic.LinearInterpolator()
    currint  = pypic.ZigZag()
    # analyzer = pypic.Analyzator()
    flt      = pyfld.Binomial2(conf.NxMesh, conf.NyMesh, conf.NzMesh)

    #enhance numerical speed of light slightly to suppress numerical Cherenkov instability
    fldprop.corr = 1.0

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
    col = coll(tile,dtf=conf.dtf,M=conf.M,K=0)

    time = lap*(conf.dtf*conf.cfl/conf.c_omp)
    for lap in range(lap, conf.Nt+1):
        debug_print(grid, "lap_start")
        #--------------------------------------------------
        #push particles
        timer.start_comp("push")
        debug_print(grid, "push")

        for tile in pytools.tiles_local(grid):
            implicit_coll(tile,col,fintp,timer)
            col.calc_residual(0)

        timer.stop_comp("push")

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

        ##################################################
        # data reduction and I/O
        cid = grid.id(0,0,0)
        c = grid.get_tile(cid)
        container = c.get_container(0)

        t.append(time)
        x.append(container.loc(0))
        y.append(container.loc(1))
        z.append(container.loc(2))
        vx.append(container.vel(0))
        vy.append(container.vel(1))
        vz.append(container.vel(2))
        xres.append(np.linalg.norm(col.Rx,axis=1))
        vres.append(np.linalg.norm(col.Rv,axis=1))

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
        time += conf.dtf*conf.cfl/conf.c_omp

    #end of loop

    timer.stop("total")
    timer.stats("total")

    output_sdc(t,x,y,z,vx,vy,vz,xres,vres,conf,'coll_py_' + conf.name + '_')

    filename = "coll_M{0}_{1}.h5".format(conf.M,conf.name)
    wp_dump(t,x,y,z,vx,vy,vz,conf,filename)

    print("")
    print("------------------------------------- END ------------------------------------")
