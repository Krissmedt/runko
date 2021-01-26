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

                    yee.ex[l,m,n] = -(conf.NxMesh/2.-iglob) * conf.einit
                    yee.ey[l,m,n] = -(conf.NyMesh/2.-jglob) * conf.einit #-beta*yee.bz[l,m,n]
                    yee.ez[l,m,n] = 2*(conf.NzMesh/2.-kglob) * conf.einit
