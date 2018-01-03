#include <vector>
#include <cmath>
#include <stdexcept>

#include "../solvers.h"
#include "../sheets.h"


namespace vlasov {



  /*! \brief Splitted Lagrangian spatial solver for Vlasov fluids
   *
   * Full solution is obtained by splitting each spatial dimension
   * x/y/z into separate solutions that are then cycled over to
   * obtain full 3D solution.
   *
   * Uses Sheets to vectorize calculations.
   *
   *
   */
  class SpatialLagrangianSolver2nd : public VlasovSpatialSolver {

    private:
      /*! \brief Volume of sheet
       * reduces sheet into volume weighted sum
       */
      Realf volWeightedSum(sheets::Sheet s) {
        Realf ret = 0.0;

        // TODO: for 2/3D we need to multiply with s.diff()
        // sheets::Sheet volWeighted = s * s.diff();
        sheets::Sheet volWeighted = s;  

        /// now reduce into sum
        ret = volWeighted.sum();

        return ret;
      };



    public:
      void solve() {

        // size_t Nintp = 4;

        // get target cell that we operate on
        // TODO fix solver interface to take cell AND node as inputs; 
        // this will avoid pointer casting?
        uint64_t cid = grid->cellId(targeti, targetj);
        Grid::CellPtr cellPtr = 
          std::dynamic_pointer_cast<VlasovCell>( grid->getCellPtr(cid));


        // Get reference to the velocity mesh we are solving
        VlasovFluid& gr0 = cellPtr->getPlasmaGrid();

        // get reference to a new mesh that we are updating into
        VlasovFluid& gr1 = cellPtr->getNewPlasmaGrid();


        // numerical flux from moving the Vlasov fluid around
        maxwell::YeeLattice& yee = cellPtr->getYee();
        yee.jx.clear();
        yee.jy.clear();
        yee.jz.clear();


        // TODO only x dir update is done here; multidimensionalize
        for (size_t dim=0; dim<1; dim++) {
          // fmt::print("Solving for dim {}\n", dim);

          // get neighbors for interpolation
          auto nindx_p1    = cellPtr->neighs(+1, 0); // i+1 neighbor
          uint64_t cid_p1  = grid->cellId( std::get<0>(nindx_p1), std::get<1>(nindx_p1) );
          Grid::CellPtr cellPtr_p1 = 
            std::dynamic_pointer_cast<VlasovCell>( grid->getCellPtr(cid_p1));

          VlasovFluid& gr_p1 = cellPtr_p1->getPlasmaGrid();


          auto nindx_m1    = cellPtr->neighs(-1, 0); // i-1 neighbor
          uint64_t cid_m1  = grid->cellId( std::get<0>(nindx_m1), std::get<1>(nindx_m1) );
          Grid::CellPtr cellPtr_m1 = 
            std::dynamic_pointer_cast<VlasovCell>( grid->getCellPtr(cid_m1));

          VlasovFluid&     gr_m1 = cellPtr_m1->getPlasmaGrid();

          double qmE, qmP;

          if (gr0.Nx == 1 && gr0.Ny == 1 && gr0.Nz == 1) {

            // No-block formatting at all
            // TODO make this into a loop

            //-------------------------------------------------- 
            // electrons
            qmE = gr0.getQ(0); // electron mass-to-charge

            vmesh::VeloMesh& v0e = gr0.electrons(0,0,0);
            vmesh::VeloMesh& v1e = gr1.electrons(0,0,0);

            vmesh::VeloMesh& vp1e = gr_p1.electrons(0, 0, 0); // xxx | 0, 1, 2, 
            vmesh::VeloMesh& vm1e = gr_m1.electrons(0, 0, 0); // 0, 1, 2, .. | xxx

            yee.jx(0,0,0) -= qmE*solve1d(v0e, vm1e, vp1e, v1e, dim, cellPtr);


            //-------------------------------------------------- 
            // positrons
            qmP = gr0.getQ(1); // positron mass-to-charge
            vmesh::VeloMesh& v0p = gr0.positrons(0,0,0);
            vmesh::VeloMesh& v1p = gr1.positrons(0,0,0);

            vmesh::VeloMesh& vp1p = gr_p1.positrons(0, 0, 0); // xxx | 0, 1, 2, 
            vmesh::VeloMesh& vm1p = gr_m1.positrons(0, 0, 0); // 0, 1, 2, .. | xxx

            yee.jx(0,0,0) -= qmP*solve1d(v0p, vm1p, vp1p, v1p, dim, cellPtr);

          } else {

            size_t first = 0;
            size_t last  = gr0.Nx-1;

            // Block structured cells
            for(size_t k = 0; k<gr0.Nz; k++) {
              for(size_t j = 0; j<gr0.Ny; j++) {

                // TODO make this into a loop
                //--------------------------------------------------
                // electrons
                qmE = gr0.getQ(0); // electron mass-to-charge

                // leftmost side blocks (-1 value from left neighbor)
                yee.jx(0,j,k) -= qmE*solve1d(
                    gr0.  electrons(first,   j,k),
                    gr_m1.electrons(last,    j,k),
                    gr0.  electrons(first+1, j,k),
                    gr1.  electrons(first,   j,k),
                    dim, cellPtr);

                // inner blocks
                for(size_t i=1; i<gr0.Nx-1; i++) {
                  
                  yee.jx(i,j,k) -= qmE*solve1d(
                      gr0.electrons(i,   j,k),
                      gr0.electrons(i-1, j,k),
                      gr0.electrons(i+1, j,k),
                      gr1.electrons(i,   j,k),
                      dim, cellPtr);

                }

                // rightmost side blocks (+1 value from right neighbor)
                yee.jx(last,j,k) -= qmE*solve1d(
                    gr0.  electrons(last,   j,k),
                    gr0.  electrons(last-1, j,k),
                    gr_p1.electrons(first,  j,k),
                    gr1.  electrons(last,   j,k),
                    dim, cellPtr);



                //--------------------------------------------------
                // positrons
                qmP = gr0.getQ(1); // positron mass-to-charge
                  
                // leftmost side blocks (-1 value from left neighbor)
                yee.jx(0,j,k) -= qmP*solve1d(
                    gr0.  positrons(first,   j,k),
                    gr_m1.positrons(last,    j,k),
                    gr0.  positrons(first+1, j,k),
                    gr1.  positrons(first,   j,k),
                    dim, cellPtr);

                // inner blocks
                for(size_t i=1; i<gr0.Nx-1; i++) {
                  
                  yee.jx(i,j,k) -= qmP*solve1d(
                      gr0.positrons(i,   j,k),
                      gr0.positrons(i-1, j,k),
                      gr0.positrons(i+1, j,k),
                      gr1.positrons(i,   j,k),
                      dim, cellPtr);

                }

                // rightmost side blocks (+1 value from right neighbor)
                yee.jx(last,j,k) -= qmP*solve1d(
                    gr0.  positrons(last,   j,k),
                    gr0.  positrons(last-1, j,k),
                    gr_p1.positrons(first,  j,k),
                    gr1.  positrons(last,   j,k),
                    dim, cellPtr);

              }
            }
          
          }


        } // end of dimension cycle

        // normalize current with momentum cell size
        for(size_t k = 0; k<gr0.Nz; k++) {
          for(size_t j = 0; j<gr0.Ny; j++) {
            for(size_t i = 0; i<gr0.Nx; i++) {

              // TODO now this assumes that electrons == positrons in grid size
              auto lens = gr0.electrons(i,j,k).lens;
              double du = lens[0];

              yee.jx(i,j,k) *= du;
            }
          }
        }
        

        return;
      }


      /*! \brief Solve 1D advection problem 
       *
       * Internally we use 2D Lagrangian interpolation.
       */

      Realf solve1d(vmesh::VeloMesh& v0,
                    vmesh::VeloMesh& vm1,
                    vmesh::VeloMesh& vp1,
                    vmesh::VeloMesh& v1,
                    size_t dim,
                    Grid::CellPtr cellPtr) {

        Realf numerical_flux = 0.0;

        // initialize cell and step size
        Realf dt = cellPtr->dt;
        Realf dx = cellPtr->dx;


        size_t Nb = v0.Nblocks[dim]; // number of slices to loop over

        // loop over every sheet in the mesh
        sheets::Sheet Up, Um, flux;
        sheets::Sheet s0, sp1, sm1;
        Realf aa;
        for(size_t i=0; i<Nb; i++) {
          sm1 = vm1.getSheet(dim, i);
          s0  = v0 .getSheet(dim, i);
          sp1 = vp1.getSheet(dim, i);

          Realf u0 = s0.sliceVal;

          // take 1D relativity into account
          // TODO in 2/3D this needs to deal with varying sheet guide grids too
          u0 = u0 / sqrt(1.0 + u0*u0);

          aa = (Realf) u0 * (dt/dx);

          // U_+1/2
          Up = (sp1 + s0)*0.5* aa   
              -(sp1 - s0)*0.5* aa*aa;

          // U_-1/2
          Um = (s0 + sm1)*0.5* aa
              -(s0 - sm1)*0.5* aa*aa;

          // dF = U_-1/2 - U_+1/2
          flux = s0 + (Um - Up);

          v1.addSheet(dim, i,  flux);

          // collect numerical volume flux that is leaking from cube to another
          // Size of elementary cell is dvx*dvy*dvz velocity cube.
          // For relativitys gamma we need sliceval vx + changing vy & vz 
          // from guide grids
          numerical_flux += u0 * volWeightedSum(flux);

        }

        return numerical_flux;
      }



    };




} // end of namespace vlasov

