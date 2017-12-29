#pragma once

#include <string>

#include "corgi/cell.h"
#include "velomesh.h"
#include "tools/Mesh.h"
#include "dataContainer.h"
#include "maxwell.h"


namespace vlasov {


/*! \brief Iterator to feed for-loops different plasma species properties
 * 
 * Holds internally pointer to an e.g., VlasovFluid and 
 * increments a counter that is used for comparison etc.
 */
template<typename M, typename T>
class PairPlasmaIterator {

  private:

    /// First species 
    size_t _beginning = 0;

    /// beyond last species 
    size_t _ending    = 2;

    /// internal pointer to the parent object/container (Mother object)
    M* ptr;

    /// Iterators own internal counting system
    size_t spcs = 0; 


  public:

    PairPlasmaIterator(M& rhs) : ptr(&rhs) {}
    PairPlasmaIterator(const M& rhs) : ptr(&rhs) {}

    PairPlasmaIterator(const PairPlasmaIterator& rhs) : ptr(rhs.ptr) {}


    /// Assignment
    PairPlasmaIterator& operator= (const PairPlasmaIterator& rhs) = default;

    /// iterate
    PairPlasmaIterator& operator++ () {
      ++this->spcs;
      return *this;
    }

    /// Referencing cell interiors
    T& operator *() {
      if(spcs == 0) return (T&) (ptr->electrons);
      else if(spcs == 1) return (T&) (ptr->positrons);
      else throw std::range_error("iterator goes beyond electrons (0) or positrons (1)");
    }


    /// iterate with steps
    // PairPlasmaIterator operator++ (int) {
    //   PairPlasmaIterator temp(*ptr);
    //   ++*this;
    //   return (temp);
    // }

    /// equal comparison done by comparing internal spcs value
    bool operator== (PairPlasmaIterator& rhs) const {
      return (ptr == rhs.ptr) && (spcs == rhs.spcs);
    }

    /// unequal comparison done by comparing internal spcs value
    bool operator!= (PairPlasmaIterator& rhs) const {
      return (ptr != rhs.ptr) || (spcs != rhs.spcs);
    }

    /// Returns an iterator pointing to the first element in the sequence
    PairPlasmaIterator begin() {
      PairPlasmaIterator temp(*ptr);
      temp.spcs = _beginning;
      return temp;
    }

    /// Returns an iterator pointing to the past-the-end element in the sequence
    PairPlasmaIterator end() {
      PairPlasmaIterator temp(*ptr);
      temp.spcs = _ending;
      return temp ;
    }


};




/*! \brief Vlasov fluid inside the cell
 *
 * Container to hold different plasma species.
 */
class VlasovFluid {

  private:
    typedef toolbox::Mesh<vmesh::VeloMesh, 0> T;

  public:

  size_t Nx;
  size_t Ny;
  size_t Nz;

  T electrons;
  T positrons;

  VlasovFluid(size_t Nx, size_t Ny, size_t Nz) : Nx(Nx), Ny(Ny), Nz(Nz),
    electrons(Nx, Ny, Nz),
    positrons(Nx, Ny, Nz) { }


  PairPlasmaIterator<VlasovFluid, T> species() {
    PairPlasmaIterator<VlasovFluid, T> ret(*this);
    return ret;
  };

};



/*! \brief Vlasov cell 
 *
 * Cell infrastructure methods are inherited from corgi::Cell
 * Maxwell field equation solver is inherited from maxwell::PlasmaCell
 */
class VlasovCell : 
                   virtual public maxwell::PlasmaCell, 
                   virtual public corgi::Cell {

  public:
    
    /// Size of the internal grid
    size_t NxGrid;
    size_t NyGrid;
    size_t NzGrid;


    /// constructor
    VlasovCell(size_t i, size_t j, 
               int o, 
               size_t NxG, size_t NyG,
               size_t NxMesh, size_t NyMesh
               ) : 
      corgi::Cell(i, j, o, NxG, NyG),
      maxwell::PlasmaCell(i,j,o,NxG, NyG, NxMesh,NyMesh),
      NxGrid(NxMesh),
      NyGrid(NyMesh),
      NzGrid(1)
      { 
        // add also initial velomeshes
        plasma.push_back( VlasovFluid(NxGrid, NyGrid, 1) );
        plasma.push_back( VlasovFluid(NxGrid, NyGrid, 1) );
      
      }

    /// destructor
    ~VlasovCell() { };

    // NOTE overwrites PlasmaCell values
    double dt = 0.01;
    double dx = 0.2;
    double dy = 1.0;
    double dz = 1.0;


    // Purely for testing class expansion
    // void bark();
    std::string bark() { return std::string("Wuff!"); };


    /// Simulation data container
    datarotators::DataContainer<VlasovFluid> plasma;


    /// Add data to the container
    /*
    void addData(vmesh::VeloMesh m) {
      vmeshes.push_back(m);
    }
    */

    // XXX defined only for python API
    // vmesh::VeloMesh getData() {
    //   return *data.get();
    // };

    VlasovFluid& getPlasmaGrid() {
      return plasma.getRef();
    };

    VlasovFluid& getNewPlasmaGrid() {
      return plasma.getNewRef();
    };


    // get pointer to the data
    // TODO: should be shared_ptr explicitly to make it memory save
    /*
    vmesh::VeloMesh* getDataPtr() {
      return vmeshes.get();
    };
    */


    /// Clip all the meshes inside cell
    void clip() {
      VlasovFluid& gr = plasma.getRef();

      for (size_t k=0; k<NzGrid; k++) {
        for (size_t j=0; j<NyGrid; j++) {
          for (size_t i=0; i<NxGrid; i++) {
            gr.electrons(i,j,k).clip();
            gr.positrons(i,j,k).clip();
          }
        }
      }

    }


    /// Cycle internal plasma container to another solution step
    void cycle() {
      plasma.cycle();
    }



};


} // end of namespace vlasov
