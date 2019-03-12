#pragma once

#include <array>
#include <mpi4cpp/mpi.h>

#include "../definitions.h"
#include "../corgi/tile.h"
#include "../corgi/corgi.h"
#include "../em-fields/tile.h"
#include "particle.h"


namespace pic {

using namespace mpi4cpp;

/*! \brief PiC tile
 *
 * Tile infrastructures are inherited from corgi::Tile
 * Maxwell field equation solver is inherited from fields::Tile
*/

template<std::size_t D>
class Tile :
  virtual public fields::Tile<D>, 
  virtual public  corgi::Tile<D> 
{

public:

  /// Size of the internal grid
  //size_t NxGrid;
  //size_t NyGrid;
  //size_t NzGrid;

  using fields::Tile<D>::mesh_lengths;


  //ParticleContainer container;
  std::vector<ParticleContainer> containers;

  //--------------------------------------------------
  // normal container methods
     
  /// get i:th container
  ParticleContainer& get_container(size_t i) { return containers[i]; };

  const ParticleContainer& get_const_container(size_t i) const { return containers[i]; };

  /// set i:th container
  void set_container(const ParticleContainer& block) { containers.push_back(block); };

  size_t Nspecies() const { return containers.size(); };



  /// constructor
  Tile(size_t nx, size_t ny, size_t nz) :
     corgi::Tile<D>(),
    fields::Tile<D>(nx,ny,nz)
  { }


  /// destructor
  ~Tile() override = default;
  //~Tile() override {
  //  std::cout << "calling pic tile destructor...\n";
  //}

  /// tile temporal and spatial scales

  using fields::Tile<D>::cfl;
  using fields::Tile<D>::dx;

  //--------------------------------------------------
  // MPI send
  virtual void
  send_data( mpi::communicator&, int orig, int tag, std::vector<mpi::request>& reqs ) override;

  /// actual tag=0 send
  void send_particle_data( mpi::communicator&, int orig, std::vector<mpi::request>& reqs );

  /// actual tag=1 send
  void send_particle_extra_data( mpi::communicator&, int orig, std::vector<mpi::request>& reqs );


  //--------------------------------------------------
  // MPI recv
  virtual void 
  recv_data(mpi::communicator&, int dest, int tag, std::vector<mpi::request>& reqs ) override;

  /// actual tag=0 recv
  void recv_particle_data(mpi::communicator&, int dest, std::vector<mpi::request>& reqs );

  /// actual tag=1 recv
  void recv_particle_extra_data(mpi::communicator&, int dest, std::vector<mpi::request>& reqs );
  //--------------------------------------------------


  /// check all particle containers for particles
  // exceeding limits
  void check_outgoing_particles();

  /// delete particles from each container that are exceeding
  // the boundaries
  void delete_transferred_particles();

  /// get particles flowing into this tile
  void get_incoming_particles(corgi::Node<D>& grid);

  /// pack all particles for MPI message
  void pack_all_particles();

  /// pack particles for MPI message
  void pack_outgoing_particles();

  /// unpack received MPI message particles
  void unpack_incoming_particles();

  /// delete all particles from each container
  void delete_all_particles();

};



} // end of namespace pic
