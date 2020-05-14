#include "verlet_loc.h"

#include <cmath>
#include <iostream>
#include "../../tools/signum.h"

using toolbox::sign;

template<size_t D, size_t V>
void pic::VerletLocPusher<D,V>::push_container(
    pic::ParticleContainer<D>& container,
    double cfl)
{
  int nparts = container.size();

  // initialize pointers to particle arrays
  real_prtcl* loc[3];
  for( int i=0; i<3; i++) loc[i] = &( container.loc(i,0) );

  real_prtcl* vel[3];
  for( int i=0; i<3; i++) vel[i] = &( container.vel(i,0) );

  real_long ex0 = 0.0, ey0 = 0.0, ez0 = 0.0;
  real_long bx0 = 0.0, by0 = 0.0, bz0 = 0.0;

  // make sure E and B tmp arrays are of correct size
  if(container.Epart.size() != (size_t)3*nparts)
    container.Epart.resize(3*nparts);
  if(container.Bpart.size() != (size_t)3*nparts)
    container.Bpart.resize(3*nparts);

  real_prtcl *ex, *ey, *ez, *bx, *by, *bz;
  ex = &( container.Epart[0*nparts] );
  ey = &( container.Epart[1*nparts] );
  ez = &( container.Epart[2*nparts] );

  bx = &( container.Bpart[0*nparts] );
  by = &( container.Bpart[1*nparts] );
  bz = &( container.Bpart[2*nparts] );

  // loop over particles
  int n1 = 0;
  int n2 = nparts;

  real_long c = cfl;
  real_long g;
  // charge (sign only)
  real_long qm = sign(container.q);

  // add division by m_s to simulate multiple species


  real_long vel0n, vel1n, vel2n;

  //TODO: SIMD
  for(int n=  n1; n<n2; n++) {
    //--------------------------------------------------
    // Velocity-Verlet Position Push

    vel0n = static_cast<real_long>( vel[0][n] );
    vel1n = static_cast<real_long>( vel[1][n] );
    vel2n = static_cast<real_long>( vel[2][n] );

    // read particle-specific fields
    ex0 = static_cast<real_long>( ex[n] );
    ey0 = static_cast<real_long>( ey[n] );
    ez0 = static_cast<real_long>( ez[n] );

    bx0 = static_cast<real_long>( bx[n] );
    by0 = static_cast<real_long>( by[n] );
    bz0 = static_cast<real_long>( bz[n] );

    // position advance;
    // NOTE: no mixed-precision calc here. Can be problematic.
    g = c / sqrt(c*c + c*c*vel0n*vel0n + c*c*vel1n*vel1n + c*c*vel2n*vel2n);

    for(size_t i=0; i<D; i++) loc[i][n] += vel[i][n]*g*c;
    loc[0][n] += c*c/2 * qm*(ex0 + g*(vel1n*bz0 - vel2n*by0));
    loc[1][n] += c*c/2 * qm*(ey0 + g*(vel2n*bx0 - vel0n*bz0));
    loc[2][n] += c*c/2 * qm*(ez0 + g*(vel0n*by0 - vel1n*bx0));
  }
}



//--------------------------------------------------
// explicit template instantiation

template class pic::VerletLocPusher<1,3>; // 1D3V
template class pic::VerletLocPusher<2,3>; // 2D3V
template class pic::VerletLocPusher<3,3>; // 3D3V
