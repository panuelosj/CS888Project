#include "GridToParticle.h"
#include <iostream>

//  ######   #######  ##    ##  ######  ######## ########
// ##    ## ##     ## ###   ## ##    ##    ##    ##     ##
// ##       ##     ## ####  ## ##          ##    ##     ##
// ##       ##     ## ## ## ##  ######     ##    ########
// ##       ##     ## ##  ####       ##    ##    ##   ##
// ##    ## ##     ## ##   ### ##    ##    ##    ##    ##
//  ######   #######  ##    ##  ######     ##    ##     ##

GridToParticle::GridToParticle() {

}

GridToParticle::GridToParticle(GridToParticleInputs in) :
  _gridSize           ( in.gridSize ),
  _gridSpacing        ( in.gridSpacing ),
  _velocityField      ( in.velocityField ),
  _velocityFieldOld   ( in.velocityFieldOld ),
  _materialField      ( in.materialField ),
  _particles          ( in.particles )
{

}

GridToParticle::~GridToParticle() {

}

// ########  ##     ## ########  ##       ####  ######
// ##     ## ##     ## ##     ## ##        ##  ##    ##
// ##     ## ##     ## ##     ## ##        ##  ##
// ########  ##     ## ########  ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##    ##
// ##         #######  ########  ######## ####  ######

// TODO: MAKE SURE TO UPDATE THE MATERIAL MATRIX YOU DORK!

void GridToParticle::transfer() {
  //

  // compute the change in velocity field (FLIP)
  _velocityFieldDelta = new MACgridVelocity( _gridSize, _gridSpacing );
  _velocityFieldDelta->copyInData(_velocityField);
  _velocityFieldDelta->subAllData(_velocityFieldOld);

  // do interpolation for every particle
  for (unsigned int idx=0; idx<_particles->nParticles(); idx++) {
    // need the particle's referencespace position (within its own cell)
    Vector2d particleReferencespacePosition = _particles->particleToReferencespace(idx);
    // need the particle's gridspace position
    Vector2i particleCell = _particles->particleToGridIndex(idx);

    // now we have a choice of interpolation methods
    // we'll make this a compile-time choice

    // interpolate velocity (PIC)
    Vector2d interpolatedPIC = _interpolate(_velocityField,
                                            particleCell,
                                            particleReferencespacePosition);

    // interpolate velocity change (FLIP)
    Vector2d interpolatedFLIP = _interpolate( _velocityFieldDelta,
                                              particleCell,
                                              particleReferencespacePosition);

    // now we can combine the PIC/FLIP contributions to get the new velocity
    // first grab old velocity
    Vector2d vOld = _particles->particleVelocity(idx);
    Vector2d vNew = PIC_FLIP_ALPHA*interpolatedPIC
                    + (1.0-PIC_FLIP_ALPHA)*(vOld + interpolatedFLIP);

    _particles->setParticleVelocity(idx, vNew);
  }

  // make sure to cleanup
  delete _velocityFieldDelta;
}


// ########  ########  #### ##     ##    ###    ######## ########
// ##     ## ##     ##  ##  ##     ##   ## ##      ##    ##
// ##     ## ##     ##  ##  ##     ##  ##   ##     ##    ##
// ########  ########   ##  ##     ## ##     ##    ##    ######
// ##        ##   ##    ##   ##   ##  #########    ##    ##
// ##        ##    ##   ##    ## ##   ##     ##    ##    ##
// ##        ##     ## ####    ###    ##     ##    ##    ########

Vector2d GridToParticle::_interpolate(MACgridVelocity *data, Vector2i index, Vector2d position) {
  Vector2d retval = Vector2d(0.0, 0.0);
#ifdef INTERPOLATE_BILINEAR
  retval = _interpolate_bilinear(data, index, position);
#endif
  return retval;
}

Vector2d GridToParticle::_interpolate_bilinear(MACgridVelocity *d, Vector2i index, Vector2d position) {
  // for caching the indices
  int i,j;
  double x,y;

  // INTERPOLATE U
  // check to see if the interpolating position is in the top or bottom half of
  // the cell and correct the position/index accordingly
  if (position.y() > 0.5) {
    // particle is in top half, we need j, j+1 as control points
      // save the indices as is
      i = index.x();
      j = index.y();
    // correct the position since the velocities are staggered
      x = position.x();
      y = position.y() - 0.5;
  } else {
    // particle is in the bottom half, we need j-1, j as control points
      // sub 1 from the y-index to get the correct numbers
      i = index.x();
      j = index.y()-1;
    // correct the position since the velocities are staggered
      x = position.x();
      y = position.y() + 0.5;
  }
  // now interpolate
  double interpolatedU = _interpolate_bilinear_scalar(  d->U(i,j),   d->U(i+1,j),
                                                        d->U(i,j+1), d->U(i+1,j+1),
                                                        x, y );

  // interpolate V
  // check to see if the interpolating position is in the left or right half of
  // the cell and correct the position/index accordingly
  if (position.x() > 0.5) {
    // particle is in the right half, we need i, i+1 as control points
      // save the indices as is
      i = index.x();
      j = index.y();
    // correct the position since the velocities are staggered
      x = position.x() - 0.5;
      y = position.y();
  } else {
    // particle is in the left half, we need i-1, i as control points
      // sub 1 from the x-index to get the correct numbers
      i = index.x()-1;
      j = index.y();
    // correct the position since the velocities are staggered
      x = position.x() + 0.5;
      y = position.y();
  }
  // now interpolate
  double interpolatedV = _interpolate_bilinear_scalar(  d->V(i,j),   d->V(i+1,j),
                                                        d->V(i,j+1), d->V(i+1,j+1),
                                                        x, y );

  return Vector2d(interpolatedU, interpolatedV);
}

double GridToParticle::_interpolate_bilinear_scalar(double c00, double c10,
                                                    double c01, double c11,
                                                    double x, double y) {
  // data is arranged as:
  // c01 -- c11
  //  |      |
  // c00 -- c10
  double x1 = c00*(1.0-x) + c10*x;
  double x2 = c01*(1.0-x) + c11*x;
  return x1*(1.0-y) + x2*y;
}
