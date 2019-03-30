#include "GridToParticle.h"

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

  // do a summation for every particle
  for (unsigned int idx=0; idx<_particles->nParticles(); idx++) {
    // need the particle's referencespace position (within its own cell)
    Vector2d particleReferencespacePosition = _particles->particleToReferencespace(idx);
    // need the particle's gridspace position
    Vector2i particleCell = _particles->particleToGridIndex(idx);

    // now we have a choice of interpolation methods
    // we'll make this a compile-time choice

  }
}


// ########  ########  #### ##     ##    ###    ######## ########
// ##     ## ##     ##  ##  ##     ##   ## ##      ##    ##
// ##     ## ##     ##  ##  ##     ##  ##   ##     ##    ##
// ########  ########   ##  ##     ## ##     ##    ##    ######
// ##        ##   ##    ##   ##   ##  #########    ##    ##
// ##        ##    ##   ##    ## ##   ##     ##    ##    ##
// ##        ##     ## ####    ###    ##     ##    ##    ########
