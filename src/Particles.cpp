#include "Particles.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

//  ######   #######  ##    ##  ######  ######## ########
// ##    ## ##     ## ###   ## ##    ##    ##    ##     ##
// ##       ##     ## ####  ## ##          ##    ##     ##
// ##       ##     ## ## ## ##  ######     ##    ########
// ##       ##     ## ##  ####       ##    ##    ##   ##
// ##    ## ##     ## ##   ### ##    ##    ##    ##    ##
//  ######   #######  ##    ##  ######     ##    ##     ##

Particles::Particles() {

}

Particles::Particles(MaterialGrid *materialField)
{
  // creates a list of particles based on the fluid cells of the provided
  // material field

  // set the number of particles we need
  _nParticles     = materialField->nFluid()*4;
  // save the material field
  _materialField  = materialField;
  // initialize the list of particles
  _positions      = MatrixXd::Zero(_nParticles, 2);
  _velocities     = MatrixXd::Zero(_nParticles, 2);
  // now seed it according to the material field
  _seedParticles();
  _jitterParticles();
}

Particles::~Particles() {

}


// ########  ##     ## ########  ##       ####  ######
// ##     ## ##     ## ##     ## ##        ##  ##    ##
// ##     ## ##     ## ##     ## ##        ##  ##
// ########  ##     ## ########  ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##    ##
// ##         #######  ########  ######## ####  ######

int Particles::nParticles() {
  return _nParticles;
}

Vector2d Particles::particleToWorldspace(unsigned int p) {
  // returns the worldspace position of a particle p
  return _positions.row(p);
}

Vector2d Particles::particleToReferencespace(unsigned int p) {
  // returns the position of particle p within its reference cell [0,1)x[0,1)
  double intpart;         // integer part of the output of modf, we just throw
                            // this away

  // grid spacings
  double dx     = _materialField->gridSpacing().x();
  double dy     = _materialField->gridSpacing().y();
  // reference position of the particles
  double xi     = modf(_positions(p,0)/dx, &intpart);
  double zeta   = modf(_positions(p,1)/dy, &intpart);
  return Vector2d(xi,zeta);
}

Vector2i Particles::particleToGridIndex(unsigned int p) {
  // returns the grid index of the cell that contains particle p

  // grid spacings
  double dx = _materialField->gridSpacing().x();
  double dy = _materialField->gridSpacing().y();
  // casting to int should give just the integer part
  int i = (int)(_positions(p,0)/dx);
  int j = (int)(_positions(p,1)/dy);
  return Vector2i(i,j);
}

Vector2d Particles::particleVelocity(unsigned int p) {
  return _velocities.row(p);
}

// ########  ########  #### ##     ##    ###    ######## ########
// ##     ## ##     ##  ##  ##     ##   ## ##      ##    ##
// ##     ## ##     ##  ##  ##     ##  ##   ##     ##    ##
// ########  ########   ##  ##     ## ##     ##    ##    ######
// ##        ##   ##    ##   ##   ##  #########    ##    ##
// ##        ##    ##   ##    ## ##   ##     ##    ##    ##
// ##        ##     ## ####    ###    ##     ##    ##    ########

void Particles::_seedParticles() {
  // seeds 2x2 jittered particles per grid cell based on the input material grid
  // note that this should only be called once on init

  unsigned int idx = 0;     // counter for the number of particles we've created
  // grid spacings
  double dx = _materialField->gridSpacing().x();
  double dy = _materialField->gridSpacing().y();

  // loop through every grid cell
  for (unsigned int i=0; i<_materialField->gridSize().x(); i++) {
    for (unsigned int j=0; j<_materialField->gridSize().y(); j++) {
      // check if the grid cell is a fluid
      if (_materialField->isFluid(i, j)) {
        // if it is, seed it with four particles
        // compute the bottom left corner of the cell
        double x0 = i*dx;
        double y0 = j*dy;

        _positions.row(idx  ) = RowVector2d(x0+0.25*dx, y0+0.25*dy);
        _positions.row(idx+1) = RowVector2d(x0+0.75*dx, y0+0.25*dy);
        _positions.row(idx+2) = RowVector2d(x0+0.25*dx, y0+0.75*dy);
        _positions.row(idx+3) = RowVector2d(x0+0.75*dx, y0+0.75*dy);

        idx += 4;
      }
    }
  }
}

void Particles::_jitterParticles() {
  // grid spacings
  double dx = _materialField->gridSpacing().x();
  double dy = _materialField->gridSpacing().y();

  for (unsigned int idx=0; idx<_nParticles; idx++) {
    _positions.row(idx) += RowVector2d(0.25*dx*_randomJitter(), 0.25*dy*_randomJitter());
  }
}

inline double Particles::_randomJitter() {
  return -_jitterFactor + ((double)rand() / (((double)RAND_MAX)/(2.0*_jitterFactor)));
}