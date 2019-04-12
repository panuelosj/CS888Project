#include "ParticleToGrid.h"

//  ######   #######  ##    ##  ######  ######## ########
// ##    ## ##     ## ###   ## ##    ##    ##    ##     ##
// ##       ##     ## ####  ## ##          ##    ##     ##
// ##       ##     ## ## ## ##  ######     ##    ########
// ##       ##     ## ##  ####       ##    ##    ##   ##
// ##    ## ##     ## ##   ### ##    ##    ##    ##    ##
//  ######   #######  ##    ##  ######     ##    ##     ##

ParticleToGrid::ParticleToGrid() {

}

ParticleToGrid::ParticleToGrid(ParticleToGridInputs in) :
  _gridSize           ( in.gridSize ),
  _gridSpacing        ( in.gridSpacing ),
  _velocityField      ( in.velocityField ),
  _velocityFieldOld   ( in.velocityFieldOld ),
  _materialField      ( in.materialField ),
  _particles          ( in.particles )
{

}

ParticleToGrid::~ParticleToGrid() {

}

// ########  ##     ## ########  ##       ####  ######
// ##     ## ##     ## ##     ## ##        ##  ##    ##
// ##     ## ##     ## ##     ## ##        ##  ##
// ########  ##     ## ########  ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##    ##
// ##         #######  ########  ######## ####  ######

void ParticleToGrid::transfer() {
  // first reset the velocity field
  _velocityField->clear();

  // we need a duplicate staggered grid to keep a running track of the sum of
  // all weights to divide by at the end
  _weights = new MACgridVelocity( _gridSize, _gridSpacing );

  // now accumulate particle data into the velocity field
  _accumulate();

  // now the particle data should be accumulated onto the _velocityField, but
  // need to be normalized by the sum of kernel weights
  _normalize();

  // now we need to extrapolate to fill the field with valid velocities
  // note that we need this to get a correct deltaV field for FLIP
  _velocityField->extrapolateVelocities();

  // now to enforce free-slip conditions, we need to set the correct boundary
  // velocities to zero:
  //_setBoundaryVelocities();

  // we need to copy the data into the _velocityFieldOld, which is used for
  // FLIP to know the initial state of the velocity field
  _saveToOldField();

  // make sure to cleanup
  delete _weights;
}


// ########  ########  #### ##     ##    ###    ######## ########
// ##     ## ##     ##  ##  ##     ##   ## ##      ##    ##
// ##     ## ##     ##  ##  ##     ##  ##   ##     ##    ##
// ########  ########   ##  ##     ## ##     ##    ##    ######
// ##        ##   ##    ##   ##   ##  #########    ##    ##
// ##        ##    ##   ##    ## ##   ##     ##    ##    ##
// ##        ##     ## ####    ###    ##     ##    ##    ########

void ParticleToGrid::_accumulate() {
  // now do a summation for every particle
  for (unsigned int idx=0; idx<_particles->nParticles(); idx++) {
    // need the particle's worldspace position
    Vector2d particlePosition = _particles->particleToWorldspace(idx);
    // need the particle's gridspace position
    Vector2i particleCell = _particles->particleToGridIndex(idx);

    // define a neighbourhood of a particle as the region within a cube of
    // length (2dx)x(2dy) centered around a particle. note that this is within
    // the index range (i-1,j-1)x(i+1,j+1) for both u and v fields.
    // we iterate through the corresponding indices
    for (int i=particleCell.x()-1; i<=particleCell.x()+1; i++) {
      for (int j=particleCell.y()-1; j<=particleCell.y()+1; j++) {
        // refer to Bridson p117, eq 7.2

        // ACCUMULATE U FIELD
        if (_velocityField->isValidUIndex(i,j)) {
          // need the face centers' worldspace position
          Vector2d uFacePosition = _velocityField->gridIndexToWorldspaceU(i, j);
          // measure the distance to the face centers
          Vector2d particleMinusU = particlePosition - uFacePosition;
          // accumulate the weighted contribution of the particle
          _velocityField->addU(i,j, _particles->particleVelocity(idx).x()*_kernel_quadraticBSpline(particleMinusU));
          // accumulate the weights so we can normalize at the end
          _weights->addU(i,j, _kernel_quadraticBSpline(particleMinusU));
        }

        // ACCUMULATE V FIELD
        if (_velocityField->isValidVIndex(i,j)) {
          // need the face centers' worldspace position
          Vector2d vFacePosition = _velocityField->gridIndexToWorldspaceV(i, j);
          // measure the distance to the face centers
          Vector2d particleMinusV = particlePosition - vFacePosition;
          // accumulate the weighted contribution of the particle
          _velocityField->addV(i,j, _particles->particleVelocity(idx).y()*_kernel_quadraticBSpline(particleMinusV));
          // accumulate the weights so we can normalize at the end
          _weights->addV(i,j, _kernel_quadraticBSpline(particleMinusV));
        }
      }
    }
  }
}

void ParticleToGrid::_normalize() {
  // Normalize U face velocities
  for (unsigned int i=0; i<_gridSize(0) + 1; i++) {
    for (unsigned int j=0; j<_gridSize(1); j++) {
      // make sure there was actually stuff saved as a weight
      if (_weights->U(i,j) >= D_EPSILON) {
        _velocityField->divideU(i,j, _weights->U(i,j));
        // set this to valid since we actually had stuff here
        _velocityField->setUValid(i,j, 1);
      } else {
        // there were no particles accumulated onto this face, so just set the
        // velocity to zero to be safe
        _velocityField->setU(i,j, 0.0);
        // and set it to invalid
        _velocityField->setUValid(i,j, 0);
      }
    }
  }

  // Normalize V face velocities
  for (unsigned int i=0; i<_gridSize(0); i++) {
    for (unsigned int j=0; j<_gridSize(1) + 1; j++) {
      // make sure there was actually stuff saved as a weight
      if (_weights->V(i,j) >= D_EPSILON) {
        _velocityField->divideV(i,j, _weights->V(i,j));
        // set this to valid since we actually had stuff here
        _velocityField->setVValid(i,j, 1);
      } else {
        // there were no particles accumulated onto this face, so just set the
        // velocity to zero to be safe
        _velocityField->setV(i,j, 0.0);
        _velocityField->setVValid(i,j, 0);
      }
    }
  }
}

void ParticleToGrid::_setBoundaryVelocities() {
  _velocityField->setBoundaryVelocities();
}

void ParticleToGrid::_saveToOldField() {
  // this function copies the velocity field data from _velocityField to
  // _velocityFieldOld, which is not changed by the velocity update and
  // pressure projection. We need to keep a copy of the velocity field at the
  // start of the timestep when using FLIP

  _velocityFieldOld->copyInData(_velocityField);
}

double ParticleToGrid::_kernel(Vector2d r) {
#ifdef KERNEL_QUADRATIC_B_SPLINE
  return _kernel_quadraticBSpline(r);
#else
  return 0.0;
#endif
}

double ParticleToGrid::_kernel_quadraticBSpline(Vector2d r) {
  // refer to Bridson p113
  double retval = _kernel_quadraticBSpline1D(r.x()/_gridSpacing.x());
  retval *= _kernel_quadraticBSpline1D(r.y()/_gridSpacing.y());

  return retval;
}

double ParticleToGrid::_kernel_quadraticBSpline1D(double r) {
  double retval = 0.0;
  if (r >= -1.5 && r < -0.5) {
    retval = 0.5*(r+1.5)*(r+1.5);
  } else if (r >= -0.5 && r < 0.5) {
    retval = 0.75 - (r*r);
  } else if (r >= 0.5 && r < 1.5) {
    retval = 0.5*(1.5-r)*(1.5-r);
  }
  return retval;
}
