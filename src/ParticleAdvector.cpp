#include "ParticleAdvector.h"

#include <iostream>
#include <float.h>

//  ######   #######  ##    ##  ######  ######## ########
// ##    ## ##     ## ###   ## ##    ##    ##    ##     ##
// ##       ##     ## ####  ## ##          ##    ##     ##
// ##       ##     ## ## ## ##  ######     ##    ########
// ##       ##     ## ##  ####       ##    ##    ##   ##
// ##    ## ##     ## ##   ### ##    ##    ##    ##    ##
//  ######   #######  ##    ##  ######     ##    ##     ##

ParticleAdvector::ParticleAdvector() {

}

ParticleAdvector::ParticleAdvector(ParticleAdvectorInputs in) :
  _dt                 ( in.dt ),
  _gridSize           ( in.gridSize ),
  _gridSpacing        ( in.gridSpacing ),
  _velocityField      ( in.velocityField ),
  _velocityFieldOld   ( in.velocityFieldOld ),
  _materialField      ( in.materialField ),
  _particles          ( in.particles ),
  _gridToParticle     ( in.gridToParticle ),
  _nParticles         ( in.particles->nParticles() )
{
  _pOld = MatrixXd::Zero(_nParticles, 2);
  _pNew = MatrixXd::Zero(_nParticles, 2);
  _vOld = MatrixXd::Zero(_nParticles, 2);
  _vNew = MatrixXd::Zero(_nParticles, 2);
}

ParticleAdvector::~ParticleAdvector() {

}

// ########  ##     ## ########  ##       ####  ######
// ##     ## ##     ## ##     ## ##        ##  ##    ##
// ##     ## ##     ## ##     ## ##        ##  ##
// ########  ##     ## ########  ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##    ##
// ##         #######  ########  ######## ####  ######

void ParticleAdvector::advect() {
  // master function that calls variations of the particle advection scheme

  // setup our substeps
  _substepTime = 0.0;
  bool finishedFullstep = false;

  unsigned int nIterations = 0;
  double dtSubstepMin = DBL_MAX;

  while (!finishedFullstep) {
    _pOld = _particles->positions();
    _vOld = _particles->velocities();

    // find the max substep size
    MatrixXd tempVelocity = MatrixXd(_nParticles, 2);
    _gridToParticle->interpolateVelocities(&_pOld, &_vOld, &tempVelocity);
    // get max speed
      // NOTE THAT WE CAN'T JUST CALL _particle->maxSpeed() HERE SINCE WE INTERPOLATED JUST NOW
      // AND WE DON'T WANT TO EDIT THE PARTICLE VELOCITIES UNTIL THE ACTUAL ADVECT STEP
    double maxSpeed = sqrt(((tempVelocity.col(0)).cwiseAbs2() + (tempVelocity.col(1)).cwiseAbs2()).maxCoeff());
    // set substep dt to the longest allowable by CFL condition
    _dtSubstep = sqrt(_gridSpacing(0)*_gridSpacing(1))/(maxSpeed+D_EPSILON);

    if (_substepTime + _dtSubstep >= _dt) {
      _dtSubstep = _dt - _substepTime;
      finishedFullstep = true;
    } else if (_substepTime + 2*_dtSubstep >= _dt) {
      _dtSubstep = 0.5*(_dt - _substepTime);
    }

    #ifdef ADVECT_RALSTONRK3
    _advectRalstonRK3();
    #endif

    _substepTime += _dtSubstep;
    nIterations++;
    if (_dtSubstep < dtSubstepMin)
      dtSubstepMin = _dtSubstep;
  }

#if defined LOG_PARTICLE_ADVECTION && !defined LOG_QUIET
  std::cout << "\t\tNumber of Substeps Taken: " << nIterations << std::endl;
  std::cout << "\t\tSmallest Substep Taken: " << dtSubstepMin << std::endl;
#endif

}


// ########  ########  #### ##     ##    ###    ######## ########
// ##     ## ##     ##  ##  ##     ##   ## ##      ##    ##
// ##     ## ##     ##  ##  ##     ##  ##   ##     ##    ##
// ########  ########   ##  ##     ## ##     ##    ##    ######
// ##        ##   ##    ##   ##   ##  #########    ##    ##
// ##        ##    ##   ##    ## ##   ##     ##    ##    ##
// ##        ##     ## ####    ###    ##     ##    ##    ########

void ParticleAdvector::_advectRalstonRK3() {
  // This is Ralson's Third-Order Runge-Kutta method
    // We use this as it is suggested by Bridson on p110

  // we would like to vectorize this for use with Eigen, so make space to save
    // the intermediate velocity states
  MatrixXd  k1 = MatrixXd(_nParticles, 2);
  MatrixXd  k2 = MatrixXd(_nParticles, 2);
  MatrixXd  k3 = MatrixXd(_nParticles, 2);


  // interpolate to get the k1 = u(xn)
    // k1 = u(xn)
  MatrixXd  p1 = _pOld;
  _gridToParticle->interpolateVelocities(&p1, &_vOld, &k1);

  // linearly advect to get xn + 0.5*dt*k1
  MatrixXd  p2 = p1 + (0.5*_dtSubstep)*k1;
  // interpolate to get the k2 = u(xn + 0.5*dt*k1)
  _gridToParticle->interpolateVelocities(&p2, &_vOld, &k2);

  // linearly advect to get xn + 0.75*dt*k2
  MatrixXd  p3 = p1 + (0.75*_dtSubstep)*k2;
  _gridToParticle->interpolateVelocities(&p3, &_vOld, &k3);

  // now we can put everything together
  _pNew = p1 + ((2.0/9.0)*_dtSubstep)*k1
             + ((3.0/9.0)*_dtSubstep)*k2
             + ((4.0/9.0)*_dtSubstep)*k3;

  _detectAndResolveSolidCollisions();

  // we have the new positions, we need to get the new velocities
    // we do one last interpolation
  MatrixXd  vNew = MatrixXd(_nParticles, 2);
  _gridToParticle->interpolateVelocities(&_pNew, &_vOld, &vNew);
  for (unsigned int idx=0; idx<_nParticles; idx++) {
    _particles->setParticleVelocity(idx, vNew.row(idx));
  }

  // now we need to update the material grid so it knows where the fluid is now
  _updateMaterialField();
}

void ParticleAdvector::_detectAndResolveSolidCollisions() {
  // now save the new positions into the particles
    // I guess we have to do this serially
  for (unsigned int idx=0; idx<_nParticles; idx++) {
    _particles->setParticlePosition(idx, _pNew.row(idx));

    // check to make sure that the particle is not stuck in a solid
    Vector2i gridIndex = _particles->particleToGridIndex(idx);
    // grab the vector pointing from the current position to the original
      // position of the particle
    Vector2d u = _pOld.row(idx) - _pNew.row(idx);

    int iterCount = 0;

    while (_materialField->isSolid(gridIndex.x(), gridIndex.y()) && iterCount < 3) {
      // here we know we are inside a solid

      // grab the positions of the bounding grid cell
      double x0 = _velocityField->gridIndexToWorldspaceU(gridIndex).x();
      double x1 = x0 + _gridSpacing(0);
      double y0 = _velocityField->gridIndexToWorldspaceV(gridIndex).y();
      double y1 = y0 + _gridSpacing(1);

      // now we need to find the scalar length that puts the particle outside
        // the solid
      double tmp  = 0.0;  // temporary variable for holding the computed length
      double k    = std::numeric_limits<double>::max();  // "best" scalar length so far

      if (abs(u.x()) > D_EPSILON) {
        // compute the distance to x=0
        tmp = (x0-_pNew.row(idx).x())/u.x();
        if (tmp > 0.0 && tmp < k) k = tmp;      // distance should be positive
        // compute the distance to x=1
        tmp = (x1-_pNew.row(idx).x())/u.x();
        if (tmp > 0.0 && tmp < k) k = tmp;      // distance should be positive
      }
      if (abs(u.y()) > D_EPSILON) {
        // compute the distance to y=0
        tmp = (y0-_pNew.row(idx).y())/u.y();
        if (tmp > 0.0 && tmp < k) k = tmp;      // distance should be positive
        // compute the distance to y=1
        tmp = (y1-_pNew.row(idx).y())/u.y();
        if (tmp > 0.0 && tmp < k) k = tmp;      // distance should be positive
      }

      assert(k > 0.0);
      // add a bit to the distance so the particle is outside the wall!
      k += D_EPSILON;

      // now push the particle outside
      _pNew.row(idx) += k*u;
      _particles->setParticlePosition(idx, _pNew.row(idx));

      // now update the gridIndex to see if the particle is still stuck in a solid
      gridIndex = _particles->particleToGridIndex(idx);

      // do a full PIC interpolation of this single particle so it doesn't keep
        // its velocity
      Vector2d newVel = _gridToParticle->interpolateOneVelocityPIC(_particles->particleToWorldspace(idx), _particles->particleVelocity(idx));
      _particles->setParticleVelocity(idx, newVel);

      iterCount++;
    }
  }
}

void ParticleAdvector::_updateMaterialField() {
  // this function sets as fluid the cells containing particles in the material field

  // first clear the fluid
  _materialField->clearFluid();

  // now find the index of every particle and set those cells to fluid
  for (unsigned int idx=0; idx<_nParticles; idx++) {
    Vector2i gridIndex = _particles->particleToGridIndex(idx);
    if (!_materialField->isSolid(gridIndex)) {
      _materialField->setMaterial(gridIndex.x(), gridIndex.y(), Material::fluid);
    }
  }

  // make sure we update the counters
  _materialField->updateCellCount();
}
