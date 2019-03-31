#include "ParticleAdvector.h"

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
  _gridToParticle     ( in.gridToParticle )
{

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

#ifdef ADVECT_RALSTONRK3
  _advectRalstonRK3();
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
  int       nParticles = _particles->nParticles();
  MatrixXd  k1 = MatrixXd(nParticles, 2);
  MatrixXd  k2 = MatrixXd(nParticles, 2);
  MatrixXd  k3 = MatrixXd(nParticles, 2);

  // now we should keep a working copy of the positions and initial velocities
  MatrixXd  p1 = _particles->positions();
  MatrixXd  v = _particles->velocities();

  // interpolate to get the k1 = u(xn)
    // k1 = u(xn)
  _gridToParticle->interpolateVelocities(&p1, &v, &k1);

  // linearly advect to get xn + 0.5*dt*k1
  MatrixXd  p2 = p1 + (0.5*_dt)*k1;
  // interpolate to get the k2 = u(xn + 0.5*dt*k1)
  _gridToParticle->interpolateVelocities(&p2, &v, &k2);

  // linearly advect to get xn + 0.75*dt*k2
  MatrixXd  p3 = p1 + (0.75*_dt)*k2;
  _gridToParticle->interpolateVelocities(&p3, &v, &k3);

  // now we can put everything together
  MatrixXd  pNew = p1 + ((2.0/9.0)*_dt)*k1
                      + ((3.0/9.0)*_dt)*k2
                      + ((4.0/9.0)*_dt)*k3;

  // now save the new positions into the particles
    // I guess we have to do this serially
  for (unsigned int idx=0; idx<_particles->nParticles(); idx++) {
    _particles->setParticlePosition(idx, pNew.row(idx));

    // check to make sure that the particle is not stuck in a solid
    Vector2i gridIndex = _particles->particleToGridIndex(idx);
    // grab the vector pointing from the current position to the original
      // position of the particle
    Vector2d u = p1.row(idx) - pNew.row(idx);

    while (_materialField->isSolid(gridIndex.x(), gridIndex.y())) {
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
        tmp = (x0-pNew.row(idx).x())/u.x();
        if (tmp > 0.0 && tmp < k) k = tmp;      // distance should be positive
        // compute the distance to x=1
        tmp = (x1-pNew.row(idx).x())/u.x();
        if (tmp > 0.0 && tmp < k) k = tmp;      // distance should be positive
      }
      if (abs(u.y()) > D_EPSILON) {
        // compute the distance to y=0
        tmp = (y0-pNew.row(idx).y())/u.y();
        if (tmp > 0.0 && tmp < k) k = tmp;      // distance should be positive
        // compute the distance to y=1
        tmp = (y1-pNew.row(idx).y())/u.y();
        if (tmp > 0.0 && tmp < k) k = tmp;      // distance should be positive
      }

      assert(k > 0.0);
      // add a bit to the distance so the particle is outside the wall!
      k += D_EPSILON;

      // now push the particle outside
      pNew.row(idx) += k*u;
      _particles->setParticlePosition(idx, pNew.row(idx));

      // now update the gridIndex to see if the particle is still stuck in a solid
      gridIndex = _particles->particleToGridIndex(idx);
    }
  }

  // we have the new positions, we need to get the new velocities
    // we do one last interpolation
  MatrixXd  vNew = MatrixXd(nParticles, 2);
  _gridToParticle->interpolateVelocities(&pNew, &v, &vNew);
  for (unsigned int idx=0; idx<_particles->nParticles(); idx++) {
    _particles->setParticleVelocity(idx, vNew.row(idx));
  }

  // now we need to update the material grid so it knows where the fluid is now
  _updateMaterialField();
}

void ParticleAdvector::_updateMaterialField() {
  // this function sets as fluid the cells containing particles in the material field

  // first clear the fluid
  _materialField->clearFluid();

  // now find the index of every particle and set those cells to fluid
  for (unsigned int idx=0; idx<_particles->nParticles(); idx++) {
    Vector2i gridIndex = _particles->particleToGridIndex(idx);
    if (!_materialField->isSolid(gridIndex)) {
      _materialField->setMaterial(gridIndex.x(), gridIndex.y(), Material::fluid);
    }
  }

  // make sure we update the counters
  _materialField->updateCellCount();
}
