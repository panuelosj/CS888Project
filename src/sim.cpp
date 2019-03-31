#include "sim.h"
#include <iostream>

//  ######   #######  ##    ##  ######  ######## ########
// ##    ## ##     ## ###   ## ##    ##    ##    ##     ##
// ##       ##     ## ####  ## ##          ##    ##     ##
// ##       ##     ## ## ## ##  ######     ##    ########
// ##       ##     ## ##  ####       ##    ##    ##   ##
// ##    ## ##     ## ##   ### ##    ##    ##    ##    ##
//  ######   #######  ##    ##  ######     ##    ##     ##

Sim::Sim() {
}

Sim::~Sim() {
}

Sim::Sim(Vector2i gridSize, Vector2d gridSpacing) :
  _gridSize( gridSize ),
  _gridSpacing( gridSpacing )
{
  // construct an empty velocity field
  _velocityField      = new MACgridVelocity( _gridSize, _gridSpacing );
  _velocityFieldOld   = new MACgridVelocity( _gridSize, _gridSpacing );
  // construct an empty material field
  _materialField      = new MaterialGrid( _gridSize, _gridSpacing );
  // construct the grid to flat vector index mapping
  _mapping            = new GridIndexMapping( _gridSize );
}

// ########  ##     ## ########  ##       ####  ######
// ##     ## ##     ## ##     ## ##        ##  ##    ##
// ##     ## ##     ## ##     ## ##        ##  ##
// ########  ##     ## ########  ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##    ##
// ##         #######  ########  ######## ####  ######

void Sim::init() {
  // wrapper function to call the internal init function
  _init();
}

void Sim::update(double timestep) {
  // wrapper function to call the internal update function
  _dt = timestep;
  _update();
}

// Wrapper functions for setting up the simulation
void Sim::setBodyAcceleration(Vector2d b) {
  _velocityField->setBodyAcceleration(b);
}
void Sim::setBlockToMaterial(unsigned int i, unsigned int j, unsigned int p, unsigned int q, Material newMaterial) {
  _materialField->setBlockToMaterial(i,j,p,q,newMaterial);
}

// =============================================================================
// ------------------------- DATA READS ----------------------------------------
// =============================================================================
int Sim::nParticles() {
  return _particles->nParticles();
}
Vector2d Sim::particlePosition(unsigned int idx) {
  return _particles->particleToWorldspace(idx);
}

// ########  ########  #### ##     ##    ###    ######## ########
// ##     ## ##     ##  ##  ##     ##   ## ##      ##    ##
// ##     ## ##     ##  ##  ##     ##  ##   ##     ##    ##
// ########  ########   ##  ##     ## ##     ##    ##    ######
// ##        ##   ##    ##   ##   ##  #########    ##    ##
// ##        ##    ##   ##    ## ##   ##     ##    ##    ##
// ##        ##     ## ####    ###    ##     ##    ##    ########

void Sim::_init() {
  // initializes the particles
  // by this time the _materialField should be set with the initial fluid conditions
  _particles = new Particles(_materialField);

  _isInitialized = true;
}

// =============================================================================
// ------------------------- MASTER TIMESTEPPING FUNCTION ----------------------
// =============================================================================

void Sim::_update() {
  assert(_isInitialized);
  _particleToGridTransfer();

  _updateVelocityField();

    std::cout << "Before Projection U: " << std::endl << _velocityField->U().transpose().colwise().reverse() << std::endl;
    std::cout << "Before Projection V: " << std::endl << _velocityField->V().transpose().colwise().reverse() << std::endl;

  _updatePressureField();

    std::cout << "After Projection U: " << std::endl << _velocityField->U().transpose().colwise().reverse() << std::endl;
    std::cout << "After Projection V: " << std::endl << _velocityField->V().transpose().colwise().reverse() << std::endl;

  //_gridToParticleTransfer();
  //std::cout << "NewParticleVelocities: " << std::endl;
  //_particles->printParticleVelocities();

  _advectParticles();
}




void Sim::_particleToGridTransfer() {
  // PARTICLE TO GRID TRANSFER
  // this function interpolates the velocity field from the particles to the
  // grid

  // create particle to grid transfer class
  _setupParticleToGridInputs();
  _particleToGrid = new ParticleToGrid(_particleToGridInputs);

  // now do it
  _particleToGrid->transfer();

  // make sure to cleanup
  delete _particleToGrid;
}

void Sim::_updateVelocityField() {
  // VELOCITY UPDATE
  _velocityField->update(_dt);
}

void Sim::_updatePressureField() {
  // PRESSURE PROJECTION
  // create pressure projection
  _setupPressureProjectionInputs();
  _pressureProjection = new PressureProjection(_pressureProjectionInputs);

  // now do it
  _pressureProjection->project();

  // make sure to cleanup
  delete _pressureProjectionInputs.fluidCells;
  delete _pressureProjection;
}

void Sim::_gridToParticleTransfer() {
  // GRID TO PARTICLE TRANSFER
  // this function interpolates date from the grid to the particles
  // this function could vary dependent on the method we're using
  // PIC: interpolate velocities directly
  // FLIP: interpolate the change in velocities

  // note that we generally would not call this class since we often need to
  // compute the transfer operations multiple times during the advection step.
  // I'm just leaving this here as a template on how to use the GridToParticle
  // class

  // create grid to particle transfer class
  _setupGridToParticleInputs();
  _gridToParticle = new GridToParticle(_gridToParticleInputs);

  // now do it
  _gridToParticle->transfer();

  // make sure to cleanup
  delete _gridToParticle;
}

void Sim::_advectParticles() {
  // PARTICLE ADVECTOR
  // this function advects the particles according to the background velocity
  // field

  // first create the grid to particle transfer operator
  _setupGridToParticleInputs();
  _gridToParticle = new GridToParticle(_gridToParticleInputs);

  // create particle advector class
  _setupParticleAdvectorInputs();
  _particleAdvector = new ParticleAdvector(_particleAdvectorInputs);

  // now do it
  _particleAdvector->advect();

  // make sure to cleanup
  delete _particleAdvector;
  delete _gridToParticle;
}

void Sim::_setupPressureProjectionInputs() {
  _pressureProjectionInputs.dt                = _dt;
  _pressureProjectionInputs.density           = _defaultDensity;
  _pressureProjectionInputs.gridSize          = _gridSize;
  _pressureProjectionInputs.gridSpacing       = _gridSpacing;
  _pressureProjectionInputs.velocityField     = _velocityField;
  _pressureProjectionInputs.materialField     = _materialField;
  _pressureProjectionInputs.mapping           = _mapping;
  _pressureProjectionInputs.fluidCells        = new GridIndices(_materialField, Material::fluid);
}

void Sim::_setupParticleToGridInputs() {
  _particleToGridInputs.gridSize              = _gridSize;
  _particleToGridInputs.gridSpacing           = _gridSpacing;
  _particleToGridInputs.velocityField         = _velocityField;
  _particleToGridInputs.velocityFieldOld      = _velocityFieldOld;
  _particleToGridInputs.materialField         = _materialField;
  _particleToGridInputs.particles             = _particles;
}

void Sim::_setupGridToParticleInputs() {
  _gridToParticleInputs.gridSize              = _gridSize;
  _gridToParticleInputs.gridSpacing           = _gridSpacing;
  _gridToParticleInputs.velocityField         = _velocityField;
  _gridToParticleInputs.velocityFieldOld      = _velocityFieldOld;
  _gridToParticleInputs.materialField         = _materialField;
  _gridToParticleInputs.particles             = _particles;
}

void Sim::_setupParticleAdvectorInputs() {
  _particleAdvectorInputs.dt                  = _dt;
  _particleAdvectorInputs.gridSize            = _gridSize;
  _particleAdvectorInputs.gridSpacing         = _gridSpacing;
  _particleAdvectorInputs.velocityField       = _velocityField;
  _particleAdvectorInputs.velocityFieldOld    = _velocityFieldOld;
  _particleAdvectorInputs.materialField       = _materialField;
  _particleAdvectorInputs.particles           = _particles;
  _particleAdvectorInputs.gridToParticle      = _gridToParticle;
}
