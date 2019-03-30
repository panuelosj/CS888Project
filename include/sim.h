#include <Eigen/Core>
#include "Config.h"
#include "Particles.h"
#include "MACgridVelocity.h"
#include "MaterialGrid.h"
#include "ParticleToGrid.h"
#include "GridToParticle.h"
#include "PressureProjection.h"

using namespace Eigen;

class Sim
{
public:
  // default constructor
  Sim();

  // construct a sim object with grid
  Sim(Vector2i gridSize, Vector2d gridSpacing);

  // default destructor
  ~Sim();

  // once the simulation is setup from main, call this
  void init();
  // update, call this to do a single timestep
  void update(double timestep);
  // setup functions for the simulation
  void setBodyAcceleration(Vector2d b);
  void setBlockToMaterial(unsigned int i, unsigned int j, unsigned int p, unsigned int q, Material newMaterial);

  // definition of grid
  Vector2i            _gridSize = Vector2i(0, 0);
  Vector2d            _gridSpacing = Vector2d(0.0, 0.0);
  // simulation classes
  Particles           *_particles;
  MACgridVelocity     *_velocityField;
  MACgridVelocity     *_velocityFieldOld;     // we need this for FLIP
  MaterialGrid        *_materialField;
  GridIndexMapping    *_mapping;
  // operators
  ParticleToGrid      *_particleToGrid;
  PressureProjection  *_pressureProjection;
  GridToParticle      *_gridToParticle;

private:
  // inits
  void _init();             // master call to setup
  // timestepper
  void _update();
  void _particleToGridTransfer();
  void _updateVelocityField();
  void _updatePressureField();
  void _gridToParticleTransfer();
  // functions that setup the inputs to the operator classes
  void _setupPressureProjectionInputs();
  void _setupParticleToGridInputs();
  void _setupGridToParticleInputs();

  // constant
  double _defaultDensity = D_DENSITY;
  // variables
  bool _isInitialized = false;
  double _timestep;
  ParticleToGridInputs      _particleToGridInputs;
  GridToParticleInputs      _gridToParticleInputs;
  PressureProjectionInputs  _pressureProjectionInputs;
};