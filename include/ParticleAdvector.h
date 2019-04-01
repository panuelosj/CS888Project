#pragma once

#include <Eigen/Core>
#include <limits>
#include "Config.h"
#include "MACgridVelocity.h"
#include "MaterialGrid.h"
#include "Particles.h"
#include "GridToParticle.h"

using namespace Eigen;

struct ParticleAdvectorInputs {
  double            dt;
  Vector2i          gridSize;
  Vector2d          gridSpacing;

  MACgridVelocity   *velocityField;
  MACgridVelocity   *velocityFieldOld;
  MaterialGrid      *materialField;
  Particles         *particles;
  GridToParticle    *gridToParticle;
};

class ParticleAdvector
{
public:
  ParticleAdvector();
  ParticleAdvector(ParticleAdvectorInputs in);
  ~ParticleAdvector();

  void advect();

private:
  void _advectRalstonRK3();
  void _detectAndResolveSolidCollisions();
  void _updateMaterialField();

  // INPUT VARIABLES
  // timestepping
  double            _dt               = 1.0;
  double            _substepTime      = 0.0;
  double            _dtSubstep        = 1.0;
  // definition of the grid
  Vector2i          _gridSize         = Vector2i(0, 0);
  Vector2d          _gridSpacing      = Vector2d(0.0, 0.0);
  unsigned int      _nParticles       = 0;
  // datastructures for our actual current state
  MACgridVelocity   *_velocityField;    // velocities
  MACgridVelocity   *_velocityFieldOld; // velocities at the start of the timestep
  MaterialGrid      *_materialField;    // voxel grid identifier
  Particles         *_particles;        // particle data
  GridToParticle    *_gridToParticle;   // grid to particle interpolation operator
  // temporary state datastructures
  MatrixXd          _pOld;              // positions at the start of a substep
  MatrixXd          _pNew;              // positions at the end of a substep
  MatrixXd          _vOld;              // velocities at the start of a substep
  MatrixXd          _vNew;              // velocities at the end of a substep
};
