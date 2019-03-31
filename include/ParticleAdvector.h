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
  void _updateMaterialField();

  // INPUT VARIABLES
  double            _dt = 1.0;
  // definition of the grid
  Vector2i          _gridSize        = Vector2i(0, 0);
  Vector2d          _gridSpacing     = Vector2d(0.0, 0.0);
  // datastructures for our actual current state
  MACgridVelocity   *_velocityField;    // velocities
  MACgridVelocity   *_velocityFieldOld; // velocities at the start of the timestep
  MaterialGrid      *_materialField;    // voxel grid identifier
  Particles         *_particles;        // particle data
  GridToParticle    *_gridToParticle;   // grid to particle interpolation operator
};
