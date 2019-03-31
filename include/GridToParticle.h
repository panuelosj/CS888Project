#pragma once

#include <Eigen/Core>
#include "Config.h"
#include "MACgridVelocity.h"
#include "MaterialGrid.h"
#include "Particles.h"

using namespace Eigen;

struct GridToParticleInputs {
  Vector2i          gridSize;
  Vector2d          gridSpacing;

  MACgridVelocity   *velocityField;
  MACgridVelocity   *velocityFieldOld;
  MaterialGrid      *materialField;
  Particles         *particles;
};

class GridToParticle
{
public:
  GridToParticle();
  GridToParticle(GridToParticleInputs in);
  ~GridToParticle();

  // GENERAL POSITION INTERPOLATION OPERATORS
  Vector2d interpolateOneVelocity(Vector2d p, Vector2d v);
  void interpolateVelocities(MatrixXd *positions, MatrixXd *velocities, MatrixXd *retval);
  // PARTICLE INDEX INTERPOLATION OPERATORS
  Vector2d interpolateOneVelocity(unsigned int idx);
  void interpolateVelocities(MatrixXd *retval);
  void transfer();

private:
  Vector2d _interpolate(MACgridVelocity *data, Vector2i index, Vector2d position);
  Vector2d _interpolate_bilinear(MACgridVelocity *d, Vector2i index, Vector2d position);
  double _interpolate_bilinear_scalar(double c00, double c10,
                                      double c01, double c11,
                                      double x, double y);

  // INPUT VARIABLES
  // definition of the grid
  Vector2i _gridSize        = Vector2i(0, 0);
  Vector2d _gridSpacing     = Vector2d(0.0, 0.0);
  // datastructures for our actual current state
  MACgridVelocity   *_velocityField;      // velocities
  MACgridVelocity   *_velocityFieldOld;   // velocities at the start of the timestep
  MACgridVelocity   *_velocityFieldDelta; // vNew - vOld
  MaterialGrid      *_materialField;      // voxel grid identifier
  Particles         *_particles;          // particle data
};
