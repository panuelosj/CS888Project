#include <Eigen/Core>
#include "Config.h"
#include "MACgridVelocity.h"
#include "MaterialGrid.h"
#include "Particles.h"

using namespace Eigen;

struct ParticleToGridInputs {
  Vector2i          gridSize;
  Vector2d          gridSpacing;

  MACgridVelocity   *velocityField;
  MACgridVelocity   *velocityFieldOld;
  MaterialGrid      *materialField;
  Particles         *particles;
};

class ParticleToGrid
{
public:
  ParticleToGrid();
  ParticleToGrid(ParticleToGridInputs in);
  ~ParticleToGrid();

  void transfer();

private:
  // operators
  void _accumulate();
  void _normalize();
  void _saveToOldField();
  // kernels
  double _kernel(Vector2d r);
  double _kernel_quadraticBSpline(Vector2d r);
  double _kernel_quadraticBSpline1D(double r);

  // INPUT VARIABLES
  // definition of the grid
  Vector2i _gridSize        = Vector2i(0, 0);
  Vector2d _gridSpacing     = Vector2d(0.0, 0.0);
  // datastructures for our actual current state
  MACgridVelocity   *_velocityField;    // velocities
  MACgridVelocity   *_velocityFieldOld; // velocities at the start of the timestep
  MACgridVelocity   *_weights;          // save the sum of weights here
  MaterialGrid      *_materialField;    // voxel grid identifier
  Particles         *_particles;        // particle data
};
