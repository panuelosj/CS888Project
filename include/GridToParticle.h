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

  void transfer();

private:
  // INPUT VARIABLES
  // definition of the grid
  Vector2i _gridSize        = Vector2i(0, 0);
  Vector2d _gridSpacing     = Vector2d(0.0, 0.0);
  // datastructures for our actual current state
  MACgridVelocity   *_velocityField;    // velocities
  MACgridVelocity   *_velocityFieldOld; // velocities at the start of the timestep
  MaterialGrid      *_materialField;    // voxel grid identifier
  Particles         *_particles;        // particle data
};
