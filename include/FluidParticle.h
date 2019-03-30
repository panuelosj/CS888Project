#include <Eigen/Core>

using namespace Eigen;

class FluidParticle
{
public:
  FluidParticle();
  ~FluidParticle();

  void init();
  void update(float timestep);

  // properties
  Vector2d position;
  Vector2d velocity;
};
