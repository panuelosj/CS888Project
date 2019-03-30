#include <Eigen/Core>

using namespace Eigen;

class MACgridPressure
{
public:
  MACgridPressure();
  MACgridPressure(Vector2i gridSize, Vector2f gridSpacing);
  ~MACgridPressure();

  // definition of grid
  Vector2i _gridSize = Vector2i(0, 0);
  Vector2f _gridSpacing = Vector2f(0.0, 0.0);
  // data
  MatrixXf _P;

private:
  void _init();
};
