#include <Eigen/Core>
#include "MaterialGrid.h"

using namespace Eigen;

// This class is responsible for keeping a list of grid indices, so as to cache
// certain cells (ex, keep a list of all valid fluid cells)
// It effectively maps a 1D index into a 2D worldspace grid index

class GridIndices
{
public:
  GridIndices();
  GridIndices(MaterialGrid *materialField, Material materialType);
  ~GridIndices();

  Vector2i index(unsigned int idx);
  unsigned int size();

private:
  VectorXi _xs;
  VectorXi _ys;
  unsigned int _nIndices = 0;
};
