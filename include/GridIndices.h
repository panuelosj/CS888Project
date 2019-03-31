#include <Eigen/Core>
#include "MaterialGrid.h"

#define GRID_INDEX_NOT_FOUND_VAL -1

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
  int find(unsigned int i, unsigned int j);

private:
  VectorXi _xs;
  VectorXi _ys;
  unsigned int _nIndices = 0;
};
