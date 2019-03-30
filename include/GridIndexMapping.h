#include <Eigen/Core>

#define GRID_TO_VECTOR_OUTSIDE_RANGE_VAL -1

using namespace Eigen;

// This class is responsible for translating a grid index into the vector index
// It effectively maps a 2D worldspace grid index to a 1D vector index

class GridIndexMapping
{
public:
  GridIndexMapping();
  GridIndexMapping(Vector2i gridSize);
  ~GridIndexMapping();

  int gridToVector(int i, int j);
  int gridToVector(Vector2i g);

private:
  inline unsigned int _flattenIndex(int i, int j);
  inline unsigned int _flattenIndex(Vector2i g);
  inline bool _isInRange(unsigned int i, unsigned int j);

  // INPUT VARIABLES
  Vector2i _gridSize = Vector2i(0, 0);
  // DERIVED VARIABLES
  int _gridVectorLength = 0;
  // CONSTRUCTED VARIABLES
  VectorXi _indices;
};
