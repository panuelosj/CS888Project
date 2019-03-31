#pragma once

#include <Eigen/Core>

using namespace Eigen;

enum class Material : int {
  empty,
  fluid,
  solid
};

class MaterialGrid
{
public:
  MaterialGrid();
  MaterialGrid(Vector2i gridSize, Vector2d gridSpacing);
  ~MaterialGrid();

  // DATA READS
  // returns the size of the grid
  Vector2i gridSize();
  Vector2d gridSpacing();
  // returns the material type at grid index (i,j)
  Material material(unsigned int i, unsigned int j);
  // returns entire material field matrix
  MatrixXi material();
  // returns true if grid index (i,j) has the corresponding material type
  bool isSolid(unsigned int i, unsigned int j);
  bool isFluid(unsigned int i, unsigned int j);
  bool isEmpty(unsigned int i, unsigned int j);
  bool isSolid(Vector2i v);
  bool isFluid(Vector2i v);
  bool isEmpty(Vector2i v);
  // returns the number of cells of each type
  unsigned int nSolid();
  unsigned int nFluid();
  unsigned int nEmpty();

  // DATA SETS
  void clearFluid();
  void setMaterial(unsigned int i, unsigned int j, Material newMaterial);
  void setBlockToMaterial(unsigned int i, unsigned int j, unsigned int p, unsigned int q, Material newMaterial);
  void updateCellCount();

private:
  // updates the count of _nSolid, _nFluid, _nEmpty
  void _updateCellCount();

  // returns true if a cell is in the valid grid range
  inline bool _isInRange(unsigned int i, unsigned int j);

  // definition of grid
  Vector2i _gridSize = Vector2i(0, 0);
  Vector2d _gridSpacing = Vector2d(0.0, 0.0);
  // data
  MatrixXi _material;
  unsigned int _nSolid = 0;
  unsigned int _nFluid = 0;
  unsigned int _nEmpty = 0;
};
