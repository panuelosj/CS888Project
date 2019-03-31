#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include "Config.h"
#include "MACgridVelocity.h"
#include "MaterialGrid.h"
#include "GridIndexMapping.h"
#include "GridIndices.h"

using namespace Eigen;

struct PressureProjectionInputs {
  double            dt;
  double            density;
  Vector2i          gridSize;
  Vector2d          gridSpacing;

  MACgridVelocity   *velocityField;
  MaterialGrid      *materialField;
  GridIndexMapping  *mapping;
  GridIndices       *fluidCells;
};

// enum for listing which element of the A matrix is being used
enum class AMatrixCell : int {
  Adiag,
  Ax,
  Ay
};

class PressureProjection
{
public:
  PressureProjection();
  PressureProjection(PressureProjectionInputs in);
  ~PressureProjection();

  void project();

private:
  void _computeNegativeDivergence();
  void _computeMatrixCoefficients();
  void _translateMatrixToSparseEigen();
  void _solve();
  void _pressureGradientUpdate();

  // CONSTANTS
  double _usolid = SOLID_VELOCITY_X;
  double _vsolid = SOLID_VELOCITY_Y;
  double _unknownVelocity = UNKNOWN_VELOCITY;

  // INPUT VARIABLES
  double _dt;                           // timestep
  double _density;                      // I'm actually not sure what this is
  // definition of the grid
  Vector2i _gridSize        = Vector2i(0, 0);
  Vector2d _gridSpacing     = Vector2d(0.0, 0.0);
  // datastructures for our actual current state
  MACgridVelocity   *_velocityField;    // velocities
  MaterialGrid      *_materialField;    // voxel grid identifier
  GridIndexMapping  *_mapping;          // grid to vector index key mapping
  GridIndices       *_fluidCells;       // list of fluid cell indices
  unsigned int      _nFluidCells;       // number of fluid cells

  // DERIVED VARIABLES
  int _gridVectorLength     = 0;
  Vector2d _invGridSpacing  = Vector2d(0.0, 0.0);

  // CONSTRUCTED VARIABLES
  VectorXd _b;                          // right hand side vector
  VectorXd _p;                          // pressure field as a flat vector
  SparseMatrix<double> _A;              // Matrix of coefficients
  MatrixXi _Acompressed;                // Compressed form of A matrix
                                          // as per Bridson p77
                                          // we store everything as an int then
                                          // multiply by scale factor later
};
