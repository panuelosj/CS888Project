#pragma once

#include <Eigen/Core>
#include "Config.h"

using namespace Eigen;

class MACgridVelocity
{
public:
  MACgridVelocity();
  MACgridVelocity(Vector2i gridSize, Vector2d gridSpacing);
  ~MACgridVelocity();

  void update(float timestep);
  void extrapolateVelocities();

  // DATA READS
  double U(unsigned int i, unsigned int j);
  double V(unsigned int i, unsigned int j);
  MatrixXd U();
  MatrixXd V();
  int uValid(unsigned int i, unsigned int j);
  int vValid(unsigned int i, unsigned int j);
  void printU();
  void printV();
  Vector2d gridIndexToWorldspaceU(unsigned int i, unsigned int j);
  Vector2d gridIndexToWorldspaceV(unsigned int i, unsigned int j);
  Vector2d gridIndexToWorldspaceU(Vector2i g);
  Vector2d gridIndexToWorldspaceV(Vector2i g);
  bool isValidUIndex(int i, int j);
  bool isValidVIndex(int i, int j);
  // DATA WRITES
  void clear();
  void copyInData(MACgridVelocity* dataIn);
  void subAllData(MACgridVelocity* dataIn);
  void setU(unsigned int i, unsigned int j, double newU);
  void setV(unsigned int i, unsigned int j, double newV);
  void setU(MatrixXd newU);
  void setV(MatrixXd newV);
  void addU(unsigned int i, unsigned int j, double addU);
  void addV(unsigned int i, unsigned int j, double addV);
  void subU(unsigned int i, unsigned int j, double subU);
  void subV(unsigned int i, unsigned int j, double subV);
  void divideU(unsigned int i, unsigned int j, double divisorU);
  void divideV(unsigned int i, unsigned int j, double divisorV);
  void setUValid(unsigned int i, unsigned int j, double newUValid);
  void setVValid(unsigned int i, unsigned int j, double newVValid);
  void setBoundaryVelocities();
  // body forces
  void setBodyAcceleration(Vector2d b);
  void addBodyAcceleration(Vector2d b);

private:
  void _update();
  void _applyBodyForces();
  void _extrapolateVelocities();
  // returns true if a cell is in the valid grid range
  inline bool _isInRangeU(unsigned int i, unsigned int j);
  inline bool _isInRangeV(unsigned int i, unsigned int j);

  // definition of grid
  Vector2i _gridSize = Vector2i(0, 0);
  Vector2d _gridSpacing = Vector2d(0.0, 0.0);
  double _outOfRangeVal = VELOCITY_OUT_OF_RANGE_VALUE;
  double _timestep = 0.0;
  Vector2d _bodyAcceleration = Vector2d(0.0, 0.0);

  // data
  MatrixXd _u;
  MatrixXd _v;
  MatrixXi _uValid;
  MatrixXi _vValid;
};
