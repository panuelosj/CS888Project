#pragma once

#include <Eigen/Core>
#include "Config.h"
#include "MaterialGrid.h"

using namespace Eigen;

class Particles
{
public:
  Particles();
  Particles(MaterialGrid *materialField);
  ~Particles();

  // DATA READS
  int nParticles();
  Vector2d particleToWorldspace(unsigned int p);
  Vector2d particleToReferencespace(unsigned int p);
  Vector2i particleToGridIndex(unsigned int p);
  Vector2d particleVelocity(unsigned int p);
  // DATA WRITES
  void setParticleVelocity(unsigned int p, Vector2d newVelocity);
  // PRINTS
  void printParticleVelocities();

private:
  void _seedParticles();                    // seeds 2x2 particles per fluid grid cell
  void _jitterParticles();
  double _randomJitter();                   // gives a random jitter based on the internal jitter factor

  // constants
  double        _jitterFactor = D_JITTER_FACTOR;
  // variables
  MatrixXd      _positions;                 // _nParticles x 2 matrix
  MatrixXd      _velocities;                // _nParticles x 2 matrix
  unsigned int  _nParticles = 0;
  MaterialGrid  *_materialField;
};
