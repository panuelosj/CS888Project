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
  // full matrix read
  MatrixXd positions();
  MatrixXd velocities();
  // single particle reads
  Vector2d particleToWorldspace(unsigned int p);
  Vector2d particleToReferencespace(unsigned int p);
  Vector2d positionToReferencespace(Vector2d p);
  Vector2i particleToGridIndex(unsigned int p);
  Vector2i positionToGridIndex(Vector2d p);
  Vector2d particleVelocity(unsigned int p);
  double particleSpeed(unsigned int p);
  double maxParticleSpeed();
  // DATA WRITES
  void setParticlePosition(unsigned int idx, Vector2d newPosition);
  void setParticleVelocity(unsigned int idx, Vector2d newVelocity);
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
