#include "MACgridVelocity.h"
#include "assert.h"
#include <iostream>

//  ######   #######  ##    ##  ######  ######## ########
// ##    ## ##     ## ###   ## ##    ##    ##    ##     ##
// ##       ##     ## ####  ## ##          ##    ##     ##
// ##       ##     ## ## ## ##  ######     ##    ########
// ##       ##     ## ##  ####       ##    ##    ##   ##
// ##    ## ##     ## ##   ### ##    ##    ##    ##    ##
//  ######   #######  ##    ##  ######     ##    ##     ##

MACgridVelocity::MACgridVelocity() {

}

MACgridVelocity::MACgridVelocity(Vector2i gridSize, Vector2d gridSpacing) :
  _gridSize( gridSize ),
  _gridSpacing( gridSpacing )
{
  // create a constant zero field velocity
  _u = MatrixXd::Zero(_gridSize(0) + 1, _gridSize(1));
  _v = MatrixXd::Zero(_gridSize(0), _gridSize(1) + 1);
  // set all fields to invalid
  _uValid = MatrixXi::Zero(_gridSize(0) + 1, _gridSize(1));
  _vValid = MatrixXi::Zero(_gridSize(0), _gridSize(1) + 1);
}

MACgridVelocity::~MACgridVelocity() {

}


// ########  ##     ## ########  ##       ####  ######
// ##     ## ##     ## ##     ## ##        ##  ##    ##
// ##     ## ##     ## ##     ## ##        ##  ##
// ########  ##     ## ########  ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##    ##
// ##         #######  ########  ######## ####  ######

void MACgridVelocity::update(float timestep) {
  // wrapper function to setup all variables and calls the private update
  _timestep = timestep;
  _update();
}

void MACgridVelocity::extrapolateVelocities() {
  // wrapper function to call the private update
  _extrapolateVelocities();
}

// =============================================================================
// ------------------------- DATA READS ----------------------------------------
// =============================================================================

// reading matrix values
double MACgridVelocity::U(unsigned int i, unsigned int j) {
  if (!_isInRangeU(i,j)) {
    return _outOfRangeVal;
  }
  else {
    return _u(i, j);
  }
}
double MACgridVelocity::V(unsigned int i, unsigned int j) {
  if (!_isInRangeV(i,j)) {
    return _outOfRangeVal;
  }
  else {
    return _v(i, j);
  }
}
// reading entire matrix
MatrixXd MACgridVelocity::U() {
  return _u;
}
MatrixXd MACgridVelocity::V() {
  return _v;
}

// printing them
void MACgridVelocity::printU() {
  std::cout << _u.transpose().colwise().reverse() << std::endl;
}
void MACgridVelocity::printV() {
  std::cout << _v.transpose().colwise().reverse() << std::endl;
}

// note that for these transformations, the direct i*dx multiplication gives the
// bottom left corner, NOT THE CELL CENTER, so the +0.5's need to be accounted
// for properly
Vector2d MACgridVelocity::gridIndexToWorldspaceU(unsigned int i, unsigned int j) {
  double x = _gridSpacing(0)*(double)i;
  double y = _gridSpacing(1)*(((double)j)+0.5);

  return Vector2d(x,y);
}
Vector2d MACgridVelocity::gridIndexToWorldspaceV(unsigned int i, unsigned int j) {
  double x = _gridSpacing(0)*(((double)i)+0.5);;
  double y = _gridSpacing(1)*(double)j;

  return Vector2d(x,y);
}

// wrapper functions to accept Vector2i inputs
Vector2d MACgridVelocity::gridIndexToWorldspaceU(Vector2i g) {
  return gridIndexToWorldspaceU(g.x(), g.y());
}
Vector2d MACgridVelocity::gridIndexToWorldspaceV(Vector2i g) {
  return gridIndexToWorldspaceV(g.x(), g.y());
}

// checks to make sure the indices are valid
bool MACgridVelocity::isValidUIndex(int i, int j) {
  bool retval = true;
  if (i<0 || i>=_gridSize(0)+1) retval = false;
  if (j<0 || j>=_gridSize(1)) retval = false;
  return retval;
}
bool MACgridVelocity::isValidVIndex(int i, int j) {
  bool retval = true;
  if (i<0 || i>=_gridSize(0)) retval = false;
  if (j<0 || j>=_gridSize(1)+1) retval = false;
  return retval;
}

// =============================================================================
// ------------------------- DATA WRITES ---------------------------------------
// =============================================================================

void MACgridVelocity::clear() {
  _u = MatrixXd::Zero(_gridSize(0) + 1, _gridSize(1));
  _v = MatrixXd::Zero(_gridSize(0), _gridSize(1) + 1);
}

void MACgridVelocity::copyInData(MACgridVelocity* dataIn) {
  _u = dataIn->U();
  _v = dataIn->V();
}

void MACgridVelocity::subAllData(MACgridVelocity* dataIn) {
  _u -= dataIn->U();
  _v -= dataIn->V();
}

// setting matrix values
void MACgridVelocity::setU(unsigned int i, unsigned int j, double newU) {
  assert(_isInRangeU(i,j));
  _u(i,j) = newU;
}
void MACgridVelocity::setV(unsigned int i, unsigned int j, double newV) {
  assert(_isInRangeV(i,j));
  _v(i,j) = newV;
}
void MACgridVelocity::addU(unsigned int i, unsigned int j, double addU) {
  assert(_isInRangeU(i,j));
  _u(i,j) += addU;
}
void MACgridVelocity::addV(unsigned int i, unsigned int j, double addV) {
  assert(_isInRangeV(i,j));
  _v(i,j) += addV;
}
void MACgridVelocity::subU(unsigned int i, unsigned int j, double subU) {
  assert(_isInRangeU(i,j));
  _u(i,j) -= subU;
}
void MACgridVelocity::subV(unsigned int i, unsigned int j, double subV) {
  assert(_isInRangeV(i,j));
  _v(i,j) -= subV;
}
void MACgridVelocity::divideU(unsigned int i, unsigned int j, double divisorU) {
  assert(_isInRangeU(i,j));
  _u(i,j) /= divisorU;
}
void MACgridVelocity::divideV(unsigned int i, unsigned int j, double divisorV) {
  assert(_isInRangeV(i,j));
  _v(i,j) /= divisorV;
}
void MACgridVelocity::setUValid(unsigned int i, unsigned int j, double newUValid) {
  assert(_isInRangeU(i,j));
  _uValid(i,j) = newUValid;
}
void MACgridVelocity::setVValid(unsigned int i, unsigned int j, double newVValid) {
  assert(_isInRangeV(i,j));
  _vValid(i,j) = newVValid;
}

void MACgridVelocity::setBodyAcceleration(Vector2d b) {
  _bodyAcceleration = b;
}
void MACgridVelocity::addBodyAcceleration(Vector2d b) {
  _bodyAcceleration += b;
}

void MACgridVelocity::setBoundaryVelocities() {
  _u.rightCols<1>() = VectorXd::Zero(_u.rows());
  _u.leftCols<1>() = VectorXd::Zero(_u.rows());
  _v.topRows<1>() = RowVectorXd::Zero(_v.cols());
  _v.bottomRows<1>() = RowVectorXd::Zero(_v.cols());
}

// ########  ########  #### ##     ##    ###    ######## ########
// ##     ## ##     ##  ##  ##     ##   ## ##      ##    ##
// ##     ## ##     ##  ##  ##     ##  ##   ##     ##    ##
// ########  ########   ##  ##     ## ##     ##    ##    ######
// ##        ##   ##    ##   ##   ##  #########    ##    ##
// ##        ##    ##   ##    ## ##   ##     ##    ##    ##
// ##        ##     ## ####    ###    ##     ##    ##    ########

void MACgridVelocity::_update() {
  _applyBodyForces();
}

void MACgridVelocity::_applyBodyForces() {
  // u = u + deltaT * bodyForce
  _u += _timestep*_bodyAcceleration(0)*MatrixXd::Ones(_u.rows(),_u.cols());
  _v += _timestep*_bodyAcceleration(1)*MatrixXd::Ones(_v.rows(),_v.cols());
}

void MACgridVelocity::_extrapolateVelocities() {
  // Function pulls the known velocities into the unknown regions
  bool finished;

  // UPDATE Us
  finished = false;
  while( !finished ) {
    finished = true;

    // loop through the values and find if any are unknown
    for (unsigned int i=0; i<_gridSize(0)+1; i++) {
      for (unsigned int j=0; j<_gridSize(1); j++) {
        // check if it is invalid
        if (_uValid(i,j) == 0) {
          // variables to keep track
          double sum = 0;       // sum of valid neighbour samples
          int nSampled = 0;     // number of valid neighbours
          // now look for neighbours which are valid and add them
            // to the running tally
          // left
          if (i-1 >= 0 && i-1 <_gridSize(0)+1) {
            if (_uValid(i-1,j) != 0) {
              sum += _u(i-1,j);
              nSampled++;
            }
          }
          // right
          if (i+1 >= 0 && i+1 <_gridSize(0)+1) {
            if (_uValid(i+1,j) != 0) {
              sum += _u(i+1,j);
              nSampled++;
            }
          }
          // down
          if (j-1 >= 0 && j-1 <_gridSize(1)) {
            if (_uValid(i,j-1) != 0) {
              sum += _u(i,j-1);
              nSampled++;
            }
          }
          // up
          if (j+1 >= 0 && j+1 <_gridSize(1)) {
            if (_uValid(i,j+1) != 0) {
              sum += _u(i,j+1);
              nSampled++;
            }
          }

          // did we find any valid samples?
          if (nSampled > 0) {
            // now set the value at the invalid point to the average of the
            // sampled points
            sum = sum / (double)nSampled;
            _u(i,j) = sum;
            // remember to set it to valid now
            _uValid(i,j) = 1;

          } else {
            // we didn't find any valid samples so this one is left unfinished
            finished = false;
          }
        }
      }
    }
  }


  // UPDATE Vs
  finished = false;
  while( !finished ) {
    finished = true;

    // loop through the values and find if any are unknown
    for (unsigned int i=0; i<_gridSize(0); i++) {
      for (unsigned int j=0; j<_gridSize(1)+1; j++) {
        // check if it is invalid
        if (_vValid(i,j) == 0) {
          // variables to keep track
          double sum = 0;       // sum of valid neighbour samples
          int nSampled = 0;     // number of valid neighbours
          // now look for neighbours which are valid and add them
            // to the running tally
          // left
          if (i-1 >= 0 && i-1 <_gridSize(0)) {
            if (_vValid(i-1,j) != 0) {
              sum += _v(i-1,j);
              nSampled++;
            }
          }
          // right
          if (i+1 >= 0 && i+1 <_gridSize(0)) {
            if (_vValid(i+1,j) != 0) {
              sum += _v(i+1,j);
              nSampled++;
            }
          }
          // down
          if (j-1 >= 0 && j-1 <_gridSize(1)+1) {
            if (_vValid(i,j-1) != 0) {
              sum += _v(i,j-1);
              nSampled++;
            }
          }
          // up
          if (j+1 >= 0 && j+1 <_gridSize(1)+1) {
            if (_vValid(i,j+1) != 0) {
              sum += _v(i,j+1);
              nSampled++;
            }
          }

          // did we find any valid samples?
          if (nSampled > 0) {
            // now set the value at the invalid point to the average of the
            // sampled points
            sum = sum / (double)nSampled;
            _v(i,j) = sum;
            // remember to set it to valid now
            _vValid(i,j) = 1;

          } else {
            // we didn't find any valid samples so this one is left unfinished
            finished = false;
          }
        }
      }
    }
  }
}


inline bool MACgridVelocity::_isInRangeU(unsigned int i, unsigned int j) {
  // checks if a grid index is a valid index
  bool retval = true;
  if (i < 0 || i >= _gridSize(0) + 1) retval = false;
  if (j < 0 || j >= _gridSize(1)) retval = false;

  return retval;
}

inline bool MACgridVelocity::_isInRangeV(unsigned int i, unsigned int j) {
  // checks if a grid index is a valid index
  bool retval = true;
  if (i < 0 || i >= _gridSize(0)) retval = false;
  if (j < 0 || j >= _gridSize(1) + 1) retval = false;

  return retval;
}
