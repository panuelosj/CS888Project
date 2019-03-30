#include "GridIndexMapping.h"

//  ######   #######  ##    ##  ######  ######## ########
// ##    ## ##     ## ###   ## ##    ##    ##    ##     ##
// ##       ##     ## ####  ## ##          ##    ##     ##
// ##       ##     ## ## ## ##  ######     ##    ########
// ##       ##     ## ##  ####       ##    ##    ##   ##
// ##    ## ##     ## ##   ### ##    ##    ##    ##    ##
//  ######   #######  ##    ##  ######     ##    ##     ##

GridIndexMapping::GridIndexMapping() {

}

GridIndexMapping::GridIndexMapping(Vector2i gridSize) :
  _gridSize( gridSize ),
  // derived variables
  _gridVectorLength( _gridSize(0) * _gridSize(1) ),
  // constructed variables
  _indices( VectorXi(_gridVectorLength) )
{
  // let's just init with the index number
  for (unsigned int idx=0; idx<_gridVectorLength; idx++) {
    _indices[idx] = idx;
  }
}

GridIndexMapping::~GridIndexMapping() {

}


// ########  ##     ## ########  ##       ####  ######
// ##     ## ##     ## ##     ## ##        ##  ##    ##
// ##     ## ##     ## ##     ## ##        ##  ##
// ########  ##     ## ########  ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##    ##
// ##         #######  ########  ######## ####  ######

int GridIndexMapping::gridToVector(int i, int j) {
  if (_isInRange(i,j)) {
    int flatidx = _flattenIndex( i, j );
    return _indices[flatidx];
  } else {
    return GRID_TO_VECTOR_OUTSIDE_RANGE_VAL;
  }
}

int GridIndexMapping::gridToVector(Vector2i g) {
  // Wrapper function to support the Eigen::Vector2i input
  return gridToVector( g(0), g(1) );
}


// ########  ########  #### ##     ##    ###    ######## ########
// ##     ## ##     ##  ##  ##     ##   ## ##      ##    ##
// ##     ## ##     ##  ##  ##     ##  ##   ##     ##    ##
// ########  ########   ##  ##     ## ##     ##    ##    ######
// ##        ##   ##    ##   ##   ##  #########    ##    ##
// ##        ##    ##   ##    ## ##   ##     ##    ##    ##
// ##        ##     ## ####    ###    ##     ##    ##    ########

inline unsigned int GridIndexMapping::_flattenIndex(int i, int j) {
  return (unsigned int)i + (unsigned int)_gridSize(0)*(unsigned int)j;
}

inline unsigned int GridIndexMapping::_flattenIndex(Vector2i g) {
  // Wrapper function to support the Eigen::Vector2i input
  return _flattenIndex( g(0), g(1) );
}

inline bool GridIndexMapping::_isInRange(unsigned int i, unsigned int j) {
  // checks if a grid index is a valid index
  bool retval = true;
  if (i < 0 || i >= _gridSize(0)) retval = false;
  if (j < 0 || j >= _gridSize(1)) retval = false;

  return retval;
}
