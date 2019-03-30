#include "MaterialGrid.h"
#include <iostream>


//  ######   #######  ##    ##  ######  ######## ########
// ##    ## ##     ## ###   ## ##    ##    ##    ##     ##
// ##       ##     ## ####  ## ##          ##    ##     ##
// ##       ##     ## ## ## ##  ######     ##    ########
// ##       ##     ## ##  ####       ##    ##    ##   ##
// ##    ## ##     ## ##   ### ##    ##    ##    ##    ##
//  ######   #######  ##    ##  ######     ##    ##     ##

MaterialGrid::MaterialGrid() {

}

MaterialGrid::MaterialGrid(Vector2i gridSize, Vector2d gridSpacing) :
  _gridSize( gridSize ),
  _gridSpacing( gridSpacing )
{
  // create empty grid with solid bounding box
  _material = MatrixXi::Ones(_gridSize(0), _gridSize(1)) * (int)Material::empty;
  _updateCellCount();
}

MaterialGrid::~MaterialGrid() {

}


// ########  ##     ## ########  ##       ####  ######
// ##     ## ##     ## ##     ## ##        ##  ##    ##
// ##     ## ##     ## ##     ## ##        ##  ##
// ########  ##     ## ########  ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##    ##
// ##         #######  ########  ######## ####  ######

// =============================================================================
// ------------------------- DATA READS ----------------------------------------
// =============================================================================

Vector2i MaterialGrid::gridSize() {
  // returns the size of the grid
  return _gridSize;
}

Vector2d MaterialGrid::gridSpacing() {
  // returns the dx and dy
  return _gridSpacing;
}

Material MaterialGrid::material(unsigned int i, unsigned int j) {
  // returns the material type at grid index (i,j)
  return (Material)_material(i,j);
}

MatrixXi MaterialGrid::material() {
  // returns the entire material field matrix
  return _material;
}

// returns true if grid index (i,j) has the corresponding material type
bool MaterialGrid::isSolid(unsigned int i, unsigned int j) {
  if (_isInRange(i,j))
    return (_material(i, j) == (int)Material::solid);
  else
    return true;              // we choose to interpret outside cells as solids
}
bool MaterialGrid::isFluid(unsigned int i, unsigned int j) {
  if (_isInRange(i,j))
    return (_material(i,j) == (int)Material::fluid);
  else
    return false;              // we choose to interpret outside cells as solids
}
bool MaterialGrid::isEmpty(unsigned int i, unsigned int j) {
  if (_isInRange(i,j))
    return (_material(i,j) == (int)Material::empty);
  else
    return false;              // we choose to interpret outside cells as solids
}

// returns number of cells of each material type
unsigned int MaterialGrid::nSolid() {
  return _nSolid;
}
unsigned int MaterialGrid::nFluid() {
  return _nFluid;
}
unsigned int MaterialGrid::nEmpty() {
  return _nEmpty;
}

// =============================================================================
// ------------------------- DATA SETS -----------------------------------------
// =============================================================================

void MaterialGrid::setBlockToMaterial(unsigned int i, unsigned int j, unsigned int p, unsigned int q, Material newMaterial) {
  // sets the material of a block starting at (i,j) of width (p,q) to a material
  _material.block(i,j,p,q) = MatrixXi::Ones(p,q) * (int)newMaterial;

  // TODO: a faster way of doing this
  _updateCellCount();
}

// ########  ########  #### ##     ##    ###    ######## ########
// ##     ## ##     ##  ##  ##     ##   ## ##      ##    ##
// ##     ## ##     ##  ##  ##     ##  ##   ##     ##    ##
// ########  ########   ##  ##     ## ##     ##    ##    ######
// ##        ##   ##    ##   ##   ##  #########    ##    ##
// ##        ##    ##   ##    ## ##   ##     ##    ##    ##
// ##        ##     ## ####    ###    ##     ##    ##    ########

void MaterialGrid::_updateCellCount() {
  _nSolid = 0;
  _nFluid = 0;
  _nEmpty = 0;

  for (unsigned int i=0; i<_gridSize(0); i++) {
    for (unsigned int j=0; j<_gridSize(1); j++) {
      if (isSolid(i,j)) _nSolid++;
      else if (isFluid(i,j)) _nFluid++;
      else if (isEmpty(i,j)) _nEmpty++;
    }
  }
}

inline bool MaterialGrid::_isInRange(unsigned int i, unsigned int j) {
  // checks if a grid index is a valid index
  bool retval = true;
  if (i < 0 || i >= _gridSize(0)) retval = false;
  if (j < 0 || j >= _gridSize(1)) retval = false;

  return retval;
}
