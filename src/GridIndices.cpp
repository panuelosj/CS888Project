#include "GridIndices.h"

//  ######   #######  ##    ##  ######  ######## ########
// ##    ## ##     ## ###   ## ##    ##    ##    ##     ##
// ##       ##     ## ####  ## ##          ##    ##     ##
// ##       ##     ## ## ## ##  ######     ##    ########
// ##       ##     ## ##  ####       ##    ##    ##   ##
// ##    ## ##     ## ##   ### ##    ##    ##    ##    ##
//  ######   #######  ##    ##  ######     ##    ##     ##

GridIndices::GridIndices() {

}

GridIndices::GridIndices(MaterialGrid *materialField, Material materialType) {
  // constructor to create a list of the correct material type based on the
  // input material field

  if (materialType == Material::solid) {
    _nIndices = materialField->nSolid();
  }
  else if (materialType == Material::fluid) {
    _nIndices = materialField->nFluid();
  }
  else if (materialType == Material::empty) {
    _nIndices = materialField->nEmpty();
  }

  // initialize vector
  _xs = VectorXi(_nIndices);
  _ys = VectorXi(_nIndices);

  // setup list
  unsigned int idx = 0;
  for (unsigned int i=0; i<materialField->gridSize().x(); i++) {
    for (unsigned int j=0; j<materialField->gridSize().y(); j++) {
      if (materialType == materialField->material(i,j)) {
        _xs(idx) = i;
        _ys(idx) = j;
        idx++;
      }
    }
  }
}

GridIndices::~GridIndices() {

}


// ########  ##     ## ########  ##       ####  ######
// ##     ## ##     ## ##     ## ##        ##  ##    ##
// ##     ## ##     ## ##     ## ##        ##  ##
// ########  ##     ## ########  ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##    ##
// ##         #######  ########  ######## ####  ######

Vector2i GridIndices::index(unsigned int idx) {
  return Vector2i(_xs(idx), _ys(idx));
}

unsigned int GridIndices::size() {
  return _nIndices;
}

int GridIndices::find(unsigned int i, unsigned int j) {
  for (unsigned int idx=0; idx<_nIndices; idx++) {
    if (_xs[idx] == i && _ys[idx] == j) return idx;
  }
  // index not found
  return GRID_INDEX_NOT_FOUND_VAL;
}

// ########  ########  #### ##     ##    ###    ######## ########
// ##     ## ##     ##  ##  ##     ##   ## ##      ##    ##
// ##     ## ##     ##  ##  ##     ##  ##   ##     ##    ##
// ########  ########   ##  ##     ## ##     ##    ##    ######
// ##        ##   ##    ##   ##   ##  #########    ##    ##
// ##        ##    ##   ##    ## ##   ##     ##    ##    ##
// ##        ##     ## ####    ###    ##     ##    ##    ########
