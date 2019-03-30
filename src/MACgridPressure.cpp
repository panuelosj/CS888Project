#include "MACgridPressure.h"

MACgridPressure::MACgridPressure() {

}

MACgridPressure::MACgridPressure(Vector2i gridSize, Vector2f gridSpacing) :
  _gridSize( gridSize ),
  _gridSpacing( gridSpacing ) {
  _init();
}

MACgridPressure::~MACgridPressure() {

}

void MACgridPressure::_init() {
  _P = MatrixXf::Zero(_gridSize(0), _gridSize(1));
}
