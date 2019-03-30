#include <iostream>
#include <Eigen/Core>
#include "Config.h"
#include "sim.h"

using namespace Eigen;

int main()
{
  // check libs
  std::cout << "Eigen version : " << EIGEN_MAJOR_VERSION << "."
      << EIGEN_MINOR_VERSION << std::endl ;

  // let's do this

  // setup problem
  Vector2i gridSize(8, 8);
  Vector2d gridLengths(1.0, 1.0);
  Vector2d gridSpacing = gridLengths.cwiseQuotient(gridSize.cast<double>());
  Sim Sim(gridSize, gridSpacing);
  Sim.setBodyAcceleration(Vector2d(GRAVITY_X, GRAVITY_Y));
  Sim.setBlockToMaterial(1,1,2,2,Material::fluid);
  Sim.setBlockToMaterial(0,0,2,1,Material::solid);
  Sim.init();

  std::cout << "Starting sim with grid size: " << std::endl << Sim._gridSize << std::endl <<
                            "and grid spacing " << std::endl << Sim._gridSpacing << std::endl;

  Sim.update(1.0);

  std::cout << "MACVelocityField U: " << std::endl << Sim._velocityField->U().transpose().colwise().reverse() << std::endl;
  std::cout << "MACVelocityField V: " << std::endl << Sim._velocityField->V().transpose().colwise().reverse() << std::endl;
  std::cout << "MaterialGrid: " << std::endl << Sim._materialField->material().transpose().colwise().reverse() << std::endl;

  return 0;
}
