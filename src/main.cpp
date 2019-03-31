#include <Eigen/Core>
#include <iostream>

#include <GLFW/glfw3.h>
#include <chrono>
#include <thread>

#include "Config.h"
#include "sim.h"

using namespace Eigen;

int main()
{
  // Check and Initialize Libs
  // Eigen
  std::cout << "Eigen version : " << EIGEN_MAJOR_VERSION << "."
      << EIGEN_MINOR_VERSION << std::endl << std::endl;

  // GLFW
  //glewExperimental = true;
  if( !glfwInit() ) {
    fprintf( stderr, "Failed to initialize GLFW\n" );
    return -1;
  } else {
    std::cout << "GLFW Started." << std::endl << std::endl;
  }

  // Start GL Window
  glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // We don't want the old OpenGL

  // Open a window and create its OpenGL context
  GLFWwindow* window; // (In the accompanying source code, this variable is global for simplicity)
  window = glfwCreateWindow( 1024, 768, "Tutorial 01", NULL, NULL);
  if( window == NULL ){
      fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
      glfwTerminate();
      return -1;
  }
  glfwMakeContextCurrent(window); // Initialize GLEW
  glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

  // Initialize our simulation
  Vector2i gridSize(8, 8);
  Vector2d gridLengths(1.0, 1.0);
  Vector2d gridSpacing = gridLengths.cwiseQuotient(gridSize.cast<double>());
  Sim Sim(gridSize, gridSpacing);
  Sim.setBodyAcceleration(Vector2d(GRAVITY_X, GRAVITY_Y));
  Sim.setBlockToMaterial(1,1,2,5,Material::fluid);
  Sim.setBlockToMaterial(0,0,8,1,Material::solid);
  Sim.init();
  std::cout << "Starting sim with grid size: " << std::endl << Sim._gridSize << std::endl <<
                            "and grid spacing " << std::endl << Sim._gridSpacing << std::endl;

  // start the simulation loop

  // setup gl points
  glPointSize(2.0);
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glColor3f(1.0, 1.0, 1.0);

  while (true) {

    // advance in time
    Sim.update(0.01);
    std::cout << "MaterialGrid: " << std::endl << Sim._materialField->material().transpose().colwise().reverse() << std::endl;

    // plot the points
    for(int i = 0; i < Sim.nParticles(); i++) {
      Vector2d p = Sim.particlePosition(i);
      glBegin(GL_POINTS);
        glVertex2f((float)p.x(), (float)p.y());
      glEnd();
      glFlush();
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(100));

  //  while ( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS ) {
  //    // wait
  //  }
  }
  std::cout << "MACVelocityField U: " << std::endl << Sim._velocityField->U().transpose().colwise().reverse() << std::endl;
  std::cout << "MACVelocityField V: " << std::endl << Sim._velocityField->V().transpose().colwise().reverse() << std::endl;
  std::cout << "MaterialGrid: " << std::endl << Sim._materialField->material().transpose().colwise().reverse() << std::endl;

  return 0;
}
