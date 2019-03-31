#include <Eigen/Core>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <chrono>
#include <thread>

#include "Config.h"
#include "sim.h"
#include "OpenGLLoadShader.h"

#define SCREEN_WIDTH 1920
#define SCREEN_HEIGHT 1080

using namespace Eigen;

int main()
{
  // Check and Initialize Libs
  // Eigen
  std::cout << "Eigen version : " << EIGEN_MAJOR_VERSION << "."
      << EIGEN_MINOR_VERSION << std::endl << std::endl;

  // GLFW
  glewExperimental = true;
  if( !glfwInit() ) {
    fprintf( stderr, "Failed to initialize GLFW\n" );
    return -1;
  } else {
    std::cout << "GLFW Started." << std::endl << std::endl;
  }

  // Initialize our simulation
  Vector2i gridSize(8, 8);
  Vector2d gridLengths(1.0, 1.0);
  Vector2d gridSpacing = gridLengths.cwiseQuotient(gridSize.cast<double>());
  Sim Sim(gridSize, gridSpacing);
  Sim.setBodyAcceleration(Vector2d(GRAVITY_X, GRAVITY_Y));
  Sim.setBlockToMaterial(1,1,2,5,Material::fluid);
  Sim.setBlockToMaterial(0,0,3,1,Material::solid);
  Sim.init();
  std::cout << "Starting sim with grid size: " << std::endl << Sim._gridSize << std::endl <<
                            "and grid spacing " << std::endl << Sim._gridSpacing << std::endl;








  // Start GL Window
  glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // We don't want the old OpenGL

  // Open a window and create its OpenGL context
  GLFWwindow* window; // (In the accompanying source code, this variable is global for simplicity)
  window = glfwCreateWindow( SCREEN_WIDTH, SCREEN_HEIGHT, "Hello World", NULL, NULL);
  if( window == NULL ){
      fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
      glfwTerminate();
      return -1;
  }
  glfwMakeContextCurrent(window); // Initialize GLEW
  glewExperimental=true; // Needed in core profile
  if (glewInit() != GLEW_OK) {
      fprintf(stderr, "Failed to initialize GLEW\n");
      return -1;
  }
  // setup position
  glViewport( 0.0f, 0.0f, SCREEN_WIDTH, SCREEN_HEIGHT );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity( );
  glOrtho( 0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0, 1 );
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity( );

  // setup GL points
  //glPointSize(1.0);
  //glClearColor(0.0, 0.0, 0.0, 1.0);
  //glClear(GL_COLOR_BUFFER_BIT);
  //glColor3f(1.0, 1.0, 1.0);
  // load GL shaders
  //GLuint programID = LoadShaders( "OpenGLvertex.vertexshader", "OpenGLfragment.fragmentshader" );
  //glUseProgram(programID);

  // start the simulation loop
  while (true) {

    // advance in time
    Sim.update(0.01);
    std::cout << "MaterialGrid: " << std::endl << Sim._materialField->material().transpose().colwise().reverse() << std::endl;

    // plot the points
    GLfloat vertex[2];
    glClearColor(1.0, 0.0, 0.0, 0.0);
    glClear( GL_COLOR_BUFFER_BIT );
    glColor3f(1.0, 1.0, 1.0);
    //glEnable( GL_POINT_SMOOTH );
    //glDisable(GL_DEPTH_TEST);
  //  glEnableClientState( GL_VERTEX_ARRAY );
    glPointSize( 10.0 );

    glBegin(GL_POINTS);
      for (int i=0; i<Sim.nParticles(); i++) {
        Vector2d p = Sim.particlePosition(i);
        //glVertex2f((float)p.x()*SCREEN_WIDTH, (float)p.y()*SCREEN_HEIGHT);
        glVertex3f(0.0f, 0.0f, 0.0f);
      }
    glEnd();
    //glDisable( GL_POINT_SMOOTH );
    glfwSwapBuffers( window );

  //  for (int i=0; i<Sim.nParticles(); i++) {
  //    vertex[0] = (GLfloat)Sim.particlePosition(i).x();
  //    vertex[1] = (GLfloat)Sim.particlePosition(i).y();
  //    glVertexPointer( 2, GL_FLOAT, 0, vertex );
  //    glDrawArrays( GL_POINTS, 0, 1 );
  //  }
  //  glDisableClientState( GL_VERTEX_ARRAY );

    //glBegin(GL_POINTS);
    //  for (int i=0; i<Sim.nParticles(); i++) {
    //    Vector2d p = Sim.particlePosition(i);
    //    glVertex2f((float)p.x(), (float)p.y());
    //  }
    //glEnd();


    std::this_thread::sleep_for(std::chrono::milliseconds(100));

  }

  glfwTerminate();
  std::cout << "MACVelocityField U: " << std::endl << Sim._velocityField->U().transpose().colwise().reverse() << std::endl;
  std::cout << "MACVelocityField V: " << std::endl << Sim._velocityField->V().transpose().colwise().reverse() << std::endl;
  std::cout << "MaterialGrid: " << std::endl << Sim._materialField->material().transpose().colwise().reverse() << std::endl;

  return 0;
}
