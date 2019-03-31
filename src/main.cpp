#include <Eigen/Core>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <chrono>
#include <thread>

#include "Config.h"
#include "sim.h"
#include "OpenGLLoadShader.h"

#define SCREEN_WIDTH 1080
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
  Vector2i gridSize(16, 16);
  Vector2d gridLengths(1.0, 1.0);
  Vector2d gridSpacing = gridLengths.cwiseQuotient(gridSize.cast<double>());
  Sim Sim(gridSize, gridSpacing);
  Sim.setBodyAcceleration(Vector2d(GRAVITY_X, GRAVITY_Y));
  Sim.setBlockToMaterial(7,0,7,10,Material::fluid);
  Sim.setBlockToMaterial(0,0,3,3,Material::solid);
  Sim.init();
  std::cout << "Starting sim with grid size: " << std::endl << Sim._gridSize << std::endl <<
                            "and grid spacing " << std::endl << Sim._gridSpacing << std::endl;








  // Start GL Window
  glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // We don't want the old OpenGL

  // Open a window and create its OpenGL context
  GLFWwindow* window; // (In the accompanying source code, this variable is global for simplicity)
  window = glfwCreateWindow( SCREEN_WIDTH, SCREEN_HEIGHT, "Hello World", NULL, NULL);
  if( window == NULL ){
      fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
      glfwTerminate();
      return -1;
  }
  glfwMakeContextCurrent(window);
  // Initialize GLEW
  glewExperimental=true; // Needed in core profile
  if (glewInit() != GLEW_OK) {
      fprintf(stderr, "Failed to initialize GLEW\n");
  		glfwTerminate();
      return -1;
  }

	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

  // create VAO
  GLuint VertexArrayID;
  glGenVertexArrays(1, &VertexArrayID);
  glBindVertexArray(VertexArrayID);


  // load GL shaders
  GLuint programID = LoadShaders( "OpenGLvertex.vertexshader", "OpenGLfragment.fragmentshader" );

  // setup position
  glViewport( 0.0f, 0.0f, SCREEN_WIDTH, SCREEN_HEIGHT );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity( );
  glOrtho( 0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0, 1 );
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity( );

  // setup GL vertex buffer


	//static const GLfloat g_vertex_buffer_data[] = {
	//	-1.0f, -1.0f, 0.0f,
	//	 1.0f, -1.0f, 0.0f,
	//	 0.0f,  1.0f, 0.0f,
	//};

  // This will identify our vertex buffer
  //GLuint vertexbuffer;
  // Generate 1 buffer, put the resulting identifier in vertexbuffer
  //glGenBuffers(1, &vertexbuffer);
  // The following commands will talk about our 'vertexbuffer' buffer
  //glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
  // Give our vertices to OpenGL.
  //glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);


  getchar();

  // start the simulation loop
  while (true) {

    // advance in time
    Sim.update(0.01);
    std::cout << "MaterialGrid: " << std::endl << Sim._materialField->material().transpose().colwise().reverse() << std::endl;



    // setup vertex buffer data
    GLfloat g_vertex_buffer_data[Sim.nParticles()*3];
    for (int i=0; i<Sim.nParticles(); i++) {
      g_vertex_buffer_data[3*i] = ((GLfloat)Sim.particlePosition(i).x()*2.0f) - 1.0f;
      g_vertex_buffer_data[3*i+1] = ((GLfloat)Sim.particlePosition(i).y()*2.0f) - 1.0f;
      g_vertex_buffer_data[3*i+2] = 0.0f;
    }

    /*static const GLfloat g_vertex_buffer_data[] = {
  		-1.0f, -1.0f, 0.0f,
  		 1.0f, -1.0f, 0.0f,
  		 0.0f,  1.0f, 0.0f,
  	};*/
  	GLuint vertexbuffer;
  	glGenBuffers(1, &vertexbuffer);
  	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_DYNAMIC_DRAW);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPointSize(10.0);
    glEnable(GL_BLEND);
    glUseProgram(programID);
    // 1st attribute buffer : vertices
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(
       0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
       3,                  // size
       GL_FLOAT,           // type
       GL_FALSE,           // normalized?
       0,                  // stride
       (void*)0            // array buffer offset
    );
    // Draw the triangle !
    glDrawArrays(GL_POINTS, 0, Sim.nParticles());
    glDisable(GL_BLEND);
  	glDisableVertexAttribArray(0);

  	// Swap buffers
  	glfwSwapBuffers(window);
/*
    // setup vertex buffer data
    GLfloat *g_vertex_buffer_data;
    g_vertex_buffer_data = new GLfloat[Sim.nParticles()*2];
    for (int i=0; i<Sim.nParticles(); i++) {
      g_vertex_buffer_data[2*i] = (GLfloat)Sim.particlePosition(i).x();
      g_vertex_buffer_data[2*i+1] = (GLfloat)Sim.particlePosition(i).y();
    }
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_DYNAMIC_DRAW);
    glClearColor(1.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glDrawArrays(GL_TRIANGLES, 0, 9);
    glDisableVertexAttribArray(0);*/
    //glDisableClientState( GL_VERTEX_ARRAY );

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
