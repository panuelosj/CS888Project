#include <Eigen/Core>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <pngwriter.h>
#include <chrono>
#include <thread>

#include "Config.h"
#include "sim.h"
#include "OpenGLLoadShader.h"

#define SCREEN_WIDTH 1080
#define SCREEN_HEIGHT 1080

using namespace Eigen;

// global variable for window context
GLFWwindow* window;
GLuint programID;
Sim* mySim;
// just make a list of colors
Vector3d colorSolid = Vector3d(0.1,0.18,0.16);
Vector3d colorFluid = Vector3d(0.55,0.6,0.68);
Vector3d colorEmpty = Vector3d(0.68,1.0,0.93);

int startOpenGL();
void endOpenGL();
void plotOpenGL();
void savePNG(unsigned int height, unsigned int width, char* filename);

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
  Vector2i gridSize(32, 64);
  Vector2d gridLengths(0.5, 1.0);
  Vector2d gridSpacing = gridLengths.cwiseQuotient(gridSize.cast<double>());
  mySim = new Sim(gridSize, gridSpacing);
  mySim->setBodyAcceleration(Vector2d(GRAVITY_X, GRAVITY_Y));

  // 16x16 full fluid
  //mySim->setBlockToMaterial(0,0,16,16,Material::fluid);
  // 16x16 half fluid
  //mySim->setBlockToMaterial(0,0,16,8,Material::fluid);
  // 16x16 half fluid drop
  //mySim->setBlockToMaterial(0,6,16,8,Material::fluid);

  // 32x32 block fluid drop
  //const char* filenameBase = "BoxDrop32x32_PICFLIP0.1/BoxDrop32x32";
  //const char* filenameBase = "BoxDrop32x32_FLIP/BoxDrop32x32";
  //const char* filenameBase = "BoxDrop32x32_PIC/BoxDrop32x32";
  //mySim->setBlockToMaterial(8,12,16,16,Material::fluid);

  // 16x16 dam break with box
  //mySim->setBlockToMaterial(4,0,8,10,Material::fluid);
  //mySim->setBlockToMaterial(0,0,3,3,Material::solid);

  // 16x16 dam break
  //mySim->setBlockToMaterial(4,0,8,15,Material::fluid);
  //mySim->setBlockToMaterial(0,0,8,15,Material::fluid);

  // 32x32 dam break
  //const char* filenameBase = "DamBreak32x32_PICFLIP0.1/DamBreak32x32";
  //const char* filenameBase = "DamBreak32x32_FLIP/DamBreak32x32";
  //const char* filenameBase = "DamBreak32x32_PIC/DamBreak32x32";
  //mySim->setBlockToMaterial(0,0,16,30,Material::fluid);

  //const char* filenameBase = "DamBreakB32x32_PICFLIP0.1/DamBreakB32x32";
  //const char* filenameBase = "DamBreakB32x32_FLIP/DamBreak32x32";
  //const char* filenameBase = "DamBreakB32x32_PIC/DamBreak32x32";
  //mySim->setBlockToMaterial(8,0,16,30,Material::fluid);

  // 16x16 plank
  //mySim->setBlockToMaterial(2,6,8,8,Material::fluid);
  //mySim->setBlockToMaterial(0,5,12,1,Material::solid);

  // 32x32 plank
  //mySim->setBlockToMaterial(4,11,16,16,Material::fluid);
  //mySim->setBlockToMaterial(0,10,24,1,Material::solid);

  // 16x16 maze
  //const char* filenameBase = "Maze16x16/Maze16x16";
  //mySim->setBlockToMaterial(1,11,10,5,Material::fluid);
  //mySim->setBlockToMaterial(0,10,12,1,Material::solid);
  //mySim->setBlockToMaterial(4,6,12,1,Material::solid);

  // 32x32 maze
  //const char* filenameBase = "Maze32x32_PICFLIP0.1/Maze32x32";
  //const char* filenameBase = "Maze32x32_FLIP/Maze32x32";
  //const char* filenameBase = "Maze32x32_PIC/Maze32x32";
  //mySim->setBlockToMaterial(2,21,20,11,Material::fluid);
  //mySim->setBlockToMaterial(0,20,24,1,Material::solid);
  //mySim->setBlockToMaterial(8,12,24,1,Material::solid);

  // 64x64 maze
  //const char* filenameBase = "Maze64x64_PICFLIP0.1/Maze64x64";
  //mySim->setBlockToMaterial(4,42,40,22,Material::fluid);
  //mySim->setBlockToMaterial(0,40,48,2,Material::solid);
  //mySim->setBlockToMaterial(16,24,48,2,Material::solid);

  // 64x64 complex
  //const char* filenameBase = "Complex64x64_PICFLIP0.1/Complex64x64";
  //mySim->setBlockToMaterial(4,42,20,20,Material::fluid);
  //mySim->setBlockToMaterial(9,13,1,1,Material::solid);
  //mySim->setBlockToMaterial(14,32,2,1,Material::solid);
  //mySim->setBlockToMaterial(18,5,7,4,Material::solid);
  //mySim->setBlockToMaterial(2,27,5,3,Material::solid);
  //mySim->setBlockToMaterial(32,0,16,1,Material::solid);
  //mySim->setBlockToMaterial(33,1,15,1,Material::solid);
  //mySim->setBlockToMaterial(34,2,14,1,Material::solid);
  //mySim->setBlockToMaterial(35,3,13,1,Material::solid);
  //mySim->setBlockToMaterial(36,4,12,1,Material::solid);
  //mySim->setBlockToMaterial(37,5,11,1,Material::solid);
  //mySim->setBlockToMaterial(38,6,10,1,Material::solid);
  //mySim->setBlockToMaterial(39,7, 9,1,Material::solid);
  //mySim->setBlockToMaterial(40,8, 8,1,Material::solid);
  //mySim->setBlockToMaterial(41,9, 7,1,Material::solid);
  //mySim->setBlockToMaterial(42,10,6,1,Material::solid);
  //mySim->setBlockToMaterial(43,11,5,1,Material::solid);
  //mySim->setBlockToMaterial(44,12,4,1,Material::solid);
  //mySim->setBlockToMaterial(45,13,3,1,Material::solid);
  //mySim->setBlockToMaterial(46,14,2,1,Material::solid);
  //mySim->setBlockToMaterial(47,15,1,1,Material::solid);

  const char* filenameBase = "Cone64x32_PICFLIP0.1/Cone64x32";
  mySim->setBlockToMaterial(6,42,20,20,Material::fluid);
  mySim->setBlockToMaterial(0,30,32,10,Material::solid);
  mySim->setBlockToMaterial(14,30,4,1,Material::empty);
  mySim->setBlockToMaterial(13,31,6,1,Material::empty);
  mySim->setBlockToMaterial(12,32,8,1,Material::empty);
  mySim->setBlockToMaterial(11,33,10,1,Material::empty);
  mySim->setBlockToMaterial(10,34,12,1,Material::empty);
  mySim->setBlockToMaterial(9,35,14,1,Material::empty);
  mySim->setBlockToMaterial(8,36,16,1,Material::empty);
  mySim->setBlockToMaterial(7,37,18,1,Material::empty);
  mySim->setBlockToMaterial(6,38,20,1,Material::empty);
  mySim->setBlockToMaterial(5,39,22,1,Material::empty);

  mySim->init();
  std::cout << "Starting sim with grid size: "
                  << std::endl << mySim->gridSize()
                  << std::endl << std::endl
            << "and grid spacing "
                  << std::endl << mySim->gridSpacing()
                  << std::endl << std::endl;

  // initialize OpenGL
  if (startOpenGL() != 0) {
    fprintf(stderr, "Failed to initialize OpenGL\n");
    return -1;
  }

  // show the initial state until user presses enter
  std::cout << std::endl;
  std::cout << "Simulation ready. Please press Enter to start:" << std::endl
            << "\tAwaiting user input..." << std::endl << std::endl;
  glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
  while( glfwGetKey(window, GLFW_KEY_ENTER) != GLFW_PRESS ) {
    plotOpenGL();
    glfwPollEvents();
  }

  unsigned int nFrame = 0;

  // start the simulation loop
  while (mySim->time() < 5.0) {
  //while (true) {
    // advance in time
    mySim->update(0.001);

    // do OpenGL
    plotOpenGL();

#ifdef LOG_SAVE_PNG
    // then save the image
    std::ostringstream filenameStream;
    filenameStream << "pngs/" << filenameBase << "_" << std::setfill('0') << std::setw(5) << nFrame << ".png";
    std::string filename = filenameStream.str();
    savePNG(SCREEN_WIDTH, SCREEN_HEIGHT, &filename[0]);
#endif

    // advance the frame counter
    nFrame++;

  }


  while( glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS ) {
    plotOpenGL();
    glfwPollEvents();
  }

  endOpenGL();

  return 0;
}


//  ######  ########    ###    ########  ########     #######  ########  ######## ##    ##  ######   ##
// ##    ##    ##      ## ##   ##     ##    ##       ##     ## ##     ## ##       ###   ## ##    ##  ##
// ##          ##     ##   ##  ##     ##    ##       ##     ## ##     ## ##       ####  ## ##        ##
//  ######     ##    ##     ## ########     ##       ##     ## ########  ######   ## ## ## ##   #### ##
//       ##    ##    ######### ##   ##      ##       ##     ## ##        ##       ##  #### ##    ##  ##
// ##    ##    ##    ##     ## ##    ##     ##       ##     ## ##        ##       ##   ### ##    ##  ##
//  ######     ##    ##     ## ##     ##    ##        #######  ##        ######## ##    ##  ######   ########

int startOpenGL() {
  // Start GL Window
  glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // We don't want the old OpenGL

  // Open a window and create its OpenGL context
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
	glClearColor(colorSolid.x(), colorSolid.y(), colorSolid.z(), 0.0f);

  // create VAO
  GLuint VertexArrayID;
  glGenVertexArrays(1, &VertexArrayID);
  glBindVertexArray(VertexArrayID);
  // load GL shaders
  programID = LoadShaders( "OpenGLvertex.vertexshader", "OpenGLfragment.fragmentshader" );
  // setup basic transformation matrices
  glViewport( 0.0f, 0.0f, SCREEN_WIDTH, SCREEN_HEIGHT );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity( );
  glOrtho( 0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0, 1 );
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity( );

  return 0;
}

void endOpenGL() {
  glfwTerminate();
}

void plotOpenGL() {
// =============================================================================
// ------------------------- PLOT GRID -----------------------------------------
// =============================================================================

  // setup vertex buffer data for the background grid
    // nGridCells * 2 triangles per cell * 3 vertices per triangle * _ number of components
  GLfloat g_vertex_buffer_data_grid[mySim->nGridCells()*2*3*2 + mySim->nGridCells()*2*3*3];
  int offsetGridVertices = 0;
  int offsetGridColors = mySim->nGridCells()*2*3*2;
  for (int i=0; i<mySim->gridSize().x(); i++) {
    for (int j=0; j<mySim->gridSize().y(); j++) {
      // get the flat index
      int idx = i + mySim->gridSize().x()*j;

      // positions
        // for a quad with vertices 0=(0,0); 1=(1,0); 2=(0,1); 3=(1,1)
        // we list vertices in order 0-1-2-1-2-3
      // make a matrix containing the vertex positions first
      MatrixXd gridVertexPositions = MatrixXd::Zero(4,2);
      gridVertexPositions.row(0) = Vector2d(i*mySim->gridSpacing().x(),j*mySim->gridSpacing().y());
      gridVertexPositions.row(1) = Vector2d((i+1)*mySim->gridSpacing().x(),j*mySim->gridSpacing().y());
      gridVertexPositions.row(2) = Vector2d(i*mySim->gridSpacing().x(),(j+1)*mySim->gridSpacing().y());
      gridVertexPositions.row(3) = Vector2d((i+1)*mySim->gridSpacing().x(),(j+1)*mySim->gridSpacing().y());

      gridVertexPositions = 2.0*gridVertexPositions - 1.0*MatrixXd::Ones(4,2);

      // now fill in the buffer by throwing the vertex coordinates in there
      g_vertex_buffer_data_grid[offsetGridVertices + 2*3*2*idx] = gridVertexPositions(0,0);
      g_vertex_buffer_data_grid[offsetGridVertices + 2*3*2*idx+1] = gridVertexPositions(0,1);
      g_vertex_buffer_data_grid[offsetGridVertices + 2*3*2*idx+2] = gridVertexPositions(1,0);
      g_vertex_buffer_data_grid[offsetGridVertices + 2*3*2*idx+3] = gridVertexPositions(1,1);
      g_vertex_buffer_data_grid[offsetGridVertices + 2*3*2*idx+4] = gridVertexPositions(2,0);
      g_vertex_buffer_data_grid[offsetGridVertices + 2*3*2*idx+5] = gridVertexPositions(2,1);
      g_vertex_buffer_data_grid[offsetGridVertices + 2*3*2*idx+6] = gridVertexPositions(1,0);
      g_vertex_buffer_data_grid[offsetGridVertices + 2*3*2*idx+7] = gridVertexPositions(1,1);
      g_vertex_buffer_data_grid[offsetGridVertices + 2*3*2*idx+8] = gridVertexPositions(2,0);
      g_vertex_buffer_data_grid[offsetGridVertices + 2*3*2*idx+9] = gridVertexPositions(2,1);
      g_vertex_buffer_data_grid[offsetGridVertices + 2*3*2*idx+10] = gridVertexPositions(3,0);
      g_vertex_buffer_data_grid[offsetGridVertices + 2*3*2*idx+11] = gridVertexPositions(3,1);

      // figure out what color it should be be
      Material cellMaterial = mySim->materialAtCell(i,j);
      Vector3d cellColor;

      if (cellMaterial == Material::solid) cellColor = colorSolid;
      else if (cellMaterial == Material::fluid) cellColor = colorFluid;
      else cellColor = colorEmpty;

      // now fill in the buffer by throwing the vertex color in there
      for (int k=0; k<6; k++) {
        g_vertex_buffer_data_grid[offsetGridColors + 2*3*3*idx + 3*k] = cellColor(0);
        g_vertex_buffer_data_grid[offsetGridColors + 2*3*3*idx + 3*k+1] = cellColor(1);
        g_vertex_buffer_data_grid[offsetGridColors + 2*3*3*idx + 3*k+2] = cellColor(2);
      }
    }
  }
  // bind vertex buffer
	GLuint vertexBufferGrid;
	glGenBuffers(1, &vertexBufferGrid);
	glBindBuffer(GL_ARRAY_BUFFER, vertexBufferGrid);
  glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data_grid), g_vertex_buffer_data_grid, GL_DYNAMIC_DRAW);

  // 1st attribute buffer : vertex positions
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, vertexBufferGrid);
  glVertexAttribPointer(
     0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
     2,                  // size
     GL_FLOAT,           // type
     GL_FALSE,           // normalized?
     0,                  // stride
     (void*)0            // array buffer offset
  );

  // 2nd attribute buffer : vertex colors
  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, vertexBufferGrid);
  glVertexAttribPointer(
    1, 3, GL_FLOAT, GL_FALSE, 0, (void*)(sizeof(GL_FLOAT)*offsetGridColors)
  );

  // Draw Points
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPointSize(10.0);
  glEnable(GL_BLEND);
  glUseProgram(programID);
  glDrawArrays(GL_TRIANGLES, 0, mySim->nGridCells()*2*3);
  glDisable(GL_BLEND);
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);


// =============================================================================
// ------------------------- PLOT PARTICLES ------------------------------------
// =============================================================================

  // setup vertex buffer data
  GLfloat g_vertex_buffer_data[mySim->nParticles()*2 + mySim->nParticles()*3];
  int offsetPositions = 0;
  int offsetColors = mySim->nParticles()*2;
  double maxParticleSpeed = mySim->maxParticleSpeed() + D_EPSILON;
  //if (maxParticleSpeed < 1.0) maxParticleSpeed = 1.0;
  for (int i=0; i<mySim->nParticles(); i++) {
    // positions
    g_vertex_buffer_data[offsetPositions + 2*i] = ((GLfloat)mySim->particlePosition(i).x()*2.0f) - 1.0f;
    g_vertex_buffer_data[offsetPositions + 2*i+1] = ((GLfloat)mySim->particlePosition(i).y()*2.0f) - 1.0f;

    // colors
    g_vertex_buffer_data[offsetColors + 3*i] = ((GLfloat)(mySim->particleSpeed(i)/maxParticleSpeed));
    g_vertex_buffer_data[offsetColors + 3*i+1] = 1.0f;
    g_vertex_buffer_data[offsetColors + 3*i+2] = 1.0f;
  }
  // bind vertex buffer
	GLuint vertexbuffer;
	glGenBuffers(1, &vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
  glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_DYNAMIC_DRAW);

  // 1st attribute buffer : vertex positions
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
  glVertexAttribPointer(
     0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
     2,                  // size
     GL_FLOAT,           // type
     GL_FALSE,           // normalized?
     0,                  // stride
     (void*)0            // array buffer offset
  );

  // 2nd attribute buffer : vertex colors
  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
  glVertexAttribPointer(
    1, 3, GL_FLOAT, GL_FALSE, 0, (void*)(sizeof(GL_FLOAT)*offsetColors)
  );

  // Draw Points
  //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPointSize(10.0);
  glEnable(GL_BLEND);
  glUseProgram(programID);
  glDrawArrays(GL_POINTS, 0, mySim->nParticles());
  glDisable(GL_BLEND);
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);

	// Swap buffers
	glfwSwapBuffers(window);
}


void savePNG(unsigned int height, unsigned int width, char* filename) {
  GLfloat* OpenGLimage = new GLfloat[height*width*3];
  glReadPixels(0, 0, width, height, GL_RGB, GL_FLOAT, OpenGLimage);
  pngwriter PNG(width, height, 1.0, filename);

  size_t x = 1;   // start the top and leftmost point of the window
  size_t y = 1;
  double R, G, B;
  for(size_t i=0; i<height*width*3; i++)
  {
        switch(i%3) //the OpenGLimage array look like [R1, G1, B1, R2, G2, B2,...]
       {
             case 2:
                   B = (double) OpenGLimage[i]; break;
             case 1:
                   G = (double) OpenGLimage[i]; break;
             case 0:
                   R = (double) OpenGLimage[i];
                   PNG.plot(x, y, R, G, B);
                   if( x == width )
                   {
                         x=1;
                         y++;
                    }
                    else
                    { x++; }
                    break;
       }
  }
  PNG.close();

  // cleanup
  delete OpenGLimage;
}
