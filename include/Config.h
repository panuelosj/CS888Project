#pragma once

//==============================================================================
//------------------- Logs -----------------------------------------------------
//==============================================================================
//#define LOG_QUIET
#define LOG_PARTICLE_ADVECTION
#define LOG_SAVE_PNG

//==============================================================================
//------------------- Particle-To-Grid Kernels ---------------------------------
//==============================================================================
// pick your poison
#define KERNEL_QUADRATIC_B_SPLINE

//==============================================================================
//------------------- Grid-To-Particle Kernels ---------------------------------
//==============================================================================
// pick your poison
#define INTERPOLATE_BILINEAR
//#define INTERPOLATE_BICUBIC   // unimplemented

// PIC/FLIP constant (0 = FLIP, 1 = PIC)
#define PIC_FLIP_ALPHA 0.1

//==============================================================================
//------------------- Particle Advection Scheme --------------------------------
//==============================================================================
// pick your poison
#define ADVECT_RALSTONRK3
//#define ADVECT_RK2            // unimplemented
//#define ADVECT_RK4            // unimplemented


#define D_DENSITY 1.0
#define D_JITTER_FACTOR 0.5

#define GRAVITY_X 0.0
#define GRAVITY_Y -9.8

#define SOLID_VELOCITY_X 0.0
#define SOLID_VELOCITY_Y 0.0
#define UNKNOWN_VELOCITY 0.0

#define VELOCITY_OUT_OF_RANGE_VALUE 0.0

//==============================================================================
//---------------------------- Epsilons ----------------------------------------
//==============================================================================
#define D_EPSILON 1.0e-10
