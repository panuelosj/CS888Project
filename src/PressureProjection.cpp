#include "PressureProjection.h"
#include <iostream>

//  ######   #######  ##    ##  ######  ######## ########
// ##    ## ##     ## ###   ## ##    ##    ##    ##     ##
// ##       ##     ## ####  ## ##          ##    ##     ##
// ##       ##     ## ## ## ##  ######     ##    ########
// ##       ##     ## ##  ####       ##    ##    ##   ##
// ##    ## ##     ## ##   ### ##    ##    ##    ##    ##
//  ######   #######  ##    ##  ######     ##    ##     ##

PressureProjection::PressureProjection() {

}

PressureProjection::PressureProjection(PressureProjectionInputs in) :
  _dt             ( in.dt ),
  _density        ( in.density ),
  _gridSize       ( in.gridSize ),
  _gridSpacing    ( in.gridSpacing ),
  _velocityField  ( in.velocityField ),
  _materialField  ( in.materialField ),
  _mapping        ( in.mapping ),
  _fluidCells     ( in.fluidCells ),
  _nFluidCells    ( _materialField->nFluid() ),
  // derived variables
  _invGridSpacing( Vector2d(1.0/_gridSpacing(0), 1.0/_gridSpacing(1)) ),
  _gridVectorLength( _gridSize(0) * _gridSize(1) ),
  // constructed variables
  _b( VectorXd(_nFluidCells) ),
  _p( VectorXd(_nFluidCells) ),
  _Acompressed( MatrixXi(_nFluidCells, 3) ),
  _A( SparseMatrix<double>(_nFluidCells,_nFluidCells) )
{

}

PressureProjection::~PressureProjection() {

}


// ########  ##     ## ########  ##       ####  ######
// ##     ## ##     ## ##     ## ##        ##  ##    ##
// ##     ## ##     ## ##     ## ##        ##  ##
// ########  ##     ## ########  ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##
// ##        ##     ## ##     ## ##        ##  ##    ##
// ##         #######  ########  ######## ####  ######

void PressureProjection::project() {
  // we will solve Ax=b
  _computeNegativeDivergence();       // setup b
  _computeMatrixCoefficients();       // setup A
  _translateMatrixToSparseEigen();    // need to convert Bridson's datastructure
                                        // into a sparse Eigen matrix so we can
                                        // use Eigen's Cholesky factorization
  _solve();                           // solve Ap = b
  _pressureGradientUpdate();          // use the computed pressure to update
                                        // the velocities to give a divergence
                                        // free field.
  _velocityField->extrapolateVelocities();  // extrapolate velocities to the unknown
                                              // regions
}


// ########  ########  #### ##     ##    ###    ######## ########
// ##     ## ##     ##  ##  ##     ##   ## ##      ##    ##
// ##     ## ##     ##  ##  ##     ##  ##   ##     ##    ##
// ########  ########   ##  ##     ## ##     ##    ##    ######
// ##        ##   ##    ##   ##   ##  #########    ##    ##
// ##        ##    ##   ##    ## ##   ##     ##    ##    ##
// ##        ##     ## ####    ###    ##     ##    ##    ########

void PressureProjection::_computeNegativeDivergence() {
  // Function is primarily based on Bridson p72 Fig 5.3 and p76 Fig 5.4
  // Constructs the right hand side vector for the linear system to
  // find the pressure correction that gives a divergence free velocity field.

  // clear out _b vector
  for (unsigned int idx=0; idx < _fluidCells->size(); idx++) {
    _b[idx] = 0.0;
  }

  for (unsigned int idx=0; idx < _fluidCells->size(); idx++) {
    // centered divergence estimate
    // Refer to Bridson p72 Fig 5.3
    int i = _fluidCells->index(idx).x();
    int j = _fluidCells->index(idx).y();

    //int v = _mapping->gridToVector(i, j);
    int v = idx;

    double rhs = 0.0;

    rhs = -_invGridSpacing(0)*(double)(_velocityField->U(i+1,j) - _velocityField->U(i,j))
          -_invGridSpacing(1)*(double)(_velocityField->V(i,j+1) - _velocityField->V(i,j));

    // modification for solid velocity
    // Refer to Bridson p76 Fig 5.4
    // only have static solid boundaries for now (use 0 boundary velocity)
    if ( _materialField->isSolid(i-1, j) ) {
      rhs -= _invGridSpacing(0)*(_velocityField->U(i, j) - _usolid);
    }
    if ( _materialField->isSolid(i+1, j) ) {
      rhs += _invGridSpacing(0)*(_velocityField->U(i+1, j) - _usolid);
    }

    if ( _materialField->isSolid(i, j-1) ) {
      rhs -= _invGridSpacing(1)*(_velocityField->V(i, j) - _vsolid);
    }
    if ( _materialField->isSolid(i, j+1) ) {
      rhs += _invGridSpacing(1)*(_velocityField->V(i, j+1) - _vsolid);
    }

    _b[v] = rhs;
  }
}

void PressureProjection::_computeMatrixCoefficients() {
  // Function is based primarily on Bridson p78 Fig 5.5.
  // Constructs the matrix A for the left hand side of the linear system to
  // find the pressure correction that gives a divergence free velocity field.
  // Note that we use Bridson's compressed A matrix datastructure here, which
  // we will later translate into Eigen's sparse matrix representation.

  // Refer to Bridson p78 Fig 5.5
  for (unsigned int idx=0; idx < _fluidCells->size(); idx++) {
    int i = _fluidCells->index(idx).x();
    int j = _fluidCells->index(idx).y();

    //int v = _mapping->gridToVector(i, j);
    int v = idx;
    // make sure to clear out the values
    _Acompressed.row(v) = Vector3i(0,0,0);

    // MODIFIED BRIDSON
    // x neighbours
    if ( _materialField->isFluid(i-1,j) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
    } else if ( _materialField->isEmpty(i-1,j) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
    }
    if ( _materialField->isFluid(i+1,j) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
      _Acompressed(v,(int)AMatrixCell::Ax)--;
    } else if ( _materialField->isEmpty(i+1,j) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
    }

    // y neighbours
    if ( _materialField->isFluid(i,j-1) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
    } else if ( _materialField->isEmpty(i,j-1) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
    }
    if ( _materialField->isFluid(i,j+1) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
      _Acompressed(v,(int)AMatrixCell::Ay)--;
    } else if ( _materialField->isEmpty(i,j+1) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
    }

    /* BRIDSON P78
    // x neighbours
    if ( _materialField->isFluid(i-1,j) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
    }
    if ( _materialField->isFluid(i+1,j) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
      _Acompressed(v,(int)AMatrixCell::Ax)--;
    } else if ( _materialField->isEmpty(i+1,j) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
    }

    // y neighbours
    if ( _materialField->isFluid(i,j-1) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
    }
    if ( _materialField->isFluid(i,j+1) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
      _Acompressed(v,(int)AMatrixCell::Ay)--;
    } else if ( _materialField->isEmpty(i,j+1) ) {
      _Acompressed(v,(int)AMatrixCell::Adiag)++;
    } */
  }
}

void PressureProjection::_translateMatrixToSparseEigen() {
  // reserve space in the sparse matrix
  // see https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
  _A.reserve(_fluidCells->size()*5);

  for (unsigned int idx=0; idx < _fluidCells->size(); idx++) {
    int i = _fluidCells->index(idx).x();
    int j = _fluidCells->index(idx).y();

    //int v = _mapping->gridToVector(i, j);
    //int vDiag = _mapping->gridToVector(i, j);
    //int vX = _mapping->gridToVector(i+1, j);
    //int vY = _mapping->gridToVector(i, j+1);
    int v = idx;
    int vDiag = idx;
    int vX = _fluidCells->find(i+1, j);
    int vY = _fluidCells->find(i, j+1);

    // compute scale factor (see Bridson p78, Fig 5.5)
    double scale = _dt / (_density*_gridSpacing(0)*_gridSpacing(1));

    // Diag
    assert(vDiag != GRID_TO_VECTOR_OUTSIDE_RANGE_VAL);
    _A.insert(vDiag,vDiag) = scale*(double)_Acompressed(v,(int)AMatrixCell::Adiag);
    // Upper
    if (vX != GRID_TO_VECTOR_OUTSIDE_RANGE_VAL)
      _A.insert(vDiag,vX) = scale*(double)_Acompressed(v,(int)AMatrixCell::Ax);
    if (vY != GRID_TO_VECTOR_OUTSIDE_RANGE_VAL)
      _A.insert(vDiag,vY) = scale*(double)_Acompressed(v,(int)AMatrixCell::Ay);
    // Lower, remember this is symmetric
    if (vX != GRID_TO_VECTOR_OUTSIDE_RANGE_VAL)
      _A.insert(vX,vDiag) = scale*(double)_Acompressed(v,(int)AMatrixCell::Ax);
    if (vY != GRID_TO_VECTOR_OUTSIDE_RANGE_VAL)
      _A.insert(vY,vDiag) = scale*(double)_Acompressed(v,(int)AMatrixCell::Ay);
  }
}

void PressureProjection::_solve() {
  //std::cout << "A: " << std::endl << _A.transpose() << std::endl;
  //std::cout << "b: " << std::endl << _b << std::endl;

  // now we just need to solve Ap = b
  ConjugateGradient<SparseMatrix<double>, Upper|Lower, IncompleteCholesky<double, Upper|Lower>> cg;
  cg.setMaxIterations(1000);
  cg.compute(_A);
  if(cg.info()!=Success) {
    std::cout << "PressureProjection: Decomposition Failed!" << std::endl;
    return;
  }

  _p = cg.solve(_b);
  if(cg.info()!=Success) {
    std::cout << "PressureProjection: Solving Failed!" << std::endl;
    return;
  }

  //std::cout << "PressureField: " << std::endl << _p << std::endl;
}

void PressureProjection::_pressureGradientUpdate() {
  // Function is based primarily on Bridson p71 Fig 5.2.
  // We can now update the velocities using the computed pressures to get a
  // divergence free field
  double scaleX = _dt / (_density*_gridSpacing(0));
  double scaleY = _dt / (_density*_gridSpacing(1));

  // UPDATE Us
  for (unsigned int i=0; i<_gridSize(0)+1; i++) {
    for (unsigned int j=0; j<_gridSize(1); j++) {
      // since we aren't recasting the _p vector into a matrix, we need the
      // grid to vector mapping index
      //int v = _mapping->gridToVector(i, j);
      //int vX = _mapping->gridToVector(i-1, j);
      //int vY = _mapping->gridToVector(i, j-1);
      int v = _fluidCells->find(i,j);
      int vX = _fluidCells->find(i-1,j);

      double pV   = 0.0;
      double pVX  = 0.0;
      if (v != GRID_INDEX_NOT_FOUND_VAL) {
        pV = _p(v);
      }
      if (vX != GRID_INDEX_NOT_FOUND_VAL) {
        pVX = _p(vX);
      }

      // update u
      if (_materialField->isFluid(i-1, j) || _materialField->isFluid(i,j)) {
        // set the velocity sample to valid
        _velocityField->setUValid(i,j,1);
        if (_materialField->isSolid(i-1, j) || _materialField->isSolid(i,j))
          _velocityField->setU(i, j, _usolid);
        else
          _velocityField->subU( i, j, scaleX*(pV-pVX) );
      } else {
        // set the velocity sample to invalid
        _velocityField->setUValid(i,j,0);
        _velocityField->setU(i, j, _unknownVelocity);
      }
    }
  }

  // UPDATE Vs
  for (unsigned int i=0; i<_gridSize(0); i++) {
    for (unsigned int j=0; j<_gridSize(1)+1; j++) {
      // since we aren't recasting the _p vector into a matrix, we need the
      // grid to vector mapping index
      //int v = _mapping->gridToVector(i, j);
      //int vX = _mapping->gridToVector(i-1, j);
      //int vY = _mapping->gridToVector(i, j-1);
      int v = _fluidCells->find(i,j);
      int vY = _fluidCells->find(i,j-1);

      double pV   = 0.0;
      double pVY  = 0.0;
      if (v != GRID_INDEX_NOT_FOUND_VAL) {
        pV = _p(v);
      }
      if (vY != GRID_INDEX_NOT_FOUND_VAL) {
        pVY = _p(vY);
      }

      // update v
      if (_materialField->isFluid(i, j-1) || _materialField->isFluid(i,j)) {
        // set the velocity sample to valid
        _velocityField->setVValid(i,j,1);
        if (_materialField->isSolid(i, j-1) || _materialField->isSolid(i,j))
          _velocityField->setV(i, j, _vsolid);
        else
            _velocityField->subV( i, j, scaleY*(pV-pVY) );
      } else {
        // set the velocity sample to invalid
        _velocityField->setVValid(i,j,0);
        _velocityField->setV(i, j, _unknownVelocity);
      }
    }
  }
}
