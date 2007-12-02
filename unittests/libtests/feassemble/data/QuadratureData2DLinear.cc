// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application quadratureapp.

#include "QuadratureData2DLinear.hh"

const int pylith::feassemble::QuadratureData2DLinear::_numVertices = 3;

const int pylith::feassemble::QuadratureData2DLinear::_spaceDim = 2;

const int pylith::feassemble::QuadratureData2DLinear::_numCells = 1;

const int pylith::feassemble::QuadratureData2DLinear::_cellDim = 2;

const int pylith::feassemble::QuadratureData2DLinear::_numBasis = 3;

const int pylith::feassemble::QuadratureData2DLinear::_numQuadPts = 1;

const double pylith::feassemble::QuadratureData2DLinear::_vertices[] = {
  2.00000000e-01, -4.00000000e-01,
  3.00000000e-01,  5.00000000e-01,
 -1.00000000e+00, -2.00000000e-01,
};

const int pylith::feassemble::QuadratureData2DLinear::_cells[] = {
       0,       1,       2,
};

const double pylith::feassemble::QuadratureData2DLinear::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
};

const double pylith::feassemble::QuadratureData2DLinear::_quadPtsRef[] = {
  3.33333333e-01,  3.33333333e-01,
};

const double pylith::feassemble::QuadratureData2DLinear::_quadWts[] = {
  5.00000000e-01,
};

const double pylith::feassemble::QuadratureData2DLinear::_quadPts[] = {
 -5.33333333e-01,  3.33333333e-01,
};

const double pylith::feassemble::QuadratureData2DLinear::_basis[] = {
 -3.33333333e-01,  6.66666667e-01,
  6.66666667e-01,};

const double pylith::feassemble::QuadratureData2DLinear::_basisDerivRef[] = {
 -5.00000000e-01, -5.00000000e-01,
  5.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  5.00000000e-01,
};

const double pylith::feassemble::QuadratureData2DLinear::_basisDeriv[] = {
  6.36363636e-01, -1.18181818e+00,
  1.81818182e-01,  1.09090909e+00,
 -8.18181818e-01,  9.09090909e-02,
};

const double pylith::feassemble::QuadratureData2DLinear::_jacobian[] = {
  5.00000000e-02, -6.00000000e-01,
  4.50000000e-01,  1.00000000e-01,
};

const double pylith::feassemble::QuadratureData2DLinear::_jacobianDet[] = {
  2.75000000e-01,
};

const double pylith::feassemble::QuadratureData2DLinear::_jacobianInv[] = {
  3.63636364e-01,  2.18181818e+00,
 -1.63636364e+00,  1.81818182e-01,
};

pylith::feassemble::QuadratureData2DLinear::QuadratureData2DLinear(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  vertices = const_cast<double*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<double*>(_verticesRef);
  quadPtsRef = const_cast<double*>(_quadPtsRef);
  quadWts = const_cast<double*>(_quadWts);
  quadPts = const_cast<double*>(_quadPts);
  basis = const_cast<double*>(_basis);
  basisDerivRef = const_cast<double*>(_basisDerivRef);
  basisDeriv = const_cast<double*>(_basisDeriv);
  jacobian = const_cast<double*>(_jacobian);
  jacobianDet = const_cast<double*>(_jacobianDet);
  jacobianInv = const_cast<double*>(_jacobianInv);
} // constructor

pylith::feassemble::QuadratureData2DLinear::~QuadratureData2DLinear(void)
{}


// End of file
