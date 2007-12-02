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

#include <portinfo>

#include "Quadrature1Din3D.hh" // implementation of class methods

#include "CellGeometry.hh" // USES CellGeometry

#include "petsc.h" // USES PetscLogFlopsNoCheck

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::Quadrature1Din3D::Quadrature1Din3D(void) : Quadrature()
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::Quadrature1Din3D::~Quadrature1Din3D(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::Quadrature1Din3D::Quadrature1Din3D(const Quadrature1Din3D& q) :
  Quadrature(q)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Compute geometric quantities for a cell at quadrature points.
void
pylith::feassemble::Quadrature1Din3D::computeGeometry(
		       const ALE::Obj<Mesh>& mesh,
		       const ALE::Obj<real_section_type>& coordinates,
		       const Mesh::point_type& cell)
{ // computeGeometry
  assert(1 == _cellDim);
  assert(3 == _spaceDim);

  _resetGeometry();

  // Get coordinates of cell's vertices
  const real_section_type::value_type* vertCoords = 
    mesh->restrict(coordinates, cell);
  assert(3 == coordinates->getFiberDimension(*mesh->depthStratum(0)->begin()));

  // Loop over quadrature points
  for (int iQuadPt=0; iQuadPt < _numQuadPts; ++iQuadPt) {
    
    // Compute coordinates of quadrature point in cell
    // x = sum[i=0,n-1] (Ni * xi)
    // y = sum[i=0,n-1] (Ni * yi)
    // z = sum[i=0,n-1] (Ni * zi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis) {
      const double basis = _basis[iQuadPt*_numBasis+iBasis];
      for (int iDim=0; iDim < _spaceDim; ++iDim)
	_quadPts[iQuadPt*_spaceDim+iDim] += 
	  basis * vertCoords[iBasis*_spaceDim+iDim];
    } // for
    
#if 0
    // Compute Jacobian at quadrature point
    // J = [dx/dp
    //      dy/dp
    //      dz/dp]
    // dx/dp = sum[i=0,n-1] (dNi/dp * xi)
    // dy/dp = sum[i=0,n-1] (dNi/dp * yi)
    // dz/dp = sum[i=0,n-1] (dNi/dp * zi)
    for (int iBasis=0; iBasis < _numBasis; ++iBasis) {
      const double deriv = _basisDerivRef[iQuadPt*_numBasis+iBasis];
      for (int iDim=0; iDim < _spaceDim; ++iDim)
	_jacobian[iQuadPt*_spaceDim+iDim] += 
	  deriv * vertCoords[iBasis*_spaceDim+iDim];
    } // for

    // Compute determinant of Jacobian at quadrature point
    // |J| = sqrt(transpose(J) J)
    double det = 0.0;
    for (int iDim=0; iDim < _spaceDim; ++iDim)
      det += _jacobian[iQuadPt*_spaceDim+iDim] * 
	_jacobian[iQuadPt*_spaceDim+iDim];
    det = sqrt(det);
    _checkJacobianDet(det, cell);
    _jacobianDet[iQuadPt] = det;
#else
    // Compute Jacobian and determinant of Jacobian at quadrature point
    assert(0 != _geometry);
    _geometry->jacobian(&_jacobian[iQuadPt*_cellDim*_spaceDim],
			&_jacobianDet[iQuadPt],
			vertCoords, &_quadPtsRef[iQuadPt*_cellDim], _spaceDim);
    _checkJacobianDet(_jacobianDet[iQuadPt], cell);
#endif

    // Compute inverse of Jacobian at quadrature point
    // Jinv = 1.0/[J]
    for (int iDim=0; iDim < _spaceDim; ++iDim)
      _jacobianInv[iQuadPt*_spaceDim+iDim] = 
	1.0 / _jacobian[iQuadPt*_spaceDim+iDim];

    // Compute derivatives of basis functions with respect to global
    // coordinates
    // dNi/dx = dNi/dp dp/dx + dNi/dq dq/dx + dNi/dr dr/dx
    for (int iBasis=0; iBasis < _numBasis; ++iBasis)
      for (int iDim=0; iDim < _spaceDim; ++iDim)
	for (int jDim=0; jDim < _cellDim; ++jDim)
	  _basisDeriv[iQuadPt*_numBasis*_spaceDim+iBasis*_spaceDim+iDim] +=
	    _basisDerivRef[iQuadPt*_numBasis*_cellDim+iBasis*_cellDim+jDim] *
	    _jacobianInv[iQuadPt*_cellDim*_spaceDim+jDim*_spaceDim+iDim];
  } // for

  PetscLogFlopsNoCheck(_numQuadPts * (1+_numBasis*_spaceDim*2 +
				      _spaceDim*1 +
				      _numBasis*_spaceDim*_cellDim*2));

} // computeGeometry


// End of file 
