// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/* Original mesh
 *
 * Cells are 0-3, vertices are 4-8.
 *
 * 4   5,6,7  8
 *
 *     ^^^^^ Face in x-y plane
 *
 * After adding cohesive elements
 *
 * Cells are 0-3,18,19, vertices are 4-17.
 *
 * 4   5,6,7  10,11,12   8
 *            13,14,15
 *     ^^^^^^^^^^^^ Cohesive element in x-y plane.
 */

#include "CohesiveKinDataTet4e.hh"

const char* pylith::faults::CohesiveKinDataTet4e::_meshFilename =
  "data/tet4e.mesh";

const int pylith::faults::CohesiveKinDataTet4e::_spaceDim = 3;

const int pylith::faults::CohesiveKinDataTet4e::_cellDim = 2;

const int pylith::faults::CohesiveKinDataTet4e::_numBasis = 3;

const int pylith::faults::CohesiveKinDataTet4e::_numQuadPts = 1;

const double pylith::faults::CohesiveKinDataTet4e::_quadPts[] = {
  -3.33333333e-01,  -3.33333333e-01,
};

const double pylith::faults::CohesiveKinDataTet4e::_quadWts[] = {
  2.0,
};

const double pylith::faults::CohesiveKinDataTet4e::_basis[] = {
  3.33333333e-01,  3.33333333e-01,
  3.33333333e-01,};

const double pylith::faults::CohesiveKinDataTet4e::_basisDeriv[] = {
 -0.50000000e+00, -0.50000000e+00,
  0.50000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.50000000e+00,
};

const double pylith::faults::CohesiveKinDataTet4e::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00,
};

const int pylith::faults::CohesiveKinDataTet4e::_id = 10;

const char* pylith::faults::CohesiveKinDataTet4e::_label = "fault";

const char* pylith::faults::CohesiveKinDataTet4e::_finalSlipFilename = 
  "data/tet4e_finalslip.spatialdb";

const char* pylith::faults::CohesiveKinDataTet4e::_slipTimeFilename = 
  "data/tet4e_sliptime.spatialdb";

const char* pylith::faults::CohesiveKinDataTet4e::_riseTimeFilename = 
  "data/tet4e_risetime.spatialdb";

const double pylith::faults::CohesiveKinDataTet4e::_fieldT[] = {
  3.1, 5.1, 7.1,
  3.2, 5.2, 7.2, // 5
  3.3, 5.3, 7.3, // 6
  3.4, 5.4, 7.4, // 7
  3.5, 5.5, 7.5, // 8
  3.6, 5.6, 7.6,
  3.7, 5.7, 7.7, // 10
  3.9, 5.9, 7.9, // 11
  4.1, 6.1, 8.1, // 12
  4.3, 6.3, 8.3, // 13
  3.8, 5.8, 7.8, // 14
  3.0, 5.0, 7.0, // 15
  4.2, 6.2, 8.2, // 16
  4.4, 6.4, 8.4, // 17
};

const double pylith::faults::CohesiveKinDataTet4e::_fieldIncr[] = {
  6.1, 7.1, 2.1,
  6.2, 7.2, 2.2, // 5
  6.3, 7.3, 2.3, // 6
  6.4, 7.4, 2.4, // 7
  6.5, 7.5, 2.5, // 8
  6.6, 7.6, 2.6,
  6.7, 7.7, 2.7, // 10
  6.9, 7.9, 2.9, // 11
  5.1, 8.1, 1.1, // 12
  5.3, 8.3, 1.3, // 13
  4.8, 8.8, 2.8, // 14
  4.0, 8.0, 2.0, // 15
  5.2, 8.2, 1.2, // 16
  5.4, 8.4, 1.4, // 17
};

const double pylith::faults::CohesiveKinDataTet4e::_jacobianLumped[] = {
  4.1, 6.1, 8.1,
  4.2, 6.2, 8.2, // 5
  4.3, 6.3, 8.3, // 6
  4.4, 6.4, 8.4, // 7
  4.5, 6.5, 8.5, // 8
  4.6, 6.6, 8.6,
  4.7, 6.7, 8.7, // 10
  4.9, 6.9, 8.9, // 11
  3.1, 5.1, 7.1, // 12
  3.3, 5.3, 7.3, // 13
  2.0/3.0, 2.0/3.0, 2.0/3.0, // 14
  1.0/3.0, 1.0/3.0, 1.0/3.0, // 15
  2.0/3.0, 2.0/3.0, 2.0/3.0, // 16
  1.0/3.0, 1.0/3.0, 1.0/3.0, // 17
};


const double pylith::faults::CohesiveKinDataTet4e::_orientation[] = {
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
  0.0, +1.0, 0.0,    0.0, 0.0, +1.0,    +1.0, 0.0, 0.0,
};

const double pylith::faults::CohesiveKinDataTet4e::_area[] = {
  2.0/3.0,
  1.0/3.0,
  2.0/3.0,
  1.0/3.0,
};

const int pylith::faults::CohesiveKinDataTet4e::_numFaultVertices = 4;
const int pylith::faults::CohesiveKinDataTet4e::_verticesFault[] = {
  3, 2, 4, 5
};
const int pylith::faults::CohesiveKinDataTet4e::_verticesLagrange[] = {
  15, 14, 16, 17
};
const int pylith::faults::CohesiveKinDataTet4e::_verticesNegative[] = {
  6, 5, 7, 8
};
const int pylith::faults::CohesiveKinDataTet4e::_verticesPositive[] = {
  11, 10, 12, 13
};

const int pylith::faults::CohesiveKinDataTet4e::_numCohesiveCells = 2;
const int pylith::faults::CohesiveKinDataTet4e::_cellMappingFault[] = {
  0, 1
};
const int pylith::faults::CohesiveKinDataTet4e::_cellMappingCohesive[] = {
  18, 19
};



const double pylith::faults::CohesiveKinDataTet4e::_residual[] = {
  0.0,  0.0,  0.0,
  -(3.8+3.0+4.2 + 3.8+4.2+4.4)/9.0,
  -(5.8+5.0+6.2 + 5.8+6.2+6.4)/9.0,
  -(7.8+7.0+8.2 + 7.8+8.2+8.4)/9.0, // 5
  -(3.8+3.0+4.2)/9.0,  -(5.8+5.0+6.2)/9.0,  -(7.8+7.0+8.2)/9.0, // 6
  -(3.8+3.0+4.2 + 3.8+4.2+4.4)/9.0,
  -(5.8+5.0+6.2 + 5.8+6.2+6.4)/9.0,
  -(7.8+7.0+8.2 + 7.8+8.2+8.4)/9.0, // 5
  -(3.8+4.2+4.4)/9.0,  -(5.8+6.2+6.4)/9.0,  -(7.8+8.2+8.4)/9.0, // 8
  0.0,  0.0,  0.0,
  +(3.8+3.0+4.2 + 3.8+4.2+4.4)/9.0,
  +(5.8+5.0+6.2 + 5.8+6.2+6.4)/9.0,
  +(7.8+7.0+8.2 + 7.8+8.2+8.4)/9.0, // 10
  +(3.8+3.0+4.2)/9.0,  +(5.8+5.0+6.2)/9.0,  +(7.8+7.0+8.2)/9.0, // 11
  +(3.8+3.0+4.2 + 3.8+4.2+4.4)/9.0,
  +(5.8+5.0+6.2 + 5.8+6.2+6.4)/9.0,
  +(7.8+7.0+8.2 + 7.8+8.2+8.4)/9.0, // 12
  +(3.8+4.2+4.4)/9.0,  +(5.8+6.2+6.4)/9.0,  +(7.8+8.2+8.4)/9.0, // 13

  (3.7-3.2 + 3.9-3.3 + 4.1-3.4)/9.0 +
  (3.7-3.2 + 4.1-3.4 + 4.3-3.5)/9.0 +
  -(0.07938069066 + 0.14140241667 + 0.18205179147)/9.0 +
  -(0.07938069066 + 0.18205179147 + 0.19904410828)/9.0,
  (5.7-5.2 + 5.9-5.3 + 6.1-5.4)/9.0 +
  (5.7-5.2 + 6.1-5.4 + 6.3-5.5)/9.0 +
  -(1.82575588523 + 1.69682900001 + 1.51709826228)/9.0 +
  -(1.82575588523 + 1.51709826228 + 1.29378670385)/9.0,
  (7.7-7.2 + 7.9-7.3 + 8.1-7.4)/9.0 +
  (7.7-7.2 + 8.1-7.4 + 8.3-7.5)/9.0 +
  -(-0.55566483464 + -0.56560966667 + -0.54615537442)/9.0 +
  -(-0.55566483464 + -0.54615537442 + -0.49761027071)/9.0, // 14

  (3.7-3.2 + 3.9-3.3 + 4.1-3.4)/9.0 +
  -(0.07938069066 + 0.14140241667 + 0.18205179147)/9.0,
  (5.7-5.2 + 5.9-5.3 + 6.1-5.4)/9.0 +
  -(1.82575588523 + 1.69682900001 + 1.51709826228)/9.0,
  (7.7-7.2 + 7.9-7.3 + 8.1-7.4)/9.0 +
  -(-0.55566483464 + -0.56560966667 + -0.54615537442)/9.0, // 15

  (3.7-3.2 + 3.9-3.3 + 4.1-3.4)/9.0 +
  (3.7-3.2 + 4.1-3.4 + 4.3-3.5)/9.0 +
  -(0.07938069066 + 0.14140241667 + 0.18205179147)/9.0 +
  -(0.07938069066 + 0.18205179147 + 0.19904410828)/9.0,
  (5.7-5.2 + 5.9-5.3 + 6.1-5.4)/9.0 +
  (5.7-5.2 + 6.1-5.4 + 6.3-5.5)/9.0 +
  -(1.82575588523 + 1.69682900001 + 1.51709826228)/9.0 +
  -(1.82575588523 + 1.51709826228 + 1.29378670385)/9.0,
  (7.7-7.2 + 7.9-7.3 + 8.1-7.4)/9.0 +
  (7.7-7.2 + 8.1-7.4 + 8.3-7.5)/9.0 +
  -(-0.55566483464 + -0.56560966667 + -0.54615537442)/9.0 +
  -(-0.55566483464 + -0.54615537442 + -0.49761027071)/9.0, // 16

  (3.7-3.2 + 4.1-3.4 + 4.3-3.5)/9.0 +
  -(0.07938069066 + 0.18205179147 + 0.19904410828)/9.0,
  (5.7-5.2 + 6.1-5.4 + 6.3-5.5)/9.0 +
  -(1.82575588523 + 1.51709826228 + 1.29378670385)/9.0,
  (7.7-7.2 + 8.1-7.4 + 8.3-7.5)/9.0 +
  -(-0.55566483464 + -0.54615537442 + -0.49761027071)/9.0, // 17

};

const double pylith::faults::CohesiveKinDataTet4e::_residualIncr[] = {
  0.0,  0.0,  0.0,
  -(3.8+3.0+4.2 + 3.8+4.2+4.4)/9.0,
  -(5.8+5.0+6.2 + 5.8+6.2+6.4)/9.0,
  -(7.8+7.0+8.2 + 7.8+8.2+8.4)/9.0, // 5
  -(3.8+3.0+4.2)/9.0,  -(5.8+5.0+6.2)/9.0,  -(7.8+7.0+8.2)/9.0, // 6
  -(3.8+3.0+4.2 + 3.8+4.2+4.4)/9.0,
  -(5.8+5.0+6.2 + 5.8+6.2+6.4)/9.0,
  -(7.8+7.0+8.2 + 7.8+8.2+8.4)/9.0, // 5
  -(3.8+4.2+4.4)/9.0,  -(5.8+6.2+6.4)/9.0,  -(7.8+8.2+8.4)/9.0, // 8
  0.0,  0.0,  0.0,
  +(3.8+3.0+4.2 + 3.8+4.2+4.4)/9.0,
  +(5.8+5.0+6.2 + 5.8+6.2+6.4)/9.0,
  +(7.8+7.0+8.2 + 7.8+8.2+8.4)/9.0, // 10
  +(3.8+3.0+4.2)/9.0,  +(5.8+5.0+6.2)/9.0,  +(7.8+7.0+8.2)/9.0, // 11
  +(3.8+3.0+4.2 + 3.8+4.2+4.4)/9.0,
  +(5.8+5.0+6.2 + 5.8+6.2+6.4)/9.0,
  +(7.8+7.0+8.2 + 7.8+8.2+8.4)/9.0, // 12
  +(3.8+4.2+4.4)/9.0,  +(5.8+6.2+6.4)/9.0,  +(7.8+8.2+8.4)/9.0, // 13

  (3.7-3.2 + 3.9-3.3 + 4.1-3.4)/9.0 +
  (3.7-3.2 + 4.1-3.4 + 4.3-3.5)/9.0 +
  -(0.07938069066 + 0.14140241667 + 0.18205179147)/9.0 +
  -(0.07938069066 + 0.18205179147 + 0.19904410828)/9.0,
  (5.7-5.2 + 5.9-5.3 + 6.1-5.4)/9.0 +
  (5.7-5.2 + 6.1-5.4 + 6.3-5.5)/9.0 +
  -(1.82575588523 + 1.69682900001 + 1.51709826228)/9.0 +
  -(1.82575588523 + 1.51709826228 + 1.29378670385)/9.0,
  (7.7-7.2 + 7.9-7.3 + 8.1-7.4)/9.0 +
  (7.7-7.2 + 8.1-7.4 + 8.3-7.5)/9.0 +
  -(-0.55566483464 + -0.56560966667 + -0.54615537442)/9.0 +
  -(-0.55566483464 + -0.54615537442 + -0.49761027071)/9.0, // 14

  (3.7-3.2 + 3.9-3.3 + 4.1-3.4)/9.0 +
  -(0.07938069066 + 0.14140241667 + 0.18205179147)/9.0,
  (5.7-5.2 + 5.9-5.3 + 6.1-5.4)/9.0 +
  -(1.82575588523 + 1.69682900001 + 1.51709826228)/9.0,
  (7.7-7.2 + 7.9-7.3 + 8.1-7.4)/9.0 +
  -(-0.55566483464 + -0.56560966667 + -0.54615537442)/9.0, // 15

  (3.7-3.2 + 3.9-3.3 + 4.1-3.4)/9.0 +
  (3.7-3.2 + 4.1-3.4 + 4.3-3.5)/9.0 +
  -(0.07938069066 + 0.14140241667 + 0.18205179147)/9.0 +
  -(0.07938069066 + 0.18205179147 + 0.19904410828)/9.0,
  (5.7-5.2 + 5.9-5.3 + 6.1-5.4)/9.0 +
  (5.7-5.2 + 6.1-5.4 + 6.3-5.5)/9.0 +
  -(1.82575588523 + 1.69682900001 + 1.51709826228)/9.0 +
  -(1.82575588523 + 1.51709826228 + 1.29378670385)/9.0,
  (7.7-7.2 + 7.9-7.3 + 8.1-7.4)/9.0 +
  (7.7-7.2 + 8.1-7.4 + 8.3-7.5)/9.0 +
  -(-0.55566483464 + -0.56560966667 + -0.54615537442)/9.0 +
  -(-0.55566483464 + -0.54615537442 + -0.49761027071)/9.0, // 16

  (3.7-3.2 + 4.1-3.4 + 4.3-3.5)/9.0 +
  -(0.07938069066 + 0.18205179147 + 0.19904410828)/9.0,
  (5.7-5.2 + 6.1-5.4 + 6.3-5.5)/9.0 +
  -(1.82575588523 + 1.51709826228 + 1.29378670385)/9.0,
  (7.7-7.2 + 8.1-7.4 + 8.3-7.5)/9.0 +
  -(-0.55566483464 + -0.54615537442 + -0.49761027071)/9.0, // 17

};

const double pylith::faults::CohesiveKinDataTet4e::_jacobian[] = {
  0.0, 0.0, 0.0, // 4x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 4y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 4z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 5x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 +2.0/3.0, 0.0, 0.0, // 14
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 5y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,+2.0/3.0, 0.0, // 14
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 5z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,+2.0/3.0, // 14
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 6x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 +1.0/3.0, 0.0, 0.0, // 15
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 6y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, +1.0/3.0, 0.0, // 15
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 6z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, +1.0/3.0, // 15
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 7x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  +2.0/3.0, 0.0, 0.0, // 16
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 7y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,+2.0/3.0, 0.0, // 16
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 7z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,+2.0/3.0, // 16
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 8x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 +1.0/3.0, 0.0, 0.0, // 17
  0.0, 0.0, 0.0, // 8y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,+1.0/3.0, 0.0, // 17
  0.0, 0.0, 0.0, // 8z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,+1.0/3.0, // 17
  0.0, 0.0, 0.0, // 9x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 9y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 9z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 10x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 -2.0/3.0, 0.0, 0.0, // 14
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 10y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,-2.0/3.0, 0.0, // 14
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 10z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,-2.0/3.0, // 14
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 11x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 -1.0/3.0, 0.0, 0.0, // 15
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 11y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,-1.0/3.0, 0.0, // 15
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 11z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,-1.0/3.0, // 15
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 12x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 -2.0/3.0, 0.0, 0.0, // 16
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 12y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,-2.0/3.0, 0.0, // 16
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 12z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,-2.0/3.0, // 16
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 13x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 -1.0/3.0, 0.0, 0.0, // 17
  0.0, 0.0, 0.0, // 13y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,-1.0/3.0, 0.0, // 17
  0.0, 0.0, 0.0, // 13z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,-1.0/3.0, // 17
  0.0, 0.0, 0.0, // 14x
 +2.0/3.0, 0.0, 0.0, // 5
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 -2.0/3.0, 0.0, 0.0, // 10
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 14y
  0.0,+2.0/3.0, 0.0, // 5
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,-2.0/3.0, 0.0, // 10
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 14z
  0.0, 0.0,+2.0/3.0, // 5
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,-2.0/3.0, // 10
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 15x
  0.0, 0.0, 0.0,
 +1.0/3.0, 0.0, 0.0, // 6
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 -1.0/3.0, 0.0, 0.0, // 11
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 15y
  0.0, 0.0, 0.0,
  0.0,+1.0/3.0, 0.0, // 6
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,-1.0/3.0, 0.0, // 11
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 15z
  0.0, 0.0, 0.0,
  0.0, 0.0,+1.0/3.0, // 6
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,-1.0/3.0, // 11
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 16x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 +2.0/3.0, 0.0, 0.0, // 7
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 -2.0/3.0, 0.0, 0.0, // 12
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 16y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,+2.0/3.0, 0.0, // 7
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,-2.0/3.0, 0.0, // 12
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 16z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,+2.0/3.0, // 7
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,-2.0/3.0, // 12
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 17x
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 +1.0/3.0, 0.0, 0.0, // 8
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
 -1.0/3.0, 0.0, 0.0, // 13
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 17y
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,+1.0/3.0, 0.0, // 8
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0,-1.0/3.0, 0.0, // 13
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, // 17z
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,+1.0/3.0, // 8
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0,-1.0/3.0, // 13
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
  0.0, 0.0, 0.0,
};

pylith::faults::CohesiveKinDataTet4e::CohesiveKinDataTet4e(void)
{ // constructor
  meshFilename = const_cast<char*>(_meshFilename);
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  quadPts = const_cast<double*>(_quadPts);
  quadWts = const_cast<double*>(_quadWts);
  basis = const_cast<double*>(_basis);
  basisDeriv = const_cast<double*>(_basisDeriv);
  verticesRef = const_cast<double*>(_verticesRef);
  id = _id;
  label = const_cast<char*>(_label);
  finalSlipFilename = const_cast<char*>(_finalSlipFilename);
  slipTimeFilename = const_cast<char*>(_slipTimeFilename);
  riseTimeFilename = const_cast<char*>(_riseTimeFilename);
  fieldT = const_cast<double*>(_fieldT);
  fieldIncr = const_cast<double*>(_fieldIncr);
  jacobianLumped = const_cast<double*>(_jacobianLumped);
  orientation = const_cast<double*>(_orientation);
  area = const_cast<double*>(_area);
  residual = const_cast<double*>(_residual);
  residualIncr = const_cast<double*>(_residualIncr);
  jacobian = const_cast<double*>(_jacobian);
  verticesFault = const_cast<int*>(_verticesFault);
  verticesLagrange = const_cast<int*>(_verticesLagrange);
  verticesNegative = const_cast<int*>(_verticesNegative);
  verticesPositive = const_cast<int*>(_verticesPositive);
  numFaultVertices = _numFaultVertices;  
  cellMappingFault = const_cast<int*>(_cellMappingFault);
  cellMappingCohesive = const_cast<int*>(_cellMappingCohesive);
  numCohesiveCells = _numCohesiveCells;  
} // constructor

pylith::faults::CohesiveKinDataTet4e::~CohesiveKinDataTet4e(void)
{}


// End of file
