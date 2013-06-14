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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application powerlawplanestrainelastic.

#include "PowerLawPlaneStrainElasticData.hh"

const int pylith::materials::PowerLawPlaneStrainElasticData::_dimension = 2;

const int pylith::materials::PowerLawPlaneStrainElasticData::_numLocs = 2;

const int pylith::materials::PowerLawPlaneStrainElasticData::_numProperties = 6;

const int pylith::materials::PowerLawPlaneStrainElasticData::_numStateVars = 3;

const int pylith::materials::PowerLawPlaneStrainElasticData::_numDBProperties = 6;

const int pylith::materials::PowerLawPlaneStrainElasticData::_numDBStateVars = 9;

const int pylith::materials::PowerLawPlaneStrainElasticData::_numPropsQuadPt = 6;

const int pylith::materials::PowerLawPlaneStrainElasticData::_numVarsQuadPt = 9;

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_lengthScale =   1.00000000e+03;

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_timeScale =   1.00000000e+00;

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_pressureScale =   2.25000000e+10;

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_densityScale =   1.00000000e+03;

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_dtStableImplicit =   4.09893495e+06;

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_dtStableExplicit =   1.92450090e-01;

const int pylith::materials::PowerLawPlaneStrainElasticData::_numPropertyValues[] = {
1,
1,
1,
1,
1,
1,
};

const int pylith::materials::PowerLawPlaneStrainElasticData::_numStateVarValues[] = {
1,
4,
4,
};

const char* pylith::materials::PowerLawPlaneStrainElasticData::_dbPropertyValues[] = {
"density",
"vs",
"vp",
"reference-strain-rate",
"reference-stress",
"power-law-exponent",
};

const char* pylith::materials::PowerLawPlaneStrainElasticData::_dbStateVarValues[] = {
"stress-zz-initial",
"viscous-strain-xx",
"viscous-strain-yy",
"viscous-strain-zz",
"viscous-strain-xy",
"stress4-xx",
"stress4-yy",
"stress4-zz",
"stress4-xy",
};

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_dbProperties[] = {
  2.50000000e+03,
  3.00000000e+03,
  5.19615242e+03,
  1.00000000e-06,
  2.00000000e+12,
  1.00000000e+00,
  2.00000000e+03,
  1.20000000e+03,
  2.07846097e+03,
  1.00000000e-06,
  1.25992105e+08,
  3.00000000e+00,
};

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_dbStateVars[] = {
  1.07540000e+04,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.57500000e+04,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
};

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_properties[] = {
  2.50000000e+03,
  2.25000000e+10,
  2.25000000e+10,
  1.00000000e-06,
  2.00000000e+12,
  1.00000000e+00,
  2.00000000e+03,
  2.88000000e+09,
  2.88000000e+09,
  1.00000000e-06,
  1.25992105e+08,
  3.00000000e+00,
};

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_stateVars[] = {
  1.07540000e+04,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.57500000e+04,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
};

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_propertiesNondim[] = {
  2.50000000e+00,
  1.00000000e+00,
  1.00000000e+00,
  1.00000000e-06,
  8.88888889e+01,
  1.00000000e+00,
  2.00000000e+00,
  1.28000000e-01,
  1.28000000e-01,
  1.00000000e-06,
  5.59964911e-03,
  3.00000000e+00,
};

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_stateVarsNondim[] = {
  4.77955556e-07,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  1.14444444e-06,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
};

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_density[] = {
  2.50000000e+03,
  2.00000000e+03,
};

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_strain[] = {
  1.10000000e-04,
  1.20000000e-04,
  1.40000000e-04,
  4.10000000e-04,
  4.20000000e-04,
  4.40000000e-04,
};

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_stress[] = {
 -1.79790000e+07,
 -1.79780000e+07,
 -8.97600000e+06,
 -2.25300000e+06,
 -2.25200000e+06,
 -1.09800000e+06,
};

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_elasticConsts[] = {
  6.75000000e+10,
  2.25000000e+10,
  0.00000000e+00,
  2.25000000e+10,
  6.75000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.50000000e+10,
  8.64000000e+09,
  2.88000000e+09,
  0.00000000e+00,
  2.88000000e+09,
  8.64000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.76000000e+09,
};

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_initialStress[] = {
  2.10000000e+04,
  2.20000000e+04,
  2.40000000e+04,
  5.10000000e+04,
  5.20000000e+04,
  5.40000000e+04,
};

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_initialStrain[] = {
  3.10000000e-04,
  3.20000000e-04,
  3.40000000e-04,
  6.10000000e-04,
  6.20000000e-04,
  6.40000000e-04,
};

const PylithScalar pylith::materials::PowerLawPlaneStrainElasticData::_stateVarsUpdated[] = {
  1.07540000e+04,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
 -1.79790000e+07,
 -1.79780000e+07,
  5.18575400e+06,
 -8.97600000e+06,
  2.57500000e+04,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
 -2.25300000e+06,
 -2.25200000e+06,
  2.41615000e+06,
 -1.09800000e+06,
};

pylith::materials::PowerLawPlaneStrainElasticData::PowerLawPlaneStrainElasticData(void)
{ // constructor
  dimension = _dimension;
  numLocs = _numLocs;
  numProperties = _numProperties;
  numStateVars = _numStateVars;
  numDBProperties = _numDBProperties;
  numDBStateVars = _numDBStateVars;
  numPropsQuadPt = _numPropsQuadPt;
  numVarsQuadPt = _numVarsQuadPt;
  lengthScale = _lengthScale;
  timeScale = _timeScale;
  pressureScale = _pressureScale;
  densityScale = _densityScale;
  dtStableImplicit = _dtStableImplicit;
  dtStableExplicit = _dtStableExplicit;
  numPropertyValues = const_cast<int*>(_numPropertyValues);
  numStateVarValues = const_cast<int*>(_numStateVarValues);
  dbPropertyValues = const_cast<char**>(_dbPropertyValues);
  dbStateVarValues = const_cast<char**>(_dbStateVarValues);
  dbProperties = const_cast<PylithScalar*>(_dbProperties);
  dbStateVars = const_cast<PylithScalar*>(_dbStateVars);
  properties = const_cast<PylithScalar*>(_properties);
  stateVars = const_cast<PylithScalar*>(_stateVars);
  propertiesNondim = const_cast<PylithScalar*>(_propertiesNondim);
  stateVarsNondim = const_cast<PylithScalar*>(_stateVarsNondim);
  density = const_cast<PylithScalar*>(_density);
  strain = const_cast<PylithScalar*>(_strain);
  stress = const_cast<PylithScalar*>(_stress);
  elasticConsts = const_cast<PylithScalar*>(_elasticConsts);
  initialStress = const_cast<PylithScalar*>(_initialStress);
  initialStrain = const_cast<PylithScalar*>(_initialStrain);
  stateVarsUpdated = const_cast<PylithScalar*>(_stateVarsUpdated);
} // constructor

pylith::materials::PowerLawPlaneStrainElasticData::~PowerLawPlaneStrainElasticData(void)
{}


// End of file
