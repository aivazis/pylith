#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#
# @file tests/pytests/topology/TestSubfield.py
#
# @brief Unit testing of Python Subfield object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.topology.Subfield import (Subfield, subfield)


class TestSubfield(TestComponent):
    """Unit testing of Subfield object.
    """
    _class = Subfield
    _factory = subfield


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSubfield))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
