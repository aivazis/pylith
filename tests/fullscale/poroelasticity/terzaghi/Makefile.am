# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------

include $(top_srcdir)/check.am

TESTS = test_pylith.py

dist_check_SCRIPTS = test_pylith.py

dist_noinst_PYTHON = \
	meshes.py \
	TestTerzaghi.py \
	terzaghi_soln.py \
	terzaghi_gendb.py \
	TestTerzaghiCompaction.py \
	terzaghi_compaction_soln.py \
	terzaghi_compaction_gendb.py	

dist_noinst_DATA = \
	geometry.jou \
	bc.jou \
	mesh_tri.jou \
	mesh_tri.exo \
	mesh_quad.jou \
	mesh_quad.exo \
	terzaghi.cfg \
	terzaghi_tri.cfg \
	terzaghi_quad.cfg \
	terzaghi_compaction.cfg \
	terzaghi_compaction_tri.cfg \
	terzaghi_compaction_quad.cfg	

noinst_TMP = \
	terzaghi_bc.spatialdb \
	terzaghi_ic.spatialdb \
	terzaghi_compaction_bc.spatialdb \
	terzaghi_compaction_ic.spatialdb	



# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/fullscale/poroelasticity/terzaghi
export-data: $(dist_noinst_PYTHON) $(dist_noinst_DATA)
	for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi


BUILT_SOURCES = export-data

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	-rm *.pyc output/*


# End of file
