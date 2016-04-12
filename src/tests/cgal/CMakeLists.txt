# (C) Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_test(
		TARGET      atlas_test_cgal_mesh_gen_from_points
		SOURCES     test_cgal_mesh_gen_from_points.cc
		CONDITION   HAVE_TESSELATION
		LIBS        atlas
        MPI         1 )

