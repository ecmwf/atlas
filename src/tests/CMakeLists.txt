# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FCTEST )

  add_fctest( TARGET fctest_atlas_mesh
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_atlas_mesh.F90
    LIBS            atlas_fortran atlas_mpl_fortran
  )

endif()

set( BOOST_UNIT_TEST_FRAMEWORK_HEADER_ONLY ON )
set( BOOST_UNIT_TEST_FRAMEWORK_LINKED OFF )

ecbuild_add_test( TARGET test_indexview
  BOOST
  SOURCES  test_indexview.cpp
  LIBS     atlas
)

ecbuild_add_test( TARGET test_gmsh
  BOOST
  SOURCES  test_gmsh.cpp
  LIBS     atlas
)

ecbuild_add_test( TARGET test_rgg
  BOOST
  SOURCES  test_rgg.cpp
  LIBS     atlas_meshgen
)

ecbuild_add_test( TARGET test_halo
  BOOST
  MPI      12
  SOURCES  test_halo.cpp
  LIBS     atlas_meshgen
)

ecbuild_add_test( TARGET test_meshgen3d
  BOOST
  SOURCES   test_meshgen3d.cpp
  LIBS      atlas_meshgen
)

if( HAVE_GRID )
  add_subdirectory( grid )
endif()

