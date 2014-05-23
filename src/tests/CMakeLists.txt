# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FCTEST )

add_fctest( TARGET fctest_atlas_mesh
  LIBS            atlas_fortran atlas_mpl_fortran
  LINKER_LANGUAGE Fortran
  SOURCES         fctest_atlas_mesh.F90
)

endif()

ecbuild_add_test( TARGET test_gmsh
  SOURCES  test_gmsh.cpp
  LIBS     atlas
  INCLUDES  ${Boost_INCLUDE_DIRS}
  CONDITION HAVE_BOOST
)

ecbuild_add_test( TARGET test_meshgen3d
  SOURCES   test_meshgen3d.cpp
  LIBS      atlas_meshgen atlas
  INCLUDES  ${Boost_INCLUDE_DIRS}
  CONDITION HAVE_BOOST
)

ecbuild_add_test( TARGET test_rgg
	SOURCES  test_rgg.cpp
	LIBS     atlas_meshgen atlas_testing
)

if( HAVE_GRID )
  add_subdirectory( grid )
endif()

