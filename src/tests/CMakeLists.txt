# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FCTEST )

  add_fctest( TARGET atlas_fctest_mesh
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_mesh.F90
    LIBS            atlas_fortran atlas_grids_fortran
  )

  add_fctest( TARGET atlas_fctest_logging
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_logging.F90
    LIBS            atlas_fortran
  )

endif()

set( BOOST_UNIT_TEST_FRAMEWORK_HEADER_ONLY ON )
set( BOOST_UNIT_TEST_FRAMEWORK_LINKED OFF )

ecbuild_add_test( TARGET atlas_test_flags
  BOOST
  SOURCES  test_flags.cc
  LIBS     ${libatlas}
)

ecbuild_add_test( TARGET atlas_test_indexview
  BOOST
  SOURCES  test_indexview.cc
  LIBS     ${libatlas}
)

ecbuild_add_test( TARGET atlas_test_haloexchange
  BOOST
  MPI      3
  SOURCES  test_haloexchange.cc
  LIBS     ${libatlas}
)

ecbuild_add_test( TARGET atlas_test_gather
  BOOST
  MPI      3
  SOURCES  test_gather.cc
  LIBS     ${libatlas}
)

ecbuild_add_test( TARGET atlas_test_functionspace
  BOOST
  SOURCES  test_functionspace.cc
  LIBS     ${libatlas} eckit_maths
)

ecbuild_add_test( TARGET atlas_test_gmsh
  BOOST
  SOURCES  test_gmsh.cc TestMeshes.h
  LIBS     ${libatlas}
)

ecbuild_add_test( TARGET atlas_test_rgg
  BOOST
  SOURCES  test_rgg.cc
  LIBS     ${libatlas}
)

ecbuild_add_test( TARGET atlas_test_parfields
  BOOST
  MPI 2
  SOURCES  test_parfields.cc
  LIBS     ${libatlas}
)

ecbuild_add_test( TARGET atlas_test_halo
  BOOST
  MPI      5
  SOURCES  test_halo.cc TestMeshes.h
  LIBS     ${libatlas}
)

ecbuild_add_test( TARGET atlas_test_distmesh
  BOOST
  MPI      5
  SOURCES  test_distmesh.cc TestMeshes.h
  LIBS     ${libatlas}
)

ecbuild_add_test( TARGET atlas_test_meshgen3d
  BOOST
  SOURCES   test_meshgen3d.cc
  LIBS      ${libatlas}
)

ecbuild_add_test( TARGET atlas_test_grids
  BOOST
  SOURCES   test_grids.cc
  LIBS      ${libatlas} atlas_grids
)



add_subdirectory( grid )

