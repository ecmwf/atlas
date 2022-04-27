# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FCTEST )

  if( NOT HAVE_TRANS )
    set( transi_HAVE_MPI 1 )
  endif()

  add_fctest( TARGET atlas_fctest_functionspace
    MPI             4
    CONDITION       eckit_HAVE_MPI AND transi_HAVE_MPI
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_functionspace.F90
    LIBS            atlas_f
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

endif()

ecbuild_add_test( TARGET atlas_test_functionspace
  SOURCES  test_functionspace.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_cellcolumns
  SOURCES  test_cellcolumns.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_cubedsphere_functionspace
  MPI        8
  CONDITION  eckit_HAVE_MPI
  SOURCES test_cubedsphere_functionspace.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET test_structuredcolumns_biperiodic
  SOURCES  test_structuredcolumns_biperiodic.cc
  LIBS     atlas
  MPI 5
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  CONDITION eckit_HAVE_MPI AND MPI_SLOTS GREATER_EQUAL 5
)

ecbuild_add_test( TARGET test_structuredcolumns_regional
  SOURCES  test_structuredcolumns_regional.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)


ecbuild_add_test( TARGET atlas_test_structuredcolumns
  SOURCES  test_structuredcolumns.cc
  LIBS     atlas
  MPI 5
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  CONDITION eckit_HAVE_MPI AND MPI_SLOTS GREATER_EQUAL 5
)

ecbuild_add_test( TARGET atlas_test_pointcloud
  SOURCES  test_pointcloud.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_reduced_halo
  SOURCES test_reduced_halo.cc
  LIBS    atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_stencil
  SOURCES  test_stencil.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_executable( TARGET atlas_test_stencil_parallel
  SOURCES  test_stencil_parallel.cc
  LIBS     atlas
  NOINSTALL
)

ecbuild_add_test( TARGET atlas_test_stencil_parallel_mpi4
  COMMAND atlas_test_stencil_parallel
  MPI 4
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  CONDITION eckit_HAVE_MPI
)

ecbuild_add_test( TARGET atlas_test_stencil_parallel_mpi16
  COMMAND atlas_test_stencil_parallel
  MPI 16
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  CONDITION eckit_HAVE_MPI AND MPI_SLOTS GREATER_EQUAL 16
)

set( _WITH_MPI )
if( eckit_HAVE_MPI )
    set( _WITH_MPI MPI 4 )
endif()

ecbuild_add_executable( TARGET atlas_test_polygons
  SOURCES  test_polygons.cc
  LIBS     atlas
  NOINSTALL
)

ecbuild_add_test( TARGET atlas_test_polygons_structuredcolumns
  COMMAND atlas_test_polygons ARGS --functionspace StructuredColumns
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  ${_WITH_MPI}
)

ecbuild_add_test( TARGET atlas_test_polygons_nodecolumns
  COMMAND atlas_test_polygons ARGS --functionspace NodeColumns
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  ${_WITH_MPI}
)

ecbuild_add_executable( TARGET atlas_test_polygons_projection
  SOURCES  test_polygons_projection.cc
  LIBS     atlas
  NOINSTALL
)

ecbuild_add_test( TARGET atlas_test_polygons_projections_structuredcolumns
  COMMAND atlas_test_polygons_projection ARGS --functionspace StructuredColumns
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  ${_WITH_MPI}
)

ecbuild_add_test( TARGET atlas_test_polygons_projections_nodecolumns
  COMMAND atlas_test_polygons_projection ARGS --functionspace NodeColumns
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  ${_WITH_MPI}
)

# Tests ATLAS-286
ecbuild_add_test( TARGET  atlas_test_structuredcolumns_haloexchange
  ${_WITH_MPI}
  SOURCES test_structuredcolumns_haloexchange.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

