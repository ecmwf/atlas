# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FCTEST )

  if( NOT HAVE_TRANS )
    set( TRANSI_HAVE_MPI 1 )
  endif()

  add_fctest( TARGET atlas_fctest_functionspace
    MPI             4
    CONDITION       ECKIT_HAVE_MPI AND TRANSI_HAVE_MPI
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

ecbuild_add_test( TARGET atlas_test_structuredcolumns
  SOURCES  test_structuredcolumns.cc
  LIBS     atlas
  MPI 5
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  CONDITION ECKIT_HAVE_MPI
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

set( exe_atlas_test_stencil_parallel $<TARGET_FILE:atlas_test_stencil_parallel> )
# When ECBUILD-432 is fixed (already fixed in develop on top of ecbuild 3.0.0), we don't need TARGET_FILE

ecbuild_add_test( TARGET atlas_test_stencil_parallel_mpi4
  COMMAND ${exe_atlas_test_stencil_parallel}
  MPI 4
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  CONDITION ECKIT_HAVE_MPI
)

ecbuild_add_test( TARGET atlas_test_stencil_parallel_mpi16
  COMMAND ${exe_atlas_test_stencil_parallel}
  MPI 16
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  CONDITION ECKIT_HAVE_MPI
)

