# (C) Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FCTEST )

  add_fctest( TARGET atlas_fctest_refcounting
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_refcounting.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_mesh
    MPI             1
    LINKER_LANGUAGE Fortran
    CONDITION       OFF
    SOURCES         fctest_mesh.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_connectivity
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_connectivity.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_elements
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_elements.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_logging
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_logging.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_error
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_error.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_state
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_state.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_functionspace
    MPI             4
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_functionspace.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_value
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_value.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_parametrisation
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_parametrisation.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_resource
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_resource.F90
    ARGS            -integer 10 -long 5000000000 -float 0.123456 -double 0.1234567890123456789 -string "hello world"
    LIBS            atlas_f
  )

  if( ATLAS_HAVE_TRANS )
    add_fctest( TARGET atlas_fctest_trans
      MPI             1
      LINKER_LANGUAGE Fortran
      SOURCES         fctest_trans.F90
      LIBS            atlas_f
    )
  endif()

endif()

#set( BOOST_UNIT_TEST_FRAMEWORK_HEADER_ONLY ON )
#set( BOOST_UNIT_TEST_FRAMEWORK_LINKED OFF )

ecbuild_add_test( TARGET atlas_test_flags
  BOOST
  MPI      1
  SOURCES  test_flags.cc
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_indexview
  BOOST
  MPI      1
  SOURCES  test_indexview.cc
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_haloexchange
  BOOST
  MPI      3
  SOURCES  test_haloexchange.cc
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_gather
  BOOST
  MPI      3
  SOURCES  test_gather.cc
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_state
  BOOST
  MPI      1
  SOURCES  test_state.cc
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_elements
  BOOST
  MPI      1
  SOURCES  test_elements.cc
  LIBS     atlas
)


ecbuild_add_test( TARGET atlas_test_accumulate_facets
  BOOST
  MPI      1
  SOURCES  test_accumulate_facets.cc
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_functionspace
  BOOST
  MPI      1
  SOURCES  test_functionspace.cc
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_gmsh
  BOOST
  MPI      1
  SOURCES  test_gmsh.cc TestMeshes.h
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_rgg
  BOOST
  MPI      1
  SOURCES  test_rgg.cc
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_ll
  BOOST
  MPI      1
  SOURCES  test_ll.cc
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_parfields
  BOOST
  MPI      2
  SOURCES  test_parfields.cc
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_halo
  BOOST
  MPI      5
  SOURCES  test_halo.cc TestMeshes.h
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_distmesh
  BOOST
  MPI      5
  SOURCES  test_distmesh.cc TestMeshes.h
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_meshgen3d
  BOOST
  MPI       1
  SOURCES   test_meshgen3d.cc
  LIBS      atlas
)

ecbuild_add_test( TARGET atlas_test_grids
  BOOST
  MPI       1
  SOURCES   test_grids.cc
  LIBS      atlas
)

ecbuild_add_test( TARGET atlas_test_pointcloud
  BOOST
  MPI       1
  SOURCES   test_pointcloud.cc
  LIBS      atlas
)

ecbuild_add_test( TARGET atlas_test_trans
  BOOST
  MPI      4
  SOURCES   test_trans.cc
  CONDITION HAVE_TRANS
  LIBS      atlas
)

ecbuild_add_test( TARGET atlas_test_trans_serial
  BOOST
  SOURCES   test_trans.cc
  CONDITION HAVE_TRANS AND NOT HAVE_MPI
  LIBS      atlas
)

add_subdirectory( cgal )
add_subdirectory( geometry )
add_subdirectory( numerics )

