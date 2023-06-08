# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FCTEST )

  add_fctest( TARGET atlas_fctest_mesh
    LINKER_LANGUAGE Fortran
    CONDITION       NOT atlas_fctest_mesh_DISABLED
    SOURCES         fctest_mesh.F90
    LIBS            atlas_f
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

  add_fctest( TARGET atlas_fctest_meshgen
    LINKER_LANGUAGE Fortran
    CONDITION       ON
    SOURCES         fctest_meshgen.F90
    LIBS            atlas_f
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

  add_fctest( TARGET atlas_fctest_connectivity
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_connectivity.F90
    LIBS            atlas_f
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

  add_fctest( TARGET atlas_fctest_elements
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_elements.F90
    LIBS            atlas_f
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

endif()

ecbuild_add_test( TARGET atlas_test_parfields
  MPI        2
  CONDITION  eckit_HAVE_MPI
  SOURCES    test_parfields.cc
  LIBS       atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_halo
  MPI        5
  CONDITION  eckit_HAVE_MPI AND MPI_SLOTS GREATER_EQUAL 5
  SOURCES    test_halo.cc ../TestMeshes.h
  LIBS       atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_distmesh
  MPI        5
  CONDITION  eckit_HAVE_MPI AND MPI_SLOTS GREATER_EQUAL 5
  SOURCES    test_distmesh.cc ../TestMeshes.h
  LIBS       atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test(
  TARGET      atlas_test_cgal_mesh_gen_from_points
  SOURCES     test_cgal_mesh_gen_from_points.cc
  CONDITION   atlas_HAVE_TESSELATION
  LIBS        atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test(
  TARGET atlas_test_mesh_build_edges
  SOURCES test_mesh_build_edges.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_mesh_node2cell
  MPI        4
  CONDITION  eckit_HAVE_MPI
  SOURCES    test_mesh_node2cell.cc
  LIBS       atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)


foreach( test connectivity stream_connectivity elements ll meshgen3d rgg )
  ecbuild_add_test( TARGET atlas_test_${test}
    SOURCES test_${test}.cc
    LIBS atlas
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )
endforeach()

ecbuild_add_test( TARGET atlas_test_healpixmeshgen
  SOURCES test_healpixmeshgen.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_cubedsphere_meshgen
  MPI        8
  CONDITION  eckit_HAVE_MPI AND MPI_SLOTS GREATER_EQUAL 8
  SOURCES test_cubedsphere_meshgen.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

if( TEST atlas_test_cubedsphere_meshgen )
    set_tests_properties( atlas_test_cubedsphere_meshgen PROPERTIES TIMEOUT 30 )
endif()

ecbuild_add_test( TARGET atlas_test_mesh_builder
  SOURCES test_mesh_builder.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_mesh_builder_parallel
  MPI        6
  CONDITION  eckit_HAVE_MPI AND MPI_SLOTS GREATER_EQUAL 6
  SOURCES test_mesh_builder_parallel.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_executable( TARGET atlas_test_mesh_reorder
  SOURCES  test_mesh_reorder.cc
  LIBS     atlas
  NOINSTALL
)

ecbuild_add_test( TARGET atlas_test_mesh_reorder_O16
  COMMAND atlas_test_mesh_reorder ARGS --grid O16
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_mesh_reorder_unstructured
  COMMAND atlas_test_mesh_reorder ARGS --mesh ${CMAKE_CURRENT_SOURCE_DIR}/test_mesh_reorder_unstructured.msh --grid "unstructured"
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)


atlas_add_cuda_test(
  TARGET atlas_test_connectivity_kernel
  SOURCES test_connectivity_kernel.cu
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_pentagon_element
  SOURCES test_pentagon_element.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)
