# (C) Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FCTEST )

  add_fctest( TARGET atlas_fctest_mesh
    MPI             1
    LINKER_LANGUAGE Fortran
    CONDITION       ON
    SOURCES         fctest_mesh.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_meshgen
    MPI             1
    LINKER_LANGUAGE Fortran
    CONDITION       ON
    SOURCES         fctest_meshgen.F90
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

endif()


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
  MPI        2
  SOURCES    test_parfields.cc
  CONDITION  HAVE_TRANS
  LIBS       atlas
)

ecbuild_add_test( TARGET atlas_test_halo
  BOOST
  MPI        5
  SOURCES    test_halo.cc ../TestMeshes.h
  CONDITION  HAVE_TRANS
  LIBS       atlas
)

ecbuild_add_test( TARGET atlas_test_distmesh
  BOOST
  MPI        5
  SOURCES    test_distmesh.cc ../TestMeshes.h
  CONDITION  HAVE_TRANS
  LIBS       atlas
)

ecbuild_add_test( TARGET atlas_test_meshgen3d
  BOOST
  MPI       1
  SOURCES   test_meshgen3d.cc
  LIBS      atlas
)

ecbuild_add_test(
  TARGET      atlas_test_cgal_mesh_gen_from_points
  SOURCES     test_cgal_mesh_gen_from_points.cc
  CONDITION   HAVE_TESSELATION
  LIBS        atlas
  MPI         1 )

