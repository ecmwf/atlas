# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

### fortranize function

function( fortranize file )
  unset( result )
  file( STRINGS "${file}" lines )
  foreach( i IN LISTS lines )
  string( REGEX REPLACE "^/\\* (#undef .*) \\*/$" "! \\1" i "${i}" )
  set( result "${result}${i}\n" )
  endforeach()
  file( WRITE "${file}" "${result}" )
endfunction()

### config headers

ecbuild_generate_config_headers( DESTINATION ${INSTALL_INCLUDE_DIR}/atlas )

configure_file( atlas_config.h.in   atlas_config.h  )
configure_file( atlas_version.h.in  atlas_version.h )
configure_file( atlas_defines.h.in  atlas_defines.h   )
configure_file( atlas_defines.h.in  atlas_defines_fortran.h )
fortranize( ${CMAKE_CURRENT_BINARY_DIR}/atlas_defines_fortran.h )

install( FILES
  ${CMAKE_CURRENT_BINARY_DIR}/atlas_defines.h
  ${CMAKE_CURRENT_BINARY_DIR}/atlas_defines_fortran.h
  ${CMAKE_CURRENT_BINARY_DIR}/atlas_config.h
  ${CMAKE_CURRENT_BINARY_DIR}/atlas_version.h
  DESTINATION
  ${INSTALL_INCLUDE_DIR}/atlas )

configure_file( atlas_version.cc.in atlas_version.cc )

### mpistubs

if( NOT HAVE_MPI )

   set( atlas_mpistubs_srcs mpistubs/mpi.h mpistubs/mpi.c )

  ecbuild_add_library( TARGET atlas_mpistubs_fortran
    AUTO_VERSION
    INSTALL_HEADERS ALL
    HEADER_DESTINATION  ${INSTALL_INCLUDE_DIR}/atlas
    CONDITION HAVE_FORTRAN AND CMAKE_Fortran_COMPILER_LOADED
    SOURCES   mpistubs/atlas_mpistubs_module.F90 )

  set( ATLAS_MPI_FORTRAN_LIBRARIES atlas_mpistubs_fortran )

else()

  set( ATLAS_MPI_FORTRAN_LIBRARIES ${MPI_Fortran_LIBRARIES} )

endif()

### sources

list( APPEND atlas_mpl_srcs
mpl/MPLArrayView.h
mpl/Checksum.h
mpl/Checksum.cc
mpl/HaloExchange.h
mpl/HaloExchange.cc
mpl/GatherScatter.h
mpl/GatherScatter.cc
mpl/MPL.h
)

list( APPEND atlas_io_srcs
io/Gmsh.h
io/Gmsh.cc
# to be abstracted into a plugin
io/Grib.cc
io/Grib.h
)

list( APPEND atlas_util_srcs
util/Array.h
util/ArrayUtil.h
util/ArrayView.h
util/ArrayView.cc
util/Checksum.h
util/Checksum.cc
util/Debug.h
util/IndexView.h
util/IndexView.cc
)

list( APPEND atlas_meshgen_srcs
meshgen/EqualAreaPartitioner.cc
meshgen/EqualAreaPartitioner.h
meshgen/RGG.cc
meshgen/RGG.h
meshgen/RGGMeshGenerator.cc
meshgen/RGGMeshGenerator.h
meshgen/T63.cc
meshgen/T95.cc
meshgen/T159.cc
meshgen/T255.cc
meshgen/T511.cc
meshgen/T1279.cc
meshgen/T2047.cc
meshgen/T3999.cc
meshgen/T7999.cc
)

list( APPEND atlas_actions_srcs
actions/BuildDualMesh.h
actions/BuildDualMesh.cc
actions/BuildEdges.h
actions/BuildEdges.cc
actions/BuildHalo.h
actions/BuildHalo.cc
actions/BuildParallelFields.h
actions/BuildParallelFields.cc
actions/BuildPeriodicBoundaries.h
actions/BuildPeriodicBoundaries.cc
actions/DistributeMesh.h
actions/DistributeMesh.cc
actions/GenerateMesh.h
actions/GenerateMesh.cc
actions/WriteLoadBalanceReport.h
actions/WriteLoadBalanceReport.cc
)


list( APPEND atlas_srcs
Field.cc
Field.h
FieldGroup.cc
FieldGroup.h
FieldSet.cc
FieldSet.h
FunctionSpace.cc
FunctionSpace.h
Grid.cc
Grid.h
GridSpec.cc
GridSpec.h
GridSpecParams.cc
GridSpecParams.h
Latitudes.cc
Latitudes.h
Mesh.cc
Mesh.h
MeshCache.cc
MeshCache.h
Metadata.cc
Metadata.h
Parameters.h
PointIndex3.cc
PointIndex3.h
PointSet.cc
PointSet.h
PolarStereoGraphic.cc
PolarStereoGraphic.h
ReducedGrid.h
ReducedGrid.cc
ReducedLatLon.cc
ReducedLatLon.h
RegularGG.cc
RegularGG.h
ReducedGG.cc
ReducedGG.h
RegularLatLon.cc
RegularLatLon.h
RotatedLatLon.cc
RotatedLatLon.h
Tesselation.cc
Tesselation.h
TriangleIntersection.cc
TriangleIntersection.h
Unstructured.cc
Unstructured.h
Util.cc
Util.h
)

### atlas c++ library

ecbuild_add_library( TARGET atlas_core
  AUTO_VERSION
  INSTALL_HEADERS ALL
  HEADER_DESTINATION  ${INSTALL_INCLUDE_DIR}/atlas
  SOURCES
  atlas.h
  atlas.cc
  atlas_mpi.h
  ${atlas_srcs}
  ${atlas_mpistubs_srcs}
  ${atlas_util_srcs}
  ${atlas_io_srcs}
  ${atlas_mpl_srcs}
  ${atlas_meshgen_srcs}
  ${atlas_actions_srcs}
  ${CMAKE_CURRENT_BINARY_DIR}/atlas_version.h
  ${CMAKE_CURRENT_BINARY_DIR}/atlas_version.cc
  ${CMAKE_CURRENT_BINARY_DIR}/atlas_config.h
  ${CMAKE_CURRENT_BINARY_DIR}/atlas_defines.h
  INCLUDES "${CGAL_INCLUDE_DIRS}"
  LIBS eckit "${ATLAS_EXTRA_LIBRARIES}" "${CGAL_LIBRARIES}" "${MPI_CXX_LIBRARIES}" )

if( HAVE_FORTRAN )
  add_subdirectory( fortran )
endif()

add_subdirectory( grids )

