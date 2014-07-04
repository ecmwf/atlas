# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

### config headers

ecbuild_generate_config_headers( DESTINATION ${INSTALL_INCLUDE_DIR}/atlas )

function( fortranize file )
  unset( result )
  file( STRINGS "${file}" lines )
  foreach( i IN LISTS lines )
    string( REGEX REPLACE "^/\\* (#undef .*) \\*/$" "! \\1" i "${i}" )
    set( result "${result}${i}\n" )
  endforeach()
  file( WRITE "${file}" "${result}" )
endfunction()

configure_file( atlas_defines.h.in  atlas_defines.h   )
configure_file( atlas_defines.h.in  atlas_defines_fortran.h )
fortranize( ${CMAKE_CURRENT_BINARY_DIR}/atlas_defines_fortran.h )

configure_file( atlas_config.h.in   atlas_config.h  )

configure_file( atlas_version.h.in  atlas_version.h )

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
                       CONDITION CMAKE_Fortran_COMPILER_LOADED
                       SOURCES   mpistubs/atlas_mpistubs_module.F90 )
  set( ATLAS_MPI_FORTRAN_LIBRARIES atlas_mpistubs_fortran )
else()
  set( ATLAS_MPI_FORTRAN_LIBRARIES ${MPI_Fortran_LIBRARIES} )
endif()

###

list( APPEND atlas_mpl_srcs
  mpl/Checksum.hpp
  mpl/Checksum.cpp
  mpl/HaloExchange.hpp
  mpl/HaloExchange.cpp
  mpl/GatherScatter.hpp
  mpl/GatherScatter.cpp
  mpl/MPL.hpp
)

list( APPEND atlas_io_srcs
  io/Gmsh.hpp
  io/Gmsh.cpp
)

list( APPEND atlas_util_srcs
  util/Array.hpp
  util/ArrayUtil.hpp
  util/ArrayView.hpp
  util/ArrayView.cpp
  util/Checksum.hpp
  util/Checksum.cpp
  util/Debug.hpp
  util/IndexView.hpp
  util/IndexView.cpp
)

list( APPEND atlas_meshgen_srcs
  meshgen/EqualAreaPartitioner.cpp
  meshgen/EqualAreaPartitioner.hpp
  meshgen/RGG.cpp
  meshgen/RGG.hpp
  meshgen/RGGMeshGenerator.cpp
  meshgen/RGGMeshGenerator.hpp
  meshgen/T63.cpp
  meshgen/T95.cpp
  meshgen/T159.cpp
  meshgen/T255.cpp
  meshgen/T511.cpp
  meshgen/T1279.cpp
  meshgen/T2047.cpp
  meshgen/T3999.cpp
  meshgen/T7999.cpp
)

list( APPEND atlas_mesh_srcs
  mesh/Field.hpp
  mesh/Field.cpp
  mesh/FieldSet.hpp
  mesh/FieldSet.cpp
  mesh/FunctionSpace.hpp
  mesh/FunctionSpace.cpp
  mesh/Mesh.hpp
  mesh/Mesh.cpp
  mesh/Metadata.hpp
  mesh/Metadata.cpp
  mesh/Parameters.hpp
  mesh/Util.hpp
  mesh/Util.cpp
)

list( APPEND atlas_actions_srcs
  actions/BuildDualMesh.hpp
  actions/BuildDualMesh.cpp
  actions/BuildEdges.hpp
  actions/BuildEdges.cpp
  actions/BuildHalo.hpp
  actions/BuildHalo.cpp
  actions/BuildParallelFields.hpp
  actions/BuildParallelFields.cpp
  actions/BuildPeriodicBoundaries.hpp
  actions/BuildPeriodicBoundaries.cpp
  actions/DistributeMesh.hpp
  actions/DistributeMesh.cpp
  actions/GenerateMesh.hpp
  actions/GenerateMesh.cpp
  actions/WriteLoadBalanceReport.hpp
  actions/WriteLoadBalanceReport.cpp
)

### atlas c++ library

ecbuild_add_library( TARGET atlas
                     AUTO_VERSION
                     SOURCES
                       ${atlas_mpistubs_srcs}
                       ${atlas_util_srcs}
                       ${atlas_io_srcs}
                       ${atlas_mesh_srcs}
                       ${atlas_mpl_srcs}
                       ${atlas_meshgen_srcs}
                       ${atlas_actions_srcs}
                       atlas_mpi.h
                     LIBS eckit ${MPI_CXX_LIBRARIES} )

add_subdirectory( fortran )

add_subdirectory( grid )

add_subdirectory( meshgen )
