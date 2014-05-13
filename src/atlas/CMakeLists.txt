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
    foreach(i IN LISTS lines)
        STRING(REGEX REPLACE "^/\\* (#undef .*) \\*/$" "! \\1" i "${i}")
        SET(result "${result}${i}\n")
    endforeach()
    file(WRITE "${file}" "${result}")
endfunction()

configure_file( atlas_defines.h.in  atlas_defines.h   )
configure_file( atlas_defines.h.in  atlas_defines.F90 )
fortranize( ${CMAKE_CURRENT_BINARY_DIR}/atlas_defines.F90 )

configure_file( atlas_config.h.in   atlas_config.h  )

configure_file( atlas_version.h.in  atlas_version.h )

install(FILES
			${CMAKE_CURRENT_BINARY_DIR}/atlas_defines.h
			${CMAKE_CURRENT_BINARY_DIR}/atlas_defines.F90
			${CMAKE_CURRENT_BINARY_DIR}/atlas_config.h
			${CMAKE_CURRENT_BINARY_DIR}/atlas_version.h
		DESTINATION
			${INSTALL_INCLUDE_DIR}/atlas )

configure_file( atlas_version.cc.in atlas_version.cc )

### mpistubs

if( NOT HAVE_MPI )
	set( atlas_mpistubs_srcs mpistubs/mpi.h mpistubs/mpi.c )
else()
	set( mpistubs_srcs "" )
endif()

###

list( APPEND atlas_mpl_srcs
		mpl/HaloExchange.hpp
		mpl/HaloExchange.cpp
		mpl/Gather.hpp
		mpl/Gather.cpp
		mpl/MPL.hpp
)

list( APPEND atlas_io_srcs
		io/Gmsh.hpp
		io/Gmsh.cpp
)

list( APPEND atlas_mesh_srcs
		mesh/Array.hpp
		mesh/ArrayView.hpp
		mesh/ArrayView.cpp
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
)

list( APPEND atlas_actions_srcs
		actions/BuildDualMesh.hpp
		actions/BuildDualMesh.cpp
		actions/BuildEdges.hpp
		actions/BuildEdges.cpp
		actions/BuildHalo.hpp
		actions/BuildHalo.cpp
		actions/BuildPeriodicBoundaries.hpp
		actions/BuildPeriodicBoundaries.cpp
)

### atlas c++ library

ecbuild_add_library( TARGET atlas
	AUTO_VERSION
	SOURCES
		${atlas_mpistubs_srcs}
		${atlas_io_srcs}
		${atlas_mesh_srcs}
		${atlas_mpl_srcs}
		${atlas_actions_srcs}
		atlas_mpi.h
  LIBS ${MPI_CXX_LIBRARIES}
) 

add_subdirectory( fortran )

add_subdirectory( grid )

add_subdirectory( meshgen )

add_subdirectory( testing )
