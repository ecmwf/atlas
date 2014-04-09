# (C) Copyright 1996-2014 ECMWF.

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
	set( mpistubs_srcs mpistubs/mpi.h mpistubs/mpi.c )
else()
	set( mpistubs_srcs "" )
endif()

### atlas c++ library

ecbuild_add_library( TARGET atlas
	AUTO_VERSION
	SOURCES
		${mpistubs_srcs}
		atlas_mpi.h
		BuildDualMesh.hpp
		BuildDualMesh.cpp
		BuildEdges.hpp
		BuildEdges.cpp
		BuildHalo.hpp
		BuildHalo.cpp
		BuildPeriodicBoundaries.hpp
		BuildPeriodicBoundaries.cpp
		Field.hpp
		Field.cpp
		FieldSet.hpp
		FieldSet.cpp
		FunctionSpace.hpp
		FunctionSpace.cpp
		Gmsh.hpp
		Gmsh.cpp
		Mesh.hpp
		Mesh.cpp
		Metadata.hpp
		Metadata.cpp
		Parameters.hpp
		Partitioner.hpp
		Partitioner.cpp
		HaloExchange.hpp
		HaloExchange.cpp
		Gather.hpp
		Gather.cpp
		MPL.hpp
  LIBS ${MPI_CXX_LIBRARIES}
) 

### fortran bindings

function(generate_fortran_bindings output filename)
	get_filename_component(base ${filename} NAME_WE)
	set(base_abs ${CMAKE_CURRENT_SOURCE_DIR}/${base})
	set(outfile ${CMAKE_CURRENT_BINARY_DIR}/${base}.f90)
	set(${output} ${${output}} ${outfile} PARENT_SCOPE)
	add_custom_command(
		OUTPUT ${outfile}
		COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/c2f.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${outfile}
		DEPENDS ${filename} )
	set_source_files_properties(${outfile} PROPERTIES GENERATED TRUE)
endfunction()

generate_fortran_bindings(FORTRAN_BINDINGS Field.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS FieldSet.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS FunctionSpace.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS Mesh.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS Metadata.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS HaloExchange.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS Gather.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS Gmsh.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS BuildPeriodicBoundaries.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS BuildEdges.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS BuildDualMesh.hpp)

### atlas fortran lib

ecbuild_add_library( TARGET atlasf
	AUTO_VERSION
	CONDITION CMAKE_Fortran_COMPILER_LOADED
	SOURCES
		${FORTRAN_BINDINGS}
		atlas_module.F90
	LIBS atlas
)

add_subdirectory( grid )

add_subdirectory( meshgen )
