# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


function(generate_fortran_bindings output filename)
    get_filename_component(base ${filename} NAME_WE)
    set(base_abs ${CMAKE_CURRENT_SOURCE_DIR}/${base})
    set(outfile ${CMAKE_CURRENT_BINARY_DIR}/${base}.f90) 
    set(${output} ${${output}} ${outfile} PARENT_SCOPE)
    add_custom_command(
        OUTPUT ${outfile}
		COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/c2f.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${outfile}
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


ecbuild_add_library( TARGET atlas
	AUTO_VERSION
	TYPE BOTH
	SOURCES
	  # C++ part
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

	  # Fortran part
		${FORTRAN_BINDINGS}
		atlas.F90

	LIBS ${TPL}
) 

