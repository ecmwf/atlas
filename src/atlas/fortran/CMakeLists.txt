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

generate_fortran_bindings(FORTRAN_BINDINGS ../mesh/Field.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../mesh/FieldSet.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../mesh/FunctionSpace.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../mesh/Mesh.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../mesh/Metadata.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../mpl/HaloExchange.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../mpl/Gather.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../io/Gmsh.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/BuildPeriodicBoundaries.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/BuildEdges.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/BuildDualMesh.hpp)


### atlas fortran lib

ecbuild_add_library( TARGET atlas_fortran
	AUTO_VERSION
	CONDITION CMAKE_Fortran_COMPILER_LOADED
	SOURCES
		${FORTRAN_BINDINGS}
		atlas_module.F90
	LIBS atlas 
)

ecbuild_add_library( TARGET atlas_mpl_fortran
	AUTO_VERSION
	CONDITION CMAKE_Fortran_COMPILER_LOADED
	SOURCES
		atlas_mpl_module.F90
	LIBS ${ATLAS_MPI_FORTRAN_LIBRARIES}
)
