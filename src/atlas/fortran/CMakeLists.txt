### fortran bindings

function(generate_fortran_bindings output filename)
	get_filename_component(base ${filename} NAME_WE)
	set(base_abs ${CMAKE_CURRENT_SOURCE_DIR}/${base})
  set(outfile ${CMAKE_CURRENT_BINARY_DIR}/${base}_c_binding.f90)
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
generate_fortran_bindings(FORTRAN_BINDINGS ../mpl/Checksum.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../io/Gmsh.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/BuildParallelFields.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/BuildPeriodicBoundaries.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/BuildHalo.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/BuildEdges.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/BuildDualMesh.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/GenerateMesh.hpp)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/WriteLoadBalanceReport.hpp)


### atlas fortran lib

ecbuild_add_library( TARGET atlas_fortran
	AUTO_VERSION
	CONDITION CMAKE_Fortran_COMPILER_LOADED
	SOURCES
		${FORTRAN_BINDINGS}
		atlas_module.F90
    atlas_C_interop.F90
  LIBS atlas
)

add_custom_target( atlas_fortran_inc SOURCES
   atlas_module_Checksum_i.f
   atlas_module_Checksum_c.f
   atlas_module_Mesh_i.f
   atlas_module_Mesh_c.f
   atlas_module_FunctionSpace_i.f
   atlas_module_FunctionSpace_c.f
   atlas_module_Field_i.f
   atlas_module_Field_c.f
   atlas_module_FieldSet_i.f
   atlas_module_FieldSet_c.f
   atlas_module_Metadata_i.f
   atlas_module_Metadata_c.f
   atlas_module_Gather_i.f
   atlas_module_Gather_c.f
   atlas_module_HaloExchange_i.f
   atlas_module_HaloExchange_c.f
)


ecbuild_add_library( TARGET atlas_mpl_fortran
	AUTO_VERSION
	CONDITION CMAKE_Fortran_COMPILER_LOADED
	SOURCES
		atlas_mpl_module.F90
	LIBS ${ATLAS_MPI_FORTRAN_LIBRARIES}
)
