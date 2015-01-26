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

generate_fortran_bindings(FORTRAN_BINDINGS ../grids/grids.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../grids/ReducedGrid.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../Field.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../FieldSet.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../FunctionSpace.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../Mesh.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../Metadata.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../mpl/HaloExchange.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../mpl/GatherScatter.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../mpl/Checksum.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../io/Gmsh.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/BuildParallelFields.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/BuildPeriodicBoundaries.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/BuildHalo.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/BuildEdges.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/BuildDualMesh.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/GenerateMesh.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../actions/WriteLoadBalanceReport.h)
generate_fortran_bindings(FORTRAN_BINDINGS atlas_logging.h)
generate_fortran_bindings(FORTRAN_BINDINGS atlas_resource.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../mpi/mpi.h)


### atlas fortran lib

ecbuild_add_library( TARGET atlas_f
    AUTO_VERSION
    CONDITION CMAKE_Fortran_COMPILER_LOADED
    SOURCES
        ${FORTRAN_BINDINGS}
        atlas_module.F90
        atlas_mpi_module.F90
        atlas_grids_module.F90
        atlas_C_interop.F90
        atlas_logging.h
        atlas_logging.cc
        atlas_resource.h
        atlas_resource.cc
  LIBS atlas ${ATLAS_MPI_FORTRAN_LIBRARIES}
)

add_custom_target( atlas_fortran_inc SOURCES
   atlas_module_Checksum_i.f
   atlas_module_Checksum_c.f
   atlas_module_Logging_i.f
   atlas_module_Logging_c.f
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
   atlas_module_GatherScatter_i.f
   atlas_module_GatherScatter_c.f
   atlas_module_Grid_i.f
   atlas_module_Grid_c.f
   atlas_module_HaloExchange_i.f
   atlas_module_HaloExchange_c.f
)

install( FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR}/atlas_module.mod
                 DESTINATION ${INSTALL_INCLUDE_DIR} )


#ecbuild_add_library( TARGET atlas_mpl_fortran
#    AUTO_VERSION
#    CONDITION CMAKE_Fortran_COMPILER_LOADED
#    SOURCES
#        atlas_mpi_module.F90
#    LIBS ${ATLAS_MPI_FORTRAN_LIBRARIES}
#)
