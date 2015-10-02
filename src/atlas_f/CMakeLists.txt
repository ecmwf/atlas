
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

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/../atlas/atlas_defines.h.in  atlas_f_defines.h )
  # fortranize( ${CMAKE_CURRENT_BINARY_DIR}/atlas_f_defines.h )

install( FILES
  ${CMAKE_CURRENT_BINARY_DIR}/atlas_f_defines.h
  ${CMAKE_CURRENT_SOURCE_DIR}/atlas_f_mpi.h
  DESTINATION
  ${INSTALL_INCLUDE_DIR}/atlas_f )


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

function(generate_fortran_bindings_custom_output output filename outfile)
    get_filename_component(base ${filename} NAME_WE)
    set(base_abs ${CMAKE_CURRENT_SOURCE_DIR}/${base})
    set(outfile ${CMAKE_CURRENT_BINARY_DIR}/${outfile})
    set(${output} ${${output}} ${outfile} PARENT_SCOPE)
    add_custom_command(
        OUTPUT ${outfile}
        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/c2f.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${outfile}
        DEPENDS ${filename} )
    set_source_files_properties(${outfile} PROPERTIES GENERATED TRUE)
endfunction()


generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/grids/grids.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/grids/ReducedGrid.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/GridDistribution.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/State.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/Field.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/FieldSet.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/FunctionSpace.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/Mesh.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/Metadata.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/Nodes.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/Config.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mpl/HaloExchange.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mpl/GatherScatter.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mpl/Checksum.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/io/Gmsh.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/actions/BuildParallelFields.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/actions/BuildPeriodicBoundaries.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/actions/BuildHalo.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/actions/BuildEdges.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/actions/BuildDualMesh.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/actions/GenerateMesh.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/actions/WriteLoadBalanceReport.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/runtime/ErrorHandling.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/atlas.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mpi/mpi.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/trans/Trans.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/functionspace/SpectralFunctionSpace.h)
generate_fortran_bindings_custom_output(FORTRAN_BINDINGS ../atlas/functionspace/NodesFunctionSpaceInterface.h NodesFunctionSpace_c_binding.f90)
generate_fortran_bindings(FORTRAN_BINDINGS atlas_logging.h)
generate_fortran_bindings(FORTRAN_BINDINGS atlas_resource.h)
generate_fortran_bindings(FORTRAN_BINDINGS atlas_value.h)


### atlas fortran lib
ecbuild_add_library( TARGET atlas_f
    AUTO_VERSION
    CONDITION CMAKE_Fortran_COMPILER_LOADED
    SOURCES
        ${FORTRAN_BINDINGS}
        atlas_object_module.F90
        atlas_refcounted_module.F90
        atlas_module.F90
        atlas_mpi_module.F90
        atlas_grids_module.F90
        atlas_C_interop.F90
        atlas_C_interop.cc
        atlas_logging.h
        atlas_logging.cc
        atlas_resource.h
        atlas_resource.cc
        atlas_value.h
        atlas_value.cc
  INCLUDES ${MPI_Fortran_INCLUDE_PATH}
  LIBS atlas ${MPI_Fortran_LIBRARIES}
)

add_custom_target( atlas_fortran_inc SOURCES
   atlas_module_Checksum_i.f
   atlas_module_Checksum_c.f
   atlas_module_Config_i.f
   atlas_module_Config_c.f
   atlas_module_Error_i.f
   atlas_module_Error_c.f
   atlas_module_Field_i.f
   atlas_module_Field_c.f
   atlas_module_FieldSet_i.f
   atlas_module_FieldSet_c.f
   atlas_module_FunctionSpace_i.f
   atlas_module_FunctionSpace_c.f
   atlas_module_GatherScatter_i.f
   atlas_module_GatherScatter_c.f
   atlas_module_Grid_i.f
   atlas_module_Grid_c.f
   atlas_module_GridDistribution_i.f
   atlas_module_GridDistribution_c.f
   atlas_module_HaloExchange_i.f
   atlas_module_HaloExchange_c.f
   atlas_module_JSON_i.f
   atlas_module_JSON_c.f
   atlas_module_Logging_i.f
   atlas_module_Logging_c.f
   atlas_module_Mesh_i.f
   atlas_module_Mesh_c.f
   atlas_module_Metadata_i.f
   atlas_module_Metadata_c.f
   atlas_module_Nodes_i.f
   atlas_module_Nodes_c.f
   atlas_module_NodesFunctionSpace_i.f
   atlas_module_NodesFunctionSpace_c.f
   atlas_module_PathName_i.f
   atlas_module_PathName_c.f
   atlas_module_SpectralFunctionSpace_i.f
   atlas_module_SpectralFunctionSpace_c.f
   atlas_module_State_i.f
   atlas_module_State_c.f
   atlas_module_Trans_i.f
   atlas_module_Trans_c.f
   atlas_module_Value_i.f
   atlas_module_Value_c.f
)

install( FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR}/atlas_module.mod
                 DESTINATION ${INSTALL_INCLUDE_DIR} )

install( FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR}/atlas_mpi_module.mod
                 DESTINATION ${INSTALL_INCLUDE_DIR} )

#ecbuild_add_library( TARGET atlas_mpl_fortran
#    AUTO_VERSION
#    CONDITION CMAKE_Fortran_COMPILER_LOADED
#    SOURCES
#        atlas_mpi_module.F90
#    LIBS ${ATLAS_MPI_FORTRAN_LIBRARIES}
#)
