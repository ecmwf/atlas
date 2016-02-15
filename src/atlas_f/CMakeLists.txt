
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

#function(generate_fortran_bindings output filename)
#    get_filename_component(base ${filename} NAME_WE)
#    set(base_abs ${CMAKE_CURRENT_SOURCE_DIR}/${base})
#    set(outfile ${CMAKE_CURRENT_BINARY_DIR}/${base}_c_binding.f90)
#    set(${output} ${${output}} ${outfile} PARENT_SCOPE)
#    add_custom_command(
#        OUTPUT ${outfile}
#        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/c2f.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${outfile}
#        DEPENDS ${filename} )
#    set_source_files_properties(${outfile} PROPERTIES GENERATED TRUE)
#endfunction()

#function(generate_fortran_bindings_custom_output output filename outfile)
#    get_filename_component(base ${filename} NAME_WE)
#    set(base_abs ${CMAKE_CURRENT_SOURCE_DIR}/${base})
#    set(outfile ${CMAKE_CURRENT_BINARY_DIR}/${outfile})
#    set(${output} ${${output}} ${outfile} PARENT_SCOPE)
#    add_custom_command(
#        OUTPUT ${outfile}
#        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/c2f.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${outfile}
#        DEPENDS ${filename} )
#    set_source_files_properties(${outfile} PROPERTIES GENERATED TRUE)
#endfunction()

function(generate_fortran_bindings output filename)

  set( options "" )
  set( single_value_args OUTPUT MODULE )
  set( multi_value_args "" )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )

  get_filename_component(base ${filename} NAME_WE)
  set(base_abs ${CMAKE_CURRENT_SOURCE_DIR}/${base})
  set(outfile ${CMAKE_CURRENT_BINARY_DIR}/${base}_c_binding.f90)

  if( _PAR_OUTPUT )
    set(outfile ${_PAR_OUTPUT})
  endif()
  set(${output} ${${output}} ${outfile} PARENT_SCOPE)

  if( _PAR_MODULE )
    add_custom_command(
        OUTPUT ${outfile}
        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/c2f.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${outfile} -m ${_PAR_MODULE}
        DEPENDS ${filename} )
  else()
    add_custom_command(
        OUTPUT ${outfile}
        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/c2f.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${outfile}
        DEPENDS ${filename} )
  endif()
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
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/Nodes.h)
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
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/functionspace/EdgeBasedFiniteVolume.h MODULE atlas_functionspace_EdgeBasedFiniteVolume_c_binding OUTPUT functionspace_EdgeBasedFiniteVolume_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/functionspace/Spectral.h MODULE atlas_functionspace_Spectral_c_binding OUTPUT functionspace_Spectral_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/functionspace/ReducedGridPoint.h MODULE atlas_functionspace_ReducedGridPoint_c_binding OUTPUT functionspace_ReducedGridPoint_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/functionspace/NodesInterface.h MODULE atlas_functionspace_Nodes_c_binding OUTPUT functionspace_Nodes_c_binding.f90)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/numerics/Nabla.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/meshgen/MeshGenerator.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/Connectivity.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/HybridElements.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/Elements.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/ElementType.h)
generate_fortran_bindings(FORTRAN_BINDINGS atlas_logging.h)
generate_fortran_bindings(FORTRAN_BINDINGS atlas_resource.h)
generate_fortran_bindings(FORTRAN_BINDINGS atlas_value.h)

if( ATLAS_HAVE_DEPRECATED_FUNCTIONSPACE )
set( DEPRECATED_FUNCTIONSPACE_SRC atlas_deprecated_functionspace_module.F90 )
endif()

### atlas fortran lib
ecbuild_add_library( TARGET atlas_f
    AUTO_VERSION
    CONDITION CMAKE_Fortran_COMPILER_LOADED
    SOURCES
        ${FORTRAN_BINDINGS}
        atlas_object_module.F90
        atlas_refcounted_module.F90
        atlas_JSON_module.F90
        atlas_Config_module.F90
        atlas_Metadata_module.F90
        atlas_FunctionSpace_module.F90
        atlas_Logging_module.F90
        atlas_Error_module.F90
        atlas_Field_module.F90
        atlas_module.F90
        atlas_mpi_module.F90
        atlas_grids_module.F90
        atlas_Connectivity_module.F90
        atlas_mesh_HybridElements_module.F90
        atlas_mesh_Elements_module.F90
        atlas_mesh_ElementType_module.F90
        atlas_GatherScatter_module.F90
        atlas_Checksum_module.F90
        atlas_HaloExchange_module.F90
        atlas_C_interop.F90
        atlas_C_interop.cc
        atlas_logging.h
        atlas_logging.cc
        atlas_resource.h
        atlas_resource.cc
        atlas_value.h
        atlas_value.cc
        ${DEPRECATED_FUNCTIONSPACE_SRC}
  INCLUDES ${MPI_Fortran_INCLUDE_PATH}
  LIBS atlas ${MPI_Fortran_LIBRARIES}
)

add_custom_target( atlas_fortran_inc SOURCES
   atlas_module_FieldSet_i.f
   atlas_module_FieldSet_c.f
   atlas_module_functionspace_EdgeBasedFiniteVolume_i.f
   atlas_module_functionspace_EdgeBasedFiniteVolume_c.f
   atlas_module_functionspace_Nodes_i.f
   atlas_module_functionspace_Nodes_c.f
   atlas_module_functionspace_ReducedGridPoint_c.f
   atlas_module_functionspace_ReducedGridPoint_i.f
   atlas_module_functionspace_Spectral_i.f
   atlas_module_functionspace_Spectral_c.f
   atlas_module_Grid_i.f
   atlas_module_Grid_c.f
   atlas_module_GridDistribution_i.f
   atlas_module_GridDistribution_c.f
   atlas_module_Mesh_i.f
   atlas_module_Mesh_c.f
   atlas_module_MeshGenerator_i.f
   atlas_module_MeshGenerator_c.f
   atlas_module_Nabla_i.f
   atlas_module_Nabla_c.f
   atlas_module_Nodes_i.f
   atlas_module_Nodes_c.f
   atlas_module_State_i.f
   atlas_module_State_c.f
   atlas_module_Trans_i.f
   atlas_module_Trans_c.f
   atlas_module_Value_i.f
   atlas_module_Value_c.f
)
