# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/atlas_f.h.in
                ${CMAKE_CURRENT_BINARY_DIR}/../atlas/atlas_f.h @ONLY )

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/atlas_f.fypp ${CMAKE_CURRENT_BINARY_DIR}/../atlas/atlas_f.fypp)

install( FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/atlas_f.fypp
  ${CMAKE_CURRENT_BINARY_DIR}/../atlas/atlas_f.h
  DESTINATION
  ${INSTALL_INCLUDE_DIR}/atlas )


### fortran bindings

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

  if( NOT _PAR_MODULE )
    set( _PAR_MODULE "atlas_${base}_c_binding" )
  endif()

  set( FINT32 "integer(c_int)"  )
  set( FINT64 "integer(c_long)" )
  set( F_GIDX ${FINT${ATLAS_BITS_GLOBAL}} )
  set( F_IDX ${FINT${ATLAS_BITS_LOCAL}} )

  add_custom_command(
    OUTPUT ${outfile}
    COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/c2f.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename}
            -o ${outfile} -m ${_PAR_MODULE}
            -t '{"idx_t":"${F_IDX}","gidx_t":"${F_GIDX}"}'
    DEPENDS ${filename} )
  set_source_files_properties(${outfile} PROPERTIES GENERATED TRUE)
endfunction()

generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/grid.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/grid/detail/grid/Structured.h
    MODULE    atlas_grid_Structured_c_binding
    OUTPUT    grid_Structured_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/grid/detail/distribution/DistributionImpl.h
    MODULE atlas_distribution_c_binding
    OUTPUT distribution_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/grid/detail/vertical/VerticalInterface.h
    MODULE atlas_vertical_c_binding
    OUTPUT vertical_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/grid/Partitioner.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/detail/MeshIntf.h
    MODULE    atlas_mesh_c_binding
    OUTPUT    atlas_mesh_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/Nodes.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/Connectivity.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/HybridElements.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/Elements.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/ElementType.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/actions/BuildParallelFields.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/actions/BuildPeriodicBoundaries.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/actions/BuildHalo.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/actions/BuildEdges.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/actions/BuildDualMesh.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/mesh/actions/WriteLoadBalanceReport.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/meshgenerator/detail/MeshGeneratorInterface.h
    MODULE atlas_meshgenerator_c_binding
    OUTPUT meshgenerator_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/output/Output.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/output/Gmsh.h
  MODULE    atlas_output_gmsh_c_binding
  OUTPUT    outout_Gmsh_c_binding.f90)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/field/State.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/field/detail/FieldInterface.h
    MODULE    atlas_field_c_binding
    OUTPUT    field_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/field/FieldSet.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/functionspace/detail/FunctionSpaceInterface.h
    MODULE atlas_functionspace_c_binding
    OUTPUT functionspace_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/functionspace/detail/SpectralInterface.h
    MODULE    atlas_functionspace_Spectral_c_binding
    OUTPUT    functionspace_Spectral_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/functionspace/detail/StructuredColumnsInterface.h
    MODULE    atlas_functionspace_StructuredColumns_c_binding
    OUTPUT    functionspace_StructuredColumns_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/functionspace/detail/NodeColumnsInterface.h
    MODULE    atlas_functionspace_NodeColumns_c_binding
    OUTPUT    functionspace_NodeColumns_c_binding.f90)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/functionspace/EdgeColumns.h
    MODULE    atlas_functionspace_EdgeColumns_c_binding
    OUTPUT    functionspace_EdgeColumns_c_binding.f90)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/numerics/Nabla.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/numerics/Nabla.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/numerics/Method.h )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/numerics/fvm/Method.h
    MODULE    atlas_fvm_method_c_binding
    OUTPUT    fvm_method_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/interpolation/Interpolation.h
    MODULE    atlas_interpolation_c_binding
    OUTPUT    interpolation_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/trans/detail/TransInterface.h
    MODULE    atlas_trans_c_binding
    OUTPUT    trans_c_binding.f90 )
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/util/Allocate.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/util/Metadata.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/util/Config.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/output/detail/GmshIO.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/parallel/HaloExchange.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/parallel/GatherScatter.h)
generate_fortran_bindings(FORTRAN_BINDINGS ../atlas/parallel/Checksum.h)
generate_fortran_bindings(FORTRAN_BINDINGS internals/atlas_read_file.h)
generate_fortran_bindings(FORTRAN_BINDINGS internals/Library.h)
generate_fortran_bindings(FORTRAN_BINDINGS runtime/atlas_trace.h MODULE atlas_trace_c_binding OUTPUT atlas_trace_c_binding.f90 )


### atlas fortran lib
ecbuild_add_library( TARGET atlas_f
    AUTO_VERSION
    SOURCES
        ${FORTRAN_BINDINGS}
        atlas_module.F90
        util/atlas_kinds_module.F90
        util/atlas_JSON_module.F90
        util/atlas_Config_module.F90
        util/atlas_Metadata_module.F90
        util/atlas_allocate_module.F90
        output/atlas_output_module.F90
        functionspace/atlas_FunctionSpace_module.F90
        functionspace/atlas_functionspace_EdgeColumns_module.F90
        functionspace/atlas_functionspace_NodeColumns_module.fypp
        functionspace/atlas_functionspace_StructuredColumns_module.F90
        functionspace/atlas_functionspace_Spectral_module.F90
        field/atlas_FieldSet_module.F90
        field/atlas_State_module.F90
        field/atlas_Field_module.fypp
        grid/atlas_Grid_module.F90
        grid/atlas_GridDistribution_module.F90
        grid/atlas_Vertical_module.F90
        grid/atlas_Partitioner_module.F90
        mesh/atlas_MeshGenerator_module.F90
        mesh/atlas_Mesh_module.F90
        mesh/atlas_mesh_Nodes_module.F90
        mesh/atlas_mesh_Edges_module.F90
        mesh/atlas_mesh_Cells_module.F90
        mesh/atlas_Connectivity_module.F90
        mesh/atlas_HybridElements_module.F90
        mesh/atlas_Elements_module.F90
        mesh/atlas_ElementType_module.F90
        mesh/atlas_mesh_actions_module.F90
        numerics/atlas_Method_module.F90
        numerics/atlas_fvm_module.F90
        numerics/atlas_Nabla_module.F90
        interpolation/atlas_Interpolation_module.F90
        parallel/atlas_GatherScatter_module.fypp
        parallel/atlas_Checksum_module.fypp
        parallel/atlas_HaloExchange_module.fypp
        trans/atlas_Trans_module.F90
        internals/atlas_read_file.h
        internals/atlas_read_file.cc
        internals/atlas_write_to_fortran_unit.F90
        internals/Library.h
        internals/Library.cc
        runtime/atlas_trace.cc
        runtime/atlas_Trace_module.F90

  PRIVATE_INCLUDES
        ${FCKIT_INCLUDE_DIRS}
)

target_link_libraries( atlas_f PRIVATE
  atlas )

target_link_libraries( atlas_f PUBLIC
  fckit )

if( ATLAS_HAVE_OMP_Fortran )
  target_link_libraries( atlas_f PRIVATE OpenMP::OpenMP_Fortran )
endif()

target_include_directories( atlas_f PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
  $<INSTALL_INTERFACE:include>
  $<INSTALL_INTERFACE:module/atlas>
)
target_include_directories( atlas_f PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
if( ECKIT_INCLUDE_DIRS )
  target_include_directories( atlas_f PRIVATE ${ECKIT_INCLUDE_DIRS} )
endif()

fckit_target_preprocess_fypp( atlas_f
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/atlas_f.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/atlas_f.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/internals/atlas_generics.fypp )
