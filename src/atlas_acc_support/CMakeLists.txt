
set( ATLAS_HAVE_ACC  0 )
set( ENABLE_ACC ON CACHE STRING "Enable OpenACC" )

if( ENABLE_ACC AND NOT CUDA_FOUND )

  ecbuild_warn( "OpenACC support not enabled as CUDA was not found" )

elseif( ENABLE_ACC )

  if( CMAKE_Fortran_COMPILER_ID MATCHES "PGI" )
    set( ACC_Fortran_FLAGS -acc -ta=tesla,cuda${CUDA_VERSION},nordc )
    set( ACC_C_FLAGS ${ACC_Fortran_FLAGS} )
    set( ACC_C_COMPILER pgcc )
    set( ATLAS_HAVE_ACC 1 )
  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
    ecbuild_warn( "OpenACC support not enabled as we did not yet implement it for Cray" )
  else()
    ecbuild_warn( "OpenACC support not enabled as we did not yet implement it for ${CMAKE_Fortran_COMPILER_ID}" )
  endif()
  if( ATLAS_HAVE_ACC )
    ecbuild_info( "OpenACC enabled" )
  endif()

endif()

if( ATLAS_HAVE_ACC )

  if( NOT (CMAKE_C_COMPILER_ID MATCHES ${CMAKE_Fortran_COMPILER_ID}) )

    foreach( _dir ${ACC_C_INCLUDE_DIRS} )
      set( ACC_C_INCLUDE ${ACC_C_INCLUDE} -I${_dir} )
    endforeach()
    add_custom_command(
        OUTPUT  ${CMAKE_BINARY_DIR}/lib/libatlas_acc_support.so ${CMAKE_CURRENT_BINARY_DIR}/atlas_acc_map_data.c.o
        COMMAND ${ACC_C_COMPILER} ${ACC_C_FLAGS} ${ACC_C_INCLUDE} -fPIC -o ${CMAKE_CURRENT_BINARY_DIR}/atlas_acc_map_data.c.o
                -c ${CMAKE_CURRENT_SOURCE_DIR}/atlas_acc_map_data.c
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/lib
        COMMAND ${ACC_C_COMPILER} ${ACC_C_FLAGS} -shared -o ${CMAKE_BINARY_DIR}/lib/libatlas_acc_support.so
                ${CMAKE_CURRENT_BINARY_DIR}/atlas_acc_map_data.c.o
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/atlas_acc_map_data.c
        COMMENT "Building atlas_acc_support with ${ACC_C_COMPILER}"
    )
    add_custom_target( build-atlas_acc_support ALL DEPENDS ${CMAKE_BINARY_DIR}/lib/libatlas_acc_support.so )
    add_library( atlas_acc_support SHARED IMPORTED GLOBAL )
    set_property( TARGET atlas_acc_support PROPERTY IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib/libatlas_acc_support.so )
    set_property( TARGET atlas_acc_support PROPERTY IMPORTED_NO_SONAME TRUE )
    set_property( TARGET atlas_acc_support PROPERTY IMPORTED_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR} )
    add_dependencies( atlas_acc_support build-atlas_acc_support )
    install( FILES ${CMAKE_BINARY_DIR}/lib/libatlas_acc_support.so DESTINATION lib/ )

  else()

    ecbuild_add_library( TARGET atlas_acc_support SOURCES atlas_acc_map_data.c )
    target_compile_options( atlas_acc_support PRIVATE ${ACC_C_FLAGS} )
    target_link_libraries( atlas_acc_support PRIVATE ${ACC_C_FLAGS} )

  endif()

endif()

set( ATLAS_HAVE_ACC ${ATLAS_HAVE_ACC} PARENT_SCOPE )
