# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( atlas_HAVE_ACC )

  if( NOT (CMAKE_C_COMPILER_ID MATCHES ${CMAKE_Fortran_COMPILER_ID}) )
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
    install( FILES ${CMAKE_BINARY_DIR}/lib/libatlas_acc_support.so DESTINATION ${INSTALL_LIB_DIR}/ )

  else()

    ecbuild_add_library( TARGET atlas_acc_support SOURCES atlas_acc_map_data.c )
    target_compile_options( atlas_acc_support PRIVATE ${ACC_C_FLAGS} )
    target_link_libraries( atlas_acc_support PRIVATE ${ACC_C_FLAGS} )
  endif()

endif()
