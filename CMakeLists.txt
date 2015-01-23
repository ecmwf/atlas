
# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

############################################################################################
# ATLAS
#
# TODO:
#	* add renumbering to improve cache performance
#	* add detection of OpenMP flags (mayeb move to ecbuild)

cmake_minimum_required( VERSION 2.8.4 FATAL_ERROR )

project( atlas CXX )

set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild/cmake")
set( ECBUILD_DEFAULT_BUILD_TYPE Release )

include( ecbuild_system )

ecbuild_requires_macro_version( 1.6 )

###############################################################################
# project

set( PERSISTENT_NAMESPACE "eckit" CACHE INTERNAL "" )

ecbuild_use_package( PROJECT eckit VERSION 0.8 REQUIRED )

if( ECKIT_HAVE_GRIB )
   message( STATUS "eckit has GRIB support")
endif()

ecbuild_declare_project()

# options & dependencies

ecbuild_find_python()

### Fortran ...

ecbuild_add_option( FEATURE FORTRAN
                    DESCRIPTION "provide Fortran bindings" )

if( HAVE_FORTRAN )

  ecbuild_enable_fortran( REQUIRED MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module )

  ecbuild_add_option( FEATURE FCTEST
                      DEFAULT ${ENABLE_TESTS}
                      DESCRIPTION "Fortran Unit Test Framework"
                      REQUIRED_PACKAGES "PROJECT fctest VERSION 0.1" )
endif()

### OMP ...

ecbuild_add_option(
	FEATURE OMP
	DESCRIPTION "support for OpenMP threaded parallelism" )
if( ENABLE_OMP )
    ecbuild_enable_omp()
else()
    ecbuild_enable_ompstubs()
endif()

### MPI ...

ecbuild_add_option( FEATURE MPI
                    DESCRIPTION "support for MPI distributed parallelism" )

if( ENABLE_MPI )

  ecbuild_enable_mpi()

  if( NOT MPI_CXX_FOUND )
    set( HAVE_MPI 0 )
  endif()

endif()

### trans_api ...

ecbuild_add_option(
	FEATURE TRANS
	DESCRIPTION "Support for trasnformations between SH and Grid space"
	REQUIRED_PACKAGES "PROJECT trans_api VERSION 0.1" )

### tesselation ...

ecbuild_add_option(
	FEATURE TESSELATION
	DESCRIPTION "support for unstructured mesh generation"
	REQUIRED_PACKAGES "CGAL QUIET" )

if( HAVE_TESSELATION )

	# CGAL needs boost

	set(Boost_USE_MULTITHREADED      ON )
	find_package( Boost 1.45.0 REQUIRED )
	set_package_properties( Boost PROPERTIES TYPE RECOMMENDED PURPOSE "CGAL support" )
	list( APPEND CGAL_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} )
	list( APPEND CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )

endif()

### Type for Global indices and unique point id's

set( ATLAS_BITS_GLOBAL 64 )

### sandbox

ecbuild_add_option( FEATURE SANDBOX
                    DEFAULT OFF
                    DESCRIPTION "build the sandbox stuff" )

############################################################################################
# export package info

list( APPEND ATLAS_TPLS      eckit )

set( ATLAS_INCLUDE_DIRS   ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src )

list( APPEND ATLAS_LIBRARIES  atlas )

if( HAVE_FORTRAN )
  list( APPEND ATLAS_INCLUDE_DIRS ${CMAKE_Fortran_MODULE_DIRECTORY} )
  set( ATLAS_LIBRARIES atlas_f ${ATLAS_LIBRARIES} )
endif()

if( ATLAS_HAVE_TRANS )
	list( APPEND ATLAS_TPLS trans_api )
endif()

if( HAVE_TESSELATION )
  list( APPEND ATLAS_TPLS CGAL )
endif()

get_directory_property( ATLAS_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${ATLAS_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND ATLAS_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND ATLAS_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND ATLAS_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()


######################################################################################
# Find  ATLAS_INCLUDE and ATLAS_LINK CMake variables, to configure in atlas_config.h

function( ecbuild_library_dependencies dependencies libraries )
  ### TODO:  MOVE THIS FUNCTION TO ECBUILD

  set( _libraries ${${libraries}} )

  foreach( _lib ${_libraries})

    unset( _location )

    if( TARGET ${_lib} ) # check if this is an existing target

      get_property( _location TARGET ${_lib} PROPERTY LOCATION )

      get_property( _config TARGET ${_lib} PROPERTY IMPORTED_CONFIGURATIONS )
      get_property( _deps TARGET ${_lib} PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES_${_config} )

      ecbuild_library_dependencies( _deps_location _deps )
      list( APPEND _location ${_deps_location} )

    else()

      set( _location ${_lib} )

    endif()

    list( APPEND _dependencies ${_location} )

  endforeach()

  list( REVERSE           _dependencies )
  list( REMOVE_DUPLICATES _dependencies )
  list( REVERSE           _dependencies )
  
  set( ${dependencies} ${_dependencies} PARENT_SCOPE )

endfunction(ecbuild_library_dependencies)

function( ecbuild_pkg_cfg_libs pkg_cfg_libs libraries )

  ### TODO:  MOVE THIS FUNCTION TO ECBUILD

  set( _libraries ${${libraries}} )

  foreach( _lib ${_libraries})

    if( ${_lib} MATCHES "-.+" )
      list( APPEND _pkg_cfg_libs ${_lib} )
    else() 

      unset( _name )
      unset( _dir  )
   

      if( TARGET ${_lib} ) # check if this is an existing target

        get_property( _location TARGET ${_lib} PROPERTY LOCATION )

        get_property( _config TARGET ${_lib} PROPERTY IMPORTED_CONFIGURATIONS )
        get_property( _deps TARGET ${_lib} PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES_${_config} )


      else()

        set( _location ${_lib} )

      endif()

      get_filename_component( _name ${_location} NAME_WE )
      get_filename_component( _dir  ${_location} PATH )

      if( NOT _name )
        set( _name ${_lib} )
      endif()

      if( NOT _dir )
        set( _dir  "${CMAKE_INSTALL_PREFIX}/lib" )
      endif()

      string( REGEX REPLACE "^lib" "" _name ${_name} )
      
      list( APPEND _pkg_cfg_libs "-L${_dir}" "-l${_name}" )

    endif()

  endforeach()

  list( REMOVE_DUPLICATES _pkg_cfg_libs )
  string( REPLACE ";" " " _pkg_cfg_libs "${_pkg_cfg_libs}" )
  
  set( ${pkg_cfg_libs} ${_pkg_cfg_libs} PARENT_SCOPE )

endfunction(ecbuild_pkg_cfg_libs)

function( ecbuild_pkg_cfg_include INCLUDE INCLUDE_DIRS )

  ### TODO:  MOVE THIS FUNCTION TO ECBUILD

  foreach( _INCDIR ${${INCLUDE_DIRS}} )
    list( APPEND _INCLUDE "-I${_INCDIR}")
  endforeach()
  list( REMOVE_DUPLICATES _INCLUDE)
  string( REPLACE ";" " " _INCLUDE "${_INCLUDE}")

  set( ${INCLUDE} ${_INCLUDE} PARENT_SCOPE )

endfunction(ecbuild_pkg_cfg_include)


############################################################################################
# pkg-config

ecbuild_library_dependencies( ATLAS_DEPS ATLAS_EXTRA_LIBRARIES )
list( APPEND ATLAS_PKG_CFG_LIBS ${ATLAS_LIBRARIES} ${ATLAS_DEPS} )
ecbuild_pkg_cfg_libs( ATLAS_PKG_CFG_LIBS ATLAS_PKG_CFG_LIBS )

ecbuild_pkg_cfg_include( ATLAS_PKG_CFG_CFLAGS ATLAS_EXTRA_INCLUDE_DIRS )
set( ATLAS_PKG_CFG_CFLAGS "-I${CMAKE_INSTALL_PREFIX}/include ${ATLAS_PKG_CFG_CFLAGS}" )
set( ATLAS_PKG_CFG_FFLAGS "-I${CMAKE_INSTALL_PREFIX}/include" )

debug_var( ATLAS_PKG_CFG_LIBS   )
debug_var( ATLAS_PKG_CFG_CFLAGS )
debug_var( ATLAS_PKG_CFG_FFLAGS )

############################################################################################
# sources

include( atlas_compiler_flags )

include_directories( ${ATLAS_INCLUDE_DIRS} ${ATLAS_EXTRA_INCLUDE_DIRS} )

add_subdirectory( src )



############################################################################################
# finalize

add_custom_target( atlas-readme SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/README.md )
ecbuild_install_project( NAME Atlas )

ecbuild_print_summary()
