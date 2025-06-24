# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( atlas_HAVE_OMP_CXX )
  set( atlas_HAVE_OMP_CXX 1 )
else()
  set( atlas_HAVE_OMP_CXX 0 )
endif()

if( atlas_HAVE_OMP_Fortran )
  set( atlas_HAVE_OMP_Fortran 1 )
else()
  set( atlas_HAVE_OMP_Fortran 0 )
endif()

if( atlas_HAVE_TESSELATION )
  set( atlas_HAVE_TESSELATION 1 )
else()
  set( atlas_HAVE_TESSELATION 0 )
endif()

if( atlas_HAVE_QHULL )
  set( atlas_HAVE_QHULL 1 )
else()
  set( atlas_HAVE_QHULL 0 )
endif()

if( atlas_HAVE_CGAL )
  set( atlas_HAVE_CGAL 1 )
else()
  set( atlas_HAVE_CGAL 0 )
endif()

if( atlas_HAVE_PROJ )
  set( atlas_HAVE_PROJ 1 )
else()
  set( atlas_HAVE_PROJ 0 )
endif()

if( atlas_HAVE_FORTRAN )
  set( atlas_HAVE_FORTRAN 1 )
else()
  set( atlas_HAVE_FORTRAN 0 )
endif()

if( atlas_HAVE_ECTRANS )
  set( atlas_HAVE_ECTRANS 1 )
  set( atlas_HAVE_ECTRANS_LAM 0 )
  if ( ectrans_HAVE_ETRANS )
    set( atlas_HAVE_ECTRANS_LAM 1 )
  endif()
else()
  set( atlas_HAVE_ECTRANS 0 )
  set( atlas_HAVE_ECTRANS_LAM 0 )
endif()

if( atlas_HAVE_PACKAGE_ECTRANS )
  set( atlas_HAVE_PACKAGE_ECTRANS 1 )
else()
  set( atlas_HAVE_PACKAGE_ECTRANS 0 )
endif()

if( atlas_HAVE_FFTW )
  set( atlas_HAVE_FFTW 1 )
else()
  set( atlas_HAVE_FFTW 0 )
endif()

if( atlas_HAVE_BOUNDSCHECKING )
  set( atlas_HAVE_BOUNDSCHECKING 1 )
else()
  set( atlas_HAVE_BOUNDSCHECKING 0 )
endif()

if( atlas_HAVE_INIT_SNAN )
  set( atlas_HAVE_INIT_SNAN 1 )
else()
  set( atlas_HAVE_INIT_SNAN 0 )
endif()

if( atlas_HAVE_CUDA )
  set( atlas_HAVE_CUDA 1 )
else()
  set( atlas_HAVE_CUDA 0 )
endif()

if( atlas_HAVE_ACC )
  set( atlas_HAVE_ACC 1 )
else()
  set( atlas_HAVE_ACC 0 )
endif()

if( atlas_HAVE_GRIDTOOLS_STORAGE )
  set( atlas_HAVE_GRIDTOOLS_STORAGE 1 )
else()
  set( atlas_HAVE_GRIDTOOLS_STORAGE 0 )
endif()

set( atlas_BUILD_TYPE_RELEASE 0 )
set( atlas_BUILD_TYPE_DEBUG   0 )
if( CMAKE_BUILD_TYPE MATCHES "Release" )
  set( atlas_BUILD_TYPE_RELEASE 1 )
endif()
if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  set( atlas_BUILD_TYPE_DEBUG 1 )
endif()

if( atlas_HAVE_ECKIT_DEVELOP )
  set( ATLAS_ECKIT_DEVELOP 1 )
else()
  set( ATLAS_ECKIT_DEVELOP 0 )
endif()

ecbuild_parse_version( ${eckit_VERSION} PREFIX ATLAS_ECKIT )
math( EXPR ATLAS_ECKIT_VERSION_INT "( 10000 * ${ATLAS_ECKIT_VERSION_MAJOR} ) + ( 100 * ${ATLAS_ECKIT_VERSION_MINOR} ) + ${ATLAS_ECKIT_VERSION_PATCH}" )

add_subdirectory( atlas_acc_support )

add_subdirectory( atlas )

if( atlas_HAVE_FORTRAN )
  if( atlas_HAVE_FCKIT_DEVELOP )
    set( ATLAS_FCKIT_DEVELOP 1 )
  else()
    set( ATLAS_FCKIT_DEVELOP 0 )
  endif()
  ecbuild_parse_version( ${fckit_VERSION} PREFIX ATLAS_FCKIT )
  math( EXPR ATLAS_FCKIT_VERSION_INT "( 10000 * ${ATLAS_FCKIT_VERSION_MAJOR} ) + ( 100 * ${ATLAS_FCKIT_VERSION_MINOR} ) + ${ATLAS_FCKIT_VERSION_PATCH}" )
 add_subdirectory( atlas_f )
endif()

add_subdirectory( apps )
add_subdirectory( tests )

if( HAVE_SANDBOX )
    add_subdirectory( sandbox )
endif()
