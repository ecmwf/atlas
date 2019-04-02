# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( ATLAS_HAVE_OMP_CXX )
  set( ATLAS_HAVE_OMP_CXX 1 )
else()
  set( ATLAS_HAVE_OMP_CXX 0 )
endif()

if( ATLAS_HAVE_OMP_Fortran )
  set( ATLAS_HAVE_OMP_Fortran 1 )
else()
  set( ATLAS_HAVE_OMP_Fortran 0 )
endif()

if( ATLAS_HAVE_TESSELATION )
  set( ATLAS_HAVE_TESSELATION 1 )
else()
  set( ATLAS_HAVE_TESSELATION 0 )
endif()

if( ATLAS_HAVE_FORTRAN )
  set( ATLAS_HAVE_FORTRAN 1 )
else()
  set( ATLAS_HAVE_FORTRAN 0 )
endif()

if( ATLAS_HAVE_TRANS )
  set( ATLAS_HAVE_TRANS 1 )
else()
  set( ATLAS_HAVE_TRANS 0 )
endif()

if( ATLAS_HAVE_FFTW )
  set( ATLAS_HAVE_FFTW 1 )
else()
  set( ATLAS_HAVE_FFTW 0 )
endif()

if( ATLAS_HAVE_BOUNDSCHECKING )
  set( ATLAS_HAVE_BOUNDSCHECKING 1 )
else()
  set( ATLAS_HAVE_BOUNDSCHECKING 0 )
endif()

if( ATLAS_HAVE_INIT_SNAN )
  set( ATLAS_HAVE_INIT_SNAN 1 )
else()
  set( ATLAS_HAVE_INIT_SNAN 0 )
endif()

if( ATLAS_HAVE_GRIDTOOLS_STORAGE )
  set( ATLAS_HAVE_GRIDTOOLS_STORAGE 1 )
else()
  set( ATLAS_HAVE_GRIDTOOLS_STORAGE 0 )
endif()

add_subdirectory( atlas_acc_support )

add_subdirectory( atlas )

if( ATLAS_HAVE_FORTRAN )
 add_subdirectory( atlas_f )
endif()

add_subdirectory( apps )
add_subdirectory( tests )

if( HAVE_SANDBOX )
    add_subdirectory( sandbox )
endif()
