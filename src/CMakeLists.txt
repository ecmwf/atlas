# (C) Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( CGAL_FOUND )
  set( CGAL_FOUND 1 )
else()
  set( CGAL_FOUND 0 )
endif()

if( ATLAS_HAVE_OMP )
  set( ATLAS_HAVE_OMP 1 )
else()
  set( ATLAS_HAVE_OMP 0 )
endif()

if( ATLAS_HAVE_MPI )
  set( ATLAS_HAVE_MPI 1 )
else()
  set( ATLAS_HAVE_MPI 0 )
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

if( ATLAS_HAVE_BOUNDSCHECKING )
  set( ATLAS_HAVE_BOUNDSCHECKING 1 )
else()
  set( ATLAS_HAVE_BOUNDSCHECKING 0 )
endif()

if( ATLAS_HAVE_MPI )

  set( ATLAS_MPI_C_LIBRARIES            ${MPI_CXX_LIBRARIES} )
  set( ATLAS_MPI_Fortran_LIBRARIES      ${MPI_Fortran_LIBRARIES} )
  ecbuild_include_mpi()

endif()

add_subdirectory( atlas )

#if( ATLAS_HAVE_FORTRAN )
#  add_subdirectory( atlas_f )
#endif()

#add_subdirectory( apps )
#add_subdirectory( tests )

#if( HAVE_SANDBOX )
#    add_subdirectory( sandbox )
#endif()
