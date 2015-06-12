# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


if( ATLAS_HAVE_MPI )

  set( ATLAS_MPI_C_LIBRARIES            ${MPI_CXX_LIBRARIES} )
  set( ATLAS_MPI_Fortran_LIBRARIES      ${MPI_Fortran_LIBRARIES} )
  ecbuild_include_mpi()

endif()

add_subdirectory( atlas )

if( ATLAS_HAVE_FORTRAN )
  add_subdirectory( atlas_f )
endif()

add_subdirectory( apps )

add_subdirectory( tests )

if( HAVE_SANDBOX )
    add_subdirectory( sandbox )
endif()
