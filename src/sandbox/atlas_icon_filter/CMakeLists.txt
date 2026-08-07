# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FORTRAN )

  ecbuild_find_package( NetCDF COMPONENTS Fortran )
  if( NOT TARGET NetCDF::NetCDF_Fortran )
    message( WARNING "NetCDF Fortran library not found, atlas-icon-filter will not be built" )
    return()
  endif()

  ecbuild_add_executable( TARGET atlas-icon-filter
    LINKER_LANGUAGE Fortran
    SOURCES atlas-icon-filter.F90
    LIBS atlas_f NetCDF::NetCDF_Fortran
    NOINSTALL
  )

endif()