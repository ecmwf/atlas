# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FORTRAN )

  ecbuild_add_executable(
    TARGET  atlas-filter-fortran
    SOURCES atlas-filter.F90
    LIBS    atlas_f
    NOINSTALL
  )

  ecbuild_add_executable(
    TARGET  atlas-interpolation-fortran
    SOURCES atlas-interpolation-fortran.F90
    LIBS    atlas_f
    LINKER_LANGUAGE Fortran
    NOINSTALL
  )

endif()
