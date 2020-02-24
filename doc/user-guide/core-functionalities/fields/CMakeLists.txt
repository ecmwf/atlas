# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_executable(
  TARGET      atlas_c-fields
  SOURCES     fields.cc
  LIBS        atlas
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_c-fields-on-grid
  SOURCES     fields-on-grid.cc
  LIBS        atlas
  NOINSTALL )

set_target_properties(atlas_c-fields
                      atlas_c-fields-on-grid
                      PROPERTIES EXCLUDE_FROM_ALL ${ENABLE_TESTS} )

if( atlas_HAVE_FORTRAN )

ecbuild_add_executable(
  TARGET      atlas_f-fields
  SOURCES     fields.F90
  LIBS        atlas_f
  LINKER_LANGUAGE Fortran
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_f-fields-on-grid
  SOURCES     fields-on-grid.F90
  LIBS        atlas_f
  LINKER_LANGUAGE Fortran
  NOINSTALL )

set_target_properties(atlas_f-fields
                      atlas_f-fields-on-grid
                      PROPERTIES EXCLUDE_FROM_ALL ${ENABLE_TESTS} )

endif()
