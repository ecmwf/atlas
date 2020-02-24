# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_executable(
  TARGET      atlas_c-meshes-Structured
  SOURCES     meshes-Structured.cc
  LIBS        atlas
  NOINSTALL )

set_target_properties(atlas_c-meshes-Structured
                      PROPERTIES EXCLUDE_FROM_ALL ${ENABLE_TESTS} )

if( atlas_HAVE_FORTRAN )

ecbuild_add_executable(
  TARGET      atlas_f-meshes-Structured
  SOURCES     meshes-Structured.F90
  LIBS        atlas_f
  LINKER_LANGUAGE Fortran
  NOINSTALL )

set_target_properties(atlas_f-meshes-Structured
                      PROPERTIES EXCLUDE_FROM_ALL ${ENABLE_TESTS} )

endif()
