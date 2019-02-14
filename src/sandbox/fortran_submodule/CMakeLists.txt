# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( CMAKE_Fortran_COMPILER_LOADED )

ecbuild_add_library( TARGET atlas_sandbox_fortran_submodule
  CONDITION OFF
  SOURCES sb_module.F90
)

endif()
