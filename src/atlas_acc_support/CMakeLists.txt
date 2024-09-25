# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( atlas_HAVE_ACC )

  if( CMAKE_CXX_COMPILER_ID MATCHES NVHPC )
    if( NOT TARGET OpenACC::OpenACC_CXX )
       ecbuild_error("ERROR: OpenACC::OpenACC_CXX TARGET not found")
    endif()
    ecbuild_add_library( TARGET atlas_acc_support SOURCES atlas_acc.cc )
    target_link_libraries( atlas_acc_support PRIVATE OpenACC::OpenACC_CXX )
  else()
    ecbuild_add_library( TARGET atlas_acc_support SOURCES atlas_acc.F90 )
    target_link_libraries( atlas_acc_support PRIVATE OpenACC::OpenACC_Fortran )
  endif()

endif()
