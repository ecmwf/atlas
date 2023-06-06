# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if( CMAKE_CXX_COMPILER_ID MATCHES Cray )


  ecbuild_add_cxx_flags("-hnomessage=3140" NAME atlas_cxx_disable_warnings ) # colon separated numbers
  ecbuild_add_fortran_flags("-hnomessage=3140" NAME atlas_fortran_disable_warnings ) # colon separated numbers

# CC-3140 crayc++: WARNING File = atlas/functionspace/NodeColumns.cc, Line = 1, Column = 1
#  The IPA optimization level was changed to "1" due to the presence of OMP
#          directives, ACC directives, or ASM intrinsics.

endif()

if( CMAKE_CXX_COMPILER_ID MATCHES NVHPC )
  ecbuild_add_cxx_flags("--diag_suppress declared_but_not_referenced --display_error_number" NAME atlas_cxx_disable_warnings )
  # For all the variables with side effects (constructor/dectructor functionality)
endif()