# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

ecbuild_add_option( FEATURE WARNING_AS_ERROR
                    DEFAULT OFF
                    DESCRIPTION "Treat compile warning as error" )

if(HAVE_WARNING_AS_ERROR)
  ecbuild_add_cxx_flags("-Werror" NO_FAIL NAME atlas_cxx_warning_as_error)
endif()

ecbuild_add_option( FEATURE WARNINGS
                    DEFAULT ON
                    DESCRIPTION "Add warnings to compiler" )

# activate warnings, ecbuild macros check the compiler recognises the options
if(HAVE_WARNINGS)

  ecbuild_add_cxx_flags("-Wall" NO_FAIL)
  ecbuild_add_cxx_flags("-Wextra" NO_FAIL)

  ecbuild_add_cxx_flags("-Wno-unused-parameter" NO_FAIL)
  ecbuild_add_cxx_flags("-Wno-sign-compare" NO_FAIL)

endif()

if( CMAKE_CXX_COMPILER_ID STREQUAL Intel )
  ecbuild_add_cxx_flags("-diag-disable=11074" NO_FAIL)   # Inline limits
  ecbuild_add_cxx_flags("-diag-disable=11076" NO_FAIL)   # Inline limits
  ecbuild_add_cxx_flags("-diag-disable=10441" NO_FAIL)   # Deprecated classic compiler
endif()

if( CMAKE_Fortran_COMPILER_ID MATCHES Intel) # Both for Intel and Intel-LLVM
  ecbuild_add_fortran_flags("-diag-disable=5462" NO_FAIL)   # Global name too long, shortened
endif()

if( CMAKE_CXX_COMPILER_ID MATCHES Cray )

  if( NOT CMAKE_CXX_COMPILER_ID MATCHES CrayClang )
    ecbuild_add_cxx_flags("-hnomessage=3140" NAME atlas_cxx_disable_warnings ) # colon separated numbers
  endif()
  ecbuild_add_fortran_flags("-hnomessage=3140" NAME atlas_fortran_disable_warnings ) # colon separated numbers

# CC-3140 crayc++: WARNING File = atlas/functionspace/NodeColumns.cc, Line = 1, Column = 1
#  The IPA optimization level was changed to "1" due to the presence of OMP
#          directives, ACC directives, or ASM intrinsics.

endif()

if( CMAKE_CXX_COMPILER_ID MATCHES NVHPC )
  ecbuild_add_cxx_flags("--diag_suppress declared_but_not_referenced --display_error_number" NAME atlas_cxx_disable_warnings )
  # For all the variables with side effects (constructor/destructor functionality)
endif()

if( CMAKE_CXX_COMPILER_ID MATCHES IntelLLVM )
  # Turn off -ffinite-math-only which gets included by some optimisation levels which assumes values can never be NaN.
  # Then results in std::isnan(value) always return false.
  ecbuild_add_cxx_flags("-fno-finite-math-only")
endif()

