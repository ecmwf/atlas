# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_test( TARGET atlas_test_trace
  SOURCES     test_trace.cc
  LIBS        atlas ${OMP_CXX}
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_TRACE_REPORT=1 
  OMP         2
)

foreach( test test_library test_library_noargs test_library_init_nofinal test_library_noinit_final )
  ecbuild_add_test( TARGET atlas_${test}
    SOURCES     ${test}.cc
    LIBS        atlas
    ENVIRONMENT ATLAS_TRACE_REPORT=1 ATLAS_DEBUG=1
  )
endforeach()

ecbuild_add_test( TARGET atlas_test_memory
  SOURCES     test_memory.cc
  LIBS        atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_TRACE_REPORT=1
)


if( HAVE_FCTEST )

add_fctest( TARGET atlas_fctest_trace
  SOURCES     fctest_trace.fypp
  LINKER_LANGUAGE Fortran
  LIBS        atlas_f ${OMP_Fortran}
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_TRACE_REPORT=1
)

endif()
