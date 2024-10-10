# (C) Copyright 2024- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

add_subdirectory(test_find_pluto)
add_subdirectory(benchmark)
add_subdirectory(sandbox)

ecbuild_add_test( TARGET pluto_test_memory_pool SOURCES pluto_test_memory_pool.cc LIBS pluto )
ecbuild_add_test( TARGET pluto_test_alignment   SOURCES pluto_test_alignment.cc   LIBS pluto )

if( HAVE_FORTRAN )
  ecbuild_add_test( TARGET pluto_test_fortran_memory_resource        SOURCES fortran/test_fortran_memory_resource.F90       LIBS pluto_f )
  ecbuild_add_test( TARGET pluto_test_fortran_allocator              SOURCES fortran/test_fortran_allocator.F90             LIBS pluto_f )
  ecbuild_add_test( TARGET pluto_test_fortran_memory_pool_resource   SOURCES fortran/test_fortran_memory_pool_resource.F90  LIBS pluto_f )
  ecbuild_add_test( TARGET pluto_test_fortran_scope                  SOURCES fortran/test_fortran_scope.F90                 LIBS pluto_f )
endif()

