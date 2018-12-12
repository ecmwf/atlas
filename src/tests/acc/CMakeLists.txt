# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


if( ATLAS_HAVE_GRIDTOOLS_STORAGE AND GRIDTOOLS_HAVE_CUDA AND ENABLE_TESTS AND HAVE_FCTEST AND ATLAS_HAVE_ACC )

  set_source_files_properties( fctest_unified_memory_with_openacc.F90 PROPERTIES COMPILE_FLAGS "-acc -ta=tesla,nordc" )

  add_fctest( 
    TARGET atlas_test_unified_memory_with_openacc
    SOURCES fctest_unified_memory_with_openacc.F90
            fctest_unified_memory_with_openacc_cxx.cc
    LIBS atlas_f
    LINKER_LANGUAGE Fortran
  )

  target_link_libraries( atlas_test_unified_memory_with_openacc "${ACC_Fortran_FLAGS}" )

endif()