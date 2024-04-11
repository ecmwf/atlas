# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


if( atlas_HAVE_GRIDTOOLS_STORAGE AND GRIDTOOLS_HAS_BACKEND_CUDA AND HAVE_TESTS AND HAVE_FCTEST AND atlas_HAVE_ACC )

  string (REPLACE ";" " " ACC_Fortran_FLAGS_STR "${ACC_Fortran_FLAGS}")


  set_source_files_properties( fctest_unified_memory_with_openacc.F90 PROPERTIES COMPILE_FLAGS ${ACC_Fortran_FLAGS_STR} )
  add_fctest( 
    TARGET atlas_test_unified_memory_with_openacc
    SOURCES fctest_unified_memory_with_openacc.F90
            fctest_unified_memory_with_openacc_cxx.cc
    LIBS atlas_f
    LINKER_LANGUAGE Fortran
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_RUN_NGPUS=1
  )
  target_link_libraries( atlas_test_unified_memory_with_openacc ${ACC_Fortran_FLAGS} CUDA::cudart )

  set_tests_properties( atlas_test_unified_memory_with_openacc PROPERTIES LABELS "gpu;acc")

  add_fctest( 
    TARGET atlas_test_connectivity_openacc
    SOURCES fctest_connectivity_openacc.F90
    LIBS atlas_f
    LINKER_LANGUAGE Fortran
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_RUN_NGPUS=1
  )
  target_link_libraries( atlas_test_connectivity_openacc ${ACC_Fortran_FLAGS} )
  set_target_properties( atlas_test_connectivity_openacc PROPERTIES COMPILE_FLAGS "${ACC_Fortran_FLAGS_STR}" )
  set_tests_properties ( atlas_test_connectivity_openacc PROPERTIES LABELS "gpu;acc")

endif()
