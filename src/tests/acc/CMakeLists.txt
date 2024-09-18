# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_CUDA AND HAVE_TESTS AND HAVE_FCTEST AND HAVE_ACC AND OpenACC_Fortran_FOUND )

  add_fctest( 
    TARGET atlas_test_unified_memory_with_openacc
    SOURCES fctest_unified_memory_with_openacc.F90
            fctest_unified_memory_with_openacc_cxx.cc
    LIBS atlas_f OpenACC::OpenACC_Fortran CUDA::cudart
    LINKER_LANGUAGE Fortran
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_RUN_NGPUS=1
  )
  target_link_libraries( atlas_test_unified_memory_with_openacc ${ACC_Fortran_FLAGS} hic )
  set_tests_properties( atlas_test_unified_memory_with_openacc PROPERTIES LABELS "gpu;acc")

  add_fctest( 
    TARGET atlas_test_connectivity_openacc
    SOURCES fctest_connectivity_openacc.F90
    LIBS atlas_f OpenACC::OpenACC_Fortran
    LINKER_LANGUAGE Fortran
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_RUN_NGPUS=1
  )
  set_tests_properties ( atlas_test_connectivity_openacc PROPERTIES LABELS "gpu;acc")

endif()
