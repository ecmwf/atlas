# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_test( TARGET atlas_test_field_missingvalue
    SOURCES     test_field_missingvalue.cc
    LIBS        atlas
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_field_foreach
  SOURCES  test_field_foreach.cc
  LIBS     atlas
  OMP      4
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

if( HAVE_FCTEST )

  add_fctest( TARGET atlas_fctest_field
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_field.F90
    LIBS            atlas_f
    ENVIRONMENT     ${ATLAS_TEST_ENVIRONMENT}
  )
  
  add_fctest( TARGET atlas_fctest_field_wrap
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_field_wrap.F90
    LIBS            atlas_f
    ENVIRONMENT     ${ATLAS_TEST_ENVIRONMENT}
  )

  add_fctest( TARGET atlas_fctest_field_view
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_field_view.F90
    LIBS            atlas_f OpenMP::OpenMP_Fortran
    OMP             3
    ENVIRONMENT     ${ATLAS_TEST_ENVIRONMENT}
    CONDITION       atlas_HAVE_OMP_Fortran
  )

  add_fctest( TARGET atlas_fctest_field_host
    CONDITION       atlas_HAVE_GRIDTOOLS_STORAGE AND ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_field_host.F90
    LIBS            atlas_f
    ENVIRONMENT     ${ATLAS_TEST_ENVIRONMENT}
  )

  add_fctest( TARGET atlas_fctest_field_device
    CONDITION       atlas_HAVE_GRIDTOOLS_STORAGE AND ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA AND atlas_HAVE_ACC
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_field_gpu.F90 external_acc_routine.F90
    LIBS            atlas_f
    ENVIRONMENT     ${ATLAS_TEST_ENVIRONMENT} ATLAS_RUN_NGPUS=1
  )

  if( TARGET atlas_fctest_field_device )
    target_compile_options( atlas_fctest_field_device PUBLIC ${ACC_Fortran_FLAGS} )
    target_link_libraries( atlas_fctest_field_device ${ACC_Fortran_FLAGS} )
  endif()

endif()

