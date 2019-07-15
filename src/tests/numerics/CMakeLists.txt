# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_test( TARGET atlas_test_fvm_nabla
  SOURCES   test_fvm_nabla.cc
  LIBS      atlas ${OMP_CXX}
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

if( HAVE_FCTEST)
  add_fctest( TARGET atlas_fctest_fvm_nabla
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_fvm_nabla.F90
    LIBS            atlas_f ${OMP_Fortran}
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )
  
  add_fctest( TARGET atlas_fctest_ifs_setup
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_ifs_setup.F90
    LIBS            atlas_f
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )
endif()
