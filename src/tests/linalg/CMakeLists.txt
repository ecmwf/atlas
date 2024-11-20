# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( atlas_HAVE_ATLAS_FUNCTIONSPACE )

ecbuild_add_test( TARGET atlas_test_linalg_sparse
    SOURCES test_linalg_sparse.cc
    LIBS    atlas
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

ecbuild_add_test( TARGET atlas_test_linalg_sparse_matrix
  SOURCES test_linalg_sparse_matrix.cc
  LIBS    atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_linalg_sparse_matrix_gpu
  SOURCES test_linalg_sparse_matrix_gpu.cc
  LIBS    atlas
  CONDITION atlas_HAVE_CUDA OR atlas_HAVE_HIP
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_linalg_dense
    SOURCES test_linalg_dense.cc
    LIBS    atlas
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

endif()