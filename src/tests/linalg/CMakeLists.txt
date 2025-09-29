# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( atlas_HAVE_ATLAS_FUNCTIONSPACE )

ecbuild_add_test( TARGET atlas_test_linalg_sparse_gpu
    SOURCES test_linalg_sparse_gpu.cc
    LIBS    atlas
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_RUN_NGPUS=1
    CONDITION HAVE_CUDA OR HAVE_HIP
)
if( TARGET atlas_test_linalg_sparse_gpu )
  set_tests_properties( atlas_test_linalg_sparse_gpu PROPERTIES LABELS "gpu")
endif()

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
  LIBS    atlas hicsparse
  CONDITION HAVE_CUDA OR HAVE_HIP
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_RUN_NGPUS=1
)
if( TARGET atlas_test_linalg_sparse_matrix_gpu )
  set_tests_properties( atlas_test_linalg_sparse_matrix_gpu PROPERTIES LABELS "gpu")
endif()


ecbuild_add_test( TARGET atlas_test_linalg_dense
    SOURCES test_linalg_dense.cc
    LIBS    atlas
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

ecbuild_add_test( TARGET atlas_test_linalg_triplet
    SOURCES test_linalg_triplet.cc
    LIBS    atlas
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

endif()