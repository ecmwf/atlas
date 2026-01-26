/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "hic/hic_dummy/dummyShouldNotBeCalled.h"

#define DUMMY_SHOULD_NOT_BE_CALLED(SYMBOL) dummyShouldNotBeCalled(#SYMBOL)
#define DUMMY_FUNCTION(SYMBOL)                                 \
    template <typename... Args>                                \
    inline dummysparseStatus_t dummy##SYMBOL(Args&&... args) { \
        DUMMY_SHOULD_NOT_BE_CALLED(hic##SYMBOL);               \
        return dummysparseStatus_t{0};                         \
    }
#define DUMMY_VALUE(SYMBOL) constexpr int dummy##SYMBOL = 0;

namespace {

using dummysparseHandle_t          = int;
using dummysparseStatus_t          = int;
using dummysparseIndexType_t       = int;
using dummysparseOrder_t           = int;
using dummysparseConstDnVecDescr_t = void*;
using dummysparseDnVecDescr_t      = void*;
using dummysparseConstDnMatDescr_t = void*;
using dummysparseDnMatDescr_t      = void*;
using dummysparseConstSpVecDescr_t = void*;
using dummysparseSpVecDescr_t      = void*;
using dummysparseConstSpMatDescr_t = void*;
using dummysparseSpMatDescr_t      = void*;
using dummysparseSpMVAlg_t         = int;
using dummysparseSpMMAlg_t         = int;

DUMMY_FUNCTION(sparseGetErrorString)
DUMMY_FUNCTION(sparseCreate)
DUMMY_FUNCTION(sparseDestroy)
DUMMY_FUNCTION(sparseCreateConstDnVec)
DUMMY_FUNCTION(sparseCreateDnVec)
DUMMY_FUNCTION(sparseDestroyDnVec)
DUMMY_FUNCTION(sparseCreateConstDnMat)
DUMMY_FUNCTION(sparseCreateDnMat)
DUMMY_FUNCTION(sparseDestroyDnMat)
DUMMY_FUNCTION(sparseCreateConstSpVec)
DUMMY_FUNCTION(sparseCreateSpVec)
DUMMY_FUNCTION(sparseDestroySpVec)
DUMMY_FUNCTION(sparseCreateConstCsr)
DUMMY_FUNCTION(sparseDestroySpMat)
DUMMY_FUNCTION(sparseSpMV_bufferSize)
DUMMY_FUNCTION(sparseSpMV_preprocess)
DUMMY_FUNCTION(sparseSpMV)
DUMMY_FUNCTION(sparseSpMM_bufferSize)
DUMMY_FUNCTION(sparseSpMM_preprocess)
DUMMY_FUNCTION(sparseSpMM)

DUMMY_VALUE(SPARSE_STATUS_SUCCESS)
DUMMY_VALUE(SPARSE_ORDER_COL)
DUMMY_VALUE(SPARSE_ORDER_ROW)
DUMMY_VALUE(SPARSE_INDEX_32I)
DUMMY_VALUE(SPARSE_INDEX_64I)
DUMMY_VALUE(SPARSE_INDEX_BASE_ZERO)
DUMMY_VALUE(SPARSE_INDEX_BASE_ONE)
DUMMY_VALUE(SPARSE_SPMV_ALG_DEFAULT)
DUMMY_VALUE(SPARSE_SPMM_ALG_DEFAULT)
DUMMY_VALUE(SPARSE_OPERATION_NON_TRANSPOSE)
DUMMY_VALUE(SPARSE_OPERATION_TRANSPOSE)
DUMMY_VALUE(SPARSE_OPERATION_CONJUGATE_TRANSPOSE)

}  // namespace

#undef DUMMY_FUNCTION
#undef DUMMY_VALUE
#undef DUMMY_SHOULD_NOT_BE_CALLED