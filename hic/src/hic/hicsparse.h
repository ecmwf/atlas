/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#pragma once

#include <sstream>
#include <stdexcept>

#include "hic/hic_config.h"
#include "hic/hic_library_types.h"
#include "hic/hic_namespace_macro.h"

#if HIC_BACKEND_CUDA
#define HICSPARSE_BACKEND_PREFIX cu
#define HICSPARSE_BACKEND_PREFIX_CAPS CU
#include <cusparse.h>
#if CUSPARSE_VERSION < 12 * 1000 + 0 * 100 + 0
// the "Const" versions only appeared with CUDA 12
#define cusparseConstDnVecDescr_t cusparseDnVecDescr_t
#define cusparseConstDnMatDescr_t cusparseDnMatDescr_t
#define cusparseConstSpVecDescr_t cusparseSpVecDescr_t
#define cusparseConstSpMatDescr_t cusparseSpMatDescr_t

cusparseStatus_t cusparseCreateConstCsr(cusparseConstSpMatDescr_t* spMatDescr, int64_t rows, int64_t cols, int64_t nnz,
                                        const void* csrRowOffsets, const void* csrColInd, const void* csrValues,
                                        cusparseIndexType_t csrRowOffsetsType, cusparseIndexType_t csrColIndType,
                                        cusparseIndexBase_t idxBase, cudaDataType valueType) {
    return cusparseCreateCsr(spMatDescr, rows, cols, nnz, const_cast<void*>(csrRowOffsets),
                             const_cast<void*>(csrColInd), const_cast<void*>(csrValues), csrRowOffsetsType,
                             csrColIndType, idxBase, valueType);
}
cusparseStatus_t cusparseCreateConstDnVec(cusparseConstDnVecDescr_t* dnVecDescr, int64_t size, const void* values,
                                          cudaDataType valueType) {
    return cusparseCreateDnVec(dnVecDescr, size, const_cast<void*>(values), valueType);
}
cusparseStatus_t cusparseCreateConstDnMat(cusparseConstDnMatDescr_t* dnMatDescr, int64_t rows, int64_t cols, int64_t ld,
                                          const void* values, cudaDataType valueType, cusparseOrder_t order) {
    return cusparseCreateDnMat(dnMatDescr, rows, cols, ld, const_cast<void*>(values), valueType, order);
}
cusparseStatus_t cusparseCreateConstSpVec(cusparseConstSpVecDescr_t* spVecDescr, int64_t size, int64_t nnz,
                                          const void* indices, const void* values, cusparseIndexType_t idxType,
                                          cusparseIndexBase_t idxBase, cudaDataType valueType) {
    return cusparseCreateSpVec(spVecDescr, size, nnz, const_cast<void*>(indices), const_cast<void*>(values), idxType,
                               idxBase, valueType);
}
#endif
#elif HIC_BACKEND_HIP
#define HICSPARSE_BACKEND_PREFIX hip
#define HICSPARSE_BACKEND_PREFIX_CAPS HIP
#include <hipsparse/hipsparse.h>
#elif HIC_BACKEND_DUMMY
#define HICSPARSE_BACKEND_PREFIX dummy
#define HICSPARSE_BACKEND_PREFIX_CAPS dummy
#include "hic/hic_dummy/hicsparse_dummy.h"
#else
#error Unsupported hic backend. Please define HIC_BACKEND_CUDA or HIC_BACKEND_HIP or HIC_BACKEND_DUMMY
#endif

#define HIC_PREFIX hic
#define HIC_PREFIX_CAPS HIC
#define HIC_CONCAT_(A, B) A##B
#define HIC_CONCAT(A, B) HIC_CONCAT_(A, B)
#define HIC_SYMBOL(API) HIC_CONCAT(HIC_PREFIX, API)
#define HIC_SYMBOL_CAPS(API) HIC_CONCAT(HIC_PREFIX_CAPS, API)

#define HIC_TYPE(TYPE) using HIC_SYMBOL(TYPE) = HIC_CONCAT(HICSPARSE_BACKEND_PREFIX, TYPE);

#define HIC_FUNCTION(FUNCTION)                                                                     \
    template <typename... Args>                                                                    \
    inline auto HIC_SYMBOL(FUNCTION)(Args && ... args)                                             \
        -> decltype(HIC_CONCAT(HICSPARSE_BACKEND_PREFIX, FUNCTION)(std::forward<Args>(args)...)) { \
        return HIC_CONCAT(HICSPARSE_BACKEND_PREFIX, FUNCTION)(std::forward<Args>(args)...);        \
    }

#define HIC_VALUE(VALUE)                                                                          \
    constexpr decltype(HIC_CONCAT(HICSPARSE_BACKEND_PREFIX_CAPS, VALUE)) HIC_SYMBOL_CAPS(VALUE) = \
        HIC_CONCAT(HICSPARSE_BACKEND_PREFIX_CAPS, VALUE);

//------------------------------------------------
HIC_NAMESPACE_BEGIN
//------------------------------------------------

HIC_TYPE(sparseHandle_t)
HIC_TYPE(sparseStatus_t)
HIC_TYPE(sparseIndexType_t)
HIC_TYPE(sparseOrder_t)
HIC_TYPE(sparseConstDnVecDescr_t)
HIC_TYPE(sparseDnVecDescr_t)
HIC_TYPE(sparseConstDnMatDescr_t)
HIC_TYPE(sparseDnMatDescr_t)
HIC_TYPE(sparseConstSpVecDescr_t)
HIC_TYPE(sparseSpVecDescr_t)
HIC_TYPE(sparseConstSpMatDescr_t)
HIC_TYPE(sparseSpMatDescr_t)
HIC_TYPE(sparseSpMVAlg_t)
HIC_TYPE(sparseSpMMAlg_t)

HIC_FUNCTION(sparseGetErrorString)
HIC_FUNCTION(sparseCreate)
HIC_FUNCTION(sparseDestroy)
HIC_FUNCTION(sparseCreateConstDnVec)
HIC_FUNCTION(sparseCreateDnVec)
HIC_FUNCTION(sparseDestroyDnVec)
HIC_FUNCTION(sparseCreateConstDnMat)
HIC_FUNCTION(sparseCreateDnMat)
HIC_FUNCTION(sparseDestroyDnMat)
HIC_FUNCTION(sparseCreateConstSpVec)
HIC_FUNCTION(sparseCreateSpVec)
HIC_FUNCTION(sparseDestroySpVec)
HIC_FUNCTION(sparseCreateCsr)
HIC_FUNCTION(sparseCreateConstCsr)
HIC_FUNCTION(sparseDestroySpMat)
HIC_FUNCTION(sparseSpMV_bufferSize)
HIC_FUNCTION(sparseSpMV_preprocess)
HIC_FUNCTION(sparseSpMV)
HIC_FUNCTION(sparseSpMM_bufferSize)
HIC_FUNCTION(sparseSpMM_preprocess)
HIC_FUNCTION(sparseSpMM)

HIC_VALUE(SPARSE_STATUS_SUCCESS)
HIC_VALUE(SPARSE_ORDER_COL)
HIC_VALUE(SPARSE_ORDER_ROW)
HIC_VALUE(SPARSE_INDEX_32I)
HIC_VALUE(SPARSE_INDEX_64I)
HIC_VALUE(SPARSE_INDEX_BASE_ZERO)
HIC_VALUE(SPARSE_INDEX_BASE_ONE)
HIC_VALUE(SPARSE_SPMV_ALG_DEFAULT)
HIC_VALUE(SPARSE_SPMM_ALG_DEFAULT)
HIC_VALUE(SPARSE_OPERATION_NON_TRANSPOSE)
HIC_VALUE(SPARSE_OPERATION_TRANSPOSE)
HIC_VALUE(SPARSE_OPERATION_CONJUGATE_TRANSPOSE)

#if HIC_BACKEND_DUMMY
#define HICSPARSE_CALL(val)
#else
#define HICSPARSE_CALL(val) hicsparse_assert((val), #val, __FILE__, __LINE__)
#endif

inline void hicsparse_assert(hicsparseStatus_t status, const char* const func, const char* const file, const int line) {
    if (status != HICSPARSE_STATUS_SUCCESS) {
        std::ostringstream msg;
        msg << "HIC Runtime Error [code=" << status << "] at: " << file << " + " << line << " : " << func << "\n";
        msg << "  Reason: " << hicsparseGetErrorString(status);
        throw std::runtime_error(msg.str());
    }
}

//------------------------------------------------
HIC_NAMESPACE_END
//------------------------------------------------

#undef HIC_FUNCTION
#undef HIC_TYPE
#undef HIC_VALUE
#undef HIC_PREFIX
#undef HIC_PREFIX_CAPS
#undef HIC_CONCAT
#undef HIC_CONCAT_
#undef HIC_SYMBOL
#undef HIC_SYMBOL_CAPS
#undef HICSPARSE_BACKEND_PREFIX
#undef HICSPARSE_BACKEND_PREFIX_CAPS
