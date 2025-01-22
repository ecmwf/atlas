/*
 * (C) Copyright 2024 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/linalg/sparse/SparseMatrixMultiply_HicSparse.h"

#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"

#include "hic/hic.h"
#include "hic/hic_library_types.h"
#include "hic/hicsparse.h"

namespace {

class HicSparseHandleRAIIWrapper {
public:
    HicSparseHandleRAIIWrapper() { hicsparseCreate(&handle_); };
    ~HicSparseHandleRAIIWrapper() { hicsparseDestroy(handle_); }
    hicsparseHandle_t value() { return handle_; }
private:
    hicsparseHandle_t handle_;
};

hicsparseHandle_t getDefaultHicSparseHandle() {
    static auto handle = HicSparseHandleRAIIWrapper();
    return handle.value();
}

template<typename T>
constexpr hicsparseIndexType_t getHicsparseIndexType() {
    using base_type = std::remove_const_t<T>;
    if constexpr (std::is_same_v<base_type, int>) {
        return HICSPARSE_INDEX_32I;
    } else {
        static_assert(std::is_same_v<base_type, long>, "Unsupported index type");
        return HICSPARSE_INDEX_64I;
    }
}

template<typename T>
constexpr auto getHicsparseValueType() {
    using base_type = std::remove_const_t<T>;
    if constexpr (std::is_same_v<base_type, float>) {
        return HIC_R_32F;
    } else {
        static_assert(std::is_same_v<base_type, double>, "Unsupported value type");\
        return HIC_R_64F;
    }
}

template<atlas::linalg::Indexing IndexLayout, typename T>
hicsparseOrder_t getHicsparseOrder(const atlas::linalg::View<T, 2>& v) {
    constexpr int row_idx = (IndexLayout == atlas::linalg::Indexing::layout_left) ? 0 : 1;
    constexpr int col_idx = (IndexLayout == atlas::linalg::Indexing::layout_left) ? 1 : 0;

    if (v.stride(row_idx) == 1) {
        return HICSPARSE_ORDER_COL;
    } else if (v.stride(col_idx) == 1) {
        return HICSPARSE_ORDER_ROW;
    } else {
        atlas::throw_Exception("Unsupported dense matrix memory order", Here());
        return HICSPARSE_ORDER_COL;
    }
}

template<typename T>
int64_t getLeadingDimension(const atlas::linalg::View<T, 2>& v) {
    if (v.stride(0) == 1) {
        return v.stride(1);
    } else if (v.stride(1) == 1) {
        return v.stride(0);
    } else {
        atlas::throw_Exception("Unsupported dense matrix memory order", Here());
        return 0;
    }
}

}

namespace atlas {
namespace linalg {
namespace sparse {

template <typename Value>
void hsSpMV(const SparseMatrixView<Value>& W, const View<const Value, 1>& src, Value beta, View<Value, 1>& tgt) {
    // Assume that W, src, and tgt are all device views

    using Index = typename SparseMatrixView<Value>::index_type;

    ATLAS_ASSERT(src.shape(0) >= W.cols());
    ATLAS_ASSERT(tgt.shape(0) >= W.rows());

    auto handle = getDefaultHicSparseHandle();

    // Create a sparse matrix descriptor
    hicsparseConstSpMatDescr_t matA;
    HICSPARSE_CALL(hicsparseCreateConstCsr(
        &matA,
        W.rows(), W.cols(), W.nnz(),
        W.outer(),    // row_offsets
        W.inner(),    // column_indices
        W.value(),    // values
        getHicsparseIndexType<Index>(),
        getHicsparseIndexType<Index>(),  
        HICSPARSE_INDEX_BASE_ZERO,
        getHicsparseValueType<Value>()));

    // Create dense matrix descriptors
    hicsparseConstDnVecDescr_t vecX;
    HICSPARSE_CALL(hicsparseCreateConstDnVec(
        &vecX,
        static_cast<int64_t>(W.cols()), 
        src.data(),
        getHicsparseValueType<Value>()));

    hicsparseDnVecDescr_t vecY;
    HICSPARSE_CALL(hicsparseCreateDnVec(
        &vecY,
        W.rows(),
        tgt.data(),
        getHicsparseValueType<Value>()));

    const Value alpha = 1;

    // Determine buffer size
    size_t bufferSize = 0;
    HICSPARSE_CALL(hicsparseSpMV_bufferSize(
        handle,
        HICSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha,
        matA,
        vecX,
        &beta,
        vecY,
        getHicsparseValueType<Value>(),
        HICSPARSE_SPMV_ALG_DEFAULT,
        &bufferSize));

    // Allocate buffer
    char* buffer;
    HIC_CALL(hicMalloc(&buffer, bufferSize));
    
    // Perform SpMV
    HICSPARSE_CALL(hicsparseSpMV(
        handle,
        HICSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha,
        matA,
        vecX,
        &beta,
        vecY,
        getHicsparseValueType<Value>(),
        HICSPARSE_SPMV_ALG_DEFAULT,
        buffer));

    HIC_CALL(hicFree(buffer));
    HICSPARSE_CALL(hicsparseDestroyDnVec(vecX));
    HICSPARSE_CALL(hicsparseDestroyDnVec(vecY));
    HICSPARSE_CALL(hicsparseDestroySpMat(matA));

    HIC_CALL(hicDeviceSynchronize());
}


template <Indexing IndexLayout, typename Value>
void hsSpMM(const SparseMatrixView<Value>& W, const View<const Value, 2>& src, Value beta, View<Value, 2>& tgt) {
    // Assume that W, src, and tgt are all device views

    using Index = typename SparseMatrixView<Value>::index_type;

    constexpr int row_idx = (IndexLayout == Indexing::layout_left) ? 0 : 1;
    constexpr int col_idx = (IndexLayout == Indexing::layout_left) ? 1 : 0;

    ATLAS_ASSERT(src.shape(row_idx) >= W.cols());
    ATLAS_ASSERT(tgt.shape(row_idx) >= W.rows());
    ATLAS_ASSERT(src.shape(col_idx) == tgt.shape(col_idx));

    auto handle = getDefaultHicSparseHandle();

    // Create a sparse matrix descriptor
    hicsparseConstSpMatDescr_t matA;
    HICSPARSE_CALL(hicsparseCreateConstCsr(
        &matA,
        W.rows(), W.cols(), W.nnz(),
        W.outer(),    // row_offsets
        W.inner(),    // column_indices
        W.value(),    // values
        getHicsparseIndexType<Index>(),
        getHicsparseIndexType<Index>(),
        HICSPARSE_INDEX_BASE_ZERO,
        getHicsparseValueType<Value>()));

    // Create dense matrix descriptors
    hicsparseConstDnMatDescr_t matB;
    HICSPARSE_CALL(hicsparseCreateConstDnMat(
        &matB,
        W.cols(), src.shape(col_idx),
        getLeadingDimension(src),
        src.data(),
        getHicsparseValueType<Value>(),
        getHicsparseOrder<IndexLayout>(src)));

    hicsparseDnMatDescr_t matC;
    HICSPARSE_CALL(hicsparseCreateDnMat(
        &matC,
        W.rows(), tgt.shape(col_idx),
        getLeadingDimension(tgt),
        tgt.data(),
        getHicsparseValueType<Value>(),
        getHicsparseOrder<IndexLayout>(tgt)));

    const Value alpha = 1;

    // Determine buffer size
    size_t bufferSize = 0;
    HICSPARSE_CALL(hicsparseSpMM_bufferSize(
        handle,
        HICSPARSE_OPERATION_NON_TRANSPOSE,
        HICSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha,
        matA,
        matB,
        &beta,
        matC,
        getHicsparseValueType<Value>(),
        HICSPARSE_SPMM_ALG_DEFAULT,
        &bufferSize));

    // Allocate buffer
    char* buffer;
    HIC_CALL(hicMalloc(&buffer, bufferSize));

    // Perform SpMM
    HICSPARSE_CALL(hicsparseSpMM(
        handle,
        HICSPARSE_OPERATION_NON_TRANSPOSE,
        HICSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha,
        matA,
        matB,
        &beta,
        matC,
        getHicsparseValueType<Value>(),
        HICSPARSE_SPMM_ALG_DEFAULT,
        buffer));

    HIC_CALL(hicFree(buffer));
    HICSPARSE_CALL(hicsparseDestroyDnMat(matC));
    HICSPARSE_CALL(hicsparseDestroyDnMat(matB));
    HICSPARSE_CALL(hicsparseDestroySpMat(matA));

    HIC_CALL(hicDeviceSynchronize());
}

template <typename ScalarValue>
void SparseMatrixMultiply<backend::hicsparse, Indexing::layout_left, 1, ScalarValue, const ScalarValue, ScalarValue>::multiply(
    const SparseMatrixView<ScalarValue>& W, const View<ScalarValue const, 1>& src, View<ScalarValue, 1>& tgt, const Configuration&) {
    ScalarValue beta = 0;
    hsSpMV(W, src, beta, tgt);
}

template <typename ScalarValue>
void SparseMatrixMultiply<backend::hicsparse, Indexing::layout_left, 1, ScalarValue, const ScalarValue, ScalarValue>::multiply_add(
    const SparseMatrixView<ScalarValue>& W, const View<ScalarValue const, 1>& src, View<ScalarValue, 1>& tgt, const Configuration&) {
    ScalarValue beta = 1;
    hsSpMV(W, src, beta, tgt);
}

template <typename ScalarValue>
void SparseMatrixMultiply<backend::hicsparse, Indexing::layout_left, 2, ScalarValue, const ScalarValue, ScalarValue>::multiply(
    const SparseMatrixView<ScalarValue>& W, const View<ScalarValue const, 2>& src, View<ScalarValue, 2>& tgt, const Configuration&) {
    ScalarValue beta = 0;
    hsSpMM<Indexing::layout_left>(W, src, beta, tgt);
}

template <typename ScalarValue>
void SparseMatrixMultiply<backend::hicsparse, Indexing::layout_left, 2, ScalarValue, const ScalarValue, ScalarValue>::multiply_add(
    const SparseMatrixView<ScalarValue>& W, const View<ScalarValue const, 2>& src, View<ScalarValue, 2>& tgt, const Configuration&) {
    ScalarValue beta = 1;
    hsSpMM<Indexing::layout_left>(W, src, beta, tgt);
}

template <typename ScalarValue>
void SparseMatrixMultiply<backend::hicsparse, Indexing::layout_right, 1, ScalarValue, const ScalarValue, ScalarValue>::multiply(
    const SparseMatrixView<ScalarValue>& W, const View<ScalarValue const, 1>& src, View<ScalarValue, 1>& tgt, const Configuration&) {
    ScalarValue beta = 0;
    hsSpMV(W, src, beta, tgt);
}

template <typename ScalarValue>
void SparseMatrixMultiply<backend::hicsparse, Indexing::layout_right, 1, ScalarValue, const ScalarValue, ScalarValue>::multiply_add(
    const SparseMatrixView<ScalarValue>& W, const View<ScalarValue const, 1>& src, View<ScalarValue, 1>& tgt, const Configuration&) {
    ScalarValue beta = 1;
    hsSpMV(W, src, beta, tgt);
}

template <typename ScalarValue>
void SparseMatrixMultiply<backend::hicsparse, Indexing::layout_right, 2, ScalarValue, const ScalarValue, ScalarValue>::multiply(
    const SparseMatrixView<ScalarValue>& W, const View<ScalarValue const, 2>& src, View<ScalarValue, 2>& tgt, const Configuration&) {
    ScalarValue beta = 0;
    hsSpMM<Indexing::layout_right>(W, src, beta, tgt);
}

template <typename ScalarValue>
void SparseMatrixMultiply<backend::hicsparse, Indexing::layout_right, 2, ScalarValue, const ScalarValue, ScalarValue>::multiply_add(
    const SparseMatrixView<ScalarValue>& W, const View<ScalarValue const, 2>& src, View<ScalarValue, 2>& tgt, const Configuration&) {
    ScalarValue beta = 1;
    hsSpMM<Indexing::layout_right>(W, src, beta, tgt);
}

#define EXPLICIT_TEMPLATE_INSTANTIATION(TYPE)                                                                 \
    template struct SparseMatrixMultiply<backend::hicsparse, Indexing::layout_left,  1, TYPE, TYPE const, TYPE>; \
    template struct SparseMatrixMultiply<backend::hicsparse, Indexing::layout_left,  2, TYPE, TYPE const, TYPE>; \
    template struct SparseMatrixMultiply<backend::hicsparse, Indexing::layout_right, 1, TYPE, TYPE const, TYPE>; \
    template struct SparseMatrixMultiply<backend::hicsparse, Indexing::layout_right, 2, TYPE, TYPE const, TYPE>;

EXPLICIT_TEMPLATE_INSTANTIATION(double);
EXPLICIT_TEMPLATE_INSTANTIATION(float);

}  // namespace sparse
}  // namespace linalg
}  // namespace atlas