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

template <typename SourceValue, typename TargetValue>
void hsSpMV(const SparseMatrix& W, const View<SourceValue, 1>& src, TargetValue beta, View<TargetValue, 1>& tgt) {
    // Assume that src and tgt are device views

    ATLAS_ASSERT(src.shape(0) >= W.cols());
    ATLAS_ASSERT(tgt.shape(0) >= W.rows());

    // Check if W is on the device and if not, copy it to the device
    if (W.deviceNeedsUpdate()) {
        W.updateDevice();
    }

    // Create sparse library handle
    // todo: use singleton class for storing hicSparse library handle.
    hicsparseHandle_t handle;
    HICSPARSE_CALL(hicsparseCreate(&handle));

    // Create a sparse matrix descriptor
    hicsparseConstSpMatDescr_t matA;
    HICSPARSE_CALL(hicsparseCreateConstCsr(
        &matA,
        W.rows(), W.cols(), W.nonZeros(),
        W.device_outer(),    // row_offsets
        W.device_inner(),    // column_indices
        W.device_data(),     // values
        getHicsparseIndexType<SparseMatrix::Index>(),
        getHicsparseIndexType<SparseMatrix::Index>(),  
        HICSPARSE_INDEX_BASE_ZERO,
        getHicsparseValueType<SparseMatrix::Scalar>()));

    // Create dense matrix descriptors
    hicsparseConstDnVecDescr_t vecX;
    HICSPARSE_CALL(hicsparseCreateConstDnVec(
        &vecX,
        static_cast<int64_t>(W.cols()), 
        src.data(),
        getHicsparseValueType<typename View<SourceValue, 1>::value_type>()));

    hicsparseDnVecDescr_t vecY;
    HICSPARSE_CALL(hicsparseCreateDnVec(
        &vecY,
        W.rows(),
        tgt.data(),
        getHicsparseValueType<typename View<TargetValue, 1>::value_type>()));

    using ComputeType = typename View<TargetValue, 1>::value_type;
    constexpr auto compute_type = getHicsparseValueType<ComputeType>();

    ComputeType alpha = 1;

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
        compute_type,
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
        compute_type,
        HICSPARSE_SPMV_ALG_DEFAULT,
        buffer));

    HIC_CALL(hicFree(buffer));
    HICSPARSE_CALL(hicsparseDestroyDnVec(vecX));
    HICSPARSE_CALL(hicsparseDestroyDnVec(vecY));
    HICSPARSE_CALL(hicsparseDestroySpMat(matA));
    HICSPARSE_CALL(hicsparseDestroy(handle));

    HIC_CALL(hicDeviceSynchronize());
}


template <Indexing IndexLayout, typename SourceValue, typename TargetValue>
void hsSpMM(const SparseMatrix& W, const View<SourceValue, 2>& src, TargetValue beta, View<TargetValue, 2>& tgt) {
    // Assume that src and tgt are device views

    constexpr int row_idx = (IndexLayout == Indexing::layout_left) ? 0 : 1;
    constexpr int col_idx = (IndexLayout == Indexing::layout_left) ? 1 : 0;

    ATLAS_ASSERT(src.shape(row_idx) >= W.cols());
    ATLAS_ASSERT(tgt.shape(row_idx) >= W.rows());
    ATLAS_ASSERT(src.shape(col_idx) == tgt.shape(col_idx));

    // Check if W is on the device and if not, copy it to the device
    if (W.deviceNeedsUpdate()) {
        W.updateDevice();
    }

    // Create sparse library handle
    // todo: use singleton class for storing hicSparse library handle.
    hicsparseHandle_t handle;
    HICSPARSE_CALL(hicsparseCreate(&handle));

    // Create a sparse matrix descriptor
    hicsparseConstSpMatDescr_t matA;
    HICSPARSE_CALL(hicsparseCreateConstCsr(
        &matA,
        W.rows(), W.cols(), W.nonZeros(),
        W.device_outer(),    // row_offsets
        W.device_inner(),    // column_indices
        W.device_data(),     // values
        getHicsparseIndexType<SparseMatrix::Index>(),
        getHicsparseIndexType<SparseMatrix::Index>(),
        HICSPARSE_INDEX_BASE_ZERO,
        getHicsparseValueType<SparseMatrix::Scalar>()));

    // Create dense matrix descriptors
    hicsparseConstDnMatDescr_t matB;
    HICSPARSE_CALL(hicsparseCreateConstDnMat(
        &matB,
        W.cols(), src.shape(col_idx),
        getLeadingDimension(src),
        src.data(),
        getHicsparseValueType<typename View<SourceValue, 2>::value_type>(),
        getHicsparseOrder<IndexLayout>(src)));

    hicsparseDnMatDescr_t matC;
    HICSPARSE_CALL(hicsparseCreateDnMat(
        &matC,
        W.rows(), tgt.shape(col_idx),
        getLeadingDimension(tgt),
        tgt.data(),
        getHicsparseValueType<typename View<TargetValue, 2>::value_type>(),
        getHicsparseOrder<IndexLayout>(tgt)));

    using ComputeType = typename View<TargetValue, 2>::value_type;
    constexpr auto compute_type = getHicsparseValueType<ComputeType>();

    ComputeType alpha = 1;

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
        compute_type,
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
        compute_type,
        HICSPARSE_SPMM_ALG_DEFAULT,
        buffer));

    HIC_CALL(hicFree(buffer));
    HICSPARSE_CALL(hicsparseDestroyDnMat(matC));
    HICSPARSE_CALL(hicsparseDestroyDnMat(matB));
    HICSPARSE_CALL(hicsparseDestroySpMat(matA));
    HICSPARSE_CALL(hicsparseDestroy(handle));

    HIC_CALL(hicDeviceSynchronize());
}

void SparseMatrixMultiply<backend::hicsparse, Indexing::layout_left, 1, double const, double>::apply(
    const SparseMatrix& W, const View<double const, 1>& src, View<double, 1>& tgt, const Configuration&) {
    double beta = 0;
    hsSpMV(W, src, beta, tgt);
}

void SparseMatrixMultiply<backend::hicsparse, Indexing::layout_left, 2, double const, double>::apply(
    const SparseMatrix& W, const View<double const, 2>& src, View<double, 2>& tgt, const Configuration&) {
    double beta = 0;
    hsSpMM<Indexing::layout_left>(W, src, beta, tgt);
}

void SparseMatrixMultiply<backend::hicsparse, Indexing::layout_right, 1, double const, double>::apply(
    const SparseMatrix& W, const View<double const, 1>& src, View<double, 1>& tgt, const Configuration&) {
    double beta = 0;
    hsSpMV(W, src, beta, tgt);
}

void SparseMatrixMultiply<backend::hicsparse, Indexing::layout_right, 2, double const, double>::apply(
    const SparseMatrix& W, const View<double const, 2>& src, View<double, 2>& tgt, const Configuration&) {
    double beta = 0;
    hsSpMM<Indexing::layout_right>(W, src, beta, tgt);
}

}  // namespace sparse
}  // namespace linalg
}  // namespace atlas