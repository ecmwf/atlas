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

#include <type_traits>

#include "eckit/linalg/SparseMatrix.h"

#include "atlas/linalg/sparse/SparseMatrixView.h"
#include "atlas/linalg/sparse/SparseMatrixStorage.h"

namespace atlas {
namespace linalg {

//----------------------------------------------------------------------------------------------------------------------

class EckitSparseMatrixNonOwningAllocator : public eckit::linalg::SparseMatrix::Allocator {
public:
    EckitSparseMatrixNonOwningAllocator(const SparseMatrixView<eckit::linalg::Scalar, eckit::linalg::Index>& view) : view_{view} {}

    eckit::linalg::SparseMatrix::Layout allocate(eckit::linalg::SparseMatrix::Shape& shape) override {
        eckit::linalg::SparseMatrix::Layout p;
        using Index      = eckit::linalg::Index;
        using OuterIndex = std::remove_pointer_t<decltype(p.outer_)>;
        using InnerIndex = std::remove_pointer_t<decltype(p.inner_)>;
        static_assert(sizeof(OuterIndex) == sizeof(Index));
        static_assert(sizeof(InnerIndex) == sizeof(Index));
        p.data_     = const_cast<eckit::linalg::Scalar*>(view_.value());
        p.outer_    = reinterpret_cast<OuterIndex*>(const_cast<Index*>(view_.outer()));
        p.inner_    = reinterpret_cast<InnerIndex*>(const_cast<Index*>(view_.inner()));
        shape.size_ = view_.nnz();
        shape.rows_ = view_.rows();
        shape.cols_ = view_.cols();
        return p;
    }

    void deallocate(eckit::linalg::SparseMatrix::Layout p, eckit::linalg::SparseMatrix::Shape) override {
        /* do nothing */
    }

    bool inSharedMemory() const override { return false; }

    void print(std::ostream& out) const override { out << "atlas::linalg::EckitSparseMatrixNonOwningAllocator[rows="<<view_.rows()<<",cols="<<view_.cols()<<",nnz="<<view_.nnz()<<"]"; }

private:
    SparseMatrixView<eckit::linalg::Scalar, eckit::linalg::Index> view_;
};

//----------------------------------------------------------------------------------------------------------------------

// Create new eckit sparsematrix from a SparseMatrixView. Warning, this matrix does not own the data!
inline eckit::linalg::SparseMatrix make_non_owning_eckit_sparse_matrix(const SparseMatrixView<eckit::linalg::Scalar, eckit::linalg::Index>& view) {
    eckit::linalg::SparseMatrix m(new EckitSparseMatrixNonOwningAllocator(view));
    return m;
}

// Create new eckit sparsematrix from a SparseMatrixStorage. Warning, this matrix does not own the data!
inline eckit::linalg::SparseMatrix make_non_owning_eckit_sparse_matrix(const SparseMatrixStorage& m) {
    return make_non_owning_eckit_sparse_matrix(
        make_host_view<eckit::linalg::Scalar, eckit::linalg::Index>(m) );
}

inline eckit::linalg::SparseMatrix make_eckit_sparse_matrix(const SparseMatrixView<eckit::linalg::Scalar, eckit::linalg::Index>& view) {
    eckit::linalg::SparseMatrix m(new EckitSparseMatrixNonOwningAllocator(view));
    return eckit::linalg::SparseMatrix{m}; // makes copy
}

// Create new eckit sparsematrix from a SparseMatrixStorage. Warning, this creates a copy!
inline eckit::linalg::SparseMatrix make_eckit_sparse_matrix(const SparseMatrixStorage& m) {
    return make_eckit_sparse_matrix(
        make_host_view<eckit::linalg::Scalar, eckit::linalg::Index>(m) );
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace linalg
}  // namespace atlas
