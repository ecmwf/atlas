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

#include "eckit/linalg/types.h"

namespace atlas {
namespace linalg {

//----------------------------------------------------------------------------------------------------------------------

template <typename Value, typename Index=eckit::linalg::Index>
class SparseMatrixView {
public:
    using value_type = Value;
    using index_type = Index;
    using size_type  = std::size_t;

    ~SparseMatrixView() = default;

    SparseMatrixView() = default;

    SparseMatrixView(
        size_type rows,
        size_type cols,
        size_type nnz,
        value_type const* value,
        index_type const* inner,
        index_type const* outer) :
        rows_(rows), cols_(cols), nnz_(nnz), outer_(outer), inner_(inner), value_(value) {
    }
    SparseMatrixView(const SparseMatrixView& other) = default;

    size_type rows() const { return rows_; }
    size_type cols() const { return cols_; }
    size_type nnz() const { return nnz_; }
    size_type outer_size() const { return rows() + 1; }
    size_type inner_size() const { return nnz(); }
    size_type value_size() const { return nnz(); }
    index_type const* outer() const { return outer_; }
    index_type const* inner() const { return inner_; }
    value_type const* value() const { return value_; }
    bool empty() const { return nnz() == 0; }

private:
    size_type rows_{0};
    size_type cols_{0};
    size_type nnz_{0};
    index_type const* outer_{nullptr};
    index_type const* inner_{nullptr};
    value_type const* value_{nullptr};
};

//----------------------------------------------------------------------------------------------------------------------
}  // namespace linalg
}  // namespace atlas
