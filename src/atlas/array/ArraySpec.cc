/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>

#include "atlas/array/ArrayDataStore.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace array {

namespace {
size_t compute_aligned_size(size_t size, size_t alignment) {
    size_t div           = size / alignment;
    size_t mod           = size % alignment;
    size_t _aligned_size = div * alignment;
    if (mod > 0) {
        _aligned_size += alignment;
    }
    return _aligned_size;
}
}  // namespace

ArraySpec::ArraySpec():
    size_(), rank_(), datatype_(DataType::KIND_REAL64), allocated_size_(), contiguous_(true), default_layout_(true) {}

ArraySpec::ArraySpec(const ArrayShape& shape): ArraySpec(shape, ArrayAlignment()) {}
ArraySpec::ArraySpec(DataType datatype, const ArrayShape& shape): ArraySpec(shape) {
    datatype_ = datatype;
}

ArraySpec::ArraySpec(const ArrayShape& shape, const ArrayAlignment& alignment): datatype_(DataType::KIND_REAL64) {
    ArrayShape aligned_shape = shape;
    aligned_shape.back()     = compute_aligned_size(aligned_shape.back(), size_t(alignment));

    rank_           = static_cast<int>(shape.size());
    size_           = 1;
    allocated_size_ = 1;
    shape_.resize(rank_);
    strides_.resize(rank_);
    layout_.resize(rank_);
    for (int j = rank_ - 1; j >= 0; --j) {
        shape_[j]   = shape[j];
        strides_[j] = allocated_size_;
        layout_[j]  = j;
        size_ *= size_t(shape_[j]);
        allocated_size_ *= size_t(aligned_shape[j]);
    }
    ATLAS_ASSERT(allocated_size_ == compute_aligned_size(size_t(shape_[0]) * size_t(strides_[0]), alignment));
    contiguous_     = (size_ == allocated_size_);
    default_layout_ = true;
    alignment_      = alignment;

#ifdef ATLAS_HAVE_FORTRAN
    allocate_fortran_specs();
#endif
};

ArraySpec::ArraySpec(DataType datatype, const ArrayShape& shape, const ArrayAlignment& alignment):
    ArraySpec(shape, alignment) {
    datatype_ = datatype;
}

ArraySpec::ArraySpec(const ArrayShape& shape, const ArrayStrides& strides):
    ArraySpec(shape, strides, ArrayAlignment()) {}

ArraySpec::ArraySpec(const ArrayShape& shape, const ArrayStrides& strides, const ArrayAlignment& alignment):
    datatype_(DataType::KIND_REAL64) {
    ATLAS_ASSERT(shape.size() == strides.size(), "dimensions of shape and stride don't match");

    rank_ = static_cast<int>(shape.size());
    size_ = 1;
    shape_.resize(rank_);
    strides_.resize(rank_);
    layout_.resize(rank_);
    for (int j = rank_ - 1; j >= 0; --j) {
        shape_[j]   = shape[j];
        strides_[j] = strides[j];
        layout_[j]  = j;
        size_ *= shape_[j];
    }
    allocated_size_ = compute_aligned_size(shape_[0] * strides_[0], alignment);
    contiguous_     = (size_ == allocated_size_);
    default_layout_ = true;

#ifdef ATLAS_HAVE_FORTRAN
    allocate_fortran_specs();
#endif
}

ArraySpec::ArraySpec(DataType datatype, const ArrayShape& shape, const ArrayStrides& strides,
                     const ArrayAlignment& alignment):
    ArraySpec(shape, strides, alignment) {
    datatype_ = datatype;
}

ArraySpec::ArraySpec(const ArrayShape& shape, const ArrayStrides& strides, const ArrayLayout& layout):
    ArraySpec(shape, strides, layout, ArrayAlignment()) {}

ArraySpec::ArraySpec(DataType datatype, const ArrayShape& shape, const ArrayStrides& strides,
                     const ArrayLayout& layout):
    ArraySpec(shape, strides, layout) {
    datatype_ = datatype;
}

ArraySpec::ArraySpec(const ArrayShape& shape, const ArrayStrides& strides, const ArrayLayout& layout,
                     const ArrayAlignment& alignment):
    datatype_(DataType::KIND_REAL64) {
    ATLAS_ASSERT(shape.size() == strides.size(), "dimensions of shape and stride don't match");

    rank_ = static_cast<int>(shape.size());
    size_ = 1;
    shape_.resize(rank_);
    strides_.resize(rank_);
    layout_.resize(rank_);
    default_layout_ = true;
    for (int j = rank_ - 1; j >= 0; --j) {
        shape_[j]   = shape[j];
        strides_[j] = strides[j];
        layout_[j]  = layout[j];
        size_ *= shape_[j];
        if (layout_[j] != idx_t(j)) {
            default_layout_ = false;
        }
    }
    allocated_size_ = compute_aligned_size(shape_[layout_[0]] * strides_[layout_[0]], alignment);
    contiguous_     = (size_ == allocated_size_);

#ifdef ATLAS_HAVE_FORTRAN
    allocate_fortran_specs();
#endif
}

ArraySpec::ArraySpec(DataType datatype, const ArrayShape& shape, const ArrayStrides& strides, const ArrayLayout& layout,
                     const ArrayAlignment& alignment):
    ArraySpec(shape, strides, layout, alignment) {
    datatype_ = datatype;
}
const std::vector<int>& ArraySpec::shapef() const {
    return shapef_;
}

const std::vector<int>& ArraySpec::stridesf() const {
    return stridesf_;
}

void ArraySpec::allocate_fortran_specs() {
    shapef_.resize(rank_);
    stridesf_.resize(rank_);
    for (idx_t j = 0; j < rank_; ++j) {
        shapef_[j] = shape_[rank_ - 1 - layout_[j]];
        stridesf_[j] = strides_[rank_ -1 - layout_[j]];
    }
}

}  // namespace array
}  // namespace atlas
