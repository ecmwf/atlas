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

#include "atlas/array/ArrayUtil.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace array {

namespace {
idx_t compute_allocated_size( idx_t size, idx_t alignment ) {
    idx_t div             = size / alignment;
    idx_t mod             = size % alignment;
    idx_t _allocated_size = div * alignment;
    if ( mod > 0 ) {
        _allocated_size += alignment;
    }
    return _allocated_size;
}
}  // namespace

ArraySpec::ArraySpec() : size_(), rank_(), allocated_size_(), contiguous_( true ), default_layout_( true ) {}

ArraySpec::ArraySpec( const ArrayShape& shape ) : ArraySpec( shape, ArrayAlignment() ) {}

ArraySpec::ArraySpec( const ArrayShape& shape, ArrayAlignment&& alignment ) {
    if ( int( alignment ) > 1 ) {
        ATLAS_NOTIMPLEMENTED;  // innermost dimension needs to be padded
    }

    rank_ = static_cast<int>( shape.size() );
    size_ = 1;
    shape_.resize( rank_ );
    strides_.resize( rank_ );
    layout_.resize( rank_ );
    for ( int j = rank_ - 1; j >= 0; --j ) {
        shape_[j]   = shape[j];
        strides_[j] = size_;
        layout_[j]  = j;
        size_ *= shape_[j];
    }
    allocated_size_ = compute_allocated_size( size_, alignment );
    contiguous_     = true;
    default_layout_ = true;

#ifdef ATLAS_HAVE_FORTRAN
    allocate_fortran_specs();
#endif
};

ArraySpec::ArraySpec( const ArrayShape& shape, const ArrayStrides& strides ) :
    ArraySpec( shape, strides, ArrayAlignment() ) {}

ArraySpec::ArraySpec( const ArrayShape& shape, const ArrayStrides& strides, ArrayAlignment&& alignment ) {
    ATLAS_ASSERT( shape.size() == strides.size(), "dimensions of shape and stride don't match" );

    rank_ = static_cast<int>( shape.size() );
    size_ = 1;
    shape_.resize( rank_ );
    strides_.resize( rank_ );
    layout_.resize( rank_ );
    for ( int j = rank_ - 1; j >= 0; --j ) {
        shape_[j]   = shape[j];
        strides_[j] = strides[j];
        layout_[j]  = j;
        size_ *= shape_[j];
    }
    allocated_size_ = compute_allocated_size( shape_[0] * strides_[0], alignment );
    contiguous_     = ( size_ == allocated_size_ );
    default_layout_ = true;

#ifdef ATLAS_HAVE_FORTRAN
    allocate_fortran_specs();
#endif
}

ArraySpec::ArraySpec( const ArrayShape& shape, const ArrayStrides& strides, const ArrayLayout& layout ) :
    ArraySpec( shape, strides, layout, ArrayAlignment() ) {}

ArraySpec::ArraySpec( const ArrayShape& shape, const ArrayStrides& strides, const ArrayLayout& layout,
                      ArrayAlignment&& alignment ) {
    ATLAS_ASSERT( shape.size() == strides.size(), "dimensions of shape and stride don't match" );

    rank_ = static_cast<int>( shape.size() );
    size_ = 1;
    shape_.resize( rank_ );
    strides_.resize( rank_ );
    layout_.resize( rank_ );
    default_layout_ = true;
    for ( int j = rank_ - 1; j >= 0; --j ) {
        shape_[j]   = shape[j];
        strides_[j] = strides[j];
        layout_[j]  = layout[j];
        size_ *= shape_[j];
        if ( layout_[j] != idx_t( j ) ) {
            default_layout_ = false;
        }
    }
    allocated_size_ = compute_allocated_size( shape_[layout_[0]] * strides_[layout_[0]], alignment );
    contiguous_     = ( size_ == allocated_size_ );

#ifdef ATLAS_HAVE_FORTRAN
    allocate_fortran_specs();
#endif
}

const std::vector<int>& ArraySpec::shapef() const {
    return shapef_;
}

const std::vector<int>& ArraySpec::stridesf() const {
    return stridesf_;
}

void ArraySpec::allocate_fortran_specs() {
    shapef_.resize( rank_ );
    for ( idx_t j = 0; j < rank_; ++j ) {
        shapef_[j] = shape_[rank_ - 1 - layout_[j]];
    }
    stridesf_.resize( strides_.size() );
    std::reverse_copy( strides_.begin(), strides_.end(), stridesf_.begin() );
}

}  // namespace array
}  // namespace atlas
