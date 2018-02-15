/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an size_tergovernmental organisation
 * nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <stddef.h>
#include <vector>

#include "atlas/array/ArrayIdx.h"
#include "atlas/array/ArrayLayout.h"
#include "atlas/array/ArrayShape.h"
#include "atlas/array/ArrayStrides.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

class ArraySpec {
private:
    size_t size_;
    size_t rank_;
    size_t allocated_size_;
    ArrayShape shape_;
    ArrayStrides strides_;
    ArrayLayout layout_;
    ArrayAlignment alignment_;
    mutable std::vector<int> shapef_;
    mutable std::vector<int> stridesf_;
    bool contiguous_;
    bool default_layout_;

public:
    ArraySpec();
    ArraySpec( const ArrayShape& );
    ArraySpec( const ArrayShape&, const ArrayStrides& );
    ArraySpec( const ArrayShape&, const ArrayStrides&, const ArrayLayout& );
    ArraySpec( const ArrayShape&, ArrayAlignment&& );
    ArraySpec( const ArrayShape&, const ArrayStrides&, ArrayAlignment&& );
    ArraySpec( const ArrayShape&, const ArrayStrides&, const ArrayLayout&, ArrayAlignment&& );
    size_t allocatedSize() const { return allocated_size_; }
    size_t size() const { return size_; }
    size_t rank() const { return rank_; }
    const ArrayShape& shape() const { return shape_; }
    const ArrayAlignment& alignment() const { return alignment_; }
    const ArrayStrides& strides() const { return strides_; }
    const ArrayLayout& layout() const { return layout_; }
    const std::vector<int>& shapef() const;
    const std::vector<int>& stridesf() const;
    bool contiguous() const { return contiguous_; }
    bool hasDefaultLayout() const { return default_layout_; }
};

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
