/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <stddef.h>
#include <vector>
#include <string_view>
#include <string>

#include "atlas/array/ArrayIdx.h"
#include "atlas/array/ArrayLayout.h"
#include "atlas/array/ArrayShape.h"
#include "atlas/array/ArrayStrides.h"
#include "atlas/array/DataType.h"
#include "atlas/library/config.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

class ArraySpec {
private:
    size_t size_;
    idx_t rank_;
    size_t allocated_size_;
    DataType datatype_;
    ArrayShape shape_;
    ArrayStrides strides_;
    ArrayStrides device_strides_;
    ArrayLayout layout_;
    ArrayAlignment alignment_;
    std::vector<int> shapef_;
    std::vector<int> stridesf_;
    std::vector<int> device_stridesf_;
    bool contiguous_;
    bool default_layout_;

public:
    ArraySpec();
    ArraySpec(const ArrayShape&);
    ArraySpec(const ArrayShape&, const ArrayStrides&);
    ArraySpec(const ArrayShape&, const ArrayStrides&, const ArrayLayout&);
    ArraySpec(const ArrayShape&, const ArrayAlignment&);
    ArraySpec(const ArrayShape&, const ArrayStrides&, const ArrayAlignment&);
    ArraySpec(const ArrayShape&, const ArrayStrides&, const ArrayLayout&, const ArrayAlignment&);
    ArraySpec(DataType, const ArrayShape&);
    ArraySpec(DataType, const ArrayShape&, const ArrayStrides&);
    ArraySpec(DataType, const ArrayShape&, const ArrayStrides&, const ArrayLayout&);
    ArraySpec(DataType, const ArrayShape&, const ArrayAlignment&);
    ArraySpec(DataType, const ArrayShape&, const ArrayStrides&, const ArrayAlignment&);
    ArraySpec(DataType, const ArrayShape&, const ArrayStrides&, const ArrayLayout&, const ArrayAlignment&);
    size_t allocatedSize() const { return allocated_size_; }
    size_t size() const { return size_; }
    idx_t rank() const { return rank_; }
    DataType datatype() const { return datatype_; }
    const ArrayShape& shape() const { return shape_; }
    const ArrayAlignment& alignment() const { return alignment_; }
    const ArrayStrides& strides() const { return strides_; }
    const ArrayStrides& device_strides() const { return device_strides_; }
    const ArrayLayout& layout() const { return layout_; }
    const std::vector<int>& shapef() const;
    const std::vector<int>& stridesf() const;
    const std::vector<int>& device_stridesf() const;
    bool contiguous() const { return contiguous_; }
    bool hasDefaultLayout() const { return default_layout_; }

private:
    void allocate_fortran_specs();
};

// --------------------------------------------------------------------------------------------

class label {
public:
    label(std::string_view s) {
        previous_ = get();
        set(s);
    }
    ~label() {
        set(previous_);
    }
    static std::string_view get();
    static void set(std::string_view);
private:
    std::string previous_;
};

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
