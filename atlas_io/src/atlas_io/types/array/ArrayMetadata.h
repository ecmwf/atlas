/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <cstdint>
#include <string>
#include <type_traits>

#include "atlas_io/Metadata.h"
#include "atlas_io/detail/DataType.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

class ArrayShape : public std::vector<size_t> {
private:
    using Base = std::vector<size_t>;

public:
    ArrayShape() {}
    ArrayShape(Base&& base): Base(std::forward<Base>(base)) {}
    template <typename idx_t>
    ArrayShape(std::initializer_list<idx_t> list): Base(list.begin(), list.end()) {}
    template <typename idx_t>
    ArrayShape(idx_t data[], size_t size): Base(data, data + size) {}
    template <typename idx_t, std::size_t N>
    ArrayShape(const std::array<idx_t, N>& list): Base(list.begin(), list.end()) {}
    template <typename idx_t>
    ArrayShape(const std::vector<idx_t>& list): Base(list.begin(), list.end()) {}
    template <typename Int1, typename = std::enable_if_t<std::is_integral_v<Int1>>>
    ArrayShape(Int1 i) {
        resize(1);
        operator[](0) = i;
    }
    template <typename Int1, typename Int2, typename = std::enable_if_t<std::is_integral_v<Int1> && std::is_integral_v<Int2>>>
    ArrayShape(Int1 i, Int2 j) {
        resize(2);
        operator[](0) = i;
        operator[](1) = j;
    }
    template <typename Int1, typename Int2, typename Int3, typename = std::enable_if_t<std::is_integral_v<Int1> && std::is_integral_v<Int2> && std::is_integral_v<Int3>>>
    ArrayShape(Int1 i, Int2 j, Int3 k) {
        resize(3);
        operator[](0) = i;
        operator[](1) = j;
        operator[](2) = k;
    }
    template <typename Int1, typename Int2, typename Int3, typename Int4, typename = std::enable_if_t<std::is_integral_v<Int1> && std::is_integral_v<Int2> && std::is_integral_v<Int3> && std::is_integral_v<Int4>>>
    ArrayShape(Int1 i, Int2 j, Int3 k, Int4 l) {
        resize(4);
        operator[](0) = i;
        operator[](1) = j;
        operator[](2) = k;
        operator[](3) = l;
    }
};

//---------------------------------------------------------------------------------------------------------------------

class ArrayMetadata {
public:
    using ArrayShape = io::ArrayShape;
    using DataType   = io::DataType;

    static std::string type() { return "array"; }

public:
    ArrayMetadata();

    explicit ArrayMetadata(const Metadata&);

    explicit ArrayMetadata(const DataType&, const ArrayShape&);

    explicit ArrayMetadata(const ArrayMetadata&);

    ArrayMetadata(ArrayMetadata&&);

    ArrayMetadata& operator=(ArrayMetadata&&);

    int rank() const { return int(shape_.size()); }

    int shape(int i) const;

    const ArrayShape& shape() const { return shape_; }

    DataType datatype() const { return datatype_; }

    size_t size() const;

    size_t bytes() const { return size() * datatype_.size(); }

    friend size_t encode_metadata(const ArrayMetadata& value, atlas::io::Metadata& out);

private:
    ArrayShape shape_;
    DataType datatype_;
};

//---------------------------------------------------------------------------------------------------------------------

size_t encode_metadata(const ArrayMetadata& value, atlas::io::Metadata& out);

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
