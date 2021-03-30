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

#include "atlas/array/ArrayShape.h"
#include "atlas/array/DataType.h"
#include "atlas/io/Metadata.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

class ArrayMetadata {
public:
    using ArrayShape = array::ArrayShape;
    using DataType   = array::DataType;

    static std::string type() { return "array"; }

public:
    ArrayMetadata();

    explicit ArrayMetadata( const Metadata& );

    explicit ArrayMetadata( const DataType&, const ArrayShape& );

    explicit ArrayMetadata( const ArrayMetadata& );

    ArrayMetadata( ArrayMetadata&& );

    ArrayMetadata& operator=( ArrayMetadata&& );

    int rank() const { return int( shape_.size() ); }

    int shape( int i ) const;

    const ArrayShape& shape() const { return shape_; }

    DataType datatype() const { return datatype_; }

    size_t size() const;

    size_t bytes() const { return size() * datatype_.size(); }

    friend size_t encode_metadata( const ArrayMetadata& value, atlas::io::Metadata& out );

private:
    ArrayShape shape_;
    DataType datatype_;
};

//---------------------------------------------------------------------------------------------------------------------

size_t encode_metadata( const ArrayMetadata& value, atlas::io::Metadata& out );

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
