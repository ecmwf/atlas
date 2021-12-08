/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "ArrayMetadata.h"

#include <algorithm>
#include <functional>
#include <numeric>

#include "atlas/runtime/Exception.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

size_t encode_metadata(const ArrayMetadata& value, atlas::io::Metadata& out) {
    out.set("type", value.type());
    out.set("shape", value.shape_);
    out.set("datatype", value.datatype_.str());
    return value.bytes();
}

//---------------------------------------------------------------------------------------------------------------------

int ArrayMetadata::shape(int i) const {
    if (i >= rank()) {
        throw_OutOfRange("shape", i, rank());
    }
    return shape_[size_t(i)];
}

//---------------------------------------------------------------------------------------------------------------------

ArrayMetadata::ArrayMetadata(const Metadata& metadata):
    datatype_(DataType::KIND_REAL64) /* circumvent absense of default constructor */ {
    std::string encoded_type;
    ATLAS_ASSERT(metadata.get("type", encoded_type), "metadata is missing 'type'");
    ATLAS_ASSERT(encoded_type == type(), "metadata has unexpected type '" + encoded_type + "'");
    ATLAS_ASSERT(metadata.get("shape", shape_), "metadata is missing 'shape'");
    std::string datatype_str;
    ATLAS_ASSERT(metadata.get("datatype", datatype_str), "metadata is missing 'datatype'");
    datatype_ = DataType(datatype_str);
}

//---------------------------------------------------------------------------------------------------------------------

size_t ArrayMetadata::size() const {
    return size_t(std::accumulate(shape_.begin(), shape_.end(), 1, std::multiplies<int>()));
}

//---------------------------------------------------------------------------------------------------------------------

ArrayMetadata::ArrayMetadata():
    shape_(), datatype_(DataType::KIND_REAL64) /* circumvent absense of default constructor */ {}

//---------------------------------------------------------------------------------------------------------------------

ArrayMetadata::ArrayMetadata(const DataType& datatype, const ArrayShape& shape): shape_(shape), datatype_(datatype) {}

//---------------------------------------------------------------------------------------------------------------------

ArrayMetadata::ArrayMetadata(const ArrayMetadata& other): ArrayMetadata{other.datatype_, other.shape_} {}

//---------------------------------------------------------------------------------------------------------------------

ArrayMetadata::ArrayMetadata(ArrayMetadata&& other): shape_{std::move(other.shape_)}, datatype_{other.datatype_} {}

//---------------------------------------------------------------------------------------------------------------------

ArrayMetadata& ArrayMetadata::operator=(ArrayMetadata&& rhs) {
    shape_    = std::move(rhs.shape_);
    datatype_ = rhs.datatype_;
    return *this;
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
