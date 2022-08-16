/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "ArrayReference.h"

#include "atlas_io/atlas_compat.h"
#include "atlas_io/detail/Assert.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

void encode_data(const ArrayReference& value, atlas::io::Data& out) {
    out = atlas::io::Data(value.data(), value.bytes());
}

//---------------------------------------------------------------------------------------------------------------------

namespace {
template <typename T>
void encode_metadata_value(const ArrayReference& value, atlas::io::Metadata& out) {
    ATLAS_IO_ASSERT(value.datatype() == make_datatype<T>());
    const T* array = reinterpret_cast<const T*>(value.data());
    std::vector<T> vector(value.size());
    std::copy(array, array + vector.size(), vector.begin());
    out.set("value", vector);
}
}  // namespace

size_t encode_metadata(const ArrayReference& value, atlas::io::Metadata& out) {
    auto bytes = encode_metadata(static_cast<const ArrayMetadata&>(value), out);
    if (value.rank() == 1 && value.size() <= 4) {
        auto kind      = value.datatype().kind();
        using DataType = ArrayReference::DataType;
        if (kind == DataType::kind<int>()) {
            encode_metadata_value<int>(value, out);
        }
        else if (kind == DataType::kind<long>()) {
            encode_metadata_value<long>(value, out);
        }
        else if (kind == DataType::kind<size_t>()) {
            encode_metadata_value<size_t>(value, out);
        }
        else if (kind == DataType::kind<float>()) {
            encode_metadata_value<float>(value, out);
        }
        else if (kind == DataType::kind<double>()) {
            encode_metadata_value<double>(value, out);
        }
    }
    return bytes;
}

//---------------------------------------------------------------------------------------------------------------------

ArrayReference::ArrayReference(const void* data, ArrayMetadata::DataType datatype,
                               const ArrayMetadata::ArrayShape& shape):
    ArrayMetadata(datatype, shape), data_(const_cast<void*>(data)) {}

//---------------------------------------------------------------------------------------------------------------------

ArrayReference::ArrayReference(ArrayReference&& other): ArrayMetadata(std::move(other)), data_(other.data_) {
    other.data_ = nullptr;
}

//---------------------------------------------------------------------------------------------------------------------

ArrayReference& ArrayReference::operator=(ArrayReference&& rhs) {
    ArrayMetadata::operator=(std::move(rhs));
    data_                  = rhs.data_;
    rhs.data_              = nullptr;
    return *this;
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
