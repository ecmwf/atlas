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

#include "atlas_io/Data.h"
#include "atlas_io/Metadata.h"
#include "atlas_io/types/array/ArrayMetadata.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

class ArrayReference : public ArrayMetadata {
public:
    ArrayReference() = default;

    ArrayReference(const void* data, DataType, const ArrayShape&);

    template <typename T>
    ArrayReference(const T* data, const ArrayShape& shape):
        ArrayMetadata(DataType::create<T>(), shape), data_(const_cast<T*>(data)) {}

    ArrayReference(ArrayReference&&);

    ArrayReference& operator=(ArrayReference&&);

    void* data() const { return data_; }

    friend void decode(const atlas::io::Metadata&, const atlas::io::Data&, ArrayReference&);

private:
    void* data_{nullptr};
};

//---------------------------------------------------------------------------------------------------------------------

size_t encode_metadata(const ArrayReference& value, atlas::io::Metadata& out);

void encode_data(const ArrayReference& value, atlas::io::Data& out);

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
