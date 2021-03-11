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

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

void encode_data( const ArrayReference& value, atlas::io::Data& out ) {
    out = atlas::io::Data( value.data(), value.bytes() );
}

//---------------------------------------------------------------------------------------------------------------------

ArrayReference::ArrayReference( const void* data, ArrayMetadata::DataType datatype,
                                const ArrayMetadata::ArrayShape& shape ) :
    ArrayMetadata( datatype, shape ), data_( const_cast<void*>( data ) ) {}

//---------------------------------------------------------------------------------------------------------------------

ArrayReference::ArrayReference( ArrayReference&& other ) : ArrayMetadata( std::move( other ) ), data_( other.data_ ) {
    other.data_ = nullptr;
}

//---------------------------------------------------------------------------------------------------------------------

ArrayReference& ArrayReference::operator=( ArrayReference&& rhs ) {
    ArrayMetadata::operator=( std::move( rhs ) );
    data_                  = rhs.data_;
    rhs.data_              = nullptr;
    return *this;
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
