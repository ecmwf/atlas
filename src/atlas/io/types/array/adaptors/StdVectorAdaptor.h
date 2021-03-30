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

#include <vector>

#include "atlas/io/Data.h"
#include "atlas/io/Exceptions.h"
#include "atlas/io/Metadata.h"
#include "atlas/io/types/array/ArrayMetadata.h"
#include "atlas/io/types/array/ArrayReference.h"
#include "atlas/runtime/Exception.h"

namespace std {

//---------------------------------------------------------------------------------------------------------------------

template <typename T>
void interprete( const std::vector<T>& vector, atlas::io::ArrayReference& out ) {
    using atlas::io::ArrayReference;
    out = ArrayReference{vector.data(), {int( vector.size() )}};
}

//---------------------------------------------------------------------------------------------------------------------

template <typename T>
void decode( const atlas::io::Metadata& m, const atlas::io::Data& encoded, std::vector<T>& out ) {
    atlas::io::ArrayMetadata array( m );
    if ( array.datatype().kind() != atlas::io::ArrayMetadata::DataType::kind<T>() ) {
        std::stringstream err;
        err << "Could not decode " << m.json() << " into std::vector<" << atlas::io::demangle<T>() << ">. "
            << "Incompatible datatype!";
        throw atlas::io::Exception( err.str(), Here() );
    }
    const T* data = static_cast<const T*>( encoded.data() );
    out.assign( data, data + array.size() );
}

//---------------------------------------------------------------------------------------------------------------------

}  // end namespace std
