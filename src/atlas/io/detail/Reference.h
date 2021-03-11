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

#include "atlas/io/Data.h"
#include "atlas/io/Metadata.h"
#include "atlas/io/detail/sfinae.h"

#include "atlas/io/Exceptions.h"
namespace atlas {
namespace io {


template <typename T>
struct Reference {
    const T* ref;
    Reference( const T& r ) : ref( &r ) {}

    friend void encode_metadata( const Reference<T>& in, atlas::io::Metadata& metadata ) {
        if ( not sfinae::encode_metadata( *in.ref, metadata ) ) {
            throw NotEncodable( *in.ref );
        }
    }

    friend void encode_data( const Reference<T>& in, atlas::io::Data& out ) {
        if ( not sfinae::encode_data( *in.ref, out ) ) {
            throw NotEncodable( *in.ref );
        }
    }
};


}  // namespace io
}  // namespace atlas
