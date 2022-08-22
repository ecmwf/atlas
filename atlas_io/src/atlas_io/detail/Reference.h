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
#include "atlas_io/detail/sfinae.h"

#include "atlas_io/Exceptions.h"
namespace atlas {
namespace io {


template <typename T>
struct Reference {
    const T* ref;
    Reference(const T& r): ref(&r) {}

    friend size_t encode_metadata(const Reference<T>& in, atlas::io::Metadata& metadata) {
        size_t size{0};
        if (not sfinae::encode_metadata(*in.ref, metadata, size)) {
            throw NotEncodable(*in.ref);
        }
        return size;
    }

    friend void encode_data(const Reference<T>& in, atlas::io::Data& out) {
        if (not sfinae::encode_data(*in.ref, out)) {
            throw NotEncodable(*in.ref);
        }
    }
};


}  // namespace io
}  // namespace atlas
