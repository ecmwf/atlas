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
#include <array>

#include "atlas_io/Data.h"
#include "atlas_io/Exceptions.h"
#include "atlas_io/Metadata.h"
#include "atlas_io/types/array/ArrayMetadata.h"
#include "atlas_io/types/array/ArrayReference.h"

namespace std {

//---------------------------------------------------------------------------------------------------------------------

template <typename T, size_t N>
void interprete(const std::vector<std::array<T, N>>& vector_of_array, atlas::io::ArrayReference& out) {
    using atlas::io::ArrayReference;
    out = ArrayReference{vector_of_array.front().data(), {vector_of_array.size(),N}};
}

//---------------------------------------------------------------------------------------------------------------------

template <typename T, size_t N>
void decode(const atlas::io::Metadata& m, const atlas::io::Data& encoded, std::vector<std::array<T, N>>& out) {
    atlas::io::ArrayMetadata array(m);
    if (array.datatype().kind() != atlas::io::ArrayMetadata::DataType::kind<T>()) {
        std::stringstream err;
        err << "Could not decode " << m.json() << " into std::vector<" << atlas::io::demangle<T>() << ">. "
            << "Incompatible datatype!";
        throw atlas::io::Exception(err.str(), Here());
    }
    if (array.rank() != 2) {
        std::stringstream err;
        err << "Could not decode " << m.json() << " into std::vector<std::array<" << atlas::io::demangle<T>() << "," << N << ">>. "
            << "Incompatible rank!";
        throw atlas::io::Exception(err.str(), Here());
    }
    if (array.shape(1) != N) {
        std::stringstream err;
        err << "Could not decode " << m.json() << " into std::vector<std::array<" << atlas::io::demangle<T>() << "," << N << ">>. "
            << "Incompatible size!";
        throw atlas::io::Exception(err.str(), Here());
    }
    const std::array<T,N>* data = static_cast<const std::array<T,N>*>(encoded.data());
    // std::copy(data, data + array.shape(0), out.begin());
    out.assign(data, data + array.shape(0));

}

//---------------------------------------------------------------------------------------------------------------------








}  // end namespace std
