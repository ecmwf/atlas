/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "AtlasIO.h"

#include <numeric>
#include <cstring>

#include "atlas/io/atlas-io.h"

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"

#include "atlas/runtime/Log.h"

namespace pluto {

template<typename T, typename Extents>
void interprete(pluto::mdspan<T,Extents> a, atlas::io::ArrayReference& out) {
    out = atlas::io::ArrayReference(a.data_handle(), atlas::io::make_datatype<int>(), atlas::io::ArrayMetadata::ArrayShape{a.extents()});
}

template<typename T, typename Extents>
void decode(const atlas::io::Metadata& metadata, const atlas::io::Data& data, pluto::mdspan<T,Extents>& out) {
    atlas::io::ArrayMetadata array(metadata);

    if (array.datatype().kind() != atlas::io::DataType::kind<T>()) {
        std::stringstream err;
        err << "Could not decode " << metadata.json() << " into mdspan with datatype " << atlas::io::DataType::str<T>() << "."
            << "Incompatible datatype!";
        throw atlas::io::Exception(err.str(), Here());
    }
    if (array.rank() != out.rank()) {
        std::stringstream err;
        err << "Could not decode " << metadata.json() << " into Array with rank " << out.rank() << "."
            << "Incompatible rank!";
        throw atlas::io::Exception(err.str(), Here());
    }

    // out.resize(array.shape());

    ATLAS_IO_ASSERT(out.size() >= array.size());

    ::memcpy(out.data_handle(), data, data.size());
}

}


namespace atlas {


void AtlasIO::read_mask(const std::string& mask_name, mdspan<int,dims<1>> mask) {
    ATLAS_TRACE();
    std::size_t size;
    atlas::io::RecordReader reader(mask_name);
    reader.read("size",size).wait();
    ATLAS_ASSERT(mask.size() == size);
    reader.read("mask", mask);
    reader.wait();
}


}
