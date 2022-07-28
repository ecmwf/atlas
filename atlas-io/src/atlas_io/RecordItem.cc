/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "RecordItem.h"

#include "eckit/filesystem/URI.h"

#include "atlas_io/atlas_compat.h"
#include "atlas_io/detail/Assert.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

RecordItem::URI::URI(const std::string& _uri) {
    eckit::URI uri{_uri};

    ATLAS_IO_ASSERT(uri.scheme() == "file");
    ATLAS_IO_ASSERT(not uri.query("key").empty());

    path   = uri.path();
    offset = 0;
    if (not uri.query("offset").empty()) {
        offset = std::stoul(uri.query("offset"));
    }
    key = uri.query("key");
}

//---------------------------------------------------------------------------------------------------------------------

RecordItem::URI::URI(const std::string& _path, uint64_t _offset, const std::string& _key):
    path(_path), offset(_offset), key(_key) {}

//---------------------------------------------------------------------------------------------------------------------

std::string RecordItem::URI::str() const {
    eckit::URI uri("file", eckit::PathName(path));
    uri.query("offset", std::to_string(offset));
    uri.query("key", key);
    return uri.asRawString();
}

//---------------------------------------------------------------------------------------------------------------------

RecordItem::RecordItem(RecordItem&& other): metadata_(std::move(other.metadata_)), data_(std::move(other.data_)) {}

//---------------------------------------------------------------------------------------------------------------------

RecordItem::RecordItem(Metadata&& metadata, Data&& data): metadata_(new Metadata(metadata)), data_(std::move(data)) {}

//---------------------------------------------------------------------------------------------------------------------

const Data& RecordItem::data() const {
    return data_;
}

//---------------------------------------------------------------------------------------------------------------------

const Metadata& RecordItem::metadata() const {
    return *metadata_;
}

//---------------------------------------------------------------------------------------------------------------------

void RecordItem::metadata(const Metadata& m) {
    *metadata_ = m;
}

//---------------------------------------------------------------------------------------------------------------------

void RecordItem::data(Data&& d) {
    data_ = std::move(d);
}

//---------------------------------------------------------------------------------------------------------------------

bool RecordItem::empty() const {
    return metadata().empty();
}

//---------------------------------------------------------------------------------------------------------------------

void RecordItem::clear() {
    data_.clear();
    metadata_.reset(new atlas::io::Metadata());
}

//---------------------------------------------------------------------------------------------------------------------

void RecordItem::decompress() {
    ATLAS_IO_ASSERT(not empty());
    if (metadata().data.compressed()) {
        data_.decompress(metadata().data.compression(), metadata().data.size());
    }
    metadata_->data.compressed(false);
}

//---------------------------------------------------------------------------------------------------------------------

void RecordItem::compress() {
    ATLAS_IO_ASSERT(not empty());
    if (not metadata().data.compressed() && metadata().data.compression() != "none") {
        data_.compress(metadata().data.compression());
        metadata_->data.compressed(true);
    }
}

//---------------------------------------------------------------------------------------------------------------------

size_t encode_metadata(const RecordItem& in, Metadata& metadata) {
    metadata.set(in.metadata());
    return in.data().size();
}

//---------------------------------------------------------------------------------------------------------------------

void encode_data(const RecordItem& in, Data& out) {
    out.assign(in.data());
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
