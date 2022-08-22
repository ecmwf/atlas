/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "RecordReader.h"

#include "atlas_io/Metadata.h"
#include "atlas_io/RecordItemReader.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

RecordReader::RecordReader(const Record::URI& ref): RecordReader(ref.path, ref.offset) {}

//---------------------------------------------------------------------------------------------------------------------

RecordReader::RecordReader(const std::string& path, uint64_t offset): session_{}, path_{path}, offset_{offset} {}

RecordReader::RecordReader(Stream stream, uint64_t offset): session_{}, stream_{stream}, path_{}, offset_{offset} {}

//---------------------------------------------------------------------------------------------------------------------

Record::URI RecordReader::uri() const {
    Record::URI _uri;
    _uri.path   = path_;
    _uri.offset = offset_;
    return _uri;
}

//---------------------------------------------------------------------------------------------------------------------

RecordItem::URI RecordReader::uri(const std::string& key) const {
    return RecordItem::URI(path_, offset_, key);
}

//---------------------------------------------------------------------------------------------------------------------

void RecordReader::trace(const std::string&, const char* file, int line, const char* func) {}

//---------------------------------------------------------------------------------------------------------------------

void RecordReader::wait(const std::string& key) {
    request(key).wait();
}

//---------------------------------------------------------------------------------------------------------------------

void RecordReader::wait() {
    // This can be optimized perhaps to overlap IO with decoding in multithreaded environment
    for (auto& pair : requests_) {
        auto& request = pair.second;
        request.wait();
    }
}

//---------------------------------------------------------------------------------------------------------------------

ReadRequest& RecordReader::request(const std::string& key) {
    return requests_.at(key);
}

//---------------------------------------------------------------------------------------------------------------------

Metadata RecordReader::metadata(const std::string& key) {
    Metadata metadata;
    RecordItemReader{uri(key)}.read(metadata);
    return metadata;
}

void RecordReader::checksum(bool b) {
    do_checksum_ = b;
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
