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

#include <cstdint>
#include <memory>
#include <string>

#include "atlas_io/RecordItem.h"
#include "atlas_io/detail/Decoder.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

class ReadRequest {
public:
    ReadRequest(ReadRequest&& other);

    template <typename T>
    ReadRequest(Stream stream, size_t offset, const std::string& key, T& value):
        ReadRequest{stream, offset, key, new Decoder(value)} {}

    template <typename T>
    ReadRequest(const std::string& URI, T& value): ReadRequest{URI, new Decoder(value)} {}

    template <typename T>
    ReadRequest(const RecordItem::URI& URI, T& value): ReadRequest{URI.str(), value} {}

    ~ReadRequest();

    void read();

    void checksum();

    void decompress();

    void decode();

    void wait();

    void checksum(bool);

private:
    ReadRequest(const std::string& URI, Decoder* decoder);
    ReadRequest(Stream, size_t offset, const std::string& key, Decoder*);

    ReadRequest()                   = delete;
    ReadRequest(const ReadRequest&) = delete;

    Stream stream_;
    size_t offset_;
    std::string key_;
    std::string uri_;
    std::unique_ptr<Decoder> decoder_;
    std::unique_ptr<RecordItem> item_;
    bool do_checksum_{true};
    bool finished_{false};
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
