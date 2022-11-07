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
#include <map>
#include <string>

#include "atlas_io/Metadata.h"
#include "atlas_io/ReadRequest.h"
#include "atlas_io/Session.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

class RecordReader {
public:
    RecordReader(const Record::URI& ref);

    RecordReader(const std::string& path, std::uint64_t offset = 0);

    RecordReader(Stream stream, std::uint64_t offset = 0);

    template <typename T>
    ReadRequest& read(const std::string& key, T& value) {
        trace("read(" + key + ")", __FILE__, __LINE__, __func__);

        if (stream_) {
            trace("stream", __FILE__, __LINE__, __func__);
            requests_.emplace(key, ReadRequest{stream_, offset_, key, value});
        }
        else {
            requests_.emplace(key, ReadRequest{uri(key), value});
        }
        if (do_checksum_ >= 0) {
            requests_.at(key).checksum(do_checksum_);
        }
        return requests_.at(key);
    }

    void wait(const std::string& key);

    void wait();

    ReadRequest& request(const std::string& key);

    Metadata metadata(const std::string& key);

    void checksum(bool);

private:
    Record::URI uri() const;

    RecordItem::URI uri(const std::string& key) const;

    void trace(const std::string&, const char* file, int line, const char* func);

private:
    Session session_;

    Stream stream_;

    std::map<std::string, ReadRequest> requests_;

    std::string path_;
    std::uint64_t offset_;

    int do_checksum_{-1};
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
