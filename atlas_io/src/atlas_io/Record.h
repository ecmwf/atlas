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
#include <vector>

#include "atlas_io/detail/Endian.h"
#include "atlas_io/detail/Time.h"
#include "atlas_io/detail/Version.h"

namespace atlas {
namespace io {

class Stream;
class ParsedRecord;
class Metadata;

//---------------------------------------------------------------------------------------------------------------------

class Record {
public:
    struct URI {
        std::string str() const;
        std::string path;
        std::uint64_t offset;
        URI() = default;
        URI(const std::string& _path, std::uint64_t _offset = 0): path(_path), offset(_offset) {}
        URI(const URI& other): path(other.path), offset(other.offset) {}
    };

private:
    std::shared_ptr<ParsedRecord> record_;

public:
    Record();

    Record(const Record&);

    bool empty() const;

    Record& read(Stream& in, bool verify_end = false);

    const Metadata& metadata(const std::string& key) const;

    Endian endian() const;

    Version version() const;

    Time time() const;

    std::uint64_t size() const;

    const std::vector<std::string>& keys() const;

    bool has(const std::string& key);

    operator const ParsedRecord&() const;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
