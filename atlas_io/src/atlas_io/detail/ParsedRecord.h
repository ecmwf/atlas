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

#include <map>
#include <string>
#include <vector>

#include "atlas_io/Metadata.h"
#include "atlas_io/detail/RecordSections.h"

namespace atlas {
namespace io {

/// Low-level Record information container.
///
/// No big data is kept here, only metadata, and information
/// on how to retrieve data at a later stage
class ParsedRecord {
public:
    RecordHead head;                                           ///< head section of parsed record
    std::vector<std::string> keys;                             ///< Keys of items encoded in parsed record
    std::map<std::string, Metadata> items;                     ///< Items encoded in parsed record
    std::vector<RecordDataIndexSection::Entry> data_sections;  ///< Description of data sections in parsed record

    /// The parse() function needs to be called during the reading of the record and
    /// completes the "items" through introspection of the "data_sections".
    /// It also computes uncompressed data size using available metadata in the "items"
    void parse();
};

}  // namespace io
}  // namespace atlas
