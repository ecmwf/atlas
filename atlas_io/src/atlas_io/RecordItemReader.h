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

#include <string>

#include "atlas_io/Record.h"
#include "atlas_io/RecordItem.h"
#include "atlas_io/Stream.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

class RecordItemReader {
public:
    RecordItemReader(Stream, size_t offset, const std::string& key);

    RecordItemReader(Stream, const std::string& key);

    RecordItemReader(const std::string& uri);

    void read(RecordItem& item);

    void read(Metadata&, bool follow_links = true);

    void read(Metadata&, Data&);

private:
    RecordItemReader(const std::string& ref, const std::string& uri);

    Stream in_;
    Record record_;

    std::string ref_{};  // directory to which relative URI's are evaluated
    RecordItem::URI uri_;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
