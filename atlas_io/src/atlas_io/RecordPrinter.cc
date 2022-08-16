/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "RecordPrinter.h"

#include <sstream>

#include "atlas_io/Exceptions.h"
#include "atlas_io/FileStream.h"
#include "atlas_io/detail/Assert.h"
#include "atlas_io/print/JSONFormat.h"
#include "atlas_io/print/TableFormat.h"

#include "atlas_io/atlas_compat.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

RecordPrinter::RecordPrinter(const eckit::PathName& path, const eckit::Configuration& config):
    RecordPrinter(path, 0, config) {}

//---------------------------------------------------------------------------------------------------------------------

RecordPrinter::RecordPrinter(const eckit::PathName& path, const std::uint64_t offset,
                             const eckit::Configuration& config):
    RecordPrinter(Record::URI{path, offset}, config) {}

//---------------------------------------------------------------------------------------------------------------------

RecordPrinter::RecordPrinter(const Record::URI& ref, const eckit::Configuration& config):
    uri_(ref), record_(Session::record(ref.path, ref.offset)) {
    if (record_.empty()) {
        auto in = InputFileStream(uri_.path);
        in.seek(uri_.offset);
        record_.read(in, true);
        ATLAS_IO_ASSERT(not record_.empty());
    }

    config.get("format", options_.format);
    config.get("details", options_.details);

    // Check if format is supported
    {
        std::vector<std::string> supported_formats{"json", "yaml", "table"};
        bool format_supported{false};
        for (auto& supported_format : supported_formats) {
            if (options_.format == supported_format) {
                format_supported = true;
                break;
            }
        }
        if (not format_supported) {
            std::stringstream s;
            s << "Format '" + options_.format + "' not supported. Supported formats:";
            for (auto& supported_format : supported_formats) {
                s << "\n  - " << supported_format;
            }
            throw Exception(s.str(), Here());
        }
    }
}

//---------------------------------------------------------------------------------------------------------------------

void RecordPrinter::print(std::ostream& out) const {
    eckit::LocalConfiguration config;
    config.set("details", options_.details);
    if (options_.format == "json") {
        JSONFormat{uri_, config}.print(out);
    }
    else if (options_.format == "yaml") {
        JSONFormat{uri_, config}.print(out);
    }
    else if (options_.format == "table") {
        TableFormat{uri_, config}.print(out);
    }
    else {
        throw Exception("Cannot print record: Unrecognized format " + options_.format + ".", Here());
    }
}

//---------------------------------------------------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const RecordPrinter& info) {
    info.print(out);
    return out;
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
