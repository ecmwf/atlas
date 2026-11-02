/*
 * (C) Copyright 2026- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "MaskIO.h"

#include <cstring>
#include <sstream>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/log/JSON.h"

#include "atlas/io/atlas-io.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"

#include "atlas/grid.h"

namespace atlas {

namespace {

constexpr char FORMAT[] = "atlas-mask";
constexpr int VERSION   = 1;

}  // namespace

void MaskIO::validate(const std::string& source, const std::vector<int>& values) {
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (values[i] != 0 && values[i] != 1) {
            ATLAS_THROW_EXCEPTION("Mask '" << source << "' has non-binary value " << values[i] << " at index " << i);
        }
    }
}

MaskData MaskIO::read(const std::string& path) {
    ATLAS_TRACE();

    io::InputFileStream stream(path);
    io::Record record;
    record.read(stream);

    const bool has_format  = record.has("format");
    const bool has_version = record.has("version");
    if (has_format != has_version) {
        ATLAS_THROW_EXCEPTION("Mask record '" << path << "' must contain both 'format' and 'version'");
    }

    MaskData data;
    std::string format;
    int version{0};
    io::RecordReader reader(path);
    if (has_format) {
        reader.read("format", format);
        reader.read("version", version);
    }
    reader.read("grid", data.grid);
    reader.read("mask", data.values);
    reader.wait();

    if (has_format) {
        if (format != FORMAT) {
            ATLAS_THROW_EXCEPTION("Mask record '" << path << "' has unsupported format '" << format << "'");
        }
        if (version != VERSION) {
            ATLAS_THROW_EXCEPTION("Mask record '" << path << "' has unsupported " << FORMAT << " version "
                                                    << version);
        }

        std::istringstream stream(data.grid);
        data.grid = Grid(util::Config(stream)).name();
    }

    if (data.grid.empty()) {
        ATLAS_THROW_EXCEPTION("Mask record '" << path << "' contains an empty grid name");
    }
    validate(path, data.values);
    return data;
}

void MaskIO::read(const std::string& path, array::mdspan<int, array::dims<1>> values) {
    const MaskData data = read(path);
    if (values.size() != data.values.size()) {
        ATLAS_THROW_EXCEPTION("Mask record '" << path << "' contains " << data.values.size() << " values, expected "
                                                << values.size());
    }
    std::memcpy(values.data_handle(), data.values.data(), data.values.size() * sizeof(int));
}

void MaskIO::write(const std::string& path, const MaskData& data) {
    if (data.grid.empty()) {
        ATLAS_THROW_EXCEPTION("Cannot write a mask with an empty grid name");
    }
    validate(path, data.values);

    eckit::LocalConfiguration compression;
    compression.set("compression", "lz4");

    io::RecordWriter writer;

    Grid grid(data.grid);
    writer.set("format", FORMAT);
    writer.set("version", VERSION);
    writer.set("grid", grid.spec().json(eckit::JSON::Formatting::compact()));
    writer.set("mask", io::ref(data.values), compression);
    writer.write(path, io::Mode::write);
}

}  // namespace atlas
