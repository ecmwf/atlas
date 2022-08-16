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

#include "eckit/filesystem/PathName.h"

#include "atlas_io/Stream.h"

namespace eckit {
class DataHandle;
}

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

enum class Mode
{
    read,
    append,
    write,
};

//---------------------------------------------------------------------------------------------------------------------

class FileStream : public Stream {
public:
    FileStream(const eckit::PathName& path, Mode openmode);
    FileStream(const eckit::PathName& path, char openmode);
    FileStream(const eckit::PathName& path, const std::string& openmode);
};

//---------------------------------------------------------------------------------------------------------------------

class InputFileStream : public FileStream {
public:
    InputFileStream(const eckit::PathName& path);
};

//---------------------------------------------------------------------------------------------------------------------

class OutputFileStream : public FileStream {
public:
    OutputFileStream(const eckit::PathName& path, Mode openmode = Mode::write);
    OutputFileStream(const eckit::PathName& path, const std::string& openmode);
    OutputFileStream(const eckit::PathName& path, char openmode);
    void close();
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
