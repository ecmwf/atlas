/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas_io/FileStream.h"

#include "eckit/io/FileHandle.h"
#include "eckit/io/PooledHandle.h"

#include "atlas_io/Session.h"
#include "atlas_io/Trace.h"

namespace atlas {
namespace io {

namespace {

//---------------------------------------------------------------------------------------------------------------------

/// DataHandle that implements file IO for reading and writing
/// Main differences with eckit::FileHandle:
///   - Automatic opening and closing of file
///   - Openmode argument:
///       * read: for reading
///       * write: for writing, will overwrite existing file
///       * append: for appending implemented via write and seek to eof.
///   - ATLAS_IO_TRACE recording
class FileHandle : public eckit::FileHandle {
public:
    FileHandle(const eckit::PathName& path, char openmode): eckit::FileHandle(path, openmode == 'a' /*overwrite*/) {
        ATLAS_IO_TRACE("FileHandle::open(" + eckit::FileHandle::path() + "," + openmode + ")");
        if (openmode == 'r') {
            openForRead();
        }
        else if (openmode == 'w' || (openmode == 'a' && not path.exists())) {
            openForWrite(0);
        }
        else if (openmode == 'a') {
            openForWrite(path.size());
            seek(eckit::Offset(path.size()));
        }
    }

    void close() override {
        if (not closed_) {
            ATLAS_IO_TRACE("FileHandle::close(" + path() + ")");
            eckit::FileHandle::close();
            closed_ = true;
        }
    }

    FileHandle(const eckit::PathName& path, Mode openmode):
        FileHandle(path, openmode == Mode::read    ? 'r'
                         : openmode == Mode::write ? 'w'
                                                   : 'a') {}

    FileHandle(const eckit::PathName& path, const std::string& openmode): FileHandle(path, openmode[0]) {}

    ~FileHandle() override { close(); }

private:
    bool closed_{false};
};

//---------------------------------------------------------------------------------------------------------------------

/// DataHandle that implements file reading only.
/// Internally there is a registry of opened files which avoids
/// opening the same file multiple times.
/// Note that close() will not actually close the file when there
/// is another PooledHandle referencing the same file.
///
/// Main difference with eckit::PooledHandle
///   - Automatic opening and closing of file
///   - ATLAS_IO_TRACE recording
class PooledHandle : public eckit::PooledHandle {
public:
    PooledHandle(const eckit::PathName& path): eckit::PooledHandle(path), path_(path) {
        ATLAS_IO_TRACE("PooledHandle::open(" + path_.baseName() + ")");
        openForRead();
    }
    ~PooledHandle() override {
        ATLAS_IO_TRACE("PooledHandle::close(" + path_.baseName() + ")");
        close();
    }
    eckit::PathName path_;
};

}  // namespace

//---------------------------------------------------------------------------------------------------------------------

FileStream::FileStream(const eckit::PathName& path, char openmode):
    Stream([&path, &openmode]() -> eckit::DataHandle* {
        eckit::DataHandle* datahandle;
        if (openmode == 'r') {
            datahandle = new PooledHandle(path);
        }
        else {
            datahandle = new FileHandle(path, openmode);
        }
        return datahandle;
    }()) {
    if (openmode == 'r') {
        // Keep the PooledHandle alive until the end of active session
        Session::store(*this);
    }
}

FileStream::FileStream(const eckit::PathName& path, Mode openmode):
    FileStream(path, openmode == Mode::read    ? 'r'
                     : openmode == Mode::write ? 'w'
                                               : 'a') {}

FileStream::FileStream(const eckit::PathName& path, const std::string& openmode): FileStream(path, openmode[0]) {}

//---------------------------------------------------------------------------------------------------------------------

InputFileStream::InputFileStream(const eckit::PathName& path): FileStream(path, Mode::read) {}

//---------------------------------------------------------------------------------------------------------------------

OutputFileStream::OutputFileStream(const eckit::PathName& path, Mode openmode): FileStream(path, openmode) {}

OutputFileStream::OutputFileStream(const eckit::PathName& path, const std::string& openmode):
    FileStream(path, openmode) {}

OutputFileStream::OutputFileStream(const eckit::PathName& path, char openmode): FileStream(path, openmode) {}

void OutputFileStream::close() {
    datahandle().close();
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
