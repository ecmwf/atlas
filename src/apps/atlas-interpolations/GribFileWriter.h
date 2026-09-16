#pragma once

#include <fstream>
#include <iostream>

#include <stdio.h>
#include "eccodes.h"

#include "eckit/filesystem/PathName.h"

#include "atlas/runtime/Exception.h"

#include "Grib.h"

// --------------------------------------------------------------------------------------------------------------

class GribFileWriter {
private:
    eckit::PathName path_;
    FILE* file_ = nullptr;
    int count_{0};
    int index_{1};

public:
    GribFileWriter(const std::string& path): path_(path) {}
    void open() {
        if (eckit::PathName(path_).exists()) {
            eckit::PathName(path_).unlink(false);
        }
        file_ = fopen(path_.asString().c_str(), "wb");
        if (!file_) {
            ATLAS_THROW_EXCEPTION("Could not open file " << path_ << " for writing.");
        }
    }
    ~GribFileWriter() {
        if (file_) {
            fclose(file_);
        }
    }
    void write(const Grib& grib) {
        /* get the coded message in a buffer */
        const void* buffer = nullptr;
        size_t size;
        if (codes_get_message(grib.handle(), &buffer, &size) != 0) {
            ATLAS_THROW_EXCEPTION("Could not get message from grib handle");
        }
        /* write the buffer in a file*/
        if (fwrite(buffer, 1, size, file_) != size) {
            ATLAS_THROW_EXCEPTION("Could write grib to file " << path_ << ".");
        }
    }
};

// --------------------------------------------------------------------------------------------------------------
