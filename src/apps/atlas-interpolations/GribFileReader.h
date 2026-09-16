#pragma once

#include <fstream>
#include <iostream>

#include <stdio.h>
#include "eccodes.h"

#include "Grib.h"
#include "eckit/filesystem/PathName.h"

// --------------------------------------------------------------------------------------------------------------

class GribFileReader {
public:
    GribFileReader(const std::string& path) {
        if (!eckit::PathName(path).exists()) {
            ATLAS_THROW_EXCEPTION("File " << path << " does not exist.");
        }
        file_ = fopen(path.c_str(), "rb");
        if (!file_) {
            ATLAS_THROW_EXCEPTION("Could not open file " << path << " for reading.");
        }
        if (codes_count_in_file(nullptr, file_, &count_)) {
            ATLAS_THROW_EXCEPTION("Could not get number of grib messages from file " << path << ".");
        };
        grib_ = Grib::from_file(file_);
    }
    ~GribFileReader() { close(); }
    const Grib& grib() const { return grib_; }
    bool next() {
        if (index_ == count_) {
            return false;
        }
        grib_ = Grib::from_file(file_);
        index_++;
        return true;
    }
    int count() const { return count_; }
    void close() {
        if (file_) {
            fclose(file_);
        }
    }

private:
    FILE* file_ = nullptr;
    Grib grib_;
    int count_{0};
    int index_{1};
};

// --------------------------------------------------------------------------------------------------------------
