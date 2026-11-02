/*
 * (C) Copyright 2026- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/array/mdspan.h"

namespace atlas {

struct MaskData {
    std::string grid;
    std::vector<int> values;
};

class MaskIO {
public:
    static MaskData read(const std::string& path);
    static void read(const std::string& path, array::mdspan<int, array::dims<1>> values);
    static void write(const std::string& path, const MaskData& data);
    static void validate(const std::string& source, const std::vector<int>& values);
};

}  // namespace atlas
