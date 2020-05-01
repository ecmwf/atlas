/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "eckit/geometry/SphereT.h"
#include "eckit/geometry/UnitSphere.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

struct DatumIFS {
    static constexpr double radius() { return 6371229.; }
};

struct DatumGRIB1 {
    static constexpr double radius() { return 6367470.; }
};

struct DatumWGS84SemiMajorAxis {
    static constexpr double radius() { return 6378137.; }
};

//------------------------------------------------------------------------------------------------------

typedef eckit::geometry::SphereT<DatumIFS> Earth;

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
