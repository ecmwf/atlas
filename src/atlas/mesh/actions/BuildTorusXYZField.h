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

#include <string>

#include "atlas/domain/Domain.h"
#include "atlas/domain/detail/RectangularDomain.h"
#include "atlas/library/config.h"

namespace atlas {
class Field;
}

namespace atlas {
class Mesh;
namespace mesh {
class Nodes;
}
}  // namespace atlas

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

/// Creates a XYZ field from the (lon,lat) field
class BuildTorusXYZField {
public:
    explicit BuildTorusXYZField(const std::string& name = "xyz");

    Field& operator()(Mesh&, const atlas::Domain&, double r0, double r1, idx_t nx, idx_t ny) const;
    Field& operator()(mesh::Nodes&, const atlas::Domain&, double r0, double r1, idx_t nx, idx_t ny) const;

private:
    std::string name_;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
