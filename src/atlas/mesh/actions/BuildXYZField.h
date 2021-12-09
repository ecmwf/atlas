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

namespace atlas {
class Field;
class Mesh;
}  // namespace atlas

namespace atlas {
namespace mesh {
class Nodes;
}
}  // namespace atlas

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

/// Creates a XYZ field from the (lon,lat) field
class BuildXYZField {
public:
    explicit BuildXYZField(const std::string& name = "xyz", bool force_recompute = false);

    Field& operator()(Mesh&) const;
    Field& operator()(mesh::Nodes&) const;

private:
    std::string name_;
    bool force_recompute_;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
